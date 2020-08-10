
#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>

#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <cmath>
#include <cfloat>
#include <cassert>
#include "ReadWriteData.h"

#include <vnl/algo/vnl_symmetric_eigensystem.h>

#include <itkArray.h>
#include <itkMatrix.h>
#include <itkImageRegionIteratorWithIndex.h>

#include "itkTDistribution.h"
#include "itkMath.h"
#include "vnl/vnl_erf.h"

namespace ants
{
// for computing F distribution
extern "C" double dbetai_(double *x, double *pin, double *qin);

extern "C" double dgamma_(double *x);

typedef struct
  {
  double statVal;
  int index;
  } StatElement;

int smallerStatElem(StatElement * elem1, StatElement * elem2)
// comparison function for sorting
{
  if( elem1->statVal > elem2->statVal )
    {
    return 1;
    }
  else if( elem1->statVal < elem2->statVal )
    {
    return -1;
    }
  else
    {
    return 0;
    }
}

template <typename TImageType>
void ReadImage(itk::SmartPointer<TImageType> & target, const char *file, bool copy)
{
  //  std::cout << " reading b " << std::string(file) << std::endl;
  typedef itk::ImageFileReader<TImageType> readertype;
  typename readertype::Pointer reader = readertype::New();
  reader->SetFileName(file);
  reader->Update();
  if( !copy )
    {
    target = (reader->GetOutput() );
    }
  else
    {
    typedef itk::ImageRegionIteratorWithIndex<TImageType> Iterator;
    Iterator vfIter2( target,  target->GetLargestPossibleRegion() );
    for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
      {
      vfIter2.Set( reader->GetOutput()->GetPixel(vfIter2.GetIndex() ) );
      }
    }
}

#define GROUPALABEL   0
#define GROUPBLABEL   1

typedef struct
  {
  int randNum;
  int index;
  } PermElement;

void generatePerm(int length, int * genPerm);

// generate a permutation of numbers from 0 to length-1
void generatePermGroup(int * groupID, int lengthGroupA, int lengthGroupB, int * genGroupID);

// generate a permutation of group assignments

int smallerPermElem(PermElement * elem1, PermElement * elem2);

static int first = 0;

void generatePermGroup(int * groupID, int lengthGroupA, int lengthGroupB,
                       int * genGroupID)
// generate a permutation of group assignments
{
  int numSubjects = lengthGroupA + lengthGroupB;

  if( !first )
    {
    first = 1; srand(time(nullptr) );
    // cout << "generatePermGroup called" << endl;
    }
  int * newPerm = new int[numSubjects];
  generatePerm(numSubjects, newPerm);
  for( int i = 0; i < numSubjects; i++ )
    {
    genGroupID[i] = groupID[newPerm[i]];
    }

  delete [] newPerm;
}

void generatePerm(int length, int * genPerm)
{
  if( !first )
    {
    first = 1; srand(time(nullptr) );
    // cout << "generatePerm called" << endl;
    }
  PermElement * newPerm = new PermElement[length];
  int           cnt;
  for( cnt = 0; cnt < length; cnt++ )
    {
    newPerm[cnt].randNum = rand();
    newPerm[cnt].index = cnt;
    }
  qsort(newPerm, length, sizeof(PermElement),
        (int (*)(const void *, const void *) )smallerPermElem);
  for( cnt = 0; cnt < length; cnt++ )
    {
    genPerm[cnt] = newPerm[cnt].index;
    }
  delete [] newPerm;
}

int smallerPermElem(PermElement * elem1, PermElement * elem2)
{
  if( elem1->randNum > elem2->randNum )
    {
    return 1;
    }
  else if( elem1->randNum < elem2->randNum )
    {
    return -1;
    }
  else
    {
    return 0;
    }
}

double computeQuantile(int numObs, double * stat, double quantile)
// computes the value in the provided statistic for the given quantile value
//
// numObs = Number of Observations
// stat = double array[numObs ] contains the test statisticss
{
  static int first = 0;

  if( !first )
    {
    first = 1;
    srand(time(nullptr) );
    }

  StatElement * sortStat = new StatElement[numObs];
  for( int perm = 0; perm < numObs; perm++ )
    {
    sortStat[perm].statVal = stat[perm];
    sortStat[perm].index = perm;
    }

  // sort, smallest first
  qsort(sortStat, numObs, sizeof(StatElement),
        (int (*)(const void *, const void *) )smallerStatElem);

  // index at value
  double quantindex = (double) numObs * quantile;
  if( quantile == 1.0 )
    {
    quantindex = numObs;
    }

  double retval = stat[sortStat[(int) quantindex].index];

  delete [] sortStat;

  return retval;
}

void computePermStatPval(int numFeatures, int numPerms,
                         double * permStat, double * permStatPval)
// computes the Pval for all permutation statistics
// the p-val is computed as the percentual ordered rank over all permutations
//
// numFeatures = Number of scalar features per Subject
// numPerms = Number of Permutations
// permStat = double array[numPerms * numFeatures] contains the test statistics
// permStatPval = double array [numPerms * numFeatures] returns the p-val of the statistics
{
  static int first = 0;

  if( !first )
    {
    first = 1;
    srand(time(nullptr) );
    }

  int feat;
  int perm;

  StatElement * sortPermStat = new StatElement[numPerms];
  for( feat = 0; feat < numFeatures; feat++ )
    {
    for( perm = 0; perm < numPerms; perm++ )
      {
      sortPermStat[perm].statVal = permStat[perm * numFeatures + feat];
      sortPermStat[perm].index = perm;
      }

    // sort, smallest first
    qsort(sortPermStat, numPerms, sizeof(StatElement),
          (int (*)(const void *, const void *) )smallerStatElem);

    double curPval = 0;
    for( perm = 0; perm < numPerms; perm++ )
      {
      // percentual rank 0..1 -> cumulative probability -> p-val
      double nextPval = 1.0 - (double) (perm + 1) / (double) numPerms;

      int curIndex = sortPermStat[perm].index;

      if( (perm == 0) || (sortPermStat[perm].statVal != sortPermStat[perm - 1].statVal) )
        {
        // current value is different from previous value (or first value),
        // thus step up p-value
        curPval = nextPval;
        }

      permStatPval[curIndex * numFeatures + feat] = curPval;
      }
    }
  delete [] sortPermStat;
}

double decode_ieee_single( unsigned char *v, int natural_order)
{
  unsigned char *data = v;
  int            s, e;
  unsigned long  src;
  long           f;
  double         value = 0.0;

  if( natural_order )
    {
    src = ( (unsigned long)data[0] << 24)
      | ( (unsigned long)data[1] << 16)
      | ( (unsigned long)data[2] << 8)
      | ( (unsigned long)data[3]);
    }
  else
    {
    src = ( (unsigned long)data[3] << 24)
      | ( (unsigned long)data[2] << 16)
      | ( (unsigned long)data[1] << 8)
      | ( (unsigned long)data[0]);
    }

  s = (src & 0x80000000UL) >> 31;
  e = (src & 0x7F800000UL) >> 23;
  f = (src & 0x007FFFFFUL);

  if( e == 255 && f != 0 )
    {
/* NaN (Not a Number) */
    value = DBL_MAX;
    }
  else if( e == 255 && f == 0 && s == 1 )
    {
/* Negative infinity */
    value = -DBL_MAX;
    }
  else if( e == 255 && f == 0 && s == 0 )
    {
/* Positive infinity */
    value = DBL_MAX;
    }
  else if( e > 0 && e < 255 )
    {
/* Normal number */
    f += 0x00800000UL;
    if( s )
      {
      f = -f;
      }
    value = ldexp(f, e - 150);
    }
  else if( e == 0 && f != 0 )
    {
/* Denormal number */
    if( s )
      {
      f = -f;
      }
    value = ldexp(f, -149);
    }
  else if( e == 0 && f == 0 && s == 1 )
    {
/* Negative zero */
    value = 0;
    }
  else if( e == 0 && f == 0 && s == 0 )
    {
/* Positive zero */
    value = 0;
    }
  else
    {
/* Never happens */
    printf("s = %d, e = %d, f = %lu\n", s, e, f);
    assert(!"Woops, unhandled case in decode_ieee_single()");
    }

  return value;
}

double factorial( double x)
{
  if( x <= 1 )
    {
    return 1;
    }
  double fac = x;
  double n = fac - 1;
  while(  n >= 1 )
    {
    fac *= n;
    n = n - 1;
    }

  return fac;
}

double betadist(  double a, double b )
{
  double numer = factorial(  a - 1) * factorial(b - 1);
  double denom = factorial( a + b - 1);

  return numer / denom;
}

double TTest(int numSubjects,   int* groupLabel, double * featureValue )
{
  int    numSubjA = 0;
  int    numSubjB = 0;
  double meanA = 0, meanB = 0;

  //  unsigned int GROUPALABEL=0;
  // unsigned int GROUPBLABEL=1;
  for( int subj = 0; subj < numSubjects; subj++ )
    {
    if( groupLabel[subj] == GROUPALABEL )
      {
      numSubjA++;
      meanA += featureValue[subj];
      }
    else if( groupLabel[subj] == GROUPBLABEL )
      {
      numSubjB++;
      meanB += featureValue[subj];
      }
    else
      {
      std::cout << " group label " << groupLabel[subj] << " does not exist" << std::endl;
      }
    }
  meanA /= (float)numSubjA;
  meanB /= (float)numSubjB;

  double varA = 0, varB = 0;
  for( int subj = 0; subj < numSubjects; subj++ )
    {
    if( groupLabel[subj] == GROUPALABEL )
      {
      varA += (featureValue[subj] - meanA) * (featureValue[subj] - meanA);
      }
    else if( groupLabel[subj] == GROUPBLABEL )
      {
      varB += (featureValue[subj] - meanB) * (featureValue[subj] - meanB);
      }
    }

  float n1 = (float) numSubjA;
  float n2 = (float) numSubjB;
  varA /= (n1); // use n1 -1  for unbiased estimator, here assume normal distribution
  varB /= (n2); // use n2 - 1 " ... "
//  float sdA=sqrt(varA);
//  float sdB=sqrt(varB);

//  float df  =  n1 + n2 - 2;
// unequal vars
  float denom = varA / n1 + varB / n2;
  // for equal vars
  //    float var =  ( (n1-1.)*newvar1 + (n2-1.)*newvar2 ) / df;
  //    denom = var*(1.0/n1+1.0/n2);
  double tt = 0;
  if( denom > 0 )
    {
    tt =  (meanA - meanB) / sqrt(denom);
    }

  return tt;
}

template <unsigned int ImageDimension>
int StudentsTestOnImages(int argc, char *argv[])
{
  typedef float                                 PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typename ImageType::Pointer mask = nullptr;
//  ReadImage<ImageType>(mask, argv[1], false);

  unsigned int numSubjectsA = std::stoi(argv[3]);
  unsigned int numSubjectsB = std::stoi(argv[4]);
  unsigned int numSubjects = numSubjectsA + numSubjectsB;
  std::string  outname = std::string(argv[2]);
  unsigned int numvals = numSubjects;
  int*         groupLabel = new int[numSubjects];
  for( unsigned int i = 0; i < numSubjectsA; i++ )
    {
    groupLabel[i] = 0;
    }
  for( unsigned int i = numSubjectsA; i < numSubjects; i++ )
    {
    groupLabel[i] = 1;
    }
  double* feature = new double[numvals];
  for( unsigned int i = 0; i < numvals; i++ )
    {
    feature[i] = 0;
    }

  std::cout << " Numvals " << numvals << std::endl;
  // Get the image dimension
  std::string               fn = std::string(argv[5]);
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(fn.c_str(), itk::ImageIOFactory::FileModeEnum::ReadMode);
  imageIO->SetFileName(fn.c_str() );
  imageIO->ReadImageInformation();
  typename ImageType::SizeType size;
  typename ImageType::SpacingType spacing;
  typename ImageType::PointType origin;
  typename ImageType::DirectionType direction;
  std::vector<double> axis;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    size[i] = imageIO->GetDimensions(i);
//       if (size[i] !=  mask->GetLargestPossibleRegion().GetSize()[i])
    // {
    //  std::cout <<  " mask not same size as data !! " << std::endl;
    //  throw std::exception();
    // }

    spacing[i] = imageIO->GetSpacing(i);
    origin[i]  = imageIO->GetOrigin(i);
    axis = imageIO->GetDirection(i);
    for( unsigned j = 0; j < ImageDimension; j++ )
      {
      if( j < imageIO->GetNumberOfDimensions() )
        {
        direction[j][i] = axis[j];
        }
      else
        {
        direction[j][i] = 0.0;
        }
      }
    }
  std::cout << " size " << size << std::endl;
  typename ImageType::RegionType region;
  region.SetSize(size );
  // ORIENTATION ALERT. the code this replaced originally didn't
  // bother setting the origins even though the directions were
  // grabbed from the ImageIO.  I'm assuming that was supposed to
  // happen, and was left out as an oversight.
  typename ImageType::Pointer StatImage = AllocImage<ImageType>(region,
                                                                spacing,
                                                                origin,
                                                                direction,
                                                                0);
  typename ImageType::Pointer PImage = AllocImage<ImageType>(region,
                                                             spacing,
                                                             origin,
                                                             direction,
                                                             0);

//   unsigned int sizeofpixel=sizeof(PixelType);

  std::vector<typename ImageType::Pointer> imagestack;
  imagestack.resize(numvals);
  for( unsigned int j = 0; j < numvals; j++ )
    {
    std::string ifn = std::string(argv[5 + j]);
    std::cout << "reading " << ifn << std::endl;
    ReadImage<ImageType>(imagestack[j], ifn.c_str(), false);
    }

  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  Iterator      vfIter(PImage, PImage->GetLargestPossibleRegion() );
  unsigned long nvox = 1;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    nvox *= PImage->GetLargestPossibleRegion().GetSize()[i];
    }

  unsigned long ct = 0;
  unsigned long prog = nvox / 20;
  std::cout << " NVals " << numvals << " NSub " << numSubjects <<  std::endl;
  for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    typename ImageType::IndexType index = vfIter.GetIndex();
    for( unsigned int subj = 0; subj < numSubjects; subj++ )
      {
      feature[subj] = imagestack[subj]->GetPixel(index);
      }
    if( ct % prog == 0 )
      {
      std::cout << " % " << (float) ct / (float) nvox << std::endl;
      }

    double stat = 0;
//      if (mask->GetPixel(index) >= 0.5)
    if( true )
      {
      stat = TTest(numSubjects, groupLabel, feature);
      }
    ct++;
    StatImage->SetPixel(  index,   stat);
    }

  typedef itk::Statistics::TDistribution DistributionType;
  typename DistributionType::Pointer distributionFunction = DistributionType::New();

  WriteImage(StatImage, outname.c_str() );

  delete [] feature;
  delete [] groupLabel;

  return 1;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int StudentsTestOnImages( std::vector<std::string> args, std::ostream* out_stream = nullptr )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "StudentsTestOnImages" );

  int     argc = args.size();
  char* * argv = new char *[args.size() + 1];
  for( unsigned int i = 0; i < args.size(); ++i )
    {
    // allocate space for the string plus a null character
    argv[i] = new char[args[i].length() + 1];
    std::strncpy( argv[i], args[i].c_str(), args[i].length() );
    // place the null character in the end
    argv[i][args[i].length()] = '\0';
    }
  argv[argc] = 0;
  // class to automatically cleanup argv upon destruction
  class Cleanup_argv
  {
public:
    Cleanup_argv( char* * argv_, int argc_plus_one_ ) : argv( argv_ ), argc_plus_one( argc_plus_one_ )
    {
    }

    ~Cleanup_argv()
    {
      for( unsigned int i = 0; i < argc_plus_one; ++i )
        {
        delete[] argv[i];
        }
      delete[] argv;
    }

private:
    char* *      argv;
    unsigned int argc_plus_one;
  };
  Cleanup_argv cleanup_argv( argv, argc + 1 );

  // antscout->set_stream( out_stream );

  std::cout <<  " df     P = 0.05  P = 0.01   P = 0.001  " << std::endl;
  std::cout << " 1             12.71     63.66     636.61  " << std::endl;
  std::cout << " 2    4.30     9.92     31.60    " << std::endl;
  std::cout << " 3    3.18     5.84     12.92" << std::endl;
  std::cout << " 4    2.78     4.60     8.61" << std::endl;
  std::cout << " 5    2.57     4.03     6.87" << std::endl;
  std::cout << " 6     2.45     3.71     5.96" << std::endl;
  std::cout << " 7    2.36     3.50     5.41" << std::endl;
  std::cout << " 8    2.31     3.36     5.04" << std::endl;
  std::cout << " 9    2.26     3.25     4.78" << std::endl;
  std::cout << " 10    2.23     3.17     4.59" << std::endl;
  std::cout << " 15    2.13     2.95     4.07" << std::endl;
  std::cout << " 20    2.09     2.85     3.85" << std::endl;
  std::cout << " 30    2.04     2.75     3.65" << std::endl;
  std::cout << " 50    2.01     2.68     3.50" << std::endl;
  std::cout << " 100    1.98     2.63     3.39  " << std::endl;

  if( argc < 6 )
    {
    std::cout << "Usage: " << argv[0] <<  " ImageDimension  OutName NGroup1 NGroup2 ControlV1*   SubjectV1*   "
             << std::endl;
    std::cout << " Assume all images the same size " << std::endl;
    std::cout << " Writes out an F-Statistic image " << std::endl;
    std::cout <<  " \n example call \n  \n ";
    std::cout << argv[0] << "  2  TEST.nii.gz 4 8 FawtJandADCcon/*SUB.nii  FawtJandADCsub/*SUB.nii  \n ";
    return 1;
    }

  switch( std::stoi(argv[1]) )
    {
    case 2:
      {
      StudentsTestOnImages<2>(argc, argv);
      }
      break;
    case 3:
      {
      StudentsTestOnImages<3>(argc, argv);
      }
      break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }

  return 0;
}
} // namespace ants
