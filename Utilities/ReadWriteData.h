#ifndef __ReadWriteData_h_
#define __ReadWriteData_h_
#include <antsAllocImage.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "itkVector.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLabeledPointSetFileReader.h"
#include "itkLabeledPointSetFileWriter.h"
#include "itkImageIntensityAndGradientToPointSetFilter.h"
#include "itkWarpImageFilter.h"
// #include "itkInverseWarpImageFilter.h"
#include "itkAffineTransform.h"
#include "itkImageRegionIterator.h"
#include "itkResampleImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkLogTensorImageFilter.h"
#include "itkExpTensorImageFilter.h"
#include "itkCastImageFilter.h"
#include <sys/stat.h>

extern bool ANTSFileExists(const std::string & strFilename);

// Nifti stores DTI values in lower tri format but itk uses upper tri
// currently, nifti io does nothing to deal with this. if this changes
// the function below should be modified/eliminated.

#if 1                           // currrently unimplemented
template <class TImageType>
void NiftiDTICheck(itk::SmartPointer<TImageType> &, const char *, bool )
{
}

#else
template <class TImageType>
void NiftiDTICheck(itk::SmartPointer<TImageType> & target, const char *file, bool makeLower)
{
  typedef typename TImageType::PixelType PixType;

  return;
  // typedef itk::ImageFileWriter<TImageType> Writer;
  // typename Writer::Pointer writer = Writer::New();
  // writer->SetInput( target );
  // writer->SetFileName( "testdt.nii" );
  // writer->Update();

  // Check for nifti file
  std::string            filename = file;
  std::string::size_type pos1 = filename.find( ".nii" );
  std::string::size_type pos2 = filename.find( ".nia" );
  if( (pos1 == std::string::npos) && (pos2 == std::string::npos) )
    {
    return;
    }

  if( PixType::Length != 6 )
    {
    return;
    }

  std::cout << "Performing lower/upper triangular format check for Nifti DTI" << std::endl;

  // swap elements 2 and 3 for lower<->upper conversion
  itk::ImageRegionIteratorWithIndex<TImageType>
  iter(target, target->GetLargestPossibleRegion() );

  unsigned int looksLikeLower = 0;
  unsigned int looksLikeUpper = 0;
  unsigned int nBadVoxels = 0;
  unsigned int count = 0;

  unsigned int el2Neg = 0;
  unsigned int el3Neg = 0;

  while( !iter.IsAtEnd() )
    {
    bool isValid = true;
    for( unsigned int i = 0; i < 6; i++ )
      {
      if( iter.Get()[i] != iter.Get()[i] )
        {
        ++nBadVoxels;
        isValid = false;
        }
      }

    double el2 = iter.Get()[2];
    double el3 = iter.Get()[3];

    if( el2 < 0 )
      {
      ++el2Neg;
      }
    if( el3 < 0 )
      {
      ++el3Neg;
      }

    if( isValid )
      {
      if( el2 > el3 )
        {
        ++looksLikeLower;
        }
      else
        {
        ++looksLikeUpper;
        }
      }

    ++count;
    ++iter;
    }

  // std::cout << "Invalid: " << nBadVoxels << "/" << count << std::endl;
  // std::cout << "Lower: " << looksLikeLower << ", Upper: " << looksLikeUpper << std::endl;
  // std::cout << "el2Neg: " << el2Neg << ", el3Neg: " << el3Neg << std::endl;

  if( ( (looksLikeUpper > looksLikeLower) && makeLower) || ( (looksLikeLower > looksLikeUpper) && !makeLower) )
    {
    std::cout << "Performing lower/upper triangular format swap for Nifti DTI" << std::endl;

    iter.GoToBegin();
    while( !iter.IsAtEnd() )
      {
      PixType pix = iter.Get();
      typename PixType::ValueType temp;
      temp = pix[2];
      pix[2] = pix[3];
      pix[3] = temp;
      iter.Set(pix);
      ++iter;
      }
    }
}

#endif

template <class TImageType>
void ReadTensorImage(itk::SmartPointer<TImageType> & target, const char *file, bool takelog = true)
{
  if( !ANTSFileExists(std::string(file) ) )
    {
    std::cerr << " file " << std::string(file) << " does not exist . " << std::endl;
    return;
    }

  typedef TImageType                      ImageType;
  typedef itk::ImageFileReader<ImageType> FileSourceType;

  typedef itk::LogTensorImageFilter<ImageType, ImageType> LogFilterType;
  typename FileSourceType::Pointer reffilter = ITK_NULLPTR;
  if( file[0] == '0' && file[1] == 'x' )
    {
    void* ptr;
    sscanf(file, "%p", (void **)&ptr);
    target = *( static_cast<typename TImageType::Pointer *>( ptr ) );
    }
  else
    {
    // Read the image files begin

    reffilter = FileSourceType::New();
    reffilter->SetFileName( file );
    try
      {
      reffilter->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cerr << "Exception caught during reference file reading " << std::endl;
      std::cerr << e << " file " << file << std::endl;
      target = ITK_NULLPTR;
      return;
      }

    target = reffilter->GetOutput();
    }

  //NiftiDTICheck<ImageType>(target, file, false);

  if( takelog )
    {
    typename LogFilterType::Pointer logFilter = LogFilterType::New();
    logFilter->SetInput( reffilter->GetOutput() );
    try
      {
      logFilter->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cerr << "Exception caught during log tensor filter " << std::endl;
      std::cerr << e << " file " << file << std::endl;
      target = ITK_NULLPTR;
      return;
      }
    target = logFilter->GetOutput();
    std::cout << "Returning Log(D) for log-euclidean math ops" << std::endl;
    }

}

template <class TImageType>
// void ReadImage(typename TImageType::Pointer target, const char *file)
bool ReadImage(itk::SmartPointer<TImageType> & target, const char *file)
{
  enum { ImageDimension = TImageType::ImageDimension };
  if( std::string(file).length() < 3 )
    {
    target = ITK_NULLPTR;
    return false;
    }

  std::string comparetype1 = std::string( "0x" );
  std::string comparetype2 = std::string( file );
  comparetype2 = comparetype2.substr( 0, 2 );
  // Read the image files begin
  if(  comparetype1 == comparetype2  )
    {
    typedef TImageType RImageType;
    void* ptr;
    sscanf(file, "%p", (void **)&ptr);
    typename RImageType::Pointer Rimage = *( static_cast<typename RImageType::Pointer *>( ptr ) );
    /** more needs to be done here to cast the pointer to an image type --- this is a work-around */
    typedef itk::CastImageFilter<RImageType, TImageType> CastFilterType;
    typename CastFilterType::Pointer caster = CastFilterType::New();
    caster->SetInput( Rimage );
    caster->UpdateLargestPossibleRegion();
    target = caster->GetOutput();
    }
  else
    {
    if( !ANTSFileExists(std::string(file) ) )
      {
      std::cerr << " file " << std::string(file) << " does not exist . " << std::endl; target = ITK_NULLPTR;
      return false;
      }
    typedef TImageType                      ImageType;
    typedef itk::ImageFileReader<ImageType> FileSourceType;

    typename FileSourceType::Pointer reffilter = FileSourceType::New();
    reffilter->SetFileName( file );
    try
      {
      reffilter->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cerr << "Exception caught during reference file reading " << std::endl;
      std::cerr << e << " file " << file << std::endl;
      target = ITK_NULLPTR;
      std::exception();
      return false;
      }

    // typename ImageType::DirectionType dir;
    // dir.SetIdentity();
    //  reffilter->GetOutput()->SetDirection(dir);

    // std::cout << " setting pointer " << std::endl;
    target = reffilter->GetOutput();
    }
  return true;
}

template <class ImageType>
typename ImageType::Pointer ReadImage(char* fn )
{
  // Read the image files begin
  typedef itk::ImageFileReader<ImageType> FileSourceType;

  typename FileSourceType::Pointer reffilter = FileSourceType::New();
  reffilter->SetFileName( fn );
  try
    {
    reffilter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception caught during image reference file reading " << std::endl;
    std::cerr << e << std::endl;
    return ITK_NULLPTR;
    }

  //typename ImageType::DirectionType dir;
  //dir.SetIdentity();
  //  reffilter->GetOutput()->SetDirection(dir);

  typename ImageType::Pointer target = reffilter->GetOutput();
  // if (reffilter->GetImageIO->GetNumberOfComponents() == 6)
  // NiftiDTICheck<ImageType>(target,fn);

  return target;
}

template <class ImageType>
typename ImageType::Pointer ReadTensorImage(char* fn, bool takelog = true )
{
  // Read the image files begin
  typedef itk::ImageFileReader<ImageType>                 FileSourceType;
  typedef itk::LogTensorImageFilter<ImageType, ImageType> LogFilterType;

  typename FileSourceType::Pointer reffilter = FileSourceType::New();
  reffilter->SetFileName( fn );
  try
    {
    reffilter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception caught during tensor image reference file reading " << std::endl;
    std::cerr << e << std::endl;
    return NULL;
    }

  typename ImageType::Pointer target = reffilter->GetOutput();

  NiftiDTICheck<ImageType>(target, fn, false);

  if( takelog )
    {
    typename LogFilterType::Pointer logFilter = LogFilterType::New();
    logFilter->SetInput( target->GetOutput() );
    logFilter->Update();
    target = logFilter->GetOutput();
    }

  return target;
}

template <class TPointSet>
// void ReadImage(typename TPointSet::Pointer target, const char *file)
bool ReadLabeledPointSet( itk::SmartPointer<TPointSet> & target, const char *file,
  bool boundaryPointsOnly = false, float samplingPercentage = 1.0 )
{
  if( std::string( file ).length() < 3 )
    {
    target = ITK_NULLPTR;
    return false;
    }

  if( !ANTSFileExists( std::string( file ) ) )
    {
    std::cerr << " file " << std::string( file ) << " does not exist . " << std::endl;
    return false;
    }

  // Read the image files begin
  typedef itk::LabeledPointSetFileReader<TPointSet>   FileSourceType;
  typename FileSourceType::Pointer reffilter = FileSourceType::New();
  reffilter->SetFileName( file );
  reffilter->SetExtractBoundaryPoints( boundaryPointsOnly );
  reffilter->SetRandomPercentage( samplingPercentage );
  try
    {
    reffilter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception caught during point set reference file reading " << std::endl;
    std::cerr << e << std::endl;
    return false;
    }

  target = reffilter->GetOutput();

  return true;
}

template <class TImage, class TMask, class TPointSet>
bool ReadImageIntensityPointSet( itk::SmartPointer<TPointSet> & target, const char *imageFile,
  const char *maskFile, std::vector<unsigned int> neighborhoodRadius, double sigma )
{
  if( std::string( imageFile ).length() < 3 )
    {
    std::cerr << " bad image file name " << std::string( imageFile ) << std::endl;
    target = ITK_NULLPTR;
    return false;
    }

  if( !ANTSFileExists( std::string( imageFile ) ) )
    {
    std::cerr << " image file " << std::string( imageFile ) << " does not exist . " << std::endl;
    target = ITK_NULLPTR;
    return false;
    }

  if( std::string( maskFile ).length() < 3 )
    {
    std::cerr << " bad mask file name " << std::string( maskFile ) << std::endl;
    target = ITK_NULLPTR;
    return false;
    }

  if( !ANTSFileExists( std::string( maskFile ) ) )
    {
    std::cerr << " mask file " << std::string( maskFile ) << " does not exist . " << std::endl;
    target = ITK_NULLPTR;
    return false;
    }

  if( neighborhoodRadius.size() != TImage::ImageDimension )
    {
    std::cerr << " size of the neighborhood radius is not equal to the image dimension." << std::endl;
    target = ITK_NULLPTR;
    return false;
    }

  typename TImage::Pointer intensityImage = ReadImage<TImage>( (char *)imageFile );
  typename TMask::Pointer maskImage = ReadImage<TMask>( (char *)maskFile );

  typedef itk::ImageIntensityAndGradientToPointSetFilter<TImage, TMask, TPointSet> FilterType;

  typename FilterType::NeighborhoodRadiusType radius;
  for( unsigned int d = 0; d < TImage::ImageDimension; d++ )
    {
    radius[d] = neighborhoodRadius[d];
    }

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1( intensityImage );
  filter->SetInput2( maskImage );
  filter->SetSigma( sigma );
  filter->SetNeighborhoodRadius( radius );
  filter->Update();

  target = filter->GetOutput();

  return true;
}

template <class TPointSet>
typename TPointSet::Pointer ReadLabeledPointSet( char* fn )
{
  if( !ANTSFileExists( std::string( fn ) ) )
    {
    std::cerr << " file " << std::string( fn ) << " does not exist . " << std::endl;
    return;
    }

  // Read the image files begin
  typedef itk::LabeledPointSetFileReader<TPointSet>   FileSourceType;
  typename FileSourceType::Pointer reffilter = FileSourceType::New();
  reffilter->SetFileName( fn );
  try
    {
    reffilter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception caught during point set reference file reading " << std::endl;
    std::cerr << e << std::endl;
    return NULL;
    }

  typename TPointSet::Pointer target = reffilter->GetOutput();

  return target;
}

template <class TPointSet>
bool WritePointSet( const TPointSet * const pointSet, const char * const file )
{
  if( std::string(file).length() < 3 )
    {
    return false;
    }

  typename itk::LabeledPointSetFileWriter<TPointSet>::Pointer writer =
    itk::LabeledPointSetFileWriter<TPointSet>::New();
  writer->SetFileName( file );
  if( !pointSet )
    {
    std::cerr << " Point set is null." << std::endl;
    std::exception();
    }
  writer->SetInput( pointSet );
  writer->Update();

  return true;
}

template <class TImageType>
bool WriteImage(const TImageType * const image, const char * const file)
{
  if( std::string(file).length() < 3 )
    {
    return false;
    }

  //  typename TImageType::DirectionType dir;
  // dir.SetIdentity();
  // image->SetDirection(dir);
  //  std::cout << " now Write direction " << image->GetOrigin() << std::endl;

  // if (writer->GetImageIO->GetNumberOfComponents() == 6)
  // NiftiDTICheck<TImageType>(image,file);

  if( file[0] == '0' && file[1] == 'x' )
    {
    void* ptr;
    sscanf(file, "%p", (void **)&ptr);
    *( static_cast<typename TImageType::Pointer *>( ptr ) ) = const_cast<TImageType *>(image);
    }
  else
    {
    typename itk::ImageFileWriter<TImageType>::Pointer writer =
      itk::ImageFileWriter<TImageType>::New();
    writer->SetFileName(file);
    if( !image )
      {
      std::cerr << "Image is null." << std::endl;
      std::exception();
      }
    writer->SetInput(image);
    writer->SetUseCompression( true );
    writer->Update();
    }
  return true;
}

template <class TImageType>
void WriteTensorImage(itk::SmartPointer<TImageType> image, const char *file, bool takeexp = true)
{
  typedef itk::ExpTensorImageFilter<TImageType, TImageType> ExpFilterType;
  typename itk::ImageFileWriter<TImageType>::Pointer writer =
    itk::ImageFileWriter<TImageType>::New();
  writer->SetFileName(file);

  typename TImageType::Pointer writeImage = image;

  if( takeexp )
    {
    typename ExpFilterType::Pointer expFilter = ExpFilterType::New();
    expFilter->SetInput( image );
    expFilter->Update();
    writeImage = expFilter->GetOutput();
    std::cout << "Taking Exp(D) before writing" << std::endl;
    }

  // convert from upper tri to lower tri
  NiftiDTICheck<TImageType>(writeImage, file, true); // BA May 30 2009 -- remove b/c ITK fixed NIFTI reader

  if( file[0] == '0' && file[1] == 'x' )
    {
    void* ptr;
    sscanf(file, "%p", (void **)&ptr);
    *( static_cast<typename TImageType::Pointer *>( ptr ) ) = writeImage;
    }
  else
    {
    writer->SetInput(writeImage);
    writer->SetUseCompression( true );
    writer->Update();
    }
}

template <class TImage, class TField>
typename TField::Pointer
ReadWarpFromFile( std::string warpfn, std::string ext)
{
  typedef TField                        FieldType;
  typedef typename FieldType::PixelType VectorType;
  enum { ImageDimension = FieldType::ImageDimension };

  typedef itk::Image<float, ImageDimension> RealImageType;
  typedef RealImageType                     ImageType;

//  typedef itk::Vector<float,itkGetStaticConstMacro(ImageDimension)>         VectorType;
//  typedef itk::Image<VectorType,itkGetStaticConstMacro(ImageDimension)>     FieldType;
// std::cout << " warp file name " << warpfn + ext << std::endl;

// First - read the vector fields
// NOTE : THE FIELD SHOULD WARP INPUT1 TO INPUT2, THUS SHOULD POINT
// FROM INPUT2 TO INPUT1
  std::string fn = warpfn + "x" + ext;
  typename RealImageType::Pointer xvec = ReadImage<ImageType>( (char *)fn.c_str() );
  //  std::cout << " done reading " << fn << std::endl;
  fn = warpfn + "y" + ext;
  typename RealImageType::Pointer yvec = ReadImage<ImageType>( (char *)fn.c_str() );
  // std::cout << " done reading " << fn << std::endl;
  fn = warpfn + "z" + ext;
  typename RealImageType::Pointer zvec = ITK_NULLPTR;
  // std::cout << " done reading " << fn << std::endl;
  if( ImageDimension == 3 )
    {
    zvec = ReadImage<ImageType>( (char *)fn.c_str() );
    }

  typename FieldType::Pointer field = AllocImage<FieldType>(xvec);

  itk::ImageRegionIteratorWithIndex<RealImageType>
  it( xvec, xvec->GetLargestPossibleRegion() );

  //  std::cout << " spacing xv " << xvec->GetSpacing()[0]
  // << " field " << field->GetSpacing()[0] << std::endl;

  unsigned int ct = 0;
  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    ct++;
    typename ImageType::IndexType index = it.GetIndex();

    VectorType disp;
    disp[0] = xvec->GetPixel(index);
    disp[1] = yvec->GetPixel(index);
    if( ImageDimension == 3 )
      {
      disp[2] = zvec->GetPixel(index);
      }

    field->SetPixel(index, disp);

    //    if (ct == 10000) std::cout << " 10000th pix " << disp << std::endl;
    }

  return field;
}



template <class TImage>
typename TImage::Pointer
MakeNewImage(typename TImage::Pointer image1, typename TImage::PixelType initval )
{
  typedef itk::ImageRegionIteratorWithIndex<TImage> Iterator;
  typename TImage::Pointer varimage = AllocImage<TImage>(image1);
  Iterator vfIter2( varimage,  varimage->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    if( initval >= 0 )
      {
      vfIter2.Set(initval);
      }
    else
      {
      vfIter2.Set(image1->GetPixel(vfIter2.GetIndex() ) );
      }
    }

  return varimage;
}

template <class TField>
void
WriteDisplacementField(const TField* const field, std::string filename)
{
  typedef TField                        FieldType;
  enum { ImageDimension = FieldType::ImageDimension };

  typedef itk::Image<float, ImageDimension> RealImageType;

  // Initialize the caster to the displacement field
  typedef itk::VectorIndexSelectionCastImageFilter<FieldType, RealImageType> IndexSelectCasterType;
  for( unsigned int dim = 0; dim < ImageDimension; dim++ )
    {
    typename IndexSelectCasterType::Pointer fieldCaster = IndexSelectCasterType::New();
    fieldCaster->SetInput( field );

    fieldCaster->SetIndex( dim );
    fieldCaster->Update();

    // Set up the output filename
    std::string outfile = filename + static_cast<char>('x' + dim) + std::string("vec.nii.gz");
    std::cout << "Writing displacements to " << outfile << " spacing "
                     << field->GetSpacing()[0] << std::endl;
    typename RealImageType::Pointer fieldcomponent = fieldCaster->GetOutput();
    fieldcomponent->SetSpacing(field->GetSpacing() );
    fieldcomponent->SetOrigin(field->GetOrigin() );
    fieldcomponent->SetDirection(field->GetDirection() );

    WriteImage<RealImageType>(fieldcomponent, outfile.c_str() );
    }
  std::cout << "...done" << std::endl;
  return;
}

template <class TField>
void
WriteDisplacementField2(const TField * const field, std::string filename, std::string app)
{
  typedef TField                        FieldType;
  enum { ImageDimension = FieldType::ImageDimension };

  typedef itk::Image<float, ImageDimension> RealImageType;

  // Initialize the caster to the displacement field
  typedef itk::VectorIndexSelectionCastImageFilter<FieldType, RealImageType> IndexSelectCasterType;
  for( unsigned int dim = 0; dim < ImageDimension; dim++ )
    {
    typename IndexSelectCasterType::Pointer fieldCaster = IndexSelectCasterType::New();
    fieldCaster->SetInput( field );

    fieldCaster->SetIndex( dim );
    fieldCaster->Update();

    // Set up the output filename
    std::string outfile = filename + static_cast<char>('x' + dim) + std::string(app);
    std::cout << "Writing displacements to " << outfile << " spacing "
                     << field->GetSpacing()[0] << std::endl;
    typename RealImageType::Pointer fieldcomponent = fieldCaster->GetOutput();
    fieldcomponent->SetSpacing(field->GetSpacing() );
    fieldcomponent->SetOrigin(field->GetOrigin() );

    WriteImage<RealImageType>(fieldcomponent, outfile.c_str() );
    }
  std::cout << "...done" << std::endl;
  return;
}

class nullBuf
: public std::streambuf
{
public:
  virtual std::streamsize xsputn( const char * itkNotUsed( s ), std::streamsize n ) ITK_OVERRIDE
    {
    return n;
    }

  virtual int overflow( int itkNotUsed( c ) ) ITK_OVERRIDE
    {
    return 1;
    }
};

class nullStream
: public std::ostream
{
  public:
    nullStream() : std::ostream( &buf ) {}
  private:
    nullBuf buf;
};

#endif
