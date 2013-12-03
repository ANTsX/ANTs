/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: LabelClustersUniquely.cxx,v $
  Language:  C++
  Date:      $Date: 2008/11/15 23:46:06 $
  Version:   $Revision: 1.18 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "antsUtilities.h"
#include <algorithm>

#include "itkDiscreteGaussianImageFilter.h"

//  RecursiveAverageImages img1  img2 weightonimg2 outputname

// We divide the 2nd input image by its mean and add it to the first
// input image with weight 1/n.
// The output overwrites the 1st img with the sum.

#include <list>
#include <vector>
#include <fstream>
#include "vnl/vnl_vector.h"

#include "itkMinimumMaximumImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkCastImageFilter.h"
#include  "ReadWriteImage.h"

namespace ants
{
template <unsigned int ImageDimension>
int  LabelUniquely(int argc, char *argv[])
{
  typedef float PixelType;
//  const unsigned int ImageDimension = AvantsImageDimension;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef typename ImageType::SizeType                                    SizeType;
  typedef typename ImageType::SpacingType                                 SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  // typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;

  typedef float                                            InternalPixelType;
  typedef int                                              ULPixelType;
  typedef itk::Image<ULPixelType, ImageDimension>          labelimagetype;
  typedef itk::CastImageFilter<ImageType, labelimagetype>  CastFilterType;
  typedef itk::CastImageFilter< labelimagetype, ImageType> CastFilterType2;

  typedef ImageType                                                          InternalImageType;
  typedef ImageType                                                          OutputImageType;
  typedef itk::ConnectedComponentImageFilter<labelimagetype, labelimagetype> FilterType;
  typedef itk::RelabelComponentImageFilter<labelimagetype, labelimagetype>   RelabelType;

  // want the average value in each cluster as defined by the mask and the value thresh and the clust thresh

  if( argc < 2 )
    {
    std::cout << "missing 1st filename" << std::endl;
    throw;
    }
  if( argc < 3 )
    {
    std::cout << "missing 2nd filename" << std::endl;
    throw;
    }
  if( argc < 4 )
    {
    std::cout << "missing cluster thresholod" << std::endl;
    throw;
    }
  std::string fn1 = std::string(argv[1]);
  float       clusterthresh = atof(argv[3]);

  typename ImageType::Pointer image1 = NULL;

  ReadImage<ImageType>(image1, fn1.c_str() );

  //  typename
  typename FilterType::Pointer filter = FilterType::New();
// typename
  typename RelabelType::Pointer relabel = RelabelType::New();

  typename CastFilterType::Pointer castInput = CastFilterType::New();
  castInput->SetInput(image1);

  filter->SetInput( castInput->GetOutput() );
  int fullyConnected = 0; // atoi( argv[5] );
  filter->SetFullyConnected( fullyConnected );
  relabel->SetInput( filter->GetOutput() );
  relabel->SetMinimumObjectSize( (unsigned int) clusterthresh );

  try
    {
    relabel->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cout << "Relabel: exception caught !" << std::endl;
    std::cout << excep << std::endl;
    }

//  float maximum=relabel->GetNumberOfObjects();
  typename CastFilterType2::Pointer castRegions = CastFilterType2::New();
  castRegions->SetInput( relabel->GetOutput() );
  castRegions->Update();
  WriteImage<ImageType>(   castRegions->GetOutput() , argv[2] );

  return 0;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int LabelClustersUniquely( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "LabelClustersUniquely" );

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

  if( argc < 3 )
    {
    std::cout << "Usage:  " << std::endl;
    std::cout << argv[0] << " ImageDimension clustersin.hdr labeledclustersout.hdr   sizethresh " << std::endl;
    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }

  switch( atoi(argv[1]) )
    {
    case 2:
      {
      LabelUniquely<2>(argc, argv + 1);
      }
      break;
    case 3:
      {
      LabelUniquely<3>(argc, argv + 1);
      }
      break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants


