/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: SmoothImage.cxx,v $
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

#include "antscout.hxx"
#include <algorithm>

#include "itkMedianImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "ReadWriteImage.h"

namespace ants
{
template <unsigned int ImageDimension>
int SmoothImage(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  std::string fn1 = std::string(argv[2]);
  float       sigma = atof(argv[3]);
  typename ImageType::Pointer image1 = NULL;
  typename ImageType::Pointer varimage = NULL;
  ReadImage<ImageType>(image1, argv[2]);

  typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> dgf;
  typedef itk::MedianImageFilter<ImageType, ImageType>           medf;
  typename dgf::Pointer filter = dgf::New();
  typename medf::Pointer filter2 = medf::New();
  bool usespacing = false;
  if( argc  >  5 )
    {
    usespacing = atoi(argv[5]);
    }
  bool usemedian = false;
  if( argc  >  6 )
    {
    usemedian = atoi(argv[6]);
    }
  if( !usespacing )
    {
    filter->SetUseImageSpacingOff();
    }
  else
    {
    filter->SetUseImageSpacingOn();
    }

  if( !usemedian )
    {
    filter->SetVariance(sigma * sigma);
    filter->SetMaximumError(.01f);
    filter->SetInput(image1);
    filter->Update();
    varimage = filter->GetOutput();
    }
  else
    {
    typename ImageType::SizeType rad;
    rad.Fill( (long unsigned int) sigma);
    filter2->SetRadius(rad);
    filter2->SetInput(image1);
    filter2->Update();
    varimage = filter2->GetOutput();
    }

  typename writertype::Pointer writer = writertype::New();
  writer->SetFileName(argv[4]);
  writer->SetInput( varimage );
  writer->Write();

  return 0;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int SmoothImage( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "SmoothImage" );

  std::remove( args.begin(), args.end(), std::string( "" ) );
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

  antscout->set_stream( out_stream );

  if( argc < 4 )
    {
    antscout << "Usage:  " << std::endl;
    antscout << argv[0]
             <<
      " ImageDimension image.ext smoothingsigma outimage.ext {sigma-is-in-spacing-coordinates-0/1} {medianfilter-0/1}"
             << std::endl;
    antscout << " if median, then sigma means radius of filtering " << std::endl;
    return 1;
    }

  switch( atoi(argv[1]) )
    {
    case 2:
      {
      SmoothImage<2>(argc, argv);
      }
      break;
    case 3:
      {
      SmoothImage<3>(argc, argv);
      }
      break;
    default:
      antscout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }

  return 0;
}
} // namespace ants
