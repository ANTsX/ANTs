/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "antsUtilities.h"
#include <algorithm>

#include "itkMedianImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "ReadWriteData.h"

namespace ants
{
template <unsigned int ImageDimension>
int SmoothImage(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;

  std::vector<float> sigmaVector = ConvertVector<float>( argv[3] );

  typename ImageType::Pointer image1 = ITK_NULLPTR;
  typename ImageType::Pointer varimage = ITK_NULLPTR;
  ReadImage<ImageType>(image1, argv[2]);

  typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> dgf;
  typedef itk::MedianImageFilter<ImageType, ImageType>           medf;
  typename dgf::Pointer filter = dgf::New();
  typename medf::Pointer filter2 = medf::New();
  bool usespacing = false;
  float gaussianmaxerror = 0.01;
  int gaussianmaxkernelwidth = 32;
  if( argc  >  5 )
    {
    usespacing = atoi(argv[5]);
    }
  bool usemedian = false;
  if( argc  >  6 )
    {
    usemedian = atoi(argv[6]);
    }
  if( argc  >  7 )
    {
    gaussianmaxerror = atof(argv[7]);
    }
  if( argc  >  8 )
    {
    gaussianmaxkernelwidth = atoi(argv[8]);
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
    if( sigmaVector.size() == 1 )
      {
      filter->SetVariance( vnl_math_sqr( sigmaVector[0] ) );
      }
    else if( sigmaVector.size() == ImageDimension )
      {
      typename dgf::ArrayType varianceArray;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        varianceArray[d] = vnl_math_sqr( sigmaVector[d] );
        }
      filter->SetVariance( varianceArray );
      }
    else
      {
      std::cerr << "Incorrect sigma vector size.  Must either be of size 1 or ImageDimension." << std::endl;
      }
    filter->SetMaximumError( gaussianmaxerror );
    filter->SetMaximumKernelWidth( gaussianmaxkernelwidth );
    filter->SetInput( image1 );
    filter->Update();
    varimage = filter->GetOutput();
    }
  else
    {
    typename ImageType::SizeType rad;
    if( sigmaVector.size() == 1 )
      {
      rad.Fill( static_cast<unsigned long>( sigmaVector[0] ) );
      }
    else if( sigmaVector.size() == ImageDimension )
      {
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        rad[d] = sigmaVector[d];
        }
      }
    else
      {
      std::cerr << "Incorrect sigma vector size.  Must either be of size 1 or ImageDimension." << std::endl;
      }
    filter2->SetRadius(rad);
    filter2->SetInput( image1 );
    filter2->Update();
    varimage = filter2->GetOutput();
    }
  WriteImage<ImageType>( varimage, argv[4] );
  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int SmoothImage( std::vector<std::string> args, std::ostream* /*out_stream = ITK_NULLPTR */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "SmoothImage" );

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
  argv[argc] = ITK_NULLPTR;
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

  if( argc < 4 )
    {
    std::cout << "Usage:  " << std::endl;
    std::cout << argv[0]
             <<
      " ImageDimension image.ext smoothingsigma outimage.ext {sigma-is-in-spacing-coordinates-0/1} {medianfilter-0/1} {GaussianSetMaximumError=0.01} {GaussianSetMaximumKernelWidth=32}"
             << std::endl;
    std::cout << " if median, then sigma means radius of filtering " << std::endl;
    std::cout << " A separate sigma can be specified for each dimension, e.g., 1.5x1x2 " << std::endl;
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
      return SmoothImage<2>(argc, argv);
      }
      break;
    case 3:
      {
      return SmoothImage<3>(argc, argv);
      }
      break;
    case 4:
      {
      return SmoothImage<4>(argc, argv);
      }
      break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
