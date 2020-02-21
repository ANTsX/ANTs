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

#include "ReadWriteData.h"

namespace ants
{
template <unsigned int ImageDimension, unsigned int NVectorComponents>
int MeasureMinMaxMean(int argc, char *argv[])
{
  typedef itk::Vector<float, NVectorComponents>                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  typename ImageType::Pointer image = nullptr;
  typename ImageType::Pointer mask = nullptr;
  PixelType     mean;  mean.Fill(0);
  PixelType     max;   max.Fill(0);
  PixelType     min;   min.Fill(9.e9);
  unsigned long ct = 0;

  ReadImage<ImageType>(image, argv[2]);
  bool takeabsval = false;
  if( argc > 4 )
    {
    takeabsval = std::stoi(argv[4]);
    }
  if( argc > 5 )
    {
    ReadImage<ImageType>(mask, argv[5]);
    }
  Iterator vfIter2( image,  image->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    bool isinside = true;
    if( mask )
      {
      if( mask->GetPixel(vfIter2.GetIndex() ).GetNorm() < 0.5 )
        {
        isinside = false;
        }
      }
    if( isinside )
      {
      PixelType val = vfIter2.Get();
      if( takeabsval )
        {
        for( unsigned int k = 0; k < NVectorComponents; k++ )
          {
          val[k] = fabs(val[k]);
          }
        }
      mean = mean + val;
      if( val.GetSquaredNorm() > max.GetSquaredNorm() )
        {
        max = val;
        }
      else if( val.GetSquaredNorm() < min.GetSquaredNorm() )
        {
        min = val;
        }
      ct++;
      }
    }

  if( ct > 0 )
    {
    mean = mean / (float)ct;
    }

  float variance = 0;
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    bool isinside = true;
    if( mask )
      {
      if( mask->GetPixel(vfIter2.GetIndex() ).GetNorm() < 0.5 )
        {
        isinside = false;
        }
      }
    if( isinside )
      {
      PixelType val = vfIter2.Get();
      if( takeabsval )
        {
        for( unsigned int k = 0; k < NVectorComponents; k++ )
          {
          val[k] = fabs(val[k]);
          }
        }
      variance += static_cast<float>( (val - mean).GetSquaredNorm() );
      }
    }

  float temp = (1.0f / (float)ct) * variance;
  std::cout <<  argv[2] << " Max : " << max << " Min : " << min << " Mean : " << mean << " Var : " <<  temp
           << " SD : " << std::sqrt(1.0f / (float)(ct - 1) * variance) << std::endl;

  if( argc > 3 )
    {
    std::ofstream logfile;
    logfile.open(argv[3], std::ofstream::app);
    if( logfile.good() )
      {
      logfile << argv[2] <<  " Max : " << max << " Min : " << min << " Mean : " << mean << std::endl;
      }
    else
      {
      std::cout << " cant open file " << argv[3] << std::endl;
      }
    logfile.close();
    }
  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int MeasureMinMaxMean( std::vector<std::string> args, std::ostream* /*out_stream = nullptr */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "MeasureMinMaxMean" );

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
  argv[argc] = nullptr;
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
    std::cout << "Basic useage ex: " << std::endl;
    std::cout << argv[0] << " ImageDimension  image.nii {log.txt} {take-absolute-value}  {mask-name} " << std::endl;
    std::cout << "  log.txt is optional  - take-abs-val reports min-max-mean of abs val image " << std::endl;
    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }
  int                       dim = std::stoi( argv[1] );
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(argv[2], itk::ImageIOFactory::FileModeEnum::ReadMode);
  imageIO->SetFileName(argv[2]);
  imageIO->ReadImageInformation();
  unsigned int ncomponents = imageIO->GetNumberOfComponents();

  int returnValue = EXIT_FAILURE;

  // Get the image dimension
  switch( std::stoi(argv[1]) )
    {
    case 2:
      {
      switch( ncomponents )
        {
        case 3:
          {
          returnValue = MeasureMinMaxMean<2, 3>(argc, argv);
          }
          break;
        case 2:
          {
          returnValue = MeasureMinMaxMean<2, 2>(argc, argv);
          }
          break;
        default:
          {
          returnValue = MeasureMinMaxMean<2, 1>(argc, argv);
          }
          break;
        }
      }
      break;
    case 3:
      {
      switch( ncomponents )
        {
        case 7:
          {
          returnValue = MeasureMinMaxMean<3, 7>(argc, argv);
          }
          break;
        case 6:
          {
          returnValue = MeasureMinMaxMean<3, 6>(argc, argv);
          }
          break;
        case 3:
          {
          returnValue = MeasureMinMaxMean<3, 3>(argc, argv);
          }
          break;
        default:
          {
          returnValue = MeasureMinMaxMean<3, 1>(argc, argv);
          }
          break;
        }
      }
      break;
    case 4:
      {
      switch( ncomponents )
        {
        case 7:
          {
          returnValue = MeasureMinMaxMean<4, 7>(argc, argv);
          }
          break;
        case 6:
          {
          returnValue = MeasureMinMaxMean<4, 6>(argc, argv);
          }
          break;
        case 4:
          {
          returnValue = MeasureMinMaxMean<4, 4>(argc, argv);
          }
          break;
        case 3:
          {
          returnValue = MeasureMinMaxMean<4, 3>(argc, argv);
          }
          break;
        case 2:
          {
          returnValue = MeasureMinMaxMean<4, 2>(argc, argv);
          }
          break;
        default:
          {
          returnValue = MeasureMinMaxMean<4, 1>(argc, argv);
          }
          break;
        }
      }
      break;
    default:
      std::cout << " not supported " << dim  << std::endl;
      return EXIT_FAILURE;
    }
  return returnValue;

}
} // namespace ants
