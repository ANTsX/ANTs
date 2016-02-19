/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "iMathFunctions.h"
#include "ReadWriteData.h"
#include "antsUtilities.h"


namespace ants
{


void WIP(int argc, char **argv)
{
  std::cout << "You have reached an umimplemented section of code by calling:" << std::endl;
  for (int i=0; i<argc; i++)
  {
    std::cout << argv[i] << " ";
  }
  std::cout << std::endl;
  exit(1);
}

template <class T>
bool from_string(T& t,
                 const std::string& s,
                 std::ios_base & (*f)(std::ios_base &) )
{
  std::istringstream iss(s);

  iss >> f >> t;

  // Check to see that there is nothing left over
  if( !iss.eof() )
    {
    return false;
    }

  return true;
}

template <class T>
std::string ants_to_string(T t)
{
  std::stringstream istream;

  istream << t;
  return istream.str();
}

std::string ANTSOptionName(const char *str)
{
  std::string            filename = str;
  std::string::size_type pos = filename.rfind( "=" );
  std::string            name = std::string( filename, 0, pos );

  return name;
}

std::string ANTSOptionValue(const char *str)
{
  std::string            filename = str;
  std::string::size_type pos = filename.rfind( "=" );
  std::string            value = std::string( filename, pos + 1, filename.length() );

  return value;
}

std::string ANTSGetFilePrefix(const char *str)
{
  const std::string      filename = str;
  std::string::size_type pos = filename.rfind( "." );
  std::string            filepre = std::string( filename, 0, pos );

#if 0 // HACK:  This does nothing useful
  if( pos != std::string::npos )
    {
    std::string extension = std::string( filename, pos, filename.length() - 1);
    if( extension == std::string(".gz") )
      {
      pos = filepre.rfind( "." );
      // extension = std::string( filepre, pos, filepre.length() - 1 );
      }
    }
#endif
  return filepre;
}

//
// iMath was a gigantic switch statement that had 3 duplicated
// lists of 'if (operation == <op>)' clauses for 2d, 3d, and 4d. I
// figured out which functions were 2D only, 3D only and 4D Only,
// which were valid for all dimensions,  which were 2d and 3d, and
// which were 3d and 4d.
// So there's a template method for each case, and they're assembled
// for each dimension in an Explicit Template Function below.
template<unsigned DIM>
int
iMathHelper2DOnly(int argc, char **argv)
{
  std::string operation = std::string(argv[3]);
  WIP(argc, argv);

  return EXIT_FAILURE;
}

template <unsigned DIM>
int
iMathHelper2DOr3D(int argc, char **argv)
{
  std::string operation = std::string(argv[3]);
  WIP(argc, argv);

  return EXIT_FAILURE;
}

template <unsigned DIM>
int
iMathHelper3DOr4D(int argc, char **argv)
{
  std::string operation = std::string(argv[3]);
  WIP(argc, argv);

  return EXIT_FAILURE;
}

template <unsigned DIM>
int
iMathHelper3DOnly(int argc, char **argv)
{
  std::string operation = std::string(argv[3]);
  WIP(argc, argv);

  return EXIT_FAILURE;
}

template <unsigned DIM>
int
iMathHelper4DOnly(int argc, char **argv)
{
  std::string operation = std::string(argv[3]);
  WIP(argc, argv);

  return EXIT_FAILURE;
}

template <unsigned DIM>
int
iMathHelperAll(int argc, char **argv)
{
  std::string operation = std::string(argv[3]);
  std::string inName = std::string(argv[4]);
  std::string outName = std::string(argv[2]);

  typedef float PixelType;
  if( operation == "BlobDetector" )
    {
    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;

    if ( argc < 6 )
    {
      std::cout << "BlobDetector: Not enough input parameters" << std::endl;
      return EXIT_FAILURE;
    }

    unsigned int nBlobs = atoi( argv[5] );

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      output = iMathBlobDetector<ImageType>(input,nBlobs);
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "BlobDetector: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }
  else if( operation == "Canny" )
    {
    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;

    if ( argc < 8 )
    {
      std::cout << "Canny: Not enough input parameters" << std::endl;
      return EXIT_FAILURE;
    }

    double sigma = atof( argv[5] );
    double lower = atof( argv[6] );
    double upper = atof( argv[7] );

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      output = iMathCanny<ImageType>(input,sigma,lower,upper);
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "Canny: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }
  else if( operation == "DistanceMap" )
    {

    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;

    bool useSpacing = iMathDistanceMapUseSpacing;

    if ( argc >= 6 )
      {
      useSpacing = atoi(argv[5]);
      }

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      output = iMathDistanceMap<ImageType>(input, useSpacing);
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "DistanceMap: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }
  else if( operation == "FillHoles" )
    {

    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;

    double holeType = iMathFillHolesHoleParam;

    if ( argc >= 6 )
      {
      holeType = atof(argv[5]);
      }

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      output = iMathFillHoles<ImageType>(input, holeType);
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "FillHoles: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }
  else if( operation == "GC" )
    {
    unsigned long radius = iMathGCRadius;

    if ( argc >= 6 )
      {
      radius = atoi(argv[5]);
      }

    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      output = iMathGC<ImageType>(input, radius);
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "GC: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }
  else if( operation == "GD" )
    {
    unsigned long radius = iMathGDRadius;

    if ( argc >= 6 )
      {
      radius = atoi(argv[5]);
      }

    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      output = iMathGD<ImageType>(input, radius);
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "GD: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }
  else if( operation == "GE" )
    {
    unsigned long radius = iMathGERadius;

    if ( argc >= 6 )
      {
      radius = atoi(argv[5]);
      }

    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      output = iMathGE<ImageType>(input, radius);
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "GE: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }
  else if( operation == "GO" )
    {
    unsigned long radius = iMathMORadius;

    if ( argc >= 6 )
      {
      radius = atoi(argv[5]);
      }

    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      output = iMathGO<ImageType>(input, radius);
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "GO: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }
  else if( operation == "GetLargestComponent" )
    {
    unsigned long minSize = iMathGetLargestComponentMinSize;
    if ( argc > 5)
    {
      minSize = atoi( argv[5] );
    }

    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      output = iMathGetLargestComponent<ImageType>(input, minSize);
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "GetLargestComponents: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }
  else if( operation == "Grad" )
    {
    double sigma = iMathGradSigma;
    bool normalize = iMathGradNormalize;

    if ( argc >= 6 )
      {
      sigma = atof(argv[5]);
      }
    if ( argc >= 7 )
      {
      normalize = (bool) atoi(argv[6]);
      }

    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      output = iMathGrad<ImageType>(input, sigma, normalize);
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "Grad: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }
  else if( operation == "HistogramEqualization" )
    {
    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;
    float alpha = 0;
    float beta  = 1;
    if ( argc >= 6 )
      {
      alpha = atof(argv[5]);
      }
    if ( argc >= 7 )
      {
      beta = atof(argv[6]);
      }

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      std::cout << " a " << alpha << " b " << beta << std::endl;
      output = iMathHistogramEqualization<ImageType>(input, alpha, beta, 1 );
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "HistogramEqualization: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }
  else if( operation == "Laplacian" )
    {
    double sigma = iMathLaplacianSigma;
    bool normalize = iMathLaplacianNormalize;

    if ( argc >= 6 )
      {
      sigma = atof(argv[5]);
      }
    if ( argc >= 7 )
      {
      normalize = (bool) atoi(argv[6]);
      }

    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      output = iMathLaplacian<ImageType>(input, sigma, normalize);
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "Laplacian: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }
  else if( operation == "MC" )
    {

    unsigned long radius = iMathMCRadius;
    PixelType value = iMathMCValue;

    if ( argc >= 6 )
      {
      radius = atoi(argv[5]);
      }
    if ( argc >= 7 )
      {
      value = (PixelType)( atof( argv[6]) );
      }

    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      output = iMathMC<ImageType>(input, radius, value);
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "MC: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }
  else if( operation == "MD" )
    {

    unsigned long radius = iMathMDRadius;
    PixelType value = iMathMDValue;

    if ( argc >= 6 )
      {
      radius = atoi(argv[5]);
      }
    if ( argc >= 7 )
      {
      value = (PixelType)( atof( argv[6]) );
      }

    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      output = iMathMD<ImageType>(input, radius, value);
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "MD: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }
  else if( operation == "ME" )
    {
    unsigned long radius = iMathMERadius;
    PixelType value = iMathMEValue;

    if ( argc >= 6 )
      {
      radius = atoi(argv[5]);
      }
    if ( argc >= 7 )
      {
      value = (PixelType)( atof( argv[6]) );
      }

    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      output = iMathME<ImageType>(input, radius, value);
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "ME: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }
  else if( operation == "MO" )
    {
    unsigned long radius = iMathMORadius;
    PixelType value = iMathMOValue;

    if ( argc >= 6 )
      {
      radius = atoi(argv[5]);
      }
    if ( argc >= 7 )
      {
      value = (PixelType)( atof( argv[6]) );
      }

    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      output = iMathMO<ImageType>(input, radius, value);
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "MO: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }
  else if( operation == "MaurerDistance" )
    {
    PixelType foreground = iMathMaurerDistanceForeground;

    if ( argc >= 6 )
      {
      foreground = (PixelType) atof(argv[5]);
      }

    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      output = iMathMaurerDistance<ImageType>(input, foreground);
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "MaurerDistance: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }
  else if( operation == "Normalize" )
    {
    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      output = iMathNormalize<ImageType>(input);
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "Normalize: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }
  else if( operation == "Pad" )
    {

    int padding = atoi(argv[5]);

    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      output = iMathPad<ImageType>(input, padding);
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "Pad: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }
  else if( operation == "PeronaMalik" )
    {
    double conductance = iMathPeronaMalikConductance;
    unsigned long nIterations = iMathPeronaMalikNIterations;

    if ( argc >= 6 )
      {
      nIterations = atoi(argv[5]);
      }
    if ( argc >= 7 )
      {
      conductance = (PixelType) atof(argv[6]);
      }

    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      output = iMathPeronaMalik<ImageType>(input, nIterations, conductance);
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "PeronaMalik: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }
  else if( operation == "Sharpen" )
    {
    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      output = iMathSharpen<ImageType>(input);
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "Sharpen: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }
  else if( operation == "TruncateIntensity" )
    {
    typedef itk::Image<float,DIM>         ImageType;
    typedef itk::Image<unsigned int,DIM>           MaskType;

    int nBins = iMathTruncateIntensityNBins;

    if ( argc < 7 )
      {
      std::cerr << "TruncateIntensity needs a lower and upper quantile" << std::endl;
      return EXIT_FAILURE;
      }

    double lowerQ = atof( argv[5] );
    double upperQ = atof( argv[6] );

    if ( argc >= 8 )
      {
      nBins= atoi(argv[7]);
      }

    typename MaskType::Pointer mask = NULL;
    if ( argc >= 9 )
      {
      ReadImage<MaskType>( mask, argv[8] );
      }

    typedef itk::Image<float,DIM> ImageType;
    typename ImageType::Pointer input = NULL;
    typename ImageType::Pointer output = NULL;

    ReadImage<ImageType>( input, inName.c_str() );
    if ( input.IsNull() )
      {
      return EXIT_FAILURE;
      }

    try
      {
      output = iMathTruncateIntensity<ImageType>(input, lowerQ, upperQ, nBins, mask);
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "TruncateIntensity: exception caught !" << std::endl;
      std::cout << excep << std::endl;
      }

    WriteImage<ImageType>( output, outName.c_str() );

    return EXIT_SUCCESS;
    }

  return EXIT_FAILURE;
}

template <unsigned DIM>
int
iMathHelper(int , char **)
{
  return 1;
}

template <>
int
iMathHelper<2>(int argc, char **argv)
{
  int returnval = iMathHelperAll<2>(argc,argv);
  if(returnval == EXIT_FAILURE)
    {
    returnval = iMathHelper2DOnly<2>(argc,argv);
    }
  if(returnval == EXIT_FAILURE)
    {
    returnval = iMathHelper2DOr3D<2>(argc,argv);
    }
  return returnval;
}

template <>
int
iMathHelper<3>(int argc, char **argv)
{
  int returnval = iMathHelperAll<3>(argc,argv);
  if(returnval == EXIT_FAILURE)
    {
    returnval = iMathHelper3DOnly<3>(argc,argv);
    }
  if(returnval == EXIT_FAILURE)
    {
    returnval = iMathHelper2DOr3D<3>(argc,argv);
    }
  if(returnval == EXIT_FAILURE)
    {
    returnval = iMathHelper3DOr4D<3>(argc,argv);
    }
  return returnval;
}

template <>
int
iMathHelper<4>(int argc, char **argv)
{
  int returnval = iMathHelperAll<4>(argc,argv);
  if(returnval == EXIT_FAILURE)
    {
    returnval = iMathHelper4DOnly<4>(argc,argv);
    }
  if(returnval == EXIT_FAILURE)
    {
    returnval = iMathHelper3DOr4D<4>(argc,argv);
    }
  return returnval;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int iMath( std::vector<std::string> args, std::ostream * itkNotUsed( out_stream ) )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "iMath" );

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

  if( argc < 5 )
    {
    std::cout << "\nUsage: " << argv[0]
              << " ImageDimension <OutputImage.ext> [operations and inputs] <Image1.ext> <Image2.ext>" << std::endl;

    std::cout << "\nUsage Information " << std::endl;
    std::cout << " ImageDimension: 2 or 3 (for 2 or 3 dimensional operations)." << std::endl;
    std::cout << " ImageDimension: 4 (for operations on 4D file, e.g. time-series data)." << std::endl;
    std::cout << " Operator: See list of valid operators below." << std::endl;


    std::cout << "Mask and Label set operations" << std::endl;
    std::cout << "-----------------------------" << std::endl;
    std::cout << "  GetLargestComponent    : Get the single largest labeled object in an image" << std::endl;
    std::cout << "    Usage                : GetLargestComponent InputImage.ext {MinObjectSize=50}" << std::endl;

    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }

  int returnvalue = EXIT_SUCCESS;

  std::string operation = std::string(argv[3]);

  unsigned int imageDimension = atoi(argv[1]);

  switch( imageDimension )
    {
    case 2:
      returnvalue = iMathHelper<2>(argc,argv);
      break;
    case 3:
      returnvalue = iMathHelper<3>(argc,argv);
      break;
    case 4:
      returnvalue = iMathHelper<4>(argc,argv);
      break;
    default:
      std::cout << " Dimension " << imageDimension << " is not supported " << std::endl;
      return EXIT_FAILURE;
    }

    if ( returnvalue == EXIT_FAILURE )
      {
      std::cout << " Operation " << operation << " not found or not supported for dimension " << imageDimension << std::endl;
      }

  return returnvalue;
}

} // namespace ants
