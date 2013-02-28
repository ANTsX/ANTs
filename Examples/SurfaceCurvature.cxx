



#include "antsUtilities.h"
#include <algorithm>

#include "itkSurfaceCurvatureBase.h"
#include "itkSurfaceImageCurvature.h"

#include "ReadWriteImage.h"

namespace ants
{
/*
void test1()
{

typedef itk::SurfaceCurvatureBase<ImageType>  ParamType;
  ParamType::Pointer Parameterizer=ParamType::New();

//  Parameterizer->TestEstimateTangentPlane(p);
  Parameterizer->FindNeighborhood();
//  Parameterizer->WeightedEstimateTangentPlane(  Parameterizer->GetOrigin() );

  Parameterizer->EstimateTangentPlane(
    Parameterizer->GetAveragePoint());
  Parameterizer->PrintFrame();


//  Parameterizer->SetOrigin(Parameterizer->GetAveragePoint());

  for(int i=0; i<3; i++){
  Parameterizer->ComputeWeightsAndDirectionalKappaAndAngles
    (Parameterizer->GetOrigin());
  Parameterizer->ComputeFrame(Parameterizer->GetOrigin());
  Parameterizer->EstimateCurvature();
  Parameterizer->PrintFrame();
  }

  Parameterizer->ComputeJoshiFrame(Parameterizer->GetOrigin());  Parameterizer->PrintFrame();
  antscout << " err 1 " << Parameterizer->ErrorEstimate(Parameterizer->GetOrigin()) <<
    " err 2 " << Parameterizer->ErrorEstimate(Parameterizer->GetOrigin(),-1) << std::endl;

}


*/

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int SurfaceCurvature( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "SurfaceCurvature" );

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

  if( argc < 3 )
    {
    antscout << " usage :  SurfaceCurvature FileNameIn FileNameOut sigma option  " << std::endl;
    antscout << " e.g  :   SurfaceCurvature    BrainIn.nii BrainOut.nii   3  0 " << std::endl;
    antscout << " option 0 means just compute mean curvature from intensity " << std::endl;
    antscout << " option 5 means characterize surface from intensity " << std::endl;
    antscout << " option 6 means compute gaussian curvature " << std::endl;
    antscout << " ... " << std::endl;
    antscout << " for surface characterization " << std::endl;
    antscout << " 1 == (+) bowl " << std::endl;
    antscout << " 2 == (-) bowl  " << std::endl;
    antscout << " 3 == (+) saddle " << std::endl;
    antscout << " 4 == (-) saddle " << std::endl;
    antscout << " 5 == (+) U " << std::endl;
    antscout << " 6 == (-) U " << std::endl;
    antscout << " 7 == flat " << std::endl;
    antscout << " 8 == a perfectly even saddle (rare) " << std::endl;
    antscout << " " << std::endl;
    antscout << " we add 128 to mean curvature results s.t. they are differentiated from background (zero) "
             << std::endl;
    return 0;
    }

  typedef itk::Image<float, 3> ImageType;
  typedef itk::Image<float, 3> floatImageType;
  enum { ImageDimension = ImageType::ImageDimension };
  typedef itk::SurfaceImageCurvature<ImageType> ParamType;
  ParamType::Pointer Parameterizer = ParamType::New();
  typedef  ImageType::PixelType PixType;

  int   opt = 0;
  float sig = 1.0;
  if( argc > 3 )
    {
    sig = atof( argv[3]);
    }
  antscout << " sigma " << sig << std::endl;
  if( argc > 4 )
    {
    opt = (int) atoi(argv[4]);
    }

  if( opt < 0 )
    {
    antscout << " error " << std::endl;
    return 0;
    }

  ImageType::Pointer input;
  ReadImage<ImageType>(input, argv[1]);
  antscout << " done reading " << std::string(argv[1]) << std::endl;

  //  float ballradius = 2.0;
  // if (argc >= 6) ballradius = (float) atof(argv[5]);
  // if (ballradius > 0 && thresh > 0) input = SegmentImage<ImageType>(input, thresh, ballradius);

  Parameterizer->SetInputImage(input);

  //  Parameterizer->ProcessLabelImage();
  Parameterizer->SetNeighborhoodRadius( 1. );
//  antscout << " set sig " ;  std::cin >> sig;
  if( sig <= 0.5 )
    {
    sig = 1.66;
    }
  antscout << " sigma " << sig << " option " << opt << std::endl;
  Parameterizer->SetSigma(sig);

  if( opt == 1 )
    {
    Parameterizer->SetUseLabel(true);
    Parameterizer->SetUseGeodesicNeighborhood(false);
    }
  else
    {
    Parameterizer->SetUseLabel(false);
    Parameterizer->SetUseGeodesicNeighborhood(false);
    float sign = 1.0;
    if( opt == 3 )
      {
      sign = -1.0;
      }
    antscout << " setting outward direction as " << sign;
    Parameterizer->SetkSign(sign);
    Parameterizer->SetThreshold(0);
    }
//  Parameterizer->ComputeSurfaceArea();
//  Parameterizer->IntegrateFunctionOverSurface();
//  Parameterizer->IntegrateFunctionOverSurface(true);

  antscout << " computing frame " << std::endl;
  if( opt != 5 && opt != 6 )
    {
    Parameterizer->ComputeFrameOverDomain( 3 );
    }
  else
    {
    Parameterizer->ComputeFrameOverDomain( opt );
    }

  //   Parameterizer->SetNeighborhoodRadius( 2 );
  //  Parameterizer->LevelSetMeanCurvature();
  //  Parameterizer->SetNeighborhoodRadius( 2.9   );
  //  Parameterizer->IntegrateFunctionOverSurface(false);
  //  Parameterizer->SetNeighborhoodRadius( 1.5  );
  //  Parameterizer->IntegrateFunctionOverSurface(true);
  //   for (int i=0; i<1; i++) Parameterizer->PostProcessGeometry();

  ImageType::Pointer output = NULL;

  //  Parameterizer->GetFunctionImage()->SetSpacing( input->GetSpacing() );
  //  Parameterizer->GetFunctionImage()->SetDirection( input->GetDirection() );
  //  Parameterizer->GetFunctionImage()->SetOrigin( input->GetOrigin() );
  //  smooth->SetSpacing(reader->GetOutput()->GetSpacing());
  // SmoothImage(Parameterizer->GetFunctionImage(),smooth,3);
  // NormalizeImage(smooth,output,mn);
  //  NormalizeImage(Parameterizer->GetFunctionImage(),output,mn);

  WriteImage<floatImageType>(Parameterizer->GetFunctionImage(), argv[2]);

  antscout << " done writing ";

  return 1;
}
} // namespace ants
