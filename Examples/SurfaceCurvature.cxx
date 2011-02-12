

// #include "curvatureapp.h"

#include "itkSurfaceCurvatureBase.h"
#include "itkSurfaceImageCurvature.h"

#include "ReadWriteImage.h"

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
  std::cout << " err 1 " << Parameterizer->ErrorEstimate(Parameterizer->GetOrigin()) <<
    " err 2 " << Parameterizer->ErrorEstimate(Parameterizer->GetOrigin(),-1) << std::endl;

}


*/

int main(int argc, char *argv[])
{
  if( argc < 3 )
    {
    std::cout << " usage :  SurfaceCurvature FileNameIn FileNameOut sigma option  " << std::endl;
    std::cout << " e.g  :   SurfaceCurvature    BrainIn.nii BrainOut.nii   3  0 " << std::endl;
    std::cout << " option 0 means just compute mean curvature from intensity " << std::endl;
    std::cout << " option 5 means characterize surface from intensity " << std::endl;
    std::cout << " option 6 means compute gaussian curvature " << std::endl;
    std::cout << " ... " << std::endl;
    std::cout << " for surface characterization " << std::endl;
    std::cout << " 1 == (+) bowl " << std::endl;
    std::cout << " 2 == (-) bowl  " << std::endl;
    std::cout << " 3 == (+) saddle " << std::endl;
    std::cout << " 4 == (-) saddle " << std::endl;
    std::cout << " 5 == (+) U " << std::endl;
    std::cout << " 6 == (-) U " << std::endl;
    std::cout << " 7 == flat " << std::endl;
    std::cout << " 8 == a perfectly even saddle (rare) " << std::endl;
    std::cout << " " << std::endl;
    std::cout << " we add 128 to mean curvature results s.t. they are differentiated from background (zero) "
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
  std::cout << " sigma " << sig << std::endl;
  if( argc > 4 )
    {
    opt = (int) atoi(argv[4]);
    }

  if( opt < 0 )
    {
    std::cout << " error " << std::endl;
    return 0;
    }

  ImageType::Pointer input;
  ReadImage<ImageType>(input, argv[1]);
  std::cout << " done reading " << std::string(argv[1]) << std::endl;

  //  float ballradius = 2.0;
  // if (argc >= 6) ballradius = (float) atof(argv[5]);
  // if (ballradius > 0 && thresh > 0) input = SegmentImage<ImageType>(input, thresh, ballradius);

  Parameterizer->SetInput(input);

  //  Parameterizer->ProcessLabelImage();
  Parameterizer->SetNeighborhoodRadius( 1. );
//  std::cout << " set sig " ;  std::cin >> sig;
  if( sig <= 0.5 )
    {
    sig = 1.66;
    }
  std::cout << " sigma " << sig << " option " << opt << std::endl;
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
    std::cout << " setting outward direction as " << sign;
    Parameterizer->SetkSign(sign);
    Parameterizer->SetThreshold(0);
    }
//  Parameterizer->ComputeSurfaceArea();
//  Parameterizer->IntegrateFunctionOverSurface();
//  Parameterizer->IntegrateFunctionOverSurface(true);

  std::cout << " computing frame " << std::endl;
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

  std::cout << " done writing ";

  return 1;
}
