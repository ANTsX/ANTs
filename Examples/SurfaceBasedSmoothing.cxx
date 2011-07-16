
#include "itkBinaryThresholdImageFilter.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "ReadWriteImage.h"
#include "itkSurfaceImageCurvature.h"

int main(int argc, char *argv[])
{
  if( argc < 3 )
    {
    std::cout << " usage :  " << argv[0] << " ImageToSmooth  sigma SurfaceImage  outname  {numrepeatsofsmoothing}"
              << std::endl;
    std::cout << " We assume the SurfaceImage has a label == 1 that defines the surface " << std::endl;
    std::cout <<   " sigma  defines the geodesic n-hood radius --- numrepeats allows one to use " << std::endl;
    std::cout << " a small geodesic n-hood repeatedly applied many times -- faster computation, same effect "
              << std::endl;
    return 0;
    }
  typedef itk::Image<float, 3> ImageType;
  enum { ImageDimension = ImageType::ImageDimension };
  typedef itk::Image<float, ImageDimension>     floatImageType;
  typedef itk::SurfaceImageCurvature<ImageType> ParamType;
  ParamType::Pointer Parameterizer = ParamType::New();

  typedef  itk::ImageFileReader<ImageType> ReaderType;
  typedef  ImageType::PixelType            PixType;

//  std::string fn="C://Data//brain15labelimage.img";
  std::string ext = ".nii";

  float opt = 0;
  float sig = 1.0;
  //  float thresh=0.0;
  if( argc > 2 )
    {
    sig = atof( argv[2]);
    }
  unsigned int numrepeats = 0;
  if( argc > 5 )
    {
    numrepeats = atoi(argv[5]);
    }

  ImageType::Pointer input;
  ReadImage<ImageType>(input, argv[1]);
  ImageType::Pointer surflabel;
  ReadImage<ImageType>(surflabel, argv[3]);

  Parameterizer->SetInput(surflabel);
  Parameterizer->SetFunctionImage(input);
  Parameterizer->SetNeighborhoodRadius( sig );
  if( sig <= 0 )
    {
    sig = 1.0;
    }
  std::cout << " sigma " << sig << " thresh " << opt << std::endl;
  Parameterizer->SetSigma(sig);
  Parameterizer->SetUseGeodesicNeighborhood(true);
  Parameterizer->SetUseLabel(true);
  Parameterizer->SetThreshold(0.5);
  //  Parameterizer->ComputeSurfaceArea();
  //  Parameterizer->IntegrateFunctionOverSurface();

  Parameterizer->SetNeighborhoodRadius( sig );
  std::cout << " begin integration NOW " << std::endl;
  Parameterizer->IntegrateFunctionOverSurface(true);
  for( unsigned int i = 0; i < numrepeats; i++ )
    {
    Parameterizer->IntegrateFunctionOverSurface(true);
    }
  std::cout << " end integration  " << std::endl;
  // Parameterizer->PostProcessGeometry();

  //  double mn=0.0;
  ImageType::Pointer      output = NULL;
  floatImageType::Pointer smooth = NULL;
  smooth = Parameterizer->GetFunctionImage();

  std::string fnname = std::string(argv[1]).substr(0, std::string(argv[1]).length() - 4);
  std::string ofn = std::string(argv[4]);
  std::cout << " writing result " << ofn <<  std::endl;
  // writer->SetFileName(ofn.c_str());
  //  writer->SetInput( smooth );
  WriteImage<ImageType>(smooth, ofn.c_str() );
  std::cout << " done writing ";

  return 1;
}
