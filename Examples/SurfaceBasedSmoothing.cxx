

#include "antsUtilities.h"
#include <algorithm>

#include "ReadWriteData.h"
#include "itkSurfaceImageCurvature.h"

namespace ants
{
// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
SurfaceBasedSmoothing(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "SurfaceBasedSmoothing");

  int     argc = args.size();
  std::vector<char *> argvStorage = ArgvFromStrings(args);
  char **             argv = argvStorage.data();

  // antscout->set_stream( out_stream );

  if (argc < 3)
  {
    std::cout << " usage :  " << argv[0] << " ImageToSmooth  sigma SurfaceImage  outname  {numrepeatsofsmoothing=0}"
              << std::endl;
    std::cout << " We assume the SurfaceImage has a label == 1 that defines the surface " << std::endl;
    std::cout << " sigma  defines the geodesic n-hood radius (default: 1.0). numrepeats allows one to use " << std::endl;
    std::cout << " a small geodesic n-hood repeatedly applied many times -- faster computation, same effect " << std::endl;
    return 0;
  }
  using ImageType = itk::Image<float, 3>;
  enum
  {
    ImageDimension = ImageType::ImageDimension
  };
  using floatImageType = itk::Image<float, ImageDimension>;
  using ParamType = itk::SurfaceImageCurvature<ImageType>;
  ParamType::Pointer Parameterizer = ParamType::New();


  //  std::string fn="C://Data//brain15labelimage.img";
  float opt = 0;
  float sig = 1.0;
  //  float thresh=0.0;
  if (argc > 2)
  {
    sig = atof(argv[2]);
  }
  unsigned int numrepeats = 0;
  if (argc > 5)
  {
    numrepeats = std::stoi(argv[5]);
  }

  ImageType::Pointer input;
  ReadImage<ImageType>(input, argv[1]);
  ImageType::Pointer surflabel;
  ReadImage<ImageType>(surflabel, argv[3]);

  Parameterizer->SetInputImage(surflabel);
  Parameterizer->SetFunctionImage(input);
  Parameterizer->SetNeighborhoodRadius(sig);
  if (sig <= 0)
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

  Parameterizer->SetNeighborhoodRadius(sig);
  std::cout << " begin integration NOW " << std::endl;
  Parameterizer->IntegrateFunctionOverSurface(true);
  for (unsigned int i = 0; i < numrepeats; i++)
  {
    Parameterizer->IntegrateFunctionOverSurface(true);
  }
  std::cout << " end integration  " << std::endl;
  // Parameterizer->PostProcessGeometry();

  //  double mn=0.0;
  ImageType::Pointer      output = nullptr;
  floatImageType::Pointer smooth = nullptr;
  smooth = Parameterizer->GetFunctionImage();

  std::string ofn = std::string(argv[4]);
  std::cout << " writing result " << ofn << std::endl;
  // writer->SetFileName(ofn.c_str());
  //  writer->SetInput( smooth );
  ANTs::WriteImage<ImageType>(smooth, ofn.c_str());
  std::cout << " done writing ";

  return 1;
}
} // namespace ants
