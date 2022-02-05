

#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>

#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "vnl/algo/vnl_determinant.h"
#include "itkWarpImageFilter.h"
#include "itkImageFileWriter.h"

#include "itkRescaleIntensityImageFilter.h"
#include "vnl/algo/vnl_determinant.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkVectorCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"

#include "ReadWriteData.h"

namespace ants
{
template <typename TField, typename TImage>
typename TImage::Pointer
GetVectorComponent(typename TField::Pointer field, unsigned int index)
{
  // Initialize the Moving to the displacement field
  using ImageType = TImage;

  typename ImageType::Pointer sfield = AllocImage<ImageType>(field);

  using Iterator = itk::ImageRegionIteratorWithIndex<TField>;
  Iterator vfIter(field, field->GetLargestPossibleRegion());
  for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
  {
    typename TField::PixelType v1 = vfIter.Get();
    sfield->SetPixel(vfIter.GetIndex(), v1[index]);
  }

  return sfield;
}

template <typename TImage>
typename TImage::Pointer
SmoothImage(typename TImage::Pointer image, float sig)
{
  // find min value
  using Iterator = itk::ImageRegionIteratorWithIndex<TImage>;
  Iterator vfIter(image, image->GetLargestPossibleRegion());
  for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
  {
    typename TImage::PixelType v1 = vfIter.Get();
    if (std::isnan(v1))
    {
      vfIter.Set(0);
    }
  }
  using dgf = itk::DiscreteGaussianImageFilter<TImage, TImage>;
  typename dgf::Pointer filter = dgf::New();
  filter->SetVariance(sig);
  filter->SetUseImageSpacing(true);
  filter->SetMaximumError(.01f);
  filter->SetInput(image);
  filter->Update();
  typename TImage::Pointer out = filter->GetOutput();

  return out;
}

template <typename TImage>
void
SmoothDeformation(typename TImage::Pointer vectorimage, float sig)
{
  using VectorType = itk::Vector<float, 3>;
  using ImageType = itk::Image<float, 3>;
  typename ImageType::Pointer subimgx = GetVectorComponent<TImage, ImageType>(vectorimage, 0);
  subimgx = SmoothImage<ImageType>(subimgx, sig);
  typename ImageType::Pointer subimgy = GetVectorComponent<TImage, ImageType>(vectorimage, 1);
  subimgy = SmoothImage<ImageType>(subimgy, sig);
  typename ImageType::Pointer subimgz = GetVectorComponent<TImage, ImageType>(vectorimage, 2);
  subimgz = SmoothImage<ImageType>(subimgz, sig);

  using IteratorType = itk::ImageRegionIteratorWithIndex<TImage>;
  IteratorType Iterator(vectorimage, vectorimage->GetLargestPossibleRegion().GetSize());
  Iterator.GoToBegin();
  while (!Iterator.IsAtEnd())
  {
    VectorType vec;
    vec[0] = subimgx->GetPixel(Iterator.GetIndex());
    vec[1] = subimgy->GetPixel(Iterator.GetIndex());
    vec[2] = subimgz->GetPixel(Iterator.GetIndex());
    Iterator.Set(vec);
    ++Iterator;
  }
}

template <typename TImage, typename TField, typename TInterp, typename TInterp2>
float
IntegrateLength(typename TImage::Pointer gmsurf,
                typename TImage::Pointer /* thickimage */,
                typename TImage::IndexType   velind,
                typename TField::Pointer     lapgrad,
                float                        itime,
                float                        starttime,
                const float                  deltaTime,
                typename TInterp::Pointer    vinterp,
                typename TImage::SpacingType spacing,
                float                        vecsign,
                float                        timesign,
                float                        gradsign)
{
  using VectorType = typename TField::PixelType;
  using DPointType = typename TField::PointType;
  using DefaultInterpolatorType = itk::VectorLinearInterpolateImageFunction<TField, float>;

  VectorType zero;
  zero.Fill(0);
  VectorType disp;
  disp.Fill(0);
  unsigned int                                          ct = 0;
  DPointType                                            pointIn1;
  DPointType                                            pointIn2;
  typename DefaultInterpolatorType::ContinuousIndexType vcontind;
  DPointType                                            pointIn3;
  enum
  {
    ImageDimension = TImage::ImageDimension
  };
  using IndexType = typename TImage::IndexType;
  unsigned int m_NumberOfTimePoints = 2;
  for (unsigned int jj = 0; jj < ImageDimension; jj++)
  {
    pointIn1[jj] = velind[jj] * lapgrad->GetSpacing()[jj];
  }
  itime = starttime;
  bool  timedone = false;
  float totalmag = 0;
  while (!timedone)
  {
    float scale = 1; // *m_DT[timeind]/m_DS[timeind];
    //     std::cout << " scale " << scale << std::endl;
    auto itimetn1 = static_cast<double>(itime - timesign * deltaTime * scale);
    auto itimetn1h = static_cast<double>(itime - timesign * deltaTime * 0.5f * scale);
    if (itimetn1h < 0)
    {
      itimetn1h = 0;
    }
    if (itimetn1h > m_NumberOfTimePoints - 1)
    {
      itimetn1h = m_NumberOfTimePoints - 1;
    }
    if (itimetn1 < 0)
    {
      itimetn1 = 0;
    }
    if (itimetn1 > m_NumberOfTimePoints - 1)
    {
      itimetn1 = m_NumberOfTimePoints - 1;
    }
    // first get current position of particle
    for (unsigned int jj = 0; jj < ImageDimension; jj++)
    {
      pointIn1[jj] = velind[jj] * lapgrad->GetSpacing()[jj];
    }
    //      std::cout << " ind " << index  << std::endl;
    // now index the time varying field at that position.
    typename DefaultInterpolatorType::OutputType f1;
    f1.Fill(0);
    typename DefaultInterpolatorType::OutputType f2;
    f2.Fill(0);
    typename DefaultInterpolatorType::OutputType f3;
    f3.Fill(0);
    typename DefaultInterpolatorType::OutputType f4;
    f4.Fill(0);

    using ContinuousIndexType = typename DefaultInterpolatorType::ContinuousIndexType;
    ContinuousIndexType Y1;
    ContinuousIndexType Y2;
    ContinuousIndexType Y3;
    ContinuousIndexType Y4;
    for (unsigned int jj = 0; jj < ImageDimension; jj++)
    {
      pointIn2[jj] = static_cast<typename DPointType::CoordRepType>(disp[jj]) + pointIn1[jj];
      vcontind[jj] = pointIn2[jj] / lapgrad->GetSpacing()[jj];
      Y1[jj] = vcontind[jj];
      Y2[jj] = vcontind[jj];
      Y3[jj] = vcontind[jj];
      Y4[jj] = vcontind[jj];
    }
    // Y1[ImageDimension]=itimetn1;
    // Y2[ImageDimension]=itimetn1h;
    // Y3[ImageDimension]=itimetn1h;
    //      Y4[ImageDimension]=itime;

    f1 = vinterp->EvaluateAtContinuousIndex(Y1);
    for (unsigned int jj = 0; jj < ImageDimension; jj++)
    {
      Y2[jj] += static_cast<typename ContinuousIndexType::CoordRepType>(static_cast<float>(f1[jj]) * deltaTime * 0.5f);
    }
    bool isinside = true;
    for (unsigned int jj = 0; jj < ImageDimension; jj++)
    {
      if (Y2[jj] < 1 || Y2[jj] > lapgrad->GetLargestPossibleRegion().GetSize()[jj] - 2)
      {
        isinside = false;
      }
    }
    if (isinside)
    {
      f2 = vinterp->EvaluateAtContinuousIndex(Y2);
    }
    for (unsigned int jj = 0; jj < ImageDimension; jj++)
    {
      Y3[jj] += static_cast<typename ContinuousIndexType::CoordRepType>(static_cast<float>(f2[jj]) * deltaTime * 0.5f);
    }
    isinside = true;
    for (unsigned int jj = 0; jj < ImageDimension; jj++)
    {
      if (Y3[jj] < 1 || Y3[jj] > lapgrad->GetLargestPossibleRegion().GetSize()[jj] - 2)
      {
        isinside = false;
      }
    }
    if (isinside)
    {
      f3 = vinterp->EvaluateAtContinuousIndex(Y3);
    }
    for (unsigned int jj = 0; jj < ImageDimension; jj++)
    {
      Y4[jj] += static_cast<typename ContinuousIndexType::CoordRepType>(static_cast<float>(f3[jj]) * deltaTime);
    }
    isinside = true;
    for (unsigned int jj = 0; jj < ImageDimension; jj++)
    {
      if (Y4[jj] < 1 || Y4[jj] > lapgrad->GetLargestPossibleRegion().GetSize()[jj] - 2)
      {
        isinside = false;
      }
    }
    if (isinside)
    {
      f4 = vinterp->EvaluateAtContinuousIndex(Y4);
    }
    using DPointCoordRepType = typename DPointType::CoordRepType;
    auto twoValue = static_cast<DPointCoordRepType>(2.0);
    for (unsigned int jj = 0; jj < ImageDimension; jj++)
    {
      pointIn3[jj] =
        pointIn2[jj] + static_cast<DPointCoordRepType>(gradsign * vecsign * deltaTime / 6.0f) *
                         (static_cast<DPointCoordRepType>(f1[jj]) + twoValue * static_cast<DPointCoordRepType>(f2[jj]) +
                          twoValue * static_cast<DPointCoordRepType>(f3[jj]) + static_cast<DPointCoordRepType>(f4[jj]));
    }

    VectorType out;
    float      mag = 0, dmag = 0, voxmag = 0;
    for (unsigned int jj = 0; jj < ImageDimension; jj++)
    {
      out[jj] = pointIn3[jj] - pointIn1[jj];
      mag += static_cast<float>(itk::Math::sqr(pointIn3[jj] - pointIn2[jj]));
      dmag += static_cast<float>(itk::Math::sqr(pointIn3[jj] - pointIn1[jj]));
      voxmag += static_cast<float>(itk::Math::sqr((pointIn3[jj] - pointIn2[jj]) / spacing[jj]));
      disp[jj] = out[jj];
    }
    voxmag = static_cast<float>(std::sqrt(voxmag));
    dmag = static_cast<float>(std::sqrt(dmag));
    totalmag += static_cast<float>(std::sqrt(mag));

    ct++;
    //      if (!propagate) //thislength=dmag;//
    //         thislength += totalmag;
    itime = itime + deltaTime * timesign;
    IndexType myind;
    for (unsigned int qq = 0; qq < ImageDimension; qq++)
    {
      myind[qq] = (unsigned long)(pointIn3[qq] / spacing[qq] + 0.5);
    }

    if (gmsurf->GetPixel(myind) < 1)
    {
      timedone = true;
    }
    if (ct > 1000)
    {
      std::cout << " stopping b/c exceed 1000 points " << voxmag << std::endl;
      timedone = true;
    }
    if (voxmag < 0.1f)
    {
      timedone = true;
    }
  }

  return totalmag;
}

template <unsigned int ImageDimension>
int
IntegrateVectorField(int argc, char * argv[])
{
  using PixelType = float;
  using VectorType = itk::Vector<float, ImageDimension>;
  using DisplacementFieldType = itk::Image<VectorType, ImageDimension>;
  using ImageType = itk::Image<PixelType, ImageDimension>;
  using SpacingType = typename ImageType::SpacingType;

  constexpr float deltaTime = 0.001;
  float           gradstep = 1.f / deltaTime; // atof(argv[3])*(-1.0);
  std::string     vectorfn = std::string(argv[1]);
  std::string     roifn = std::string(argv[2]);
  int             argct = 3;
  argct++;
  std::string lenoutname = std::string("");
  if (argc > argct)
  {
    lenoutname = std::string(argv[argct]);
  }
  argct++;
  if (argc > argct)
  {
    gradstep *= static_cast<float>(atof(argv[argct]));
  }
  argct++;

  typename ImageType::Pointer ROIimage;
  ReadImage<ImageType>(ROIimage, roifn.c_str());
  typename ImageType::Pointer thickimage;
  ReadImage<ImageType>(thickimage, roifn.c_str());
  thickimage->FillBuffer(0);
  typename DisplacementFieldType::Pointer VECimage;
  ReadImage<DisplacementFieldType>(VECimage, vectorfn.c_str());
  SpacingType spacing = ROIimage->GetSpacing();
  using IteratorType = itk::ImageRegionIteratorWithIndex<ImageType>;
  IteratorType Iterator(ROIimage, ROIimage->GetLargestPossibleRegion().GetSize());

  double timezero = 0;         // 1
  double timeone = 1;          // (s[ImageDimension]-1-timezero);
  float  starttime = timezero; // timezero;
  float  finishtime = timeone; // s[ImageDimension]-1;//timeone;

  typename DisplacementFieldType::IndexType velind;
  float                                     timesign = 1.0;
  if (starttime > finishtime)
  {
    timesign = -1.0;
  }
  using TimeVaryingVelocityFieldType = DisplacementFieldType;
  // UNUSED: typedef typename DisplacementFieldType::PointType                                      DPointType;
  using DefaultInterpolatorType = itk::VectorLinearInterpolateImageFunction<TimeVaryingVelocityFieldType, float>;
  typename DefaultInterpolatorType::Pointer vinterp = DefaultInterpolatorType::New();
  using ScalarInterpolatorType = itk::LinearInterpolateImageFunction<ImageType, float>;
  VectorType zero;
  zero.Fill(0);

  using VIteratorType = itk::ImageRegionIteratorWithIndex<DisplacementFieldType>;
  VIteratorType VIterator(VECimage, VECimage->GetLargestPossibleRegion().GetSize());
  VIterator.GoToBegin();
  while (!VIterator.IsAtEnd())
  {
    VectorType vec = VIterator.Get();
    float      mag = 0;
    for (unsigned int qq = 0; qq < ImageDimension; qq++)
    {
      mag += vec[qq] * vec[qq];
    }
    mag = sqrt(mag);
    if (mag > 0)
    {
      vec = vec / mag;
    }
    VIterator.Set(vec * gradstep);
    ++VIterator;
  }

  Iterator.GoToBegin();
  while (!Iterator.IsAtEnd())
  {
    velind = Iterator.GetIndex();
    float      itime = starttime;
    VectorType disp;
    disp.Fill(0.0);
    if (itk::Math::FloatAlmostEqual(ROIimage->GetPixel(velind), static_cast<PixelType>(2)))
    {
      vinterp->SetInputImage(VECimage);
      float  gradsign = -1.0;
      double vecsign = -1.0;
      float  len1 =
        IntegrateLength<ImageType, DisplacementFieldType, DefaultInterpolatorType, ScalarInterpolatorType>(ROIimage,
                                                                                                           thickimage,
                                                                                                           velind,
                                                                                                           VECimage,
                                                                                                           itime,
                                                                                                           starttime,
                                                                                                           deltaTime,
                                                                                                           vinterp,
                                                                                                           spacing,
                                                                                                           vecsign,
                                                                                                           gradsign,
                                                                                                           timesign);

      gradsign = 1.0;
      vecsign = 1;
      const float len2 =
        IntegrateLength<ImageType, DisplacementFieldType, DefaultInterpolatorType, ScalarInterpolatorType>(ROIimage,
                                                                                                           thickimage,
                                                                                                           velind,
                                                                                                           VECimage,
                                                                                                           itime,
                                                                                                           starttime,
                                                                                                           deltaTime,
                                                                                                           vinterp,
                                                                                                           spacing,
                                                                                                           vecsign,
                                                                                                           gradsign,
                                                                                                           timesign);

      float totalength = len1 + len2;
      thickimage->SetPixel(velind, totalength);
      if ((totalength) > 0)
      {
        std::cout << " len1 " << len1 << " len2 " << len2 << " ind " << velind << std::endl;
      }
    }
    ++Iterator;
  }

  ANTs::WriteImage<ImageType>(thickimage, lenoutname.c_str());

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
ANTSIntegrateVectorField(std::vector<std::string> args, std::ostream * /*out_stream = nullptr*/)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "ANTSIntegrateVectorField");
  int     argc = args.size();
  char ** argv = new char *[args.size() + 1];
  for (unsigned int i = 0; i < args.size(); ++i)
  {
    // allocate space for the string plus a null character
    argv[i] = new char[args[i].length() + 1];
    std::strncpy(argv[i], args[i].c_str(), args[i].length());
    // place the null character in the end
    argv[i][args[i].length()] = '\0';
  }
  argv[argc] = nullptr;
  // class to automatically cleanup argv upon destruction
  class Cleanup_argv
  {
  public:
    Cleanup_argv(char ** argv_, int argc_plus_one_)
      : argv(argv_)
      , argc_plus_one(argc_plus_one_)
    {}

    ~Cleanup_argv()
    {
      for (unsigned int i = 0; i < argc_plus_one; ++i)
      {
        delete[] argv[i];
      }
      delete[] argv;
    }

  private:
    char **      argv;
    unsigned int argc_plus_one;
  };
  Cleanup_argv cleanup_argv(argv, argc + 1);

  // antscout->set_stream( out_stream );

  if (argc < 4)
  {
    std::cout << "Usage:   " << argv[0]
              << "  VecImageIN.nii.gz ROIMaskIN.nii.gz FibersOUT.vtk  LengthImageOUT.nii.gz   " << std::endl;
    std::cout << " The vector field should have vectors as voxels , the ROI is an integer image, fibers out will be "
                 "vtk text files .... "
              << std::endl;
    std::cout << "  ROI-Mask controls where the integration is performed and the start point region ... " << std::endl;
    std::cout << " e.g. the brain will have value 1 , the ROI has value 2 , then all starting seed points "
              << std::endl;
    std::cout << " for the integration will start in the region labeled 2 and be constrained to the region labeled 1. "
              << std::endl;
    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  std::string               ifn = std::string(argv[1]);
  itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(ifn.c_str(), itk::IOFileModeEnum::ReadMode);
  imageIO->SetFileName(ifn.c_str());
  imageIO->ReadImageInformation();
  unsigned int dim = imageIO->GetNumberOfDimensions();

  switch (dim)
  {
    case 2:
    {
      IntegrateVectorField<2>(argc, argv);
    }
    break;
    case 3:
    {
      IntegrateVectorField<3>(argc, argv);
    }
    break;
    case 4:
    {
      IntegrateVectorField<4>(argc, argv);
    }
    break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
} // namespace ants
