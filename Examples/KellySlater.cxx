

#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>

#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "vnl/algo/vnl_determinant.h"

#include "itkWarpImageFilter.h"
#include "itkWarpImageMultiTransformFilter.h"
#include "itkDisplacementFieldFromMultiTransformFilter.h"
#include "itkFastMarchingUpwindGradientImageFilter.h"
#include "itkFastMarchingUpwindGradientImageFilter.h"
#include "itkImageFileWriter.h"

#include "itkANTSImageRegistrationOptimizer.h"
#include "vnl/algo/vnl_determinant.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkCovariantVector.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkVectorCurvatureAnisotropicDiffusionImageFilter.h"

#include "ReadWriteData.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"

#include "itkGradientRecursiveGaussianImageFilter.h"

#include "itkCentralDifferenceImageFunction.h"
#include "itkSurfaceCurvatureBase.h"
#include "itkSurfaceImageCurvature.h"

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
MaurerDistanceMap(typename TImage::PixelType pixlo, typename TImage::PixelType pixhi, typename TImage::Pointer input)
{
  // std::cout << " DDMap " << std::endl;

  using ImageType = TImage;

  using FilterType = itk::SignedMaurerDistanceMapImageFilter<ImageType, ImageType>;

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetSquaredDistance(false);
  //  filter->InputIsBinaryOn();
  filter->SetUseImageSpacing(true);
  filter->SetBackgroundValue(0);
  filter->SetInput(BinaryThreshold<TImage>(pixlo, pixhi, pixhi, input));
  filter->Update();

  //  ANTs::WriteImage<ImageType>(filter->GetOutput(),"temp1.nii");
  return filter->GetOutput();
}

template <typename TImage>
typename TImage::Pointer
SmoothImage(typename TImage::Pointer image, double sig)
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
SmoothDeformation(typename TImage::Pointer vectorimage, double sig)
{
  enum
  {
    ImageDimension = TImage::ImageDimension
  };
  using RealType = typename TImage::PixelType;
  using VectorType = itk::Vector<RealType, ImageDimension>;
  using ImageType = itk::Image<RealType, ImageDimension>;
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

// this has to have never been called because it doesn't actually
// copy anything
template <typename TImage, typename TDisplacementField>
typename TImage::Pointer
CopyImage(TDisplacementField * field)
{
  enum
  {
    ImageDimension = TImage::ImageDimension
  };
  //  unsigned int row=0;
  // unsigned int col=0;
  using PixelType = typename TImage::PixelType;
  using RealImageType = itk::Image<PixelType, ImageDimension>;

  typename RealImageType::Pointer m_RealImage = nullptr;
  m_RealImage = AllocImage<RealImageType>(field, 0);

  return m_RealImage;
}

template <typename TImage>
typename TImage::Pointer
LabelSurface(typename TImage::PixelType foreground,
             typename TImage::PixelType newval,
             typename TImage::Pointer   input,
             double                     distthresh)
{
  std::cout << " Label Surf " << std::endl;
  using ImageType = TImage;
  enum
  {
    ImageDimension = ImageType::ImageDimension
  };
  // ORIENTATION ALERT: Original code set origin & spacing from
  // examplar without also setting directions.
  typename ImageType::Pointer Image = AllocImage<ImageType>(input, 0.0);

  using iteratorType = itk::NeighborhoodIterator<ImageType>;

  typename iteratorType::RadiusType rad;
  for (int j = 0; j < ImageDimension; j++)
  {
    rad[j] = (unsigned int)(distthresh + 0.5);
  }
  iteratorType GHood(rad, input, input->GetLargestPossibleRegion());

  GHood.GoToBegin();

  //  std::cout << " foreg " << (int) foreground;
  while (!GHood.IsAtEnd())
  {
    typename TImage::PixelType p = GHood.GetCenterPixel();
    typename TImage::IndexType ind = GHood.GetIndex();
    typename TImage::IndexType ind2;
    if (itk::Math::FloatAlmostEqual(p, foreground))
    {
      bool atedge = false;
      for (unsigned int i = 0; i < GHood.Size(); i++)
      {
        ind2 = GHood.GetIndex(i);
        double dist = 0.0;
        for (int j = 0; j < ImageDimension; j++)
        {
          dist += (double)(ind[j] - ind2[j]) * (double)(ind[j] - ind2[j]);
        }
        dist = sqrt(dist);
        if (!itk::Math::FloatAlmostEqual(GHood.GetPixel(i), foreground) && dist < distthresh)
        {
          atedge = true;
        }
      }
      if (atedge && itk::Math::FloatAlmostEqual(p, foreground))
      {
        Image->SetPixel(ind, newval);
      }
      else
      {
        Image->SetPixel(ind, 0);
      }
    }
    ++GHood;
  }

  return Image;
}

template <typename TImage, typename TField>
typename TField::Pointer
LaplacianGrad(typename TImage::Pointer wm, typename TImage::Pointer gm, double sig)
{
  using IndexType = typename TImage::IndexType;
  IndexType ind;
  using ImageType = TImage;
  using GradientImageType = TField;
  using GradientImageFilterType = itk::GradientRecursiveGaussianImageFilter<ImageType, GradientImageType>;
  using GradientImageFilterPointer = typename GradientImageFilterType::Pointer;

  typename TField::Pointer sfield = AllocImage<TField>(wm);

  typename TImage::Pointer laplacian = SmoothImage<TImage>(wm, 3);
  laplacian->FillBuffer(0);
  using IteratorType = itk::ImageRegionIteratorWithIndex<TImage>;
  IteratorType Iterator(wm, wm->GetLargestPossibleRegion().GetSize());
  Iterator.GoToBegin();

  // initialize L(wm)=1, L(gm)=0.5, else 0
  while (!Iterator.IsAtEnd())
  {
    ind = Iterator.GetIndex();
    if (!itk::Math::FloatAlmostEqual(wm->GetPixel(ind), itk::NumericTraits<typename ImageType::PixelType>::ZeroValue()))
    {
      laplacian->SetPixel(ind, 1);
    }
    else
    {
      laplacian->SetPixel(ind, 0.);
    }
    ++Iterator;
  }
  // smooth and then reset the values
  for (unsigned int iterations = 0; iterations < 100; iterations++)
  {
    laplacian = SmoothImage<TImage>(laplacian, sqrt(sig));
    Iterator.GoToBegin();
    while (!Iterator.IsAtEnd())
    {
      ind = Iterator.GetIndex();
      if (!itk::Math::FloatAlmostEqual(wm->GetPixel(ind),
                                       itk::NumericTraits<typename ImageType::PixelType>::ZeroValue()))
      {
        laplacian->SetPixel(ind, itk::NumericTraits<typename ImageType::PixelType>::OneValue());
      }
      else if (itk::Math::FloatAlmostEqual(gm->GetPixel(ind),
                                           itk::NumericTraits<typename ImageType::PixelType>::ZeroValue()) &&
               itk::Math::FloatAlmostEqual(wm->GetPixel(ind),
                                           itk::NumericTraits<typename ImageType::PixelType>::ZeroValue()))
      {
        laplacian->SetPixel(ind, itk::NumericTraits<typename ImageType::PixelType>::ZeroValue());
      }
      ++Iterator;
    }
  }

  GradientImageFilterPointer filter = GradientImageFilterType::New();
  filter->SetInput(laplacian);
  filter->SetSigma(sig);
  filter->Update();
  return filter->GetOutput();
}

template <typename TImage, typename TField>
typename TField::Pointer
ExpDiffMap(typename TField::Pointer velofield, typename TImage::Pointer wm, double sign, unsigned int numtimepoints)
{
  using ImageType = TImage;
  using DisplacementFieldType = TField;
  using PixelType = typename TField::PixelType;
  typename TField::PixelType zero, disp;
  enum
  {
    ImageDimension = TImage::ImageDimension
  };
  disp.Fill(0);
  zero.Fill(0);

  typename DisplacementFieldType::Pointer incrfield = AllocImage<DisplacementFieldType>(velofield, zero);

  using IteratorType = itk::ImageRegionIteratorWithIndex<ImageType>;
  IteratorType Iterator(wm, wm->GetLargestPossibleRegion().GetSize());
  Iterator.GoToBegin();
  while (!Iterator.IsAtEnd())
  {
    incrfield->SetPixel(Iterator.GetIndex(), velofield->GetPixel(Iterator.GetIndex()) * sign);
    ++Iterator;
  }

  // generate phi
  using AffineTransformType = itk::MatrixOffsetTransformBase<PixelType, ImageDimension, ImageDimension>;
  using WarperType = itk::DisplacementFieldFromMultiTransformFilter<TField, TField, AffineTransformType>;
  typename WarperType::Pointer warper = WarperType::New();
  warper->SetOutputParametersFromImage(velofield);
  warper->DetermineFirstDeformNoInterp();

  unsigned int ttiter = 0;
  while (ttiter < numtimepoints) // 10 time integration points
  {
    ttiter++;
    warper->PushBackDisplacementFieldTransform(incrfield);
  }

  warper->Update();
  return warper->GetOutput();
}

template <typename TImage, typename TField>
typename TField::Pointer
DiReCTCompose(typename TField::Pointer velofield, typename TField::Pointer diffmap)
{
  using PixelType = typename TField::PixelType;
  typename TField::PixelType zero, disp;
  enum
  {
    ImageDimension = TImage::ImageDimension
  };
  disp.Fill(0);
  zero.Fill(0);

  using AffineTransformType = itk::MatrixOffsetTransformBase<PixelType, ImageDimension, ImageDimension>;
  using WarperType = itk::DisplacementFieldFromMultiTransformFilter<TField, TField, AffineTransformType>;
  typename WarperType::Pointer warper = WarperType::New();
  warper->SetOutputParametersFromImage(velofield);
  warper->DetermineFirstDeformNoInterp();
  warper->PushBackDisplacementFieldTransform(diffmap);
  warper->PushBackDisplacementFieldTransform(velofield);
  warper->Update();
  return warper->GetOutput();
}

template <typename TImage, typename TField>
void
InvertField(typename TField::Pointer field,
            typename TField::Pointer inverseFieldIN,
            double                   weight = 1.0,
            double                   toler = 0.1,
            int                      maxiter = 20,
            bool /* print */ = false)
{
  enum
  {
    ImageDimension = TImage::ImageDimension
  };
  using DisplacementFieldType = TField;
  using DisplacementFieldPointer = typename TField::Pointer;
  using VectorType = typename TField::PixelType;
  using ImageType = TImage;
  using ImagePointer = typename TImage::Pointer;
  double       mytoler = toler;
  unsigned int mymaxiter = maxiter;

  VectorType zero;
  zero.Fill(0);
  //  if (this->GetElapsedIterations() < 2 ) maxiter=10;

  ImagePointer realImage = AllocImage<ImageType>(field);

  using VectorType = typename DisplacementFieldType::PixelType;
  using IndexType = typename DisplacementFieldType::IndexType;
  using Iterator = itk::ImageRegionIteratorWithIndex<DisplacementFieldType>;

  using ROType = itk::ANTSImageRegistrationOptimizer<ImageDimension, double>;
  typename ROType::Pointer m_MFR = ROType::New();

  DisplacementFieldPointer inverseField = AllocImage<DisplacementFieldType>(field, zero);

  DisplacementFieldPointer lagrangianInitCond = AllocImage<DisplacementFieldType>(field);

  DisplacementFieldPointer eulerianInitCond = AllocImage<DisplacementFieldType>(field);

  using SizeType = typename DisplacementFieldType::SizeType;
  SizeType size = field->GetLargestPossibleRegion().GetSize();

  typename ImageType::SpacingType spacing = field->GetSpacing();
  unsigned long                   npix = 1;
  for (int j = 0; j < ImageDimension; j++) // only use in-plane spacing
  {
    npix *= field->GetLargestPossibleRegion().GetSize()[j];
  }

  double   max = 0;
  Iterator iter(field, field->GetLargestPossibleRegion());
  for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
  {
    IndexType  index = iter.GetIndex();
    VectorType vec1 = iter.Get();
    VectorType newvec = vec1 * weight;
    lagrangianInitCond->SetPixel(index, newvec);
    inverseField->SetPixel(index, inverseFieldIN->GetPixel(index));

    double mag = 0;
    for (unsigned int jj = 0; jj < ImageDimension; jj++)
    {
      mag += newvec[jj] * newvec[jj];
    }
    mag = sqrt(mag);
    if (mag > max)
    {
      max = mag;
    }
  }

  eulerianInitCond->FillBuffer(zero);

  double scale = (1.) / max;
  if (scale > 1.)
  {
    scale = 1.0;
  }
  //    double initscale=scale;
  Iterator vfIter(inverseField, inverseField->GetLargestPossibleRegion());

  //  int num=10;
  //  for (int its=0; its<num; its++)
  double       difmag = 10.0;
  unsigned int ct = 0;
  double       meandif = 1.e8;
  //    int badct=0;
  //  while (difmag > subpix && meandif > subpix*0.1 && badct < 2 )//&& ct < 20 && denergy > 0)
  //    double length=0.0;
  double stepl = 2.;

  double epsilon = (double)size[0] / 256;
  if (epsilon > 1)
  {
    epsilon = 1;
  }

  while (difmag > mytoler && ct < mymaxiter && meandif > 0.001)
  {
    meandif = 0.0;

    // this field says what position the eulerian field should contain in the E domain
    m_MFR->ComposeDiffs(inverseField, lagrangianInitCond, eulerianInitCond, 1);
    difmag = 0.0;
    for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
    {
      IndexType  index = vfIter.GetIndex();
      VectorType update = eulerianInitCond->GetPixel(index);
      double     mag = 0;
      for (int j = 0; j < ImageDimension; j++)
      {
        update[j] *= (-1.0);
        mag += (update[j] / spacing[j]) * (update[j] / spacing[j]);
      }
      mag = sqrt(mag);
      meandif += mag;
      if (mag > difmag)
      {
        difmag = mag;
      }
      //      if (mag < 1.e-2) update.Fill(0);

      eulerianInitCond->SetPixel(index, update);
      realImage->SetPixel(index, mag);
    }
    meandif /= (double)npix;
    if (ct == 0)
    {
      epsilon = 0.75;
    }
    else
    {
      epsilon = 0.5;
    }
    stepl = difmag * epsilon;
    for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
    {
      double     val = realImage->GetPixel(vfIter.GetIndex());
      VectorType update = eulerianInitCond->GetPixel(vfIter.GetIndex());
      if (val > stepl)
      {
        update = update * (stepl / val);
      }
      VectorType upd = vfIter.Get() + update * (epsilon);
      vfIter.Set(upd);
    }
    ct++;
  }
  for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
  {
    IndexType index = iter.GetIndex();
    inverseFieldIN->SetPixel(index, inverseField->GetPixel(index));
  }

  //    std::cout <<" difmag " << difmag << ": its " << ct <<  " len " << m_MFR->MeasureDeformation(inverseField ) <<
  // std::endl;
}

template <unsigned int ImageDimension>
int
LaplacianThicknessExpDiff2(int argc, char * argv[])
{
  int         argct = 2;
  std::string segfn = std::string(argv[argct]);
  argct++;

  std::string wfn = std::string(argv[argct]);
  argct++;
  std::string gfn = std::string(argv[argct]);
  argct++;
  std::string outname = std::string(argv[argct]);
  argct++;
  unsigned int numtimepoints = 10;

  using RealType = double;
  RealType gradstep = (RealType)(-1.0) * 0.5; // (ImageDimension-1);
  if (argc > argct)
  {
    gradstep = atof(argv[argct]) * (-1.0);
  }
  gradstep *= 1.0 / (RealType)numtimepoints * 10;
  argct++;
  unsigned int alltheits = 50;
  if (argc > argct)
  {
    alltheits = std::stoi(argv[argct]);
  }
  argct++;
  RealType thickprior = 6.0;
  if (argc > argct)
  {
    thickprior = atof(argv[argct]);
  }
  argct++;
  // bool useCurvaturePrior = false;
  // if( argc > argct )
  //   {
  //   useCurvaturePrior = std::stoi(argv[argct]);
  //   }
  argct++;
  RealType smoothingsigma = 1.5;
  if (argc > argct)
  {
    smoothingsigma = atof(argv[argct]);
  }
  argct++;
  // bool useEuclidean = true;
  // if( argc > argct )
  //   {
  //   useEuclidean = std::stoi(argv[argct]);
  //   }
  argct++;
  std::cout << " smooth " << smoothingsigma << " thp " << thickprior << " gs " << gradstep << std::endl;
  using PixelType = RealType;
  using VectorType = itk::Vector<RealType, ImageDimension>;
  using DisplacementFieldType = itk::Image<VectorType, ImageDimension>;
  using ImageType = itk::Image<PixelType, ImageDimension>;
  using IndexType = typename ImageType::IndexType;
  using ROType = itk::ANTSImageRegistrationOptimizer<ImageDimension, RealType>;
  typename ROType::Pointer m_MFR = ROType::New();

  typename ImageType::Pointer segmentationimage;
  ReadImage<ImageType>(segmentationimage, segfn.c_str());
  typename ImageType::Pointer wm;
  ReadImage<ImageType>(wm, wfn.c_str());
  typename ImageType::DirectionType omat = wm->GetDirection();
  typename ImageType::DirectionType fmat = wm->GetDirection();
  fmat.SetIdentity();
  std::cout << " Setting Identity Direction  " << fmat << std::endl;
  wm->SetDirection(fmat);
  typename ImageType::Pointer totalimage;
  ReadImage<ImageType>(totalimage, wfn.c_str());
  totalimage->SetDirection(fmat);
  typename ImageType::Pointer hitimage;
  ReadImage<ImageType>(hitimage, wfn.c_str());
  hitimage->SetDirection(fmat);
  typename ImageType::Pointer gm;
  ReadImage<ImageType>(gm, gfn.c_str());
  gm->SetDirection(fmat);
  wm->SetDirection(fmat);
  segmentationimage->SetDirection(fmat);

  typename DisplacementFieldType::Pointer lapgrad;
  typename ImageType::Pointer             gmb = BinaryThreshold<ImageType>(2, 2, 1, segmentationimage); // fixme
  typename ImageType::Pointer             wmb = BinaryThreshold<ImageType>(3, 3, 1, segmentationimage); // fixme
  typename ImageType::Pointer             laplacian = SmoothImage<ImageType>(wm, smoothingsigma);
  lapgrad = LaplacianGrad<ImageType, DisplacementFieldType>(wmb, gmb, 1);

  VectorType zero;
  zero.Fill(0);
  typename DisplacementFieldType::Pointer corrfield = AllocImage<DisplacementFieldType>(wm, zero);

  typename DisplacementFieldType::Pointer incrfield = AllocImage<DisplacementFieldType>(wm, zero);

  typename DisplacementFieldType::Pointer invfield = AllocImage<DisplacementFieldType>(wm, zero);

  typename DisplacementFieldType::Pointer incrinvfield = AllocImage<DisplacementFieldType>(wm, zero);

  typename DisplacementFieldType::Pointer velofield = AllocImage<DisplacementFieldType>(wm, zero);

  //  LabelSurface(typename TImage::PixelType foreground,
  //       typename TImage::PixelType newval, typename TImage::Pointer input, RealType distthresh )
  RealType                    distthresh = 1.5;
  typename ImageType::Pointer wmgrow = Morphological<ImageType>(wmb, 0, 1, 1);
  typename ImageType::Pointer bsurf = LabelSurface<ImageType>(1, 1, wmgrow, distthresh); // or wmb ?
  typename ImageType::Pointer speedprior = nullptr;
  ANTs::WriteImage<ImageType>(bsurf, "surf.nii.gz");
  //    typename RealTypeImageType::Pointer distfromboundary =
  //  typename ImageType::Pointer surf=MaurerDistanceMap<ImageType>(0.5,1.e9,bsurf);
  // surf= SmoothImage<ImageType>(surf,3);
  typename ImageType::Pointer finalthickimage = BinaryThreshold<ImageType>(3, 3, 1, segmentationimage); // fixme

  gmb = BinaryThreshold<ImageType>(2, 3, 1, segmentationimage); // fixme
  typename ImageType::Pointer gmgrow = Morphological<ImageType>(gmb, 1, 1, 1);
  typename ImageType::Pointer gmsurf = LabelSurface<ImageType>(1, 1, gmgrow, distthresh); // or wmb ?
  //  ANTs::WriteImage<ImageType>(gmsurf,"surfdefgm.nii.gz");
  //  ANTs::WriteImage<ImageType>(bsurf,"surfdefwm.nii.gz");

  using TimeVaryingVelocityFieldType = DisplacementFieldType;
  using DefaultInterpolatorType = itk::VectorLinearInterpolateImageFunction<TimeVaryingVelocityFieldType, RealType>;
  typename DefaultInterpolatorType::Pointer vinterp = DefaultInterpolatorType::New();
  vinterp->SetInputImage(lapgrad);
  using ScalarInterpolatorType = itk::LinearInterpolateImageFunction<ImageType, RealType>;
  typename ScalarInterpolatorType::Pointer ginterp = ScalarInterpolatorType::New();
  typename ScalarInterpolatorType::Pointer winterp = ScalarInterpolatorType::New();
  winterp->SetInputImage(wm);
  ginterp->SetInputImage(gm);

  using IteratorType = itk::ImageRegionIteratorWithIndex<ImageType>;
  using VIteratorType = itk::ImageRegionIteratorWithIndex<DisplacementFieldType>;
  VIteratorType VIterator(lapgrad, lapgrad->GetLargestPossibleRegion().GetSize());
  VIterator.GoToBegin();
  while (!VIterator.IsAtEnd())
  {
    // the velocity field solution value
    VectorType vec = VIterator.Get();
    RealType   mag = 0;
    for (unsigned dd = 0; dd < ImageDimension; dd++)
    {
      mag += vec[dd] * vec[dd];
    }
    mag = sqrt(mag);
    //      if (mag > 0) vec=vec/mag;
    VIterator.Set((vec)*gradstep);
    ++VIterator;
  }

  //  m_MFR->SmoothDisplacementFieldGauss(lapgrad,1.7);
  std::cout << " Scaling done " << std::endl;

  typename ImageType::Pointer thickimage = laplacian;
  VectorType                  disp;
  VectorType                  incdisp;
  disp.Fill(0.0);
  incdisp.Fill(0.0);
  IteratorType Iterator(wm, wm->GetLargestPossibleRegion().GetSize());
  RealType     totalerr = 1.e8, lasterr = 1.e10;
  unsigned     its = 0;
  wmgrow->FillBuffer(0);
  RealType      dmag = 0;
  RealType      thicknesserror = 0;
  unsigned long thickerrct = 0;
  unsigned int  badct = 0;
  RealType      thickoffset = 0;
  bool          checknans = true;

  while (its < alltheits && badct < 4)
  {
    its++;
    if (totalerr > lasterr)
    {
      badct++;
      std::cout << " badct " << badct << std::endl;
    }
    else
    {
      badct = 0;
    }
    lasterr = totalerr;
    // Optimization Error initialized for this iteration
    totalerr = 0;
    incrfield->FillBuffer(zero);
    incrfield->FillBuffer(zero);
    incrinvfield->FillBuffer(zero);

    // generate phi
    //      corrfield->FillBuffer(zero);
    invfield->FillBuffer(zero);
    unsigned int ttiter = 0;

    thickimage->FillBuffer(0);
    hitimage->FillBuffer(0);
    totalimage->FillBuffer(0);
    thicknesserror = 0;
    thickerrct = 1;
    bool                        debug = false;
    bool                        spatprior = false;
    typename ImageType::Pointer priorim = nullptr;
    if (speedprior)
    {
      spatprior = true;
      priorim = speedprior;
    }
    typename ImageType::Pointer wpriorim = nullptr;
    RealType                    origthickprior = thickprior;

    while (ttiter < numtimepoints) // N time integration points
    {
      //      m_MFR->Compose(incrinvfield,invfield,nullptr);
      m_MFR->ComposeDiffs(invfield, incrinvfield, invfield, 1);

      if (debug)
      {
        std::cout << " exp " << std::endl;
      }
      // Integrate the negative velocity field to generate diffeomorphism corrfield step 3(a)
      //      corrfield=ExpDiffMap<ImageType,DisplacementFieldType>( velofield,  wm, -1, numtimepoints-ttiter);
      //      std::cout  << " corrf len " << m_MFR->MeasureDeformation( corrfield ) << std::endl;
      if (debug)
      {
        std::cout << " gmdef " << std::endl;
      }
      typename ImageType::Pointer gmdef = gm; // m_MFR->WarpImageBackward(gm,corrfield);
      totalerr = 0;

      typename ImageType::Pointer surfdef = m_MFR->WarpImageBackward(wm, invfield);
      if (debug)
      {
        std::cout << " thkdef " << std::endl;
      }
      typename ImageType::Pointer thkdef = m_MFR->WarpImageBackward(thickimage, invfield);
      if (debug)
      {
        std::cout << " thindef " << std::endl;
      }
      typename ImageType::Pointer thindef = m_MFR->WarpImageBackward(bsurf, invfield);
      if (spatprior)
      {
        wpriorim = m_MFR->WarpImageBackward(priorim, invfield);
      }

      using GradientImageType = DisplacementFieldType;
      using GradientImageFilterType = itk::GradientRecursiveGaussianImageFilter<ImageType, GradientImageType>;
      using GradientImageFilterPointer = typename GradientImageFilterType::Pointer;
      GradientImageFilterPointer gfilter = GradientImageFilterType::New();
      gfilter->SetInput(surfdef);
      gfilter->SetSigma(smoothingsigma);
      gfilter->Update();
      typename DisplacementFieldType::Pointer lapgrad2 = gfilter->GetOutput();

      // this is the "speed" image
      typename ImageType::Pointer speed_image = CopyImage<ImageType, DisplacementFieldType>(invfield);
      IteratorType                xxIterator(speed_image, speed_image->GetLargestPossibleRegion().GetSize());
      xxIterator.GoToBegin();
      RealType maxlapgrad2mag = 0;
      while (!xxIterator.IsAtEnd())
      {
        typename ImageType::IndexType speedindex = xxIterator.GetIndex();
        if (itk::Math::FloatAlmostEqual(segmentationimage->GetPixel(speedindex),
                                        static_cast<typename ImageType::PixelType>(2.0))) // fixme
        {
          thickprior = origthickprior;
          VectorType wgradval = lapgrad2->GetPixel(speedindex);
          RealType   wmag = 0;
          for (unsigned kq = 0; kq < ImageDimension; kq++)
          {
            wmag += wgradval[kq] * wgradval[kq];
          }
          if (fabs(wmag) < 1.e-6)
          {
            wmag = 0;
          }
          wmag = sqrt(wmag);
          if (checknans)
          {
            if (std::isnan(wmag) || std::isinf(wmag) ||
                itk::Math::FloatAlmostEqual(wmag, itk::NumericTraits<RealType>::ZeroValue()))
            {
              wgradval.Fill(0);
              lapgrad2->SetPixel(speedindex, wgradval);
              wmag = 0;
            }
            else
            {
              lapgrad2->SetPixel(speedindex, wgradval / wmag);
            }
          }
          totalerr += fabs(surfdef->GetPixel(speedindex) - gmdef->GetPixel(speedindex));
          //          RealType thkval=thkdef->GetPixel(speedindex);
          //          RealType thkval=finalthickimage->GetPixel(speedindex);
          //          RealType fval=1; //(thickprior-thkval);
          //          if ( fval > 0 ) fval=1; else fval=-1;
          // speed function here IMPORTANT!!
          RealType dd = (surfdef->GetPixel(speedindex) - gmdef->GetPixel(speedindex)) * gradstep;
          dd *= gm->GetPixel(speedindex);
          if (checknans)
          {
            if (std::isnan(dd) || std::isinf(dd))
            {
              dd = 0;
            }
          }
          speed_image->SetPixel(speedindex, dd);
          if (wmag * dd > maxlapgrad2mag)
          {
            maxlapgrad2mag = wmag * dd;
          }
        }
        else
        {
          speed_image->SetPixel(speedindex, 0);
        }
        ++xxIterator;
      }

      if (maxlapgrad2mag < 1.e-4)
      {
        maxlapgrad2mag = 1.e9;
      }
      if (ttiter == numtimepoints - 1)
      {
        if (ImageDimension == 2)
        {
          ANTs::WriteImage<ImageType>(surfdef, "surfdef.nii.gz");
        }
        if (ImageDimension == 2)
        {
          ANTs::WriteImage<ImageType>(thindef, "thindef.nii.gz");
        }
        if (ImageDimension == 2)
        {
          ANTs::WriteImage<ImageType>(gmdef, "gmdef.nii.gz");
        }
        if (ImageDimension == 2)
        {
          ANTs::WriteImage<ImageType>(thkdef, "thick2.nii.gz");
        }
      }
      /* Now that we have the gradient image, we need to visit each voxel and compute objective function */
      //      std::cout << " maxlapgrad2mag " << maxlapgrad2mag << std::endl;
      Iterator.GoToBegin();
      while (!Iterator.IsAtEnd())
      {
        typename DisplacementFieldType::IndexType velind = Iterator.GetIndex();
        VectorType wgradval = lapgrad2->GetPixel(velind); // *5.0/(maxlapgrad2mag*(RealType)numtimepoints);
        disp = wgradval * speed_image->GetPixel(velind);
        incrfield->SetPixel(velind, incrfield->GetPixel(velind) + disp);

        if (ttiter == 0) // make euclidean distance image
        {
          dmag = 0;
          disp = corrfield->GetPixel(velind);
          for (unsigned int jj = 0; jj < ImageDimension; jj++)
          {
            dmag += disp[jj] * disp[jj];
          }
          RealType bval = bsurf->GetPixel(velind);
          if (checknans)
          {
            if (std::isnan(dmag) || std::isinf(dmag))
            {
              dmag = 0;
            }
            if (std::isnan(bval) || std::isinf(bval))
            {
              bval = 0;
            }
          }
          /** Change 2-26-2010 = incoporate gm prob in length ... */
          dmag = sqrt(dmag) * bval; // *gm->GetPixel(velind);
          thickimage->SetPixel(velind, dmag);
          totalimage->SetPixel(velind, dmag);
          hitimage->SetPixel(velind, bval);
        }
        else if (itk::Math::FloatAlmostEqual(segmentationimage->GetPixel(velind),
                                             static_cast<typename ImageType::PixelType>(2.0))) // fixme
        {
          RealType thkval = thkdef->GetPixel(velind);
          RealType putval = thindef->GetPixel(velind);
          hitimage->SetPixel(velind, hitimage->GetPixel(velind) + putval);
          totalimage->SetPixel(velind, totalimage->GetPixel(velind) + thkval);
        }

        ++Iterator;
      }

      Iterator.GoToBegin();
      while (!Iterator.IsAtEnd())
      {
        IndexType velind = Iterator.GetIndex();
        bool      shouldbezero = false;
        if (itk::Math::FloatAlmostEqual(segmentationimage->GetPixel(velind),
                                        itk::NumericTraits<typename ImageType::PixelType>::ZeroValue())) // fixme
        {
          shouldbezero = true;
        }
        if (!shouldbezero)
        {
          if (itk::Math::FloatAlmostEqual(bsurf->GetPixel(velind),
                                          itk::NumericTraits<typename ImageType::PixelType>::ZeroValue()) &&
              itk::Math::FloatAlmostEqual(gmsurf->GetPixel(velind),
                                          itk::NumericTraits<typename ImageType::PixelType>::ZeroValue()) &&
              !itk::Math::FloatAlmostEqual(segmentationimage->GetPixel(velind),
                                           static_cast<typename ImageType::PixelType>(2.0)))
          {
            shouldbezero = true;
          }
        }
        if (shouldbezero)
        {
          velofield->SetPixel(velind, zero);
          corrfield->SetPixel(velind, zero);
          invfield->SetPixel(velind, zero);
        }
        incrinvfield->SetPixel(velind, velofield->GetPixel(velind));
        ++Iterator;
      }

      if (ttiter == 0)
      {
        corrfield->FillBuffer(zero);
      }

      InvertField<ImageType, DisplacementFieldType>(invfield, corrfield, 1.0, 0.1, 20, true);
      InvertField<ImageType, DisplacementFieldType>(corrfield, invfield, 1.0, 0.1, 20, true);
      //      InvertField<ImageType,DisplacementFieldType>( invfield, corrfield, 1.0,0.1,20,true);
      //  InvertField<ImageType,DisplacementFieldType>( corrfield, invfield, 1.0,0.1,20,true);
      ttiter++;
    }

    Iterator.GoToBegin();
    RealType maxth = 0;
    while (!Iterator.IsAtEnd())
    {
      typename DisplacementFieldType::IndexType velind = Iterator.GetIndex();
      // increment velocity field at every voxel v = v + u, step 4
      velofield->SetPixel(Iterator.GetIndex(),
                          velofield->GetPixel(Iterator.GetIndex()) + incrfield->GetPixel(Iterator.GetIndex()));
      RealType hitval = hitimage->GetPixel(velind);
      RealType thkval = 0;
      if (hitval > 0.001) /** potential source of problem 2 -- this value could be smaller ... */
      {
        thkval = totalimage->GetPixel(velind) / hitval - thickoffset;
      }
      if (thkval > 10)
      {
        std::cout << "thkval " << thkval << " hitval " << hitval << " total " << totalimage->GetPixel(velind)
                  << std::endl;
      }
      if (thkval < 0)
      {
        thkval = 0;
      }
      if (itk::Math::FloatAlmostEqual(segmentationimage->GetPixel(velind),
                                      static_cast<typename ImageType::PixelType>(2.0)))
      {
        finalthickimage->SetPixel(velind, thkval);
      }
      else
      {
        finalthickimage->SetPixel(velind, 0);
      }
      if (thkval > maxth)
      {
        maxth = thkval;
      }
      if (finalthickimage->GetPixel(velind) > thickprior)
      {
        finalthickimage->SetPixel(velind, thickprior);
      }
      ++Iterator;
    }

    if (debug)
    {
      std::cout << " now smooth " << std::endl;
    }
    m_MFR->SmoothDisplacementFieldGauss(velofield, smoothingsigma);
    ANTs::WriteImage<DisplacementFieldType>(corrfield, "corrfield.nii.gz");
    ANTs::WriteImage<DisplacementFieldType>(invfield, "invfield.nii.gz");

    //    std::string velofieldname = outname + "velofield";
    // WriteDisplacementField<DisplacementFieldType>(velofield,velofieldname.c_str());
    // std::string incrfieldname = outname + "incrfield";
    // WriteDisplacementField<DisplacementFieldType>(incrfield,incrfieldname.c_str());

    // std::string tname = outname + "dork1.nii.gz";
    // ANTs::WriteImage<ImageType>(hitimage,tname.c_str());
    // tname = outname + "dork2.nii.gz";
    // ANTs::WriteImage<ImageType>(totalimage,tname.c_str());
    if (thickerrct == 0)
    {
      thickerrct = 1;
    }
    std::cout << " error " << totalerr << " at it " << its << " th-err " << thicknesserror / (RealType)thickerrct
              << " max thick " << maxth << std::endl;
    //    std::string sulcthickname =outname + "sulcthick.nii";
    //    if (ImageDimension==2) WriteJpg<ImageType>(finalthickimage,"thick.jpg");
    //    std::string velofieldname = outname + "velofield";
    // WriteDisplacementField<DisplacementFieldType>(velofield,velofieldname.c_str());
    if (debug)
    {
      std::cout << "outside it " << its << std::endl;
    }
    // std::cin.get();
  }

  finalthickimage->SetDirection(omat);
  ANTs::WriteImage<ImageType>(finalthickimage, outname.c_str());
  finalthickimage->SetDirection(fmat);

  return 0;

  thickimage->FillBuffer(0);
  typename ImageType::Pointer thkdef = m_MFR->WarpImageBackward(finalthickimage, invfield);
  Iterator.GoToBegin();
  while (!Iterator.IsAtEnd())
  {
    RealType tt1 = finalthickimage->GetPixel(Iterator.GetIndex());
    RealType tt = thkdef->GetPixel(Iterator.GetIndex());
    if (tt1 > tt)
    {
      tt = tt1;
    }
    thickimage->SetPixel(Iterator.GetIndex(), tt);
    ++Iterator;
  }

  std::string thickname = outname;
  thickimage->SetDirection(omat);
  ANTs::WriteImage<ImageType>(thickimage, thickname.c_str());

  return 0;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
KellySlater(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "KellySlater");

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

  if (argc < 6)
  {
    std::cout << "Usage:   " << argv[0]
              << " ImageDimension Segmentation.nii.gz WMProb.nii.gz GMProb.nii.gz   Out.nii {GradStep-1-2D,2-3D}   "
                 "{#Its-~50}  {ThickPriorValue-6} {Bool-use-curvature-prior} {smoothing} {BoolUseEuclidean?}"
              << std::endl;
    std::cout << " this is a kind of binary image registration thing with diffeomorphisms " << std::endl;
    std::cout << " Segmentation.nii.gz -- should contain the value 3 where WM exists and the value 2 where GM exists "
              << std::endl;
    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  unsigned int dim = std::stoi(argv[1]);
  std::cout << " dim " << dim << std::endl;

  switch (dim)
  {
    case 2:
    {
      LaplacianThicknessExpDiff2<2>(argc, argv);
    }
    break;
    case 3:
    {
      LaplacianThicknessExpDiff2<3>(argc, argv);
    }
    break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
} // namespace ants
