/*=========================================================================

 Program:   Advanced Normalization Tools

 Copyright (c) ConsortiumOfANTS. All rights reserved.
 See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/
#ifndef __itkANTSImageRegistrationOptimizer_h
#define __itkANTSImageRegistrationOptimizer_h

#include "antsAllocImage.h"
#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkVectorGaussianInterpolateImageFunction.h"
#include "antsCommandLineParser.h"
#include "itkShiftScaleImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkImage.h"
#include "itkMacro.h"
#include "ReadWriteData.h"
#include "itkCenteredEuler3DTransform.h"
#include "itkQuaternionRigidTransform.h"
#include "itkANTSAffine3DTransform.h"
#include "itkANTSCenteredAffine2DTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkFiniteDifferenceFunction.h"
#include "itkFixedArray.h"
#include "itkANTSSimilarityMetric.h"
#include "itkPDEDeformableRegistrationFilter.h"
#include "itkResampleImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkWarpImageMultiTransformFilter.h"
#include "itkDisplacementFieldFromMultiTransformFilter.h"
#include "itkWarpImageWAffineFilter.h"
#include "itkPointSet.h"
#include "itkVector.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkGeneralToBSplineDisplacementFieldFilter.h"
#include "itkBSplineControlPointImageFunction.h"
#include "ANTS_affine_registration2.h"
#include "itkVectorFieldGradientImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"

namespace itk
{
template <unsigned int TDimension = 3, typename TReal = float>
class ANTSImageRegistrationOptimizer final : public Object
{
public:
  /** Standard class typedefs. */
  typedef ANTSImageRegistrationOptimizer Self;
  typedef Object                         Superclass;
  typedef SmartPointer<Self>             Pointer;
  typedef SmartPointer<const Self>       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(ANTSImageRegistrationOptimizer);
  static constexpr unsigned int Dimension = TDimension;
  static constexpr unsigned int ImageDimension = TDimension;

  typedef double TComp;
  typedef TReal  RealType;

  typedef Image<RealType, Self::Dimension> ImageType;
  typedef typename ImageType::Pointer                        ImagePointer;

  typedef itk::MatrixOffsetTransformBase<TComp, ImageDimension, ImageDimension> TransformType;

  /** Point Types  for landmarks and labeled point-sets */
  typedef itk::ANTSLabeledPointSet<Dimension>        LabeledPointSetType;
  typedef typename LabeledPointSetType::Pointer      LabeledPointSetPointer;
  typedef typename LabeledPointSetType::PointSetType PointSetType;
  typedef typename PointSetType::Pointer             PointSetPointer;
  typedef typename PointSetType::PointType           PointType;
  typedef typename PointSetType::PixelType           PointDataType;
  typedef typename ImageType::PointType              ImagePointType;

  typedef TransformType                             AffineTransformType;
  typedef typename AffineTransformType::Pointer     AffineTransformPointer;
  typedef OptAffine<AffineTransformType, ImageType> OptAffineType;

  typedef itk::Vector<TReal, ImageDimension>      VectorType;
  typedef itk::Image<VectorType, ImageDimension>  DisplacementFieldType;
  typedef typename DisplacementFieldType::Pointer DisplacementFieldPointer;

  typedef itk::Image<VectorType, ImageDimension + 1>     TimeVaryingVelocityFieldType;
  typedef typename TimeVaryingVelocityFieldType::Pointer TimeVaryingVelocityFieldPointer;
  typedef itk::VectorLinearInterpolateImageFunction<TimeVaryingVelocityFieldType, TReal> VelocityFieldInterpolatorType;
  typedef itk::VectorGaussianInterpolateImageFunction<TimeVaryingVelocityFieldType, TReal>
                                                    VelocityFieldInterpolatorType2;
  typedef typename DisplacementFieldType::IndexType IndexType;

  typedef ants::CommandLineParser         ParserType;
  typedef typename ParserType::OptionType OptionType;

  typedef GeneralToBSplineDisplacementFieldFilter<DisplacementFieldType> BSplineFilterType;
  typedef FixedArray<RealType, Self::ImageDimension>   ArrayType;

  /** Typedefs for similarity metrics */
  typedef ANTSSimilarityMetric<Self::Dimension, TReal> SimilarityMetricType;
  typedef typename SimilarityMetricType::Pointer                         SimilarityMetricPointer;
  typedef std::vector<SimilarityMetricPointer>                           SimilarityMetricListType;

  /** FiniteDifferenceFunction type. */
  typedef FiniteDifferenceFunction<DisplacementFieldType>     FiniteDifferenceFunctionType;
  typedef typename FiniteDifferenceFunctionType::TimeStepType TimeStepType;
  typedef typename FiniteDifferenceFunctionType::Pointer      FiniteDifferenceFunctionPointer;
  typedef AvantsPDEDeformableRegistrationFunction<ImageType, ImageType, DisplacementFieldType> MetricBaseType;
  typedef typename MetricBaseType::Pointer                                                     MetricBaseTypePointer;

  /* Jacobian and other calculations */
  typedef itk::VectorFieldGradientImageFunction<DisplacementFieldType> JacobianFunctionType;

  /** Set functions */
  void
  SetAffineTransform(AffineTransformPointer A)
  {
    this->m_AffineTransform = A;
  }

  void
  SetDisplacementField(DisplacementFieldPointer A)
  {
    this->m_DisplacementField = A;
  }

  void
  SetTimeVaryingVelocityField(TimeVaryingVelocityFieldPointer A)
  {
    this->m_TimeVaryingVelocity = A;
    this->m_VelocityFieldInterpolator->SetInputImage(this->m_TimeVaryingVelocity);
  }

  void
  SetInverseDisplacementField(DisplacementFieldPointer A)
  {
    this->m_InverseDisplacementField = A;
  }

  void
  SetMaskImage(ImagePointer m)
  {
    this->m_MaskImage = m;
  }

  void
  SetReferenceSpaceImage(ImagePointer m)
  {
    this->m_ReferenceSpaceImage = m;
  }

  void
  SetFixedImageAffineTransform(AffineTransformPointer A)
  {
    this->m_FixedImageAffineTransform = A;
  }

  AffineTransformPointer
  GetFixedImageAffineTransform()
  {
    return this->m_FixedImageAffineTransform;
  }

  /** Get functions */
  AffineTransformPointer
  GetAffineTransform()
  {
    return this->m_AffineTransform;
  }

  DisplacementFieldPointer
  GetDisplacementField()
  {
    return this->m_DisplacementField;
  }

  DisplacementFieldPointer
  GetInverseDisplacementField()
  {
    return this->m_InverseDisplacementField;
  }

  /** Initialize all parameters */

  void
  SetNumberOfLevels(unsigned int i)
  {
    this->m_NumberOfLevels = i;
  }

  void
  SetParser(typename ParserType::Pointer P)
  {
    this->m_Parser = P;
  }

  /** Basic operations */
  DisplacementFieldPointer
  CopyDisplacementField(DisplacementFieldPointer input);

  std::string
  localANTSGetFilePrefix(const char * str)
  {
    std::string            filename = str;
    std::string::size_type pos = filename.rfind(".");
    std::string            filepre = std::string(filename, 0, pos);

#if 0 // HACK: This looks like an error in logic.  It has no effect
    //Perhaps this shoudl be put into a common function, or use ITK code to accomplish the desired task.
    //as is done in BRAINS.
    if( pos != std::string::npos )
      {
      std::string extension = std::string( filename, pos, filename.length() - 1);
      if( extension == std::string(".gz") )
        {
        pos = filepre.rfind( "." );
        extension = std::string( filepre, pos, filepre.length() - 1 );
        }
      }
#endif
    return filepre;
  }

  void
  SmoothDisplacementField(DisplacementFieldPointer field, bool TrueEqualsGradElseTotal)
  {
    typename ParserType::OptionType::Pointer regularizationOption = this->m_Parser->GetOption("regularization");

    if ((regularizationOption->GetFunction(0)->GetName()).find("DMFFD") != std::string::npos)
    {
      if ((!TrueEqualsGradElseTotal &&
           itk::Math::FloatAlmostEqual(this->m_TotalSmoothingparam, itk::NumericTraits<TReal>::ZeroValue())) ||
          (TrueEqualsGradElseTotal &&
           itk::Math::FloatAlmostEqual(this->m_GradSmoothingparam, itk::NumericTraits<TReal>::ZeroValue())))
      {
        return;
      }
      ArrayType    meshSize;
      unsigned int splineOrder = this->m_BSplineFieldOrder;
      TReal        bsplineKernelVariance = static_cast<TReal>(splineOrder + 1) / static_cast<TReal>(12.0);
      unsigned int numberOfLevels = 1;

      if (TrueEqualsGradElseTotal)
      {
        if (this->m_GradSmoothingparam < itk::NumericTraits<TReal>::ZeroValue())
        {
          meshSize = this->m_GradSmoothingMeshSize;
          for (unsigned int d = 0; d < ImageDimension; d++)
          {
            meshSize[d] *= static_cast<unsigned int>(std::pow(2.0, static_cast<int>(this->m_CurrentLevel)));
          }
        }
        else
        {
          TReal spanLength = std::sqrt(this->m_GradSmoothingparam / bsplineKernelVariance);
          for (unsigned int d = 0; d < ImageDimension; d++)
          {
            meshSize[d] = static_cast<unsigned int>(static_cast<TReal>(field->GetLargestPossibleRegion().GetSize()[d]) /
                                                      spanLength +
                                                    static_cast<TReal>(0.5));
          }
        }
        this->SmoothDisplacementFieldBSpline(field, meshSize, splineOrder, numberOfLevels);
      }
      else
      {
        if (this->m_TotalSmoothingparam < itk::NumericTraits<TReal>::ZeroValue())
        {
          meshSize = this->m_TotalSmoothingMeshSize;
          for (unsigned int d = 0; d < ImageDimension; d++)
          {
            meshSize[d] *= static_cast<unsigned int>(std::pow(2.0, static_cast<int>(this->m_CurrentLevel)));
          }
        }
        else
        {
          TReal spanLength = std::sqrt(this->m_TotalSmoothingparam / bsplineKernelVariance);
          for (unsigned int d = 0; d < ImageDimension; d++)
          {
            meshSize[d] = static_cast<unsigned int>(field->GetLargestPossibleRegion().GetSize()[d] / spanLength +
                                                    static_cast<TReal>(0.5));
          }
        }

        RealType maxMagnitude = 0.0;

        ImageRegionIterator<DisplacementFieldType> It(field, field->GetLargestPossibleRegion());
        for (It.GoToBegin(); !It.IsAtEnd(); ++It)
        {
          RealType magnitude = (It.Get()).GetNorm();
          if (magnitude > maxMagnitude)
          {
            maxMagnitude = magnitude;
          }
        }
        this->SmoothDisplacementFieldBSpline(field, meshSize, splineOrder, numberOfLevels);

        if (maxMagnitude > itk::NumericTraits<RealType>::ZeroValue())
        {
          for (It.GoToBegin(); !It.IsAtEnd(); ++It)
          {
            It.Set(It.Get() / maxMagnitude);
          }
        }
      }
    }
    else // Gaussian
    {
      TReal sig = 0;
      if (TrueEqualsGradElseTotal)
      {
        sig = this->m_GradSmoothingparam;
      }
      else
      {
        sig = this->m_TotalSmoothingparam;
      }
      this->SmoothDisplacementFieldGauss(field, sig);
    }
  }

  void
  SmoothDisplacementFieldGauss(DisplacementFieldPointer field = nullptr,
                               TReal                    sig = 0.0,
                               bool                     useparamimage = false,
                               unsigned int             lodim = ImageDimension);

  //  TReal = smoothingparam, int = maxdim to smooth
  void
  SmoothVelocityGauss(TimeVaryingVelocityFieldPointer field, TReal, unsigned int);

  void
  SmoothDisplacementFieldBSpline(DisplacementFieldPointer field,
                                 ArrayType                meshSize,
                                 unsigned int             splineorder,
                                 unsigned int             numberoflevels);

  DisplacementFieldPointer
  ComputeUpdateFieldAlternatingMin(DisplacementFieldPointer fixedwarp,
                                   DisplacementFieldPointer movingwarp,
                                   PointSetPointer          fpoints = nullptr,
                                   PointSetPointer          wpoints = nullptr,
                                   DisplacementFieldPointer updateFieldInv = nullptr,
                                   bool                     updateenergy = true);

  DisplacementFieldPointer
  ComputeUpdateField(DisplacementFieldPointer fixedwarp,
                     DisplacementFieldPointer movingwarp,
                     PointSetPointer          fpoints = nullptr,
                     PointSetPointer          wpoints = nullptr,
                     DisplacementFieldPointer updateFieldInv = nullptr,
                     bool                     updateenergy = true);

  TimeVaryingVelocityFieldPointer
  ExpandVelocity()
  {
    using VelocityFieldSpacingType = typename TimeVaryingVelocityFieldType::SpacingType;
    using VelocityFieldSizeType = typename TimeVaryingVelocityFieldType::SizeType;

    VelocityFieldSpacingType outputSpacing;
    VelocityFieldSizeType    outputSize;

    VelocityFieldSpacingType inputSpacing = this->m_TimeVaryingVelocity->GetSpacing();
    VelocityFieldSizeType    inputSize = this->m_TimeVaryingVelocity->GetLargestPossibleRegion().GetSize();

    for (unsigned int d = 0; d < ImageDimension; d++)
    {
      outputSize[d] = static_cast<typename VelocityFieldSizeType::SizeValueType>(this->m_CurrentDomainSize[d]);
      outputSpacing[d] = inputSpacing[d] * static_cast<double>(inputSize[d]) / static_cast<double>(outputSize[d]);
    }

    using ResamplerType = ResampleImageFilter<TimeVaryingVelocityFieldType, TimeVaryingVelocityFieldType>;
    typename ResamplerType::Pointer resampler = ResamplerType::New();
    resampler->SetInput(this->m_TimeVaryingVelocity);
    resampler->SetOutputOrigin(this->m_TimeVaryingVelocity->GetOrigin());
    resampler->SetOutputDirection(this->m_TimeVaryingVelocity->GetDirection());
    resampler->SetOutputSpacing(outputSpacing);
    resampler->SetSize(outputSize);
    resampler->Update();

    typename TimeVaryingVelocityFieldType::Pointer expandedField = resampler->GetOutput();
    expandedField->DisconnectPipeline();

    return (expandedField);
  }

  DisplacementFieldPointer
  ExpandField(DisplacementFieldPointer field, typename ImageType::SpacingType targetSpacing)
  {
    using DisplacementFieldSizeType = typename DisplacementFieldType::SizeType;

    DisplacementFieldSizeType outputSize;
    for (unsigned int d = 0; d < ImageDimension; d++)
    {
      outputSize[d] = static_cast<typename DisplacementFieldSizeType::SizeValueType>(this->m_CurrentDomainSize[d]);
    }

    using ResamplerType = ResampleImageFilter<DisplacementFieldType, DisplacementFieldType>;
    typename ResamplerType::Pointer resampler = ResamplerType::New();
    resampler->SetInput(field);
    resampler->SetOutputOrigin(field->GetOrigin());
    resampler->SetOutputDirection(field->GetDirection());
    resampler->SetOutputSpacing(targetSpacing);
    resampler->SetSize(outputSize);
    resampler->Update();

    typename DisplacementFieldType::Pointer expandedField = resampler->GetOutput();
    expandedField->DisconnectPipeline();

    return (expandedField);
  }

  ImagePointer
  GetVectorComponent(DisplacementFieldPointer field, unsigned int index)
  {
    // Initialize the Moving to the displacement field
    typedef DisplacementFieldType FieldType;

    typename ImageType::Pointer sfield = AllocImage<ImageType>(field);

    typedef itk::ImageRegionIteratorWithIndex<FieldType> Iterator;
    Iterator                                             vfIter(field, field->GetLargestPossibleRegion());
    for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
    {
      VectorType v1 = vfIter.Get();
      sfield->SetPixel(vfIter.GetIndex(), v1[index]);
    }

    return sfield;
  }

  ImagePointer
  SubsampleImage(ImagePointer,
                 RealType,
                 typename ImageType::PointType     outputOrigin,
                 typename ImageType::DirectionType outputDirection,
                 AffineTransformPointer            aff = nullptr);

  DisplacementFieldPointer
  SubsampleField(DisplacementFieldPointer        field,
                 typename ImageType::SizeType    targetSize,
                 typename ImageType::SpacingType targetSpacing)
  {
    std::cout << "FIXME -- NOT DONE CORRECTLY " << std::endl;
    std::cout << "FIXME -- NOT DONE CORRECTLY " << std::endl;
    std::cout << "FIXME -- NOT DONE CORRECTLY " << std::endl;
    std::cout << "FIXME -- NOT DONE CORRECTLY " << std::endl;
    std::cout << " SUBSAM FIELD SUBSAM FIELD SUBSAM FIELD " << std::endl;

    typename DisplacementFieldType::Pointer sfield;
    for (unsigned int i = 0; i < ImageDimension; i++)
    {
      typename ImageType::Pointer precomp = this->GetVectorComponent(field, i);
      typename ImageType::Pointer comp = this->SubsampleImage(precomp, targetSize, targetSpacing);
      if (i == 0)
      {
        sfield = AllocImage<DisplacementFieldType>(comp);
      }

      typedef itk::ImageRegionIteratorWithIndex<DisplacementFieldType> Iterator;
      typedef typename DisplacementFieldType::PixelType                DispVectorType;
      DispVectorType                                                   v1;
      DispVectorType                                                   zero;
      zero.Fill(0.0);
      Iterator vfIter(sfield, sfield->GetLargestPossibleRegion());
      for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
      {
        v1 = vfIter.Get();
        v1[i] = comp->GetPixel(vfIter.GetIndex());
        vfIter.Set(v1);
      }
    }

    return sfield;
  }

  PointSetPointer
  WarpMultiTransform(ImagePointer             referenceimage,
                     ImagePointer             movingImage,
                     PointSetPointer          movingpoints,
                     AffineTransformPointer   aff,
                     DisplacementFieldPointer totalField,
                     bool                     doinverse,
                     AffineTransformPointer   fixedaff)
  {
    if (!movingpoints)
    {
      std::cout << " NULL POINTS " << std::endl;
      return nullptr;
    }

    AffineTransformPointer affinverse = nullptr;
    if (aff)
    {
      affinverse = AffineTransformType::New();
      aff->GetInverse(affinverse);
    }
    AffineTransformPointer fixedaffinverse = nullptr;
    if (fixedaff)
    {
      fixedaffinverse = AffineTransformType::New();
      fixedaff->GetInverse(fixedaffinverse);
    }

    typedef itk::WarpImageMultiTransformFilter<ImageType, ImageType, DisplacementFieldType, TransformType> WarperType;
    typename WarperType::Pointer warper = WarperType::New();
    warper->SetInput(movingImage);
    warper->SetEdgePaddingValue(0);
    warper->SetSmoothScale(1);
    if (!doinverse)
    {
      if (totalField)
      {
        warper->PushBackDisplacementFieldTransform(totalField);
      }
      if (fixedaff)
      {
        warper->PushBackAffineTransform(fixedaff);
      }
      else if (aff)
      {
        warper->PushBackAffineTransform(aff);
      }
    }
    else
    {
      if (aff)
      {
        warper->PushBackAffineTransform(affinverse);
      }
      else if (fixedaff)
      {
        warper->PushBackAffineTransform(fixedaffinverse);
      }
      if (totalField)
      {
        warper->PushBackDisplacementFieldTransform(totalField);
      }
    }

    warper->SetOutputOrigin(referenceimage->GetOrigin());
    typename ImageType::SizeType size = referenceimage->GetLargestPossibleRegion().GetSize();
    if (totalField)
    {
      size = totalField->GetLargestPossibleRegion().GetSize();
    }
    warper->SetOutputSize(size);
    typename ImageType::SpacingType spacing = referenceimage->GetSpacing();
    if (totalField)
    {
      spacing = totalField->GetSpacing();
    }
    warper->SetOutputSpacing(spacing);
    warper->SetOutputDirection(referenceimage->GetDirection());
    totalField->SetOrigin(referenceimage->GetOrigin());
    totalField->SetDirection(referenceimage->GetDirection());

    // warper->Update();
    //      std::cout << " updated in point warp " << std::endl;
    PointSetPointer     outputMesh = PointSetType::New();
    unsigned long       count = 0;
    const unsigned long sz1 = movingpoints->GetNumberOfPoints();
    if (this->m_Debug)
    {
      std::cout << " BEFORE #points " << sz1 << std::endl;
    }
    for (unsigned long ii = 0; ii < sz1; ii++)
    {
      PointType  point;
      const bool isValid = movingpoints->GetPoint(ii, &point);
      if (!isValid)
      {
        itkExceptionMacro("Invalid moving point at : " << ii);
      }
      else
      {
        PointDataType label = 0;
        movingpoints->GetPointData(ii, &label);
        // convert pointtype to imagepointtype
        ImagePointType pt;
        for (unsigned int jj = 0; jj < ImageDimension; jj++)
        {
          pt[jj] = point[jj];
        }
        ImagePointType wpt;
        const bool     bisinside = warper->MultiTransformSinglePoint(pt, wpt);
        if (bisinside)
        {
          PointType wpoint;
          for (unsigned int jj = 0; jj < ImageDimension; jj++)
          {
            wpoint[jj] = wpt[jj];
          }
          outputMesh->SetPointData(count, label);
          outputMesh->SetPoint(count, wpoint);
          //      if (ii % 100 == 0) std::cout << " pt " << pt << " wpt " << wpt << std::endl;
          count++;
        }
      }
    }
    if (this->m_Debug)
    {
      std::cout << " AFTER #points " << count << std::endl;
    }
    //      if (count != sz1 ) std::cout << " WARNING:  POINTS ARE MAPPING OUT OF IMAGE DOMAIN " << 1.0 - (TReal)
    // count/(TReal)(sz1+1) << std::endl;
    return outputMesh;
  }

  ImagePointer
  WarpMultiTransform(ImagePointer             referenceimage,
                     ImagePointer             movingImage,
                     AffineTransformPointer   aff,
                     DisplacementFieldPointer totalField,
                     bool                     doinverse,
                     AffineTransformPointer   fixedaff)
  {
    typedef typename ImageType::DirectionType DirectionType;
    // NOT USED: DirectionType rdirection = referenceimage->GetDirection();
    // NOT USED: DirectionType mdirection = movingImage->GetDirection();

    AffineTransformPointer affinverse = nullptr;
    if (aff)
    {
      affinverse = AffineTransformType::New();
      aff->GetInverse(affinverse);
    }
    AffineTransformPointer fixedaffinverse = nullptr;
    if (fixedaff)
    {
      fixedaffinverse = AffineTransformType::New();
      fixedaff->GetInverse(fixedaffinverse);
    }

    DirectionType iddir;
    iddir.Fill(0);
    for (unsigned int i = 0; i < ImageDimension; i++)
    {
      iddir[i][i] = 1;
    }

    typedef itk::LinearInterpolateImageFunction<ImageType, TComp>          InterpolatorType1;
    typedef itk::NearestNeighborInterpolateImageFunction<ImageType, TComp> InterpolatorType2;
    typedef itk::BSplineInterpolateImageFunction<ImageType, TComp>         InterpolatorType3;
    typename InterpolatorType1::Pointer                                    interp1 = InterpolatorType1::New();
    typename InterpolatorType2::Pointer                                    interpnn = InterpolatorType2::New();
    typename InterpolatorType3::Pointer                                    interpcu = InterpolatorType3::New();

    typedef itk::WarpImageMultiTransformFilter<ImageType, ImageType, DisplacementFieldType, TransformType> WarperType;
    typename WarperType::Pointer warper = WarperType::New();
    warper->SetInput(movingImage);
    warper->SetEdgePaddingValue(0);
    warper->SetSmoothScale(1);
    warper->SetInterpolator(interp1);
    if (this->m_UseNN)
    {
      warper->SetInterpolator(interpnn);
    }
    if (!doinverse)
    {
      if (totalField)
      {
        warper->PushBackDisplacementFieldTransform(totalField);
      }
      if (fixedaff)
      {
        warper->PushBackAffineTransform(fixedaff);
      }
      else if (aff)
      {
        warper->PushBackAffineTransform(aff);
      }
    }
    else
    {
      if (aff)
      {
        warper->PushBackAffineTransform(affinverse);
      }
      else if (fixedaff)
      {
        warper->PushBackAffineTransform(fixedaffinverse);
      }
      if (totalField)
      {
        warper->PushBackDisplacementFieldTransform(totalField);
      }
    }

    warper->SetOutputOrigin(referenceimage->GetOrigin());
    typename ImageType::SizeType size = referenceimage->GetLargestPossibleRegion().GetSize();
    if (totalField)
    {
      size = totalField->GetLargestPossibleRegion().GetSize();
    }
    warper->SetOutputSize(size);
    typename ImageType::SpacingType spacing = referenceimage->GetSpacing();
    if (totalField)
    {
      spacing = totalField->GetSpacing();
    }
    warper->SetOutputSpacing(spacing);
    warper->SetOutputDirection(referenceimage->GetDirection());
    totalField->SetOrigin(referenceimage->GetOrigin());
    totalField->SetDirection(referenceimage->GetDirection());

    warper->Update();
    if (this->m_Debug)
    {
      std::cout << " updated ok -- warped image output size "
                << warper->GetOutput()->GetLargestPossibleRegion().GetSize() << " requested size "
                << totalField->GetLargestPossibleRegion().GetSize() << std::endl;
    }

    typename ImageType::Pointer outimg = warper->GetOutput();

    return outimg;
  }

  ImagePointer
  SmoothImageToScale(ImagePointer image, TReal scalingFactor)
  {
    typename ImageType::SpacingType          inputSpacing = image->GetSpacing();
    typename ImageType::RegionType::SizeType inputSize = image->GetRequestedRegion().GetSize();

    typename ImageType::SpacingType          outputSpacing;
    typename ImageType::RegionType::SizeType outputSize;

    RealType minimumSpacing = inputSpacing.GetVnlVector().min_value();
    //    RealType maximumSpacing = inputSpacing.GetVnlVector().max_value();
    for (unsigned int d = 0; d < Dimension; d++)
    {
      RealType scaling = std::min(scalingFactor * minimumSpacing / static_cast<RealType>(inputSpacing[d]),
                                  static_cast<RealType>(inputSize[d]) / static_cast<RealType>(32.0));
      outputSpacing[d] = inputSpacing[d] * static_cast<double>(scaling);
      outputSize[d] =
        static_cast<unsigned long>(static_cast<RealType>(inputSpacing[d]) * static_cast<RealType>(inputSize[d]) /
                                     static_cast<RealType>(outputSpacing[d]) +
                                   static_cast<RealType>(0.5));

      typedef RecursiveGaussianImageFilter<ImageType, ImageType> GaussianFilterType;
      typename GaussianFilterType::Pointer                       smoother = GaussianFilterType::New();
      smoother->SetInputImage(image);
      smoother->SetDirection(d);
      smoother->SetNormalizeAcrossScale(false);
      TReal sig = (outputSpacing[d] / inputSpacing[d] - 1.0) * 0.2; // /(TReal)ImageDimension;
      smoother->SetSigma(sig);

      if (smoother->GetSigma() > 0.0)
      {
        smoother->Update();
        image = smoother->GetOutput();
      }
    }

    image = this->NormalizeImage(image);

    return image;
  }

  ImagePointer
  GaussianSmoothImage(ImagePointer image, TReal sigma)
  {
    typedef DiscreteGaussianImageFilter<ImageType, ImageType> SmootherType;
    typename SmootherType::Pointer                            smoother = SmootherType::New();
    smoother->SetVariance(itk::Math::sqr(sigma));
    smoother->SetMaximumError(0.01);
    smoother->SetInput(image);

    ImagePointer smoothImage = smoother->GetOutput();
    smoothImage->Update();
    smoothImage->DisconnectPipeline();

    smoothImage = this->NormalizeImage(smoothImage);

    return smoothImage;
  }

  typename ANTSImageRegistrationOptimizer<TDimension, TReal>::DisplacementFieldPointer
  IntegrateConstantVelocity(DisplacementFieldPointer totalField, unsigned int ntimesteps, TReal timeweight);

  /** Base optimization functions */
  // AffineTransformPointer AffineOptimization(AffineTransformPointer &aff_init, OptAffine &affine_opt); // {return
  // nullptr;}
  AffineTransformPointer
  AffineOptimization(OptAffineType & affine_opt); // {return nullptr;}

  std::string
  GetTransformationModel()
  {
    return this->m_TransformationModel;
  }

  void
  SetTransformationModel(std::string s)
  {
    this->m_TransformationModel = s;
    std::cout << " Requested Transformation Model:  " << this->m_TransformationModel << " : Using " << std::endl;
    if (this->m_TransformationModel == std::string("Elast"))
    {
      std::cout << "Elastic model for transformation. " << std::endl;
    }
    else if (this->m_TransformationModel == std::string("SyN"))
    {
      std::cout << "SyN diffeomorphic model for transformation. " << std::endl;
    }
    else if (this->m_TransformationModel == std::string("GreedyExp"))
    {
      std::cout
        << "Greedy Exp Diff model for transformation.   Similar to Diffeomorphic Demons.  Params same as Exp model. "
        << std::endl;
      this->m_TransformationModel = std::string("GreedyExp");
    }
    else
    {
      std::cout << "Exp Diff model for transformation. " << std::endl;
      this->m_TransformationModel = std::string("Exp");
    }
  }

  void
  SetUpParameters()
  {
    /** Univariate Deformable Mapping */

    // set up parameters for deformation restriction

    typename OptionType::Pointer restrictDeformationOption = this->m_Parser->GetOption("Restrict-Deformation");
    if (restrictDeformationOption && restrictDeformationOption->GetNumberOfFunctions())
    {
      std::string temp = restrictDeformationOption->GetFunction()->GetName();

      this->m_RestrictDeformation = this->m_Parser->template ConvertVector<TReal>(temp);
      if (this->m_RestrictDeformation.size() != ImageDimension)
      {
        std::cout << " You input a vector of size :  " << this->m_RestrictDeformation.size()
                  << " for --Restrict-Deformation.  The vector length does not match the image dimension.  Ignoring.  "
                  << std::endl;
        for (unsigned int jj = 0; jj < this->m_RestrictDeformation.size(); jj++)
        {
          this->m_RestrictDeformation[jj] = 0;
        }
      }
    }

    // set up max iterations per level

    typename OptionType::Pointer numberOfIterationsOption = this->m_Parser->GetOption("number-of-iterations");
    if (numberOfIterationsOption && numberOfIterationsOption->GetNumberOfFunctions())
    {
      std::string temp = this->m_Parser->GetOption("number-of-iterations")->GetFunction()->GetName();
      this->m_Iterations = this->m_Parser->template ConvertVector<unsigned int>(temp);
      this->SetNumberOfLevels(this->m_Iterations.size());
    }
    this->m_UseROI = false;
    typename OptionType::Pointer roiOption = this->m_Parser->GetOption("roi");
    if (roiOption && roiOption->GetNumberOfFunctions())
    {
      std::string temp = roiOption->GetFunction()->GetName();
      this->m_RoiNumbers = this->m_Parser->template ConvertVector<TReal>(temp);
      if (temp.length() > 3)
      {
        this->m_UseROI = true;
      }
    }

    typename ParserType::OptionType::Pointer oOption = this->m_Parser->GetOption("output-naming");
    if (oOption->GetNumberOfFunctions())
    {
      this->m_OutputNamingConvention = oOption->GetFunction(0)->GetName();
    }

    typename ParserType::OptionType::Pointer thicknessOption = this->m_Parser->GetOption("geodesic");
    if (thicknessOption->GetNumberOfFunctions() &&
        (thicknessOption->GetFunction(0)->GetName() == "true" || thicknessOption->GetFunction(0)->GetName() == "1"))
    {
      // asymm forces
      this->m_ComputeThickness = 1;
      this->m_SyNFullTime = 2;
    }
    else if (thicknessOption->GetNumberOfFunctions() && thicknessOption->GetFunction(0)->GetName() == "2")
    {
      // symmetric forces
      this->m_ComputeThickness = 1;
      this->m_SyNFullTime = 1;
    }
    else
    {
      this->m_ComputeThickness = 0; // not full time varying stuff
    }
    /**
     * Get transformation model and associated parameters
     */
    typename ParserType::OptionType::Pointer transformOption = this->m_Parser->GetOption("transformation-model");
    this->SetTransformationModel(transformOption->GetFunction(0)->GetName());
    if (transformOption->GetFunction(0)->GetNumberOfParameters() >= 1)
    {
      std::string parameter = transformOption->GetFunction(0)->GetParameter(0);
      TReal       _temp = this->m_Parser->template Convert<TReal>(parameter);
      this->m_Gradstep = _temp;
      this->m_GradstepAltered = _temp;
    }
    else
    {
      this->m_Gradstep = 0.5;
      this->m_GradstepAltered = 0.5;
    }
    if (transformOption->GetFunction(0)->GetNumberOfParameters() >= 2)
    {
      std::string parameter = transformOption->GetFunction(0)->GetParameter(1);
      this->m_NTimeSteps = this->m_Parser->template Convert<unsigned int>(parameter);
    }
    else
    {
      this->m_NTimeSteps = 1;
    }
    if (transformOption->GetFunction(0)->GetNumberOfParameters() >= 3)
    {
      std::string parameter = transformOption->GetFunction(0)->GetParameter(2);
      this->m_DeltaTime = this->m_Parser->template Convert<TReal>(parameter);
      if (this->m_DeltaTime > 1)
      {
        this->m_DeltaTime = 1;
      }
      if (this->m_DeltaTime <= 0)
      {
        this->m_DeltaTime = 0.001;
      }
      std::cout << " set DT " << this->m_DeltaTime << std::endl;
      this->m_SyNType = 1;
    }
    else
    {
      this->m_DeltaTime = 0.1;
    }
    //    if ( transformOption->GetFunction( 0 )->GetNumberOfParameters() >= 3 )
    //      {
    //      std::string parameter = transformOption->GetFunction( 0 )->GetParameter( 2 );
    //      this->m_SymmetryType
    //        = this->m_Parser->template Convert<unsigned int>( parameter );
    //      }

    /**
     * Get regularization and associated parameters
     */
    this->m_GradSmoothingparam = -1;
    this->m_TotalSmoothingparam = -1;
    this->m_GradSmoothingMeshSize.Fill(0);
    this->m_TotalSmoothingMeshSize.Fill(0);

    typename ParserType::OptionType::Pointer regularizationOption = this->m_Parser->GetOption("regularization");
    if (regularizationOption->GetFunction(0)->GetName() == "Gauss")
    {
      if (regularizationOption->GetFunction(0)->GetNumberOfParameters() >= 1)
      {
        std::string parameter = regularizationOption->GetFunction(0)->GetParameter(0);
        this->m_GradSmoothingparam = this->m_Parser->template Convert<TReal>(parameter);
      }
      else
      {
        this->m_GradSmoothingparam = 3;
      }
      if (regularizationOption->GetFunction(0)->GetNumberOfParameters() >= 2)
      {
        std::string parameter = regularizationOption->GetFunction(0)->GetParameter(1);
        this->m_TotalSmoothingparam = this->m_Parser->template Convert<TReal>(parameter);
      }
      else
      {
        this->m_TotalSmoothingparam = 0.5;
      }
      if (regularizationOption->GetFunction(0)->GetNumberOfParameters() >= 3)
      {
        std::string parameter = regularizationOption->GetFunction(0)->GetParameter(2);
        this->m_GaussianTruncation = this->m_Parser->template Convert<TReal>(parameter);
      }
      else
      {
        this->m_GaussianTruncation = 256;
      }
      std::cout << "  Grad Step " << this->m_Gradstep << " total-smoothing " << this->m_TotalSmoothingparam
                << " gradient-smoothing " << this->m_GradSmoothingparam << std::endl;
    }
    else if ((regularizationOption->GetFunction(0)->GetName()).find("DMFFD") != std::string::npos)
    {
      if (regularizationOption->GetFunction(0)->GetNumberOfParameters() >= 1)
      {
        std::string parameter = regularizationOption->GetFunction(0)->GetParameter(0);
        if (parameter.find("x") != std::string::npos)
        {
          std::vector<unsigned int> gradMeshSize = this->m_Parser->template ConvertVector<unsigned int>(parameter);
          for (unsigned int d = 0; d < ImageDimension; d++)
          {
            this->m_GradSmoothingMeshSize[d] = gradMeshSize[d];
          }
        }
        else
        {
          this->m_GradSmoothingparam = this->m_Parser->template Convert<TReal>(parameter);
        }
      }
      else
      {
        this->m_GradSmoothingparam = 3.0;
      }
      if (regularizationOption->GetFunction(0)->GetNumberOfParameters() >= 2)
      {
        std::string parameter = regularizationOption->GetFunction(0)->GetParameter(1);
        if (parameter.find("x") != std::string::npos)
        {
          std::vector<unsigned int> totalMeshSize = this->m_Parser->template ConvertVector<unsigned int>(parameter);
          for (unsigned int d = 0; d < ImageDimension; d++)
          {
            this->m_TotalSmoothingMeshSize[d] = totalMeshSize[d];
          }
        }
        else
        {
          this->m_TotalSmoothingparam = this->m_Parser->template Convert<TReal>(parameter);
        }
      }
      else
      {
        this->m_TotalSmoothingparam = 0.5;
      }
      if (regularizationOption->GetFunction(0)->GetNumberOfParameters() >= 3)
      {
        std::string parameter = regularizationOption->GetFunction(0)->GetParameter(2);
        this->m_BSplineFieldOrder = this->m_Parser->template Convert<unsigned int>(parameter);
      }
      else
      {
        this->m_BSplineFieldOrder = 3;
      }
      std::cout << "  Grad Step " << this->m_Gradstep << " total-smoothing " << this->m_TotalSmoothingparam
                << " gradient-smoothing " << this->m_GradSmoothingparam << " bspline-field-order "
                << this->m_BSplineFieldOrder << std::endl;
    }
    else
    {
      this->m_GradSmoothingparam = 3;
      this->m_TotalSmoothingparam = 0.5;
      std::cout << " Default Regularization is Gaussian smoothing with : " << this->m_GradSmoothingparam << " & "
                << this->m_TotalSmoothingparam << std::endl;
      //      itkExceptionMacro( "Invalid regularization: " << regularizationOption->GetFunction( 0 )->GetName() );
    }
  }

  void
  ComputeMultiResolutionParameters(ImagePointer fixedImage)
  {
    VectorType zero;

    zero.Fill(0);
    /** Compute scale factors */
    this->m_FullDomainSpacing = fixedImage->GetSpacing();
    this->m_FullDomainSize = fixedImage->GetRequestedRegion().GetSize();
    this->m_CurrentDomainSpacing = fixedImage->GetSpacing();
    this->m_CurrentDomainSize = fixedImage->GetRequestedRegion().GetSize();
    this->m_CurrentDomainDirection = fixedImage->GetDirection();
    this->m_FullDomainOrigin.Fill(0);
    this->m_CurrentDomainOrigin.Fill(0);
    /** alter the input size based on information gained from the ROI information - if available */
    if (this->m_UseROI)
    {
      for (unsigned int ii = 0; ii < ImageDimension; ii++)
      {
        this->m_FullDomainSize[ii] =
          (typename ImageType::SizeType::SizeValueType)this->m_RoiNumbers[ii + ImageDimension];
        this->m_FullDomainOrigin[ii] = this->m_RoiNumbers[ii];
      }
      std::cout << " ROI #s : size " << this->m_FullDomainSize << " orig " << this->m_FullDomainOrigin << std::endl;
    }

    RealType minimumSpacing = this->m_FullDomainSpacing.GetVnlVector().min_value();
    //            RealType maximumSpacing = this->m_FullDomainSpacing.GetVnlVector().max_value();
    for (unsigned int d = 0; d < Dimension; d++)
    {
      RealType scaling = this->m_ScaleFactor;

      // If the user doesn't specify subsampling factors, we go to the
      // default behavior

      if (this->m_SubsamplingFactors.size() == 0)
      {
        scaling = std::min(this->m_ScaleFactor * minimumSpacing / static_cast<RealType>(this->m_FullDomainSpacing[d]),
                           static_cast<RealType>(this->m_FullDomainSize[d]) / static_cast<RealType>(32.0));
      }
      if (scaling < itk::NumericTraits<RealType>::OneValue())
      {
        scaling = itk::NumericTraits<RealType>::OneValue();
      }
      this->m_CurrentDomainSpacing[d] = this->m_FullDomainSpacing[d] * static_cast<double>(scaling);
      this->m_CurrentDomainSize[d] =
        static_cast<unsigned long>(this->m_FullDomainSpacing[d] * static_cast<double>(this->m_FullDomainSize[d]) /
                                     this->m_CurrentDomainSpacing[d] +
                                   0.5);
      this->m_CurrentDomainOrigin[d] =
        static_cast<unsigned long>(this->m_FullDomainSpacing[d] * static_cast<double>(this->m_FullDomainOrigin[d]) /
                                     this->m_CurrentDomainSpacing[d] +
                                   0.5);
    }

    //            this->m_Debug=true;
    if (this->m_Debug)
    {
      std::cout << " outsize " << this->m_CurrentDomainSize << " curspc " << this->m_CurrentDomainSpacing << " fullspc "
                << this->m_FullDomainSpacing << " fullsz " << this->m_FullDomainSize << std::endl;
    }
    //            this->m_Debug=false;

    if (!this->m_DisplacementField)
    { /*FIXME -- need initial deformation strategy */
      typename ImageType::RegionType region;
      region.SetSize(this->m_CurrentDomainSize);
      this->m_DisplacementField = AllocImage<DisplacementFieldType>(
        region, this->m_CurrentDomainSpacing, fixedImage->GetOrigin(), fixedImage->GetDirection(), zero);
      std::cout << " allocated def field " << this->m_DisplacementField->GetDirection() << std::endl;
      // std::exception();
    }
    else
    {
      this->m_DisplacementField = this->ExpandField(this->m_DisplacementField, this->m_CurrentDomainSpacing);
      if (this->m_TimeVaryingVelocity)
      {
        this->ExpandVelocity();
      }
    }
  }

  ImagePointer
  NormalizeImage(ImagePointer image)
  {
    typedef itk::MinimumMaximumImageFilter<ImageType> MinMaxFilterType;
    typename MinMaxFilterType::Pointer                minMaxFilter = MinMaxFilterType::New();

    minMaxFilter->SetInput(image);
    minMaxFilter->Update();

    TReal min = minMaxFilter->GetMinimum();
    TReal shift = -min;
    TReal scale = static_cast<TReal>(minMaxFilter->GetMaximum());
    scale += shift;
    scale = itk::NumericTraits<TReal>::OneValue() / scale;

    typedef itk::ShiftScaleImageFilter<ImageType, ImageType> FilterType;
    typename FilterType::Pointer                             filter = FilterType::New();

    filter->SetInput(image);
    filter->SetShift(shift);
    filter->SetScale(scale);
    filter->Update();

    return filter->GetOutput();
  }

  void
  DeformableOptimization()
  {
    DisplacementFieldPointer updateField = nullptr;

    this->SetUpParameters();
    VectorType zero;
    zero.Fill(0);
    std::cout << " setting N-TimeSteps = " << this->m_NTimeSteps << " trunc " << this->m_GaussianTruncation
              << std::endl;

    // Get subsample factors and gaussian smoothing sigmas if specified
    // by the user.
    typename OptionType::Pointer subsamplingOption = this->m_Parser->GetOption("subsampling-factors");
    if (subsamplingOption && subsamplingOption->GetNumberOfFunctions())
    {
      std::string subsamplingfactors = subsamplingOption->GetFunction()->GetName();
      if (subsamplingfactors.size() > 0)
      {
        std::vector<float> factors = this->m_Parser->template ConvertVector<float>(subsamplingfactors);
        if (factors.size() == this->m_NumberOfLevels)
        {
          this->m_SubsamplingFactors.SetSize(this->m_NumberOfLevels);
          for (unsigned int d = 0; d < this->m_NumberOfLevels; d++)
          {
            this->m_SubsamplingFactors[d] = factors[d];
          }
        }
        //         else
        //           {
        //           itkWarningMacro( "The number of levels does not match the size of factors."
        //                            << "  Using default settings." );
        //           }
      }
    }

    typename OptionType::Pointer gaussianSmoothingSigmasOption = this->m_Parser->GetOption("gaussian-smoothing-sigmas");

    if (gaussianSmoothingSigmasOption && gaussianSmoothingSigmasOption->GetNumberOfFunctions())
    {
      std::string        gaussiansmoothingsigmas = gaussianSmoothingSigmasOption->GetFunction()->GetName();
      std::vector<float> sigmas = this->m_Parser->template ConvertVector<float>(gaussiansmoothingsigmas);
      if (sigmas.size() == this->m_NumberOfLevels)
      {
        this->m_GaussianSmoothingSigmas.SetSize(this->m_NumberOfLevels);
        for (unsigned int d = 0; d < this->m_NumberOfLevels; d++)
        {
          this->m_GaussianSmoothingSigmas[d] = sigmas[d];
        }
      }
      //       else
      //         {
      //         itkWarningMacro( "The number of levels does not match the size of sigmas."
      //                          << "  Using default settings." );
      //         }
    }

    unsigned int maxits = 0;
    for (unsigned int currentLevel = 0; currentLevel < this->m_NumberOfLevels; currentLevel++)
    {
      if (this->m_Iterations[currentLevel] > maxits)
      {
        maxits = this->m_Iterations[currentLevel];
      }
    }
    if (maxits == 0)
    {
      this->m_DisplacementField = nullptr;
      this->m_InverseDisplacementField = nullptr;
      //    this->ComputeMultiResolutionParameters(this->m_SimilarityMetrics[0]->GetFixedImage());
      return;
    }

    /* this is a hack to force univariate mappings   in the future,
       we will re-cast this framework  s.t. multivariate images can be used */
    unsigned int numberOfMetrics = this->m_SimilarityMetrics.size();
    for (unsigned int metricCount = 1; metricCount < numberOfMetrics; metricCount++)
    {
      this->m_SimilarityMetrics[metricCount]->GetFixedImage()->SetOrigin(
        this->m_SimilarityMetrics[0]->GetFixedImage()->GetOrigin());
      this->m_SimilarityMetrics[metricCount]->GetFixedImage()->SetDirection(
        this->m_SimilarityMetrics[0]->GetFixedImage()->GetDirection());
      this->m_SimilarityMetrics[metricCount]->GetMovingImage()->SetOrigin(
        this->m_SimilarityMetrics[0]->GetMovingImage()->GetOrigin());
      this->m_SimilarityMetrics[metricCount]->GetMovingImage()->SetDirection(
        this->m_SimilarityMetrics[0]->GetMovingImage()->GetDirection());
    }
    /* here, we assign all point set pointers to any single
       non-null point-set pointer */
    for (unsigned int metricCount = 0; metricCount < numberOfMetrics; metricCount++)
    {
      for (unsigned int metricCount2 = 0; metricCount2 < numberOfMetrics; metricCount2++)
      {
        if (this->m_SimilarityMetrics[metricCount]->GetFixedPointSet())
        {
          this->m_SimilarityMetrics[metricCount2]->SetFixedPointSet(
            this->m_SimilarityMetrics[metricCount]->GetFixedPointSet());
        }
        if (this->m_SimilarityMetrics[metricCount]->GetMovingPointSet())
        {
          this->m_SimilarityMetrics[metricCount2]->SetMovingPointSet(
            this->m_SimilarityMetrics[metricCount]->GetMovingPointSet());
        }
      }
    }
    this->m_SmoothFixedImages.resize(numberOfMetrics, nullptr);
    this->m_SmoothMovingImages.resize(numberOfMetrics, nullptr);
    for (unsigned int currentLevel = 0; currentLevel < this->m_NumberOfLevels; currentLevel++)
    {
      this->m_CurrentLevel = currentLevel;
      typedef Vector<TReal, 1>                      ProfilePointDataType;
      typedef Image<ProfilePointDataType, 1>        CurveType;
      typedef PointSet<ProfilePointDataType, 1>     EnergyProfileType;
      typedef typename EnergyProfileType::PointType ProfilePointType;
      typedef typename EnergyProfileType::Pointer   EnergyProfilePointer;
      std::vector<EnergyProfilePointer>             energyProfiles;
      energyProfiles.resize(numberOfMetrics);
      for (unsigned int qq = 0; qq < numberOfMetrics; qq++)
      {
        energyProfiles[qq] = EnergyProfileType::New();
        energyProfiles[qq]->Initialize();
      }

      ImagePointer fixedImage;
      ImagePointer movingImage;
      this->m_GradstepAltered = this->m_Gradstep;

      if (this->m_SubsamplingFactors.Size() == 0)
      {
        this->m_ScaleFactor = std::pow(2.0, (int)static_cast<RealType>(this->m_NumberOfLevels - currentLevel - 1));
      }
      else
      {
        this->m_ScaleFactor = this->m_SubsamplingFactors[currentLevel];
      }
      std::cout << " ScaleFactor " << this->m_ScaleFactor;
      if (this->m_GaussianSmoothingSigmas.size() > 0)
      {
        std::cout << " smoothing sigma " << this->m_GaussianSmoothingSigmas[currentLevel];
      }
      std::cout << " nlev " << this->m_NumberOfLevels << " curl " << currentLevel << std::endl;

      /** FIXME -- here we assume the metrics all have the same image */
      fixedImage = this->m_SimilarityMetrics[0]->GetFixedImage();
      movingImage = this->m_SimilarityMetrics[0]->GetMovingImage();
      if (!this->m_ReferenceSpaceImage)
      {
        this->m_ReferenceSpaceImage = fixedImage;
      }
      this->ComputeMultiResolutionParameters(this->m_ReferenceSpaceImage);
      std::cout << " Its at this level " << this->m_Iterations[currentLevel] << std::endl;
      /*  generate smoothed images for all metrics */
      for (unsigned int metricCount = 0; metricCount < numberOfMetrics; metricCount++)
      {
        if (this->m_GaussianSmoothingSigmas.size() == 0)
        {
          this->m_SmoothFixedImages[metricCount] =
            this->SmoothImageToScale(this->m_SimilarityMetrics[metricCount]->GetFixedImage(), this->m_ScaleFactor);
          this->m_SmoothMovingImages[metricCount] =
            this->SmoothImageToScale(this->m_SimilarityMetrics[metricCount]->GetMovingImage(), this->m_ScaleFactor);
        }
        else
        {
          this->m_SmoothFixedImages[metricCount] = this->GaussianSmoothImage(
            this->m_SimilarityMetrics[metricCount]->GetFixedImage(), this->m_GaussianSmoothingSigmas[currentLevel]);
          this->m_SmoothMovingImages[metricCount] = this->GaussianSmoothImage(
            this->m_SimilarityMetrics[metricCount]->GetMovingImage(), this->m_GaussianSmoothingSigmas[currentLevel]);
        }
      }
      fixedImage = this->m_SmoothFixedImages[0];
      movingImage = this->m_SmoothMovingImages[0];

      unsigned int nmet = this->m_SimilarityMetrics.size();
      this->m_LastEnergy.resize(nmet, 1.e12);
      this->m_Energy.resize(nmet, 1.e9);
      this->m_EnergyBad.resize(nmet, 0);
      bool converged = false;
      this->m_CurrentIteration = 0;

      if (this->GetTransformationModel() != std::string("SyN"))
      {
        this->m_FixedImageAffineTransform = nullptr;
      }
      while (!converged)
      {
        for (unsigned int metricCount = 0; metricCount < numberOfMetrics; metricCount++)
        {
          this->m_SimilarityMetrics[metricCount]->GetModifiableMetric()->SetIterations(this->m_CurrentIteration);
        }

        if (this->GetTransformationModel() == std::string("Elast"))
        {
          if (this->m_Iterations[currentLevel] > 0)
          {
            this->ElasticRegistrationUpdate(fixedImage, movingImage);
          }
        }
        else if (this->GetTransformationModel() == std::string("SyN"))
        {
          if (currentLevel > 0)
          {
            this->m_SyNF = this->ExpandField(this->m_SyNF, this->m_CurrentDomainSpacing);
            this->m_SyNFInv = this->ExpandField(this->m_SyNFInv, this->m_CurrentDomainSpacing);
            this->m_SyNM = this->ExpandField(this->m_SyNM, this->m_CurrentDomainSpacing);
            this->m_SyNMInv = this->ExpandField(this->m_SyNMInv, this->m_CurrentDomainSpacing);
          }
          if (this->m_Iterations[currentLevel] > 0)
          {
            if (this->m_SyNType && this->m_ComputeThickness)
            {
              this->DiReCTUpdate(fixedImage,
                                 movingImage,
                                 this->m_SimilarityMetrics[0]->GetFixedPointSet(),
                                 this->m_SimilarityMetrics[0]->GetMovingPointSet());
            }
            else if (this->m_SyNType)
            {
              this->SyNTVRegistrationUpdate(fixedImage,
                                            movingImage,
                                            this->m_SimilarityMetrics[0]->GetFixedPointSet(),
                                            this->m_SimilarityMetrics[0]->GetMovingPointSet());
            }
            else
            {
              this->SyNRegistrationUpdate(fixedImage,
                                          movingImage,
                                          this->m_SimilarityMetrics[0]->GetFixedPointSet(),
                                          this->m_SimilarityMetrics[0]->GetMovingPointSet());
            }
          }
          else if (this->m_SyNType)
          {
            this->UpdateTimeVaryingVelocityFieldWithSyNFandSyNM();
          }
          //            this->CopyOrAddToVelocityField( this->m_SyNF,  0 , false);
        }
        else if (this->GetTransformationModel() == std::string("Exp"))
        {
          if (this->m_Iterations[currentLevel] > 0)
          {
            this->DiffeomorphicExpRegistrationUpdate(fixedImage,
                                                     movingImage,
                                                     this->m_SimilarityMetrics[0]->GetFixedPointSet(),
                                                     this->m_SimilarityMetrics[0]->GetMovingPointSet());
          }
        }
        else if (this->GetTransformationModel() == std::string("GreedyExp"))
        {
          if (this->m_Iterations[currentLevel] > 0)
          {
            this->GreedyExpRegistrationUpdate(fixedImage,
                                              movingImage,
                                              this->m_SimilarityMetrics[0]->GetFixedPointSet(),
                                              this->m_SimilarityMetrics[0]->GetMovingPointSet());
          }
        }

        this->m_CurrentIteration++;
        /**
         * This is where we track the energy profile to check for convergence.
         */
        for (unsigned int qq = 0; qq < numberOfMetrics; qq++)
        {
          ProfilePointType point;
          point[0] = this->m_CurrentIteration - 1;

          ProfilePointDataType energy;
          energy[0] = this->m_Energy[qq];

          energyProfiles[qq]->SetPoint(this->m_CurrentIteration - 1, point);
          energyProfiles[qq]->SetPointData(this->m_CurrentIteration - 1, energy);
        }

        /**
         * If there are a sufficent number of iterations, fit a quadratic
         * single B-spline span to the number of energy profile points
         * in the first metric. To test convergence, evaluate the derivative
         * at the end of the profile to determine if >= 0.  To change to a
         * window of the energy profile, simply change the origin (assuming that
         * the desired window will start at the user-specified origin and
         * end at the current iteration).
         */

        // Added by Paul: allow weighted metric to be used here
        typename ParserType::OptionType::Pointer regularizationOption =
          this->m_Parser->GetOption("use-all-metrics-for-convergence");
        bool use_all_metrics = (atoi(regularizationOption->GetFunction()->GetName().c_str()) > 0);

        unsigned int domtar = 12;
        if (this->m_CurrentIteration > domtar)
        {
          typedef BSplineScatteredDataPointSetToImageFilter<EnergyProfileType, CurveType> BSplinerType;
          typename BSplinerType::Pointer bspliner = BSplinerType::New();

          typename CurveType::PointType origin;
          unsigned int                  domainorigin = 0;
          unsigned int                  domainsize = this->m_CurrentIteration - domainorigin;
          if (this->m_CurrentIteration > domtar)
          {
            domainsize = domtar;
            domainorigin = this->m_CurrentIteration - domainsize;
          }
          origin.Fill(domainorigin);
          typename CurveType::SizeType size;
          size.Fill(domainsize);
          typename CurveType::SpacingType _spacing;
          _spacing.Fill(1);

          typename EnergyProfileType::Pointer energyProfileWindow = EnergyProfileType::New();
          energyProfileWindow->Initialize();

          unsigned int windowBegin = static_cast<unsigned int>(origin[0]);
          TReal        totale = 0;
          for (unsigned int qq = windowBegin; qq < this->m_CurrentIteration; qq++)
          {
            ProfilePointType point;
            point[0] = qq;
            ProfilePointDataType energy;
            energy.Fill(0);
            if (use_all_metrics)
            {
              for (unsigned int im = 0; im < numberOfMetrics; im++)
              {
                ProfilePointType pm;
                pm[0] = qq;
                ProfilePointDataType em;
                em.Fill(0);
                energyProfiles[im]->GetPointData(qq, &em);
                RealType weight = this->m_SimilarityMetrics[im]->GetWeightScalar();
                energy[0] += weight * em[0];
              }
            }
            else
            {
              energyProfiles[0]->GetPointData(qq, &energy);
            }
            totale += energy[0];
            energyProfileWindow->SetPoint(qq - windowBegin, point);
            energyProfileWindow->SetPointData(qq - windowBegin, energy);
          }
          //      std::cout <<" totale " << totale << std::endl;
          if (totale > 0)
          {
            totale *= -itk::NumericTraits<TReal>::OneValue();
          }
          for (unsigned int qq = windowBegin; qq < this->m_CurrentIteration; qq++)
          {
            ProfilePointDataType energy;
            energy.Fill(0);
            if (use_all_metrics)
            {
              for (unsigned int im = 0; im < numberOfMetrics; im++)
              {
                ProfilePointType pm;
                pm[0] = qq;
                ProfilePointDataType em;
                em.Fill(0);
                energyProfiles[im]->GetPointData(qq, &em);
                RealType weight = this->m_SimilarityMetrics[im]->GetWeightScalar();
                energy[0] += weight * em[0];
              }
            }
            else
            {
              energyProfiles[0]->GetPointData(qq, &energy);
            }
            energyProfileWindow->SetPointData(qq - windowBegin, energy / totale);
          }

          bspliner->SetInput(energyProfileWindow);
          bspliner->SetOrigin(origin);
          bspliner->SetSpacing(_spacing);
          bspliner->SetSize(size);
          bspliner->SetNumberOfLevels(1);
          unsigned int order = 1;
          bspliner->SetSplineOrder(order);
          typename BSplinerType::ArrayType ncps;
          ncps.Fill(order + 1); // single span, order = 2
          bspliner->SetNumberOfControlPoints(ncps);
          bspliner->Update();

          typedef BSplineControlPointImageFunction<CurveType> BSplinerType2;
          typename BSplinerType2::Pointer                     bspliner2 = BSplinerType2::New();
          bspliner2->SetOrigin(origin);
          bspliner2->SetSpacing(_spacing);
          bspliner2->SetSize(size);
          bspliner2->SetSplineOrder(1);
          bspliner2->SetInputImage(bspliner->GetPhiLattice());

          ProfilePointType endPoint;
          endPoint[0] = static_cast<TReal>(this->m_CurrentIteration - domainsize * 0.5);
          typename BSplinerType2::GradientType gradient = bspliner2->EvaluateGradientAtParametricPoint(endPoint);
          this->m_ESlope = gradient[0][0];
          if (this->m_ESlope < static_cast<TReal>(0.0001) && this->m_CurrentIteration > domtar)
          {
            converged = true;
          }
          std::cout << " E-Slope " << this->m_ESlope; // << std::endl;
        }
        for (unsigned int qq = 0; qq < this->m_Energy.size(); qq++)
        {
          if (qq == 0)
          {
            std::cout << " iteration " << this->m_CurrentIteration;
          }

          std::cout << " energy " << qq << " : " << this->m_Energy[qq]; //  << " Last " <<
                                                                        // this->m_LastEnergy[qq];
          if (this->m_LastEnergy[qq] < this->m_Energy[qq])
          {
            this->m_EnergyBad[qq]++;
          }
        }
        for (unsigned int qq = 0; qq < this->m_Energy.size(); qq++)
        {
          if (this->m_CurrentIteration <= 1)
          {
            this->m_EnergyBad[qq] = 0;
          }
        }

        // if ( this->m_EnergyBad[0] > 2)
        //         {
        //          this->m_GradstepAltered*=0.8;
        //        std::cout <<" reducing gradstep " <<  this->m_GradstepAltered;
        // this->m_EnergyBad[this->m_Energy.size()-1]=0;
        // }
        std::cout << std::endl;

        if (this->m_CurrentIteration >= this->m_Iterations[currentLevel])
        {
          converged = true;
        }
        //        || this->m_EnergyBad[0] >= 6 )
        //
        if (converged && this->m_CurrentIteration >= this->m_Iterations[currentLevel])
        {
          std::cout << " tired convergence: reached max iterations " << std::endl;
        }
        else if (converged)
        {
          std::cout << " Converged due to oscillation in optimization ";
          for (unsigned int qq = 0; qq < this->m_Energy.size(); qq++)
          {
            std::cout << " metric " << qq << " bad " << this->m_EnergyBad[qq] << "  ";
          }
          std::cout << std::endl;
        }
      }
    }

    if (this->GetTransformationModel() == std::string("SyN"))
    {
      //         TReal timestep=1.0/(TReal)this->m_NTimeSteps;
      //         unsigned int nts=this->m_NTimeSteps;
      if (this->m_SyNType)
      {
        //         this->m_SyNFInv = this->IntegrateConstantVelocity(this->m_SyNF, nts, timestep*(-1.));
        //         this->m_SyNMInv = this->IntegrateConstantVelocity(this->m_SyNM, nts, timestep*(-1.));
        //         this->m_SyNF= this->IntegrateConstantVelocity(this->m_SyNF, nts, timestep);
        //         this->m_SyNM= this->IntegrateConstantVelocity(this->m_SyNM,
        //         nts, timestep);
        //     DisplacementFieldPointer fdiffmap = this->IntegrateVelocity(0,0.5);
        //     this->m_SyNFInv = this->IntegrateVelocity(0.5,0);
        //     DisplacementFieldPointer mdiffmap = this->IntegrateVelocity(0.5,1);
        //     this->m_SyNMInv = this->IntegrateVelocity(1,0.5);
        //     this->m_SyNM=this->CopyDisplacementField(mdiffmap);
        //     this->m_SyNF=this->CopyDisplacementField(fdiffmap);
        this->m_DisplacementField = this->IntegrateVelocity(0, 1);
        //     ImagePointer wmimage= this->WarpMultiTransform(
        //  this->m_SmoothFixedImages[0],this->m_SmoothMovingImages[0], this->m_AffineTransform,
        // this->m_DisplacementField, false , this->m_ScaleFactor );
        this->m_InverseDisplacementField = this->IntegrateVelocity(1, 0);
      }
      else
      {
        this->m_InverseDisplacementField = this->CopyDisplacementField(this->m_SyNM);
        this->ComposeDiffs(this->m_SyNF, this->m_SyNMInv, this->m_DisplacementField, 1);
        this->ComposeDiffs(this->m_SyNM, this->m_SyNFInv, this->m_InverseDisplacementField, 1);
      }
    }
    else if (this->GetTransformationModel() == std::string("Exp"))
    {
      DisplacementFieldPointer diffmap =
        this->IntegrateConstantVelocity(this->m_DisplacementField, (unsigned int)this->m_NTimeSteps, 1);
      // 1.0/ (TReal)this->m_NTimeSteps);
      DisplacementFieldPointer invdiffmap =
        this->IntegrateConstantVelocity(this->m_DisplacementField, (unsigned int)this->m_NTimeSteps, -1);
      // -1.0/(TReal)this->m_NTimeSteps);
      this->m_InverseDisplacementField = invdiffmap;
      this->m_DisplacementField = diffmap;
      AffineTransformPointer invaff = nullptr;
      if (this->m_AffineTransform)
      {
        invaff = AffineTransformType::New();
        this->m_AffineTransform->GetInverse(invaff);
        if (this->m_Debug)
        {
          std::cout << " ??????invaff " << this->m_AffineTransform << std::endl << std::endl;
        }
        if (this->m_Debug)
        {
          std::cout << " invaff?????? " << invaff << std::endl << std::endl;
        }
      }
    }
    else if (this->GetTransformationModel() == std::string("GreedyExp"))
    {
      DisplacementFieldPointer diffmap = this->m_DisplacementField;
      this->m_InverseDisplacementField = nullptr;
      this->m_DisplacementField = diffmap;
      AffineTransformPointer invaff = nullptr;
      if (this->m_AffineTransform)
      {
        invaff = AffineTransformType::New();
        this->m_AffineTransform->GetInverse(invaff);
        if (this->m_Debug)
        {
          std::cout << " ??????invaff " << this->m_AffineTransform << std::endl << std::endl;
        }
        if (this->m_Debug)
        {
          std::cout << " invaff?????? " << invaff << std::endl << std::endl;
        }
      }
    }

    this->m_DisplacementField->SetOrigin(this->m_ReferenceSpaceImage->GetOrigin());
    this->m_DisplacementField->SetDirection(this->m_ReferenceSpaceImage->GetDirection());
    if (this->m_InverseDisplacementField)
    {
      this->m_InverseDisplacementField->SetOrigin(this->m_ReferenceSpaceImage->GetOrigin());
      this->m_InverseDisplacementField->SetDirection(this->m_ReferenceSpaceImage->GetDirection());
    }

    if (this->m_TimeVaryingVelocity)
    {
      std::string outname =
        localANTSGetFilePrefix(this->m_OutputNamingConvention.c_str()) + std::string("velocity.mhd");
      typename itk::ImageFileWriter<TimeVaryingVelocityFieldType>::Pointer writer =
        itk::ImageFileWriter<TimeVaryingVelocityFieldType>::New();
      writer->SetFileName(outname.c_str());
      writer->SetInput(this->m_TimeVaryingVelocity);
      writer->UpdateLargestPossibleRegion();
      //    writer->Write();
      std::cout << " write tv field " << outname << std::endl;
      //        ANTs::WriteImage<TimeVaryingVelocityFieldType>( this->m_TimeVaryingVelocity , outname.c_str());
    }
  }

  void
  DiffeomorphicExpRegistrationUpdate(ImagePointer    fixedImage,
                                     ImagePointer    movingImage,
                                     PointSetPointer fpoints = nullptr,
                                     PointSetPointer mpoints = nullptr);

  void
  GreedyExpRegistrationUpdate(ImagePointer    fixedImage,
                              ImagePointer    movingImage,
                              PointSetPointer fpoints = nullptr,
                              PointSetPointer mpoints = nullptr);

  void
  SyNRegistrationUpdate(ImagePointer    fixedImage,
                        ImagePointer    movingImage,
                        PointSetPointer fpoints = nullptr,
                        PointSetPointer mpoints = nullptr);

  void
  SyNExpRegistrationUpdate(ImagePointer    fixedImage,
                           ImagePointer    movingImage,
                           PointSetPointer fpoints = nullptr,
                           PointSetPointer mpoints = nullptr);

  void
  SyNTVRegistrationUpdate(ImagePointer    fixedImage,
                          ImagePointer    movingImage,
                          PointSetPointer fpoints = nullptr,
                          PointSetPointer mpoints = nullptr);

  void
  DiReCTUpdate(ImagePointer    fixedImage,
               ImagePointer    movingImage,
               PointSetPointer fpoints = nullptr,
               PointSetPointer mpoints = nullptr);

  /** allows one to copy or add a field to a time index within the velocity
   * field
   */
  void
  UpdateTimeVaryingVelocityFieldWithSyNFandSyNM();

  void
  CopyOrAddToVelocityField(TimeVaryingVelocityFieldPointer velocity,
                           DisplacementFieldPointer        update1,
                           DisplacementFieldPointer        update2,
                           TReal                           timept);

  // void CopyOrAddToVelocityField( DisplacementFieldPointer update,  unsigned int timeindex,  bool
  // CopyIsTrueOtherwiseAdd);

  void ElasticRegistrationUpdate(ImagePointer /* fixedImage */, ImagePointer /* xxxxmovingImage */)
  {
    VectorType zero;
    zero.Fill(0);
    DisplacementFieldPointer updateField;

    updateField = this->ComputeUpdateField(this->m_DisplacementField, nullptr, nullptr, nullptr, nullptr);

    typedef ImageRegionIteratorWithIndex<DisplacementFieldType> Iterator;
    Iterator dIter(this->m_DisplacementField, this->m_DisplacementField->GetLargestPossibleRegion());
    for (dIter.GoToBegin(); !dIter.IsAtEnd(); ++dIter)
    {
      typename ImageType::IndexType index = dIter.GetIndex();
      VectorType                    vec = updateField->GetPixel(index);
      dIter.Set(dIter.Get() + vec * this->m_Gradstep);
    }

    if (this->m_Debug)
    {
      std::cout << " updated elast "
                << " up-sz " << updateField->GetLargestPossibleRegion() << std::endl;
      std::cout << " t-sz " << this->m_DisplacementField->GetLargestPossibleRegion() << std::endl;
    }
    this->SmoothDisplacementField(this->m_DisplacementField, false);

    return;
  }

  ImagePointer
  WarpImageBackward(ImagePointer image, DisplacementFieldPointer field)
  {
    typedef WarpImageFilter<ImageType, ImageType, DisplacementFieldType> WarperType;
    typename WarperType::Pointer                                         warper = WarperType::New();
    warper->SetInput(image);
    warper->SetDisplacementField(field);
    warper->SetEdgePaddingValue(0);
    warper->SetOutputSpacing(field->GetSpacing());
    warper->SetOutputOrigin(field->GetOrigin());
    warper->Update();
    return warper->GetOutput();
  }

  void
  ComposeDiffs(DisplacementFieldPointer fieldtowarpby,
               DisplacementFieldPointer field,
               DisplacementFieldPointer fieldout,
               TReal                    sign);

  void
  SetSimilarityMetrics(SimilarityMetricListType S)
  {
    this->m_SimilarityMetrics = S;
  }

  void
  SetFixedPointSet(PointSetPointer p)
  {
    this->m_FixedPointSet = p;
  }

  void
  SetMovingPointSet(PointSetPointer p)
  {
    this->m_MovingPointSet = p;
  }

  void
  SetDeltaTime(TReal t)
  {
    this->m_DeltaTime = t;
  }

  TReal
  InvertField(DisplacementFieldPointer field,
              DisplacementFieldPointer inverseField,
              TReal                    weight = 1.0,
              TReal                    toler = 0.1,
              int                      maxiter = 20,
              bool /* print */ = false)
  {
    TReal        mytoler = toler;
    unsigned int mymaxiter = maxiter;

    typename ParserType::OptionType::Pointer thicknessOption = this->m_Parser->GetOption("go-faster");
    if (thicknessOption->GetFunction(0)->GetName() == "true" || thicknessOption->GetFunction(0)->GetName() == "1")
    {
      mytoler = 0.5;
      maxiter = 12;
    }

    VectorType zero;
    zero.Fill(0);
    //  if (this->GetElapsedIterations() < 2 ) maxiter=10;

    ImagePointer TRealImage = AllocImage<ImageType>(field);

    typedef typename DisplacementFieldType::PixelType           DispVectorType;
    typedef typename DisplacementFieldType::IndexType           DispIndexType;
    typedef ImageRegionIteratorWithIndex<DisplacementFieldType> Iterator;

    DisplacementFieldPointer lagrangianInitCond = AllocImage<DisplacementFieldType>(field);
    DisplacementFieldPointer eulerianInitCond = AllocImage<DisplacementFieldType>(field);

    typedef typename DisplacementFieldType::SizeType SizeType;
    SizeType                                         size = field->GetLargestPossibleRegion().GetSize();

    typename ImageType::SpacingType spacing = field->GetSpacing();

    unsigned long npix = 1;
    for (unsigned int j = 0; j < ImageDimension; j++) // only use in-plane spacing
    {
      npix *= field->GetLargestPossibleRegion().GetSize()[j];
    }

    TReal    max = 0;
    Iterator iter(field, field->GetLargestPossibleRegion());
    for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
    {
      DispIndexType  index = iter.GetIndex();
      DispVectorType vec1 = iter.Get();
      DispVectorType newvec = vec1 * weight;
      lagrangianInitCond->SetPixel(index, newvec);
      TReal mag = 0;
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

    TReal scale = itk::NumericTraits<TReal>::OneValue() / max;
    if (scale > itk::NumericTraits<TReal>::OneValue())
    {
      scale = 1.0;
    }
    //    TReal initscale=scale;
    Iterator vfIter(inverseField, inverseField->GetLargestPossibleRegion());

    //  int num=10;
    //  for (int its=0; its<num; its++)
    TReal        difmag = 10.0;
    unsigned int ct = 0;

    TReal meandif = 1.e8;
    //    int badct=0;
    //  while (difmag > subpix && meandif > subpix*0.1 && badct < 2 )//&& ct < 20 && denergy > 0)
    //    TReal length=0.0;
    TReal stepl = 2.;

    TReal epsilon = (TReal)size[0] / 256;
    if (epsilon > 1)
    {
      epsilon = 1;
    }

    while (difmag > mytoler && ct < mymaxiter && meandif > static_cast<TReal>(0.001))
    {
      meandif = 0.0;

      // this field says what position the eulerian field should contain in the E domain
      this->ComposeDiffs(inverseField, lagrangianInitCond, eulerianInitCond, 1);
      difmag = 0.0;
      for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
      {
        DispIndexType  index = vfIter.GetIndex();
        DispVectorType update = eulerianInitCond->GetPixel(index);
        TReal          mag = 0;
        for (unsigned int j = 0; j < ImageDimension; j++)
        {
          update[j] *= -itk::NumericTraits<TReal>::OneValue();
          mag += static_cast<TReal>(itk::Math::sqr((update[j] / static_cast<TReal>(spacing[j]))));
        }
        mag = sqrt(mag);
        meandif += mag;
        if (mag > difmag)
        {
          difmag = mag;
        }
        //      if (mag < 1.e-2) update.Fill(0);

        eulerianInitCond->SetPixel(index, update);
        TRealImage->SetPixel(index, mag);
      }
      meandif /= (TReal)npix;
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
        TReal          val = TRealImage->GetPixel(vfIter.GetIndex());
        DispVectorType update = eulerianInitCond->GetPixel(vfIter.GetIndex());
        if (val > stepl)
        {
          update = update * (stepl / val);
        }
        DispVectorType upd = vfIter.Get() + update * (epsilon);
        vfIter.Set(upd);
      }
      ct++;
    }

    // std::cout <<" difmag " << difmag << ": its " << ct <<  std::endl;

    return difmag;
  }

  void
  SetUseNearestNeighborInterpolation(bool useNN)
  {
    this->m_UseNN = useNN;
  }

  void
  SetUseBSplineInterpolation(bool useNN)
  {
    this->m_UseBSplineInterpolation = useNN;
  }

  VectorType
  IntegratePointVelocity(TReal starttimein, TReal finishtimein, IndexType startPoint);

protected:
  DisplacementFieldPointer IntegrateVelocity(TReal, TReal);
  DisplacementFieldPointer
  IntegrateLandmarkSetVelocity(TReal, TReal, PointSetPointer movingpoints, ImagePointer referenceimage);

  ImagePointer
  MakeSubImage(ImagePointer bigimage)
  {
    typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;

    typename ImageType::RegionType region;

    region.SetSize(this->m_CurrentDomainSize);
    typename ImageType::IndexType index;
    index.Fill(0);
    region.SetIndex(index);
    ImagePointer varimage =
      AllocImage<ImageType>(region, this->m_CurrentDomainSpacing, bigimage->GetOrigin(), bigimage->GetDirection(), 0);

    typename ImageType::IndexType cornerind;
    cornerind.Fill(0);
    for (unsigned int ii = 0; ii < ImageDimension; ii++)
    {
      TReal diff = (TReal)this->m_CurrentDomainOrigin[ii] - (TReal)this->m_CurrentDomainSize[ii] / 2;
      if (diff < 0)
      {
        diff = 0;
      }
      cornerind[ii] = (unsigned long)diff;
    }
    //  std::cout << " corner index " << cornerind << std::endl;
    Iterator vfIter2(bigimage, bigimage->GetLargestPossibleRegion());
    for (vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2)
    {
      typename ImageType::IndexType origindex = vfIter2.GetIndex();
      typename ImageType::IndexType _index = vfIter2.GetIndex();
      bool                          oktosample = true;
      for (unsigned int ii = 0; ii < ImageDimension; ii++)
      {
        TReal centerbasedcoord = (origindex[ii] - this->m_CurrentDomainOrigin[ii]);
        //                TReal diff =
        //                    index[ii]=origindex[ii]-cornerind[ii];
        if (fabs(centerbasedcoord) > (this->m_CurrentDomainSize[ii] / 2 - 1))
        {
          oktosample = false;
        }
      }
      if (oktosample)
      {
        //      std::cout << " index " << index <<  " origindex " << origindex << " ok? " << oktosample <<
        // std::endl;
        varimage->SetPixel(_index, bigimage->GetPixel(origindex));
      }
    }
    // std::cout << " sizes " << varimage->GetLargestPossibleRegion().GetSize() << " bigimage " <<
    // bigimage->GetLargestPossibleRegion().GetSize() << std::endl;
    return varimage;
  }

  TReal
  MeasureDeformation(DisplacementFieldPointer field, int /* option */ = 0)
  {
    typedef typename DisplacementFieldType::PixelType           DispVectorType;
    typedef typename DisplacementFieldType::IndexType           DispIndexType;
    typedef typename DisplacementFieldType::SizeType            SizeType;
    typedef ImageRegionIteratorWithIndex<DisplacementFieldType> Iterator;
    // all we have to do here is add the local field to the global field.
    Iterator      vfIter(field, field->GetLargestPossibleRegion());
    SizeType      size = field->GetLargestPossibleRegion().GetSize();
    // unsigned long ct = 1;
    TReal         totalmag = 0;
    TReal         maxstep = 0;
    //  this->m_EuclideanNorm=0;

    typename ImageType::SpacingType myspacing = field->GetSpacing();
    for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
    {
      DispIndexType  index = vfIter.GetIndex();
      DispIndexType  rindex = vfIter.GetIndex();
      DispIndexType  lindex = vfIter.GetIndex();
      DispVectorType update = vfIter.Get();
      TReal          mag = 0.0;
      TReal          stepl = 0.0;
      for (unsigned int i = 0; i < ImageDimension; i++)
      {
        rindex = index;
        lindex = index;
        if ((int)index[i] < (int)(size[i] - 2))
        {
          rindex[i] = rindex[i] + 1;
        }
        if (index[i] > 2)
        {
          lindex[i] = lindex[i] - 1;
        }
        DispVectorType rupdate = field->GetPixel(rindex);
        DispVectorType lupdate = field->GetPixel(lindex);
        DispVectorType dif = rupdate - lupdate;
        for (unsigned int tt = 0; tt < ImageDimension; tt++)
        {
          stepl += static_cast<TReal>(itk::Math::sqr(update[tt] / static_cast<TReal>(myspacing[tt])));
          mag += static_cast<TReal>(itk::Math::sqr(dif[tt] / static_cast<TReal>(myspacing[tt])));
        }
      }
      stepl = sqrt(stepl);
      mag = sqrt(mag);
      if (stepl > maxstep)
      {
        maxstep = stepl;
      }
      // ct++;
      totalmag += mag;
      //    this->m_EuclideanNorm+=stepl;
    }
    // this->m_EuclideanNorm/=ct;
    // this->m_ElasticPathLength = totalmag/ct;
    // this->m_LinftyNorm = maxstep;
    //  std::cout << " Elast path length " << this->m_ElasticPathLength <<  " L inf norm " << this->m_LinftyNorm <<
    // std::endl;
    // if (this->m_ElasticPathLength >= this->m_ArcLengthGoal)
    //  if (maxstep >= this->m_ArcLengthGoal)
    {
      //  this->StopRegistration();
      // scale the field to the right length
      //  TReal scale=this->m_ArcLengthGoal/this->m_ElasticPathLength;
      //  for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )vfIter.Set(vfIter.Get()*scale);
    }
    // if (this->m_LinftyNorm <= 0) this->m_LinftyNorm=1;
    // if (this->m_ElasticPathLength <= 0) this->m_ElasticPathLength=0;
    // if (this->m_EuclideanNorm <= 0) this->m_EuclideanNorm=0;

    // if (option==0) return this->m_ElasticPathLength;
    // else if (option==2) return this->m_EuclideanNorm;
    // else
    return maxstep;
  }

  ANTSImageRegistrationOptimizer();
  ~ANTSImageRegistrationOptimizer() override = default;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  ANTSImageRegistrationOptimizer(const Self &) = delete;
  void
  operator=(const Self &) = delete;

  typename VelocityFieldInterpolatorType::Pointer m_VelocityFieldInterpolator;

  typename ImageType::SizeType      m_CurrentDomainSize;
  typename ImageType::PointType     m_CurrentDomainOrigin;
  typename ImageType::SpacingType   m_CurrentDomainSpacing;
  typename ImageType::DirectionType m_CurrentDomainDirection;
  typename ImageType::SizeType      m_FullDomainSize;
  typename ImageType::PointType     m_FullDomainOrigin;
  typename ImageType::SpacingType   m_FullDomainSpacing;

  AffineTransformPointer   m_AffineTransform;
  AffineTransformPointer   m_FixedImageAffineTransform;
  DisplacementFieldPointer m_DisplacementField;
  DisplacementFieldPointer m_InverseDisplacementField;

  std::vector<TReal>        m_GradientDescentParameters;
  std::vector<TReal>        m_MetricScalarWeights;
  std::vector<ImagePointer> m_SmoothFixedImages;
  std::vector<ImagePointer> m_SmoothMovingImages;

  bool                         m_Debug;
  unsigned int                 m_NumberOfLevels;
  typename ParserType::Pointer m_Parser;
  SimilarityMetricListType     m_SimilarityMetrics;
  ImagePointer                 m_MaskImage;
  ImagePointer                 m_ReferenceSpaceImage;
  TReal                        m_ScaleFactor;
  bool                         m_UseMulti;
  bool                         m_UseROI;
  bool                         m_UseNN;
  bool                         m_UseBSplineInterpolation;
  unsigned int                 m_CurrentIteration;
  unsigned int                 m_CurrentLevel;
  std::string                  m_TransformationModel;
  std::string                  m_OutputNamingConvention;
  PointSetPointer              m_FixedPointSet;
  PointSetPointer              m_MovingPointSet;
  std::vector<unsigned int>    m_Iterations;
  std::vector<TReal>           m_RestrictDeformation;
  std::vector<TReal>           m_RoiNumbers;
  TReal                        m_GradSmoothingparam;
  TReal                        m_TotalSmoothingparam;
  TReal                        m_Gradstep;
  TReal                        m_GradstepAltered;
  TReal                        m_NTimeSteps;
  TReal                        m_GaussianTruncation;
  TReal                        m_DeltaTime;
  TReal                        m_ESlope;

  /** energy stuff */
  std::vector<TReal>        m_Energy;
  std::vector<TReal>        m_LastEnergy;
  std::vector<unsigned int> m_EnergyBad;

  /** for SyN only */
  DisplacementFieldPointer        m_SyNF;
  DisplacementFieldPointer        m_SyNFInv;
  DisplacementFieldPointer        m_SyNM;
  DisplacementFieldPointer        m_SyNMInv;
  TimeVaryingVelocityFieldPointer m_TimeVaryingVelocity;
  TimeVaryingVelocityFieldPointer m_LastTimeVaryingVelocity;
  TimeVaryingVelocityFieldPointer m_LastTimeVaryingUpdate;
  unsigned int                    m_SyNType;

  /** for BSpline stuff */
  unsigned int m_BSplineFieldOrder;
  ArrayType    m_GradSmoothingMeshSize;
  ArrayType    m_TotalSmoothingMeshSize;

  /** For thickness calculation */
  ImagePointer m_HitImage;
  ImagePointer m_ThickImage;
  unsigned int m_ComputeThickness;
  unsigned int m_SyNFullTime;

  /** Subsampling/Gaussian smoothing parameters */
  Array<float> m_GaussianSmoothingSigmas;
  Array<float> m_SubsamplingFactors;
};
} // namespace itk
// end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkANTSImageRegistrationOptimizer.cxx"
#endif

#endif
