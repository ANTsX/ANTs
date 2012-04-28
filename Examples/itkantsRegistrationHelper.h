/*=========================================================================
*
*  Copyright Insight Software Consortium
*
*  Licensed under the Apache License, Version 2.0 (the "License");
*  you may not use this file except in compliance with the License.
*  You may obtain a copy of the License at
*
*         http://www.apache.org/licenses/LICENSE-2.0.txt
*
*  Unless required by applicable law or agreed to in writing, software
*  distributed under the License is distributed on an "AS IS" BASIS,
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*  See the License for the specific language governing permissions and
*  limitations under the License.
*
*=========================================================================*/
#ifndef __antsRegistrationHelper_h
#define __antsRegistrationHelper_h
#include <iostream>
#include <deque>
#include "itkObject.h"
#include "itkWeakPointer.h"
#include "itkDisplacementFieldTransform.h"
#include "itkTimeVaryingVelocityFieldTransform.h"
#include "itkCompositeTransform.h"
#include "itkImageToImageMetricv4.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImage.h"

namespace itk
{
namespace ants
{
template <unsigned VImageDimension>
class RegistrationHelper : public Object
{
public:
  /** Standard class typedefs */
  typedef RegistrationHelper       Self;
  typedef Object                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  typedef WeakPointer<const Self>  ConstWeakPointer;

  typedef double                            RealType;
  typedef float                             PixelType;
  typedef Image<PixelType, VImageDimension> ImageType;

  typedef itk::Transform<double, VImageDimension, VImageDimension>     TransformType;
  typedef CompositeTransform<RealType, VImageDimension>                CompositeTransformType;
  typedef DisplacementFieldTransform<RealType, VImageDimension>        DisplacementFieldTransformType;
  typedef TimeVaryingVelocityFieldTransform<RealType, VImageDimension> TimeVaryingVelocityFieldTransformType;
  typedef ImageToImageMetricv4<ImageType, ImageType>                   MetricType;
  typedef ImageMaskSpatialObject<VImageDimension>                      ImageMaskSpatialObjectType;
  typedef typename ImageMaskSpatialObjectType::ImageType               MaskImageType;

  enum MetricEnumeration
    {
    CC = 0,
    MI = 1,
    Mattes = 2,
    MeanSquares = 3,
    Demons = 4,
    GC = 5,
    IllegalMetric = 6
    };
  enum SamplingStrategy
    {
    none = 0,
    regular = 1,
    random = 2
    };
  class Metric
  {
public:
    Metric(MetricEnumeration metricType, typename ImageType::Pointer & fixedImage,
           typename ImageType::Pointer & movingImage, double weighting, SamplingStrategy samplingStrategy,
           int numberOfBins,
           unsigned int radius,
           double samplingPercentage) :
      m_MetricType(metricType),
      m_FixedImage(fixedImage),
      m_MovingImage(movingImage),
      m_Weighting(weighting),
      m_SamplingStrategy(samplingStrategy),
      m_NumberOfBins(numberOfBins),
      m_Radius(radius),
      m_SamplingPercentage(samplingPercentage)
    {
    }

    const std::string GetMetricAsString() const
    {
      switch( this->m_MetricType )
        {
        case CC:
      { return "CC"; }
        case MI:
      { return "MI"; }
        case Mattes:
      { return "Mattes";; }
        case MeanSquares:
      { return "MeanSquares"; }
        case Demons:
      { return "Demons"; }
        case GC:
      { return "GC"; }
        default:
          {
          }
          break;
        }
      return "";
    }

    MetricEnumeration m_MetricType;
    typename ImageType::Pointer m_FixedImage;
    typename ImageType::Pointer m_MovingImage;
    double           m_Weighting;
    SamplingStrategy m_SamplingStrategy;
    int              m_NumberOfBins;
    unsigned int     m_Radius;
    double           m_SamplingPercentage;
  };

  typedef std::deque<Metric> MetricListType;

  enum XfrmMethod
    {
    Rigid = 0,
    Affine = 1,
    CompositeAffine = 2,
    Similarity = 3,
    Translation = 4,
    BSpline = 5,
    GaussianDisplacementField = 6,
    BSplineDisplacementField = 7,
    TimeVaryingVelocityField = 8,
    TimeVaryingBSplineVelocityField = 9,
    SyN = 10,
    BSplineSyN = 11,
    UnknownXfrm = 12
    };

  class TransformMethod
  {
public:
    TransformMethod() : m_XfrmMethod(Rigid),
      m_GradientStep(0),
      m_UpdateFieldVarianceInVarianceSpace(0.0),
      m_TotalFieldVarianceInVarianceSpace(0.0),
      m_SplineOrder(3),
      m_UpdateFieldTimeSigma(0.0),
      m_TotalFieldTimeSigma(0.0),
      m_NumberOfTimeIndices(0),
      m_NumberOfTimePointSamples(4)
    {
    }

    std::string XfrmMethodAsString() const
    {
      switch( this->m_XfrmMethod )
        {
        case Rigid:
      { return std::string("Rigid"); }
        case Affine:
      { return std::string("Affine"); }
        case CompositeAffine:
      { return std::string("CompositeAffine"); }
        case Similarity:
      { return std::string("Similarity"); }
        case Translation:
      { return std::string("Translation"); }
        case BSpline:
      { return std::string("BSpline"); }
        case GaussianDisplacementField:
      { return std::string("GaussianDisplacementField"); }
        case BSplineDisplacementField:
      { return std::string("BSplineDisplacementField"); }
        case TimeVaryingVelocityField:
      { return std::string("TimeVaryingVelocityField"); }
        case TimeVaryingBSplineVelocityField:
      { return std::string("TimeVaryingBSplineVelocityField"); }
        case SyN:
      { return std::string("SyN"); }
        case BSplineSyN:
      { return std::string("BSplineSyN"); }
        case UnknownXfrm: return std::string("UnknownXfrm");
        }
      return std::string("Impossible");
    }

    XfrmMethod m_XfrmMethod;
    // all transforms
    double m_GradientStep;
    // BSPline
    std::vector<unsigned int> m_MeshSizeAtBaseLevel;
    // GaussianDisplacementField
    double m_UpdateFieldVarianceInVarianceSpace;
    double m_TotalFieldVarianceInVarianceSpace;
    // BSplineDisplacementField
    std::vector<unsigned int> m_TotalFieldMeshSizeAtBaseLevel;
    std::vector<unsigned int> m_UpdateFieldMeshSizeAtBaseLevel;
    unsigned int              m_SplineOrder; // also TimeVaryingBSplineVelocityField
    // TimeVaryingVelocityField
    double       m_UpdateFieldTimeSigma;
    double       m_TotalFieldTimeSigma;
    unsigned int m_NumberOfTimeIndices;
    // TimeVaryingBSplineVelocityField
    std::vector<unsigned int> m_VelocityFieldMeshSize;
    unsigned int              m_NumberOfTimePointSamples;
  };
  typedef std::deque<TransformMethod> TransformMethodListType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(RegistrationHelper, Object);

  /** Dimension of the image.  This constant is used by functions that are
   * templated over image type (as opposed to being templated over pixel type
   * and dimension) when they need compile time access to the dimension of
   * the image. */
  itkStaticConstMacro(ImageDimension, unsigned int, VImageDimension);

  /**
   * add a metric, corresponding to the registration stage
   */
  void AddMetric(MetricEnumeration metricType,
                 typename ImageType::Pointer & fixedImage,
                 typename ImageType::Pointer & movingImage,
                 double weighting,
                 SamplingStrategy samplingStrategy,
                 int numberOfBins,
                 unsigned int radius,
                 double samplingPercentage);

  /**
   * return the enumerated Metric type based on a string representation
   */
  MetricEnumeration StringToMetricType(const std::string & str) const;

  /**
   * return the enumerated transform method specified by str
   */
  XfrmMethod StringToXfrmMethod(const std::string & str) const;

  /**
   * set the fixed initial transform.
   */
  void SetFixedInitialTransform(const TransformType *initialTransform);

  /**
   * set the moving initial transform.
   */
  void SetMovingInitialTransform(const TransformType *initialTransform);

  /**
   * add a rigid transform
   */
  void AddRigidTransform(double GradientStep);

  /**
   * add an affine transform
   */
  void AddAffineTransform(double GradientStep);

  /**
   * add a composite affine transform
   */
  void AddCompositeAffineTransform(double GradientStep);

  /**
   * add a similarity transform
   */
  void AddSimilarityTransform(double GradientStep);

  /**
   * add a translation transform
   */
  void AddTranslationTransform(double GradientStep);

  /**
   * add a spline transform
   */
  void AddBSplineTransform(double GradientStep, std::vector<unsigned int> & MeshSizeAtBaseLevel);

  /**
   * add gaussian displacement transform
   */
  void AddGaussianDisplacementFieldTransform(double GradientStep, double UpdateFieldVarianceInVarianceSpace,
                                             double TotalFieldVarianceInVarianceSpace);

  /**
   * add bspline displacement transform
   */
  void AddBSplineDisplacementFieldTransform(double GradientStep,
                                            std::vector<unsigned int> & UpdateFieldMeshSizeAtBaseLevel,
                                            std::vector<unsigned int> & TotalFieldMeshSizeAtBaseLevel,
                                            unsigned int SplineOrder);

  /**
   * add a time varying velocity field transform
   */
  void AddTimeVaryingVelocityFieldTransform(double GradientStep, unsigned int NumberOfTimeIndices,
                                            double UpdateFieldVarianceInVarianceSpace, double UpdateFieldTimeSigma,
                                            double TotalFieldVarianceInVarianceSpace, double TotalFieldTimeSigma);

  /**
   * add a time varying b spline velocity field transform
   */
  void AddTimeVaryingBSplineVelocityFieldTransform(double GradientStep, std::vector<unsigned int> VelocityFieldMeshSize,
                                                   unsigned int NumberOfTimePointSamples, unsigned int SplineOrder);

  /**
   * add a SyN transform
   */
  void AddSyNTransform(double GradientStep, double UpdateFieldVarianceInVarianceSpace,
                       double TotalFieldVarianceInVarianceSpace);

  /**
   * add a B-spline SyN transform
   */
  void AddBSplineSyNTransform(double GradientStep, std::vector<unsigned int> & UpdateFieldMeshSizeAtBaseLevel,
                              std::vector<unsigned int> & TotalFieldMeshSizeAtBaseLevel, unsigned int SplineOrder);

  /**
   * Add the collected iterations list
   */
  void SetIterations(const std::vector<std::vector<unsigned int> > & Iterations);

  /**
   * Add the collected convergence thresholds
   */
  void SetConvergenceThresholds(const std::vector<double> & thresholds);

  /**
   * Add the collected convergence window sizes
   */
  void SetConvergenceWindowSizes(const std::vector<unsigned int> & windowSizes);

  /**
   * Add the collected smoothing sigmas list
   */
  void SetSmoothingSigmas(const std::vector<std::vector<float> > & SmoothingSigmas);

  /**
   * Add the collected shrink factors list
   */
  void SetShrinkFactors(const std::vector<std::vector<unsigned int> > & ShrinkFactors);

  /**
   * turn on histogram matching of the input images
   */
  itkSetMacro(UseHistogramMatching, bool);
  itkGetMacro(UseHistogramMatching, bool);

  /**
   * turn on the option that lets you estimate the learning rate step size only at the beginning of each level.
   * useful as a second stage of fine-scale registration.
   */
  itkSetMacro(DoEstimateLearningRateOnce, bool);
  itkGetMacro(DoEstimateLearningRateOnce, bool);

  /**
   * turn on winsorize image intensity normalization
   */
  void SetWinsorizeImageIntensities(bool Winsorize, float LowerQuantile = 0.0, float UpperQuantile = 1.0);

  itkGetObjectMacro(CompositeTransform, CompositeTransformType);

  /**
   *  Get the Warped Image & Inverse Warped Images
   */
  typename ImageType::Pointer GetWarpedImage() const;

  typename ImageType::Pointer GetInverseWarpedImage() const;

  /**
   * Set the fixed/moving image masks with a spatial object
   */
  void SetFixedImageMask(typename ImageMaskSpatialObjectType::Pointer & fixedImageMask)
  {
    this->m_FixedImageMask = fixedImageMask;
  }

  void SetMovingImageMask(typename ImageMaskSpatialObjectType::Pointer & movingImageMask)
  {
    this->m_MovingImageMask = movingImageMask;
  }

  /**
   * Set the fixed/moving mask image. this will be used to instantiate
   * ImageMaskSpatialObject masks.
   */
  void SetFixedImageMask(typename MaskImageType::Pointer & fixedImageMask);
  void SetMovingImageMask(typename MaskImageType::Pointer & movingImageMask);

  /**
   * Do the registration. Will return EXIT_FAILURE if there is any
   * problem completing the registration.
   */
  int DoRegistration();

  /**
   * print out the internal registration helper state
   */
  void PrintState() const;

  void SetLogStream(std::ostream & logStream)
  {
    this->m_LogStream = &logStream;
  }

protected:
  RegistrationHelper();
  virtual ~RegistrationHelper();
private:
  int ValidateParameters();

  std::ostream & Logger() const
  {
    return *m_LogStream;
  }

  typename CompositeTransformType::Pointer m_CompositeTransform;
  typename CompositeTransformType::Pointer m_FixedInitialTransform;
  typename ImageMaskSpatialObjectType::Pointer     m_FixedImageMask;
  typename ImageMaskSpatialObjectType::Pointer     m_MovingImageMask;
  unsigned int                            m_NumberOfStages;
  MetricListType                          m_Metrics;
  TransformMethodListType                 m_TransformMethods;
  std::vector<std::vector<unsigned int> > m_Iterations;
  std::vector<double>                     m_ConvergenceThresholds;
  std::vector<unsigned int>               m_ConvergenceWindowSizes;
  std::vector<std::vector<float> >        m_SmoothingSigmas;
  std::vector<std::vector<unsigned int> > m_ShrinkFactors;
  bool                                    m_UseHistogramMatching;
  bool                                    m_WinsorizeImageIntensities;
  bool                                    m_DoEstimateLearningRateOnce;
  double                                  m_LowerQuantile;
  double                                  m_UpperQuantile;
  std::ostream *                          m_LogStream;
};
} // namespace ants
} // namespace itk

#include "itkantsRegistrationHelper.hxx"

#endif
