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
#include "itkCompositeTransform.h"
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

  typedef double                                        RealType;
  typedef float                                         PixelType;
  typedef Image<PixelType, VImageDimension>             ImageType;
  typedef CompositeTransform<RealType, VImageDimension> CompositeTransformType;

  class InitialTransform
  {
public:
    InitialTransform(const std::string & filename, bool useInverse) :
      m_Filename(filename), m_UseInverse(useInverse)
    {
    }

    std::string m_Filename;
    bool        m_UseInverse;
  };
  typedef std::deque<InitialTransform> InitialTransformListType;

  enum MetricEnumeration
    {
    CC = 0,
    MI = 1,
    Mattes = 2,
    MeanSquares = 3,
    GC = 4,
    IllegalMetric = 5
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
    UnknownXfrm = 11
    };

  class TransformMethod
  {
public:
    TransformMethod() : m_XfrmMethod(Rigid),
      m_GradientStep(0),
      m_UpdateFieldSigmaInPhysicalSpace(0.0),
      m_TotalFieldSigmaInPhysicalSpace(0.0),
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
    double m_UpdateFieldSigmaInPhysicalSpace;
    double m_TotalFieldSigmaInPhysicalSpace;
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

  itkSetMacro(WriteOutputs, bool);
  itkGetMacro(WriteOutputs, bool);

  itkSetStringMacro(OutputTransformPrefix);
  itkGetStringMacro(OutputTransformPrefix);

  itkSetStringMacro(OutputWarpedImageName);
  itkGetStringMacro(OutputWarpedImageName);

  itkSetStringMacro(OutputInverseWarpedImageName);
  itkGetStringMacro(OutputInverseWarpedImageName);

  void AddMetric(MetricEnumeration metricType,
                 typename ImageType::Pointer & fixedImage,
                 typename ImageType::Pointer & movingImage,
                 double weighting,
                 SamplingStrategy samplingStrategy,
                 int numberOfBins,
                 unsigned int radius,
                 double samplingPercentage);

  MetricEnumeration StringToMetricType(const std::string & str) const;

  XfrmMethod StringToXfrmMethod(const std::string & str) const;

  void AddInitialTransform(const std::string & filename, bool useInverse);

  void AddRigidTransform(double GradientStep);

  void AddAffineTransform(double GradientStep);

  void AddCompositeAffineTransform(double GradientStep);

  void AddSimilarityTransform(double GradientStep);

  void AddTranslationTransform(double GradientStep);

  void AddBSplineTransform(double GradientStep, std::vector<unsigned int> & MeshSizeAtBaseLevel);

  void AddGaussianDisplacementFieldTransform(double GradientStep, double UpdateFieldSigmaInPhysicalSpace,
                                             double TotalFieldSigmaInPhysicalSpace);

  void AddBSplineDisplacementFieldTransform(double GradientStep,
                                            std::vector<unsigned int> & UpdateFieldMeshSizeAtBaseLevel,
                                            std::vector<unsigned int> & TotalFieldMeshSizeAtBaseLevel,
                                            unsigned int SplineOrder);

  void AddTimeVaryingVelocityFieldTransform(double GradientStep, unsigned int NumberOfTimeIndices,
                                            double UpdateFieldSigmaInPhysicalSpace, double UpdateFieldTimeSigma,
                                            double TotalFieldSigmaInPhysicalSpace, double TotalFieldTimeSigma);

  void AddTimeVaryingBSplineVelocityFieldTransform(double GradientStep, std::vector<unsigned int> VelocityFieldMeshSize,
                                                   unsigned int NumberOfTimePointSamples, unsigned int SplineOrder);

  void AddSyNTransform(double GradientStep, double UpdateFieldSigmaInPhysicalSpace,
                       double TotalFieldSigmaInPhysicalSpace);

  void SetIterations(const std::vector<std::vector<unsigned int> > & Iterations);

  void SetSmoothingSigmas(const std::vector<std::vector<float> > & SmoothingSigmas);

  void SetShrinkFactors(const std::vector<std::vector<unsigned int> > & ShrinkFactors);

  itkSetMacro(UseHistogramMatching, bool);
  itkGetMacro(UseHistogramMatching, bool);

  void SetWinsorizeImageIntensities(bool Winsorize, float LowerQuantile = 0.0, float UpperQuantile = 1.0);

  itkGetObjectMacro(CompositeTransform, CompositeTransformType);

  itkGetObjectMacro(WarpedImage, ImageType);
  itkGetObjectMacro(InverseWarpedImage, ImageType);

  int DoRegistration();

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

  int SetupInitialTransform(typename CompositeTransformType::Pointer & compositeTransform);
  std::ostream & Logger() const
  {
    return *m_LogStream;
  }

  typename CompositeTransformType::Pointer m_CompositeTransform;
  typename ImageType::Pointer              m_WarpedImage;
  typename ImageType::Pointer              m_InverseWarpedImage;
  bool                                    m_WriteOutputs;
  unsigned int                            m_NumberOfStages;
  std::string                             m_OutputTransformPrefix;
  std::string                             m_OutputWarpedImageName;
  std::string                             m_OutputInverseWarpedImageName;
  InitialTransformListType                m_InitialTransforms;
  MetricListType                          m_Metrics;
  TransformMethodListType                 m_TransformMethods;
  std::vector<std::vector<unsigned int> > m_Iterations;
  std::vector<std::vector<float> >        m_SmoothingSigmas;
  std::vector<std::vector<unsigned int> > m_ShrinkFactors;
  bool                                    m_UseHistogramMatching;
  bool                                    m_WinsorizeImageIntensities;
  double                                  m_LowerQuantile;
  double                                  m_UpperQuantile;
  std::ostream *                          m_LogStream;
};
} // namespace ants
} // namespace itk

#include "itkantsRegistrationHelper.hxx"

#endif
