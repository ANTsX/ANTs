/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/ANTsX/ANTs/blob/main/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDiReCTImageFilter_h
#define __itkDiReCTImageFilter_h

#include "ITKGetterSetterMacroShims.h"
#include "itkImageToImageFilter.h"

#include "itkVector.h"
#include <vnl/vnl_sparse_matrix.h>


namespace itk
{
/** \class DiReCTImageFilter
 * \brief Diffeomorphic Registration-based Cortical Thickness measurement.
 *
 * To estimate the cortical thickness, the following inputs are required:
 *   - Segmentation image in which the csf, grey matter, and white matter voxels
 *     are all labeled with values of 1, 2, and 3, respectively.
 *   - Corresponding grey matter and white matter probability maps.
 *
 * \author Nick Tustison, Brian Avants
 *
 * \par REFERENCE
 * S. R. Das, B. B. Avants, M. Grossman, and J. C. Gee, "Registration based
 * cortical thickness measurement," Neuroimage 2009, 45:867--879.
 *
 * Nicholas J. Tustison, Philip A. Cook, Arno Klein, Gang Song, Sandhitsu
 * R. Das, Jeffrey T. Duda, Benjamin M. Kandel, Niels van Strien, James R.
 * Stone, James C. Gee, and Brian B. Avants. Large-Scale Evaluation of ANTs
 * and FreeSurfer Cortical Thickness Measurements. NeuroImage, 99:166-179,
 * Oct 2014.
 *
 */

template <typename TInputImage, typename TOutputImage>
class DiReCTImageFilter final : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  using Self = DiReCTImageFilter;
  using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  using DataObjectPointer = DataObject::Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Extract dimension from input and output image. */
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;

  /** Convenient typedefs for simplifying declarations. */
  using InputImageType = TInputImage;
  using InputImagePointer = typename InputImageType::Pointer;
  using InputPixelType = typename InputImageType::PixelType;
  using RegionType = typename InputImageType::RegionType;
  using IndexType = typename InputImageType::IndexType;

  using LabelType = InputPixelType;

  using OutputImageType = TOutputImage;
  using RealType = typename OutputImageType::PixelType;
  using RealImageType = TOutputImage;
  using RealImagePointer = typename RealImageType::Pointer;

  using VectorType = Vector<RealType, ImageDimension>;
  using DisplacementFieldType = Image<VectorType, ImageDimension>;
  using DisplacementFieldPointer = typename DisplacementFieldType::Pointer;
  using CumulativeVelocityFieldType = Image<VectorType, ImageDimension+1>;
  using CumulativeVelocityFieldPointer = typename CumulativeVelocityFieldType::Pointer;
  using VectorValueType = typename VectorType::ValueType;
  using PointType = typename DisplacementFieldType::PointType;
  using SparseMatrixType = vnl_sparse_matrix<RealType>;

  /**
   * Set the segmentation image.  The segmentation image is a labeled image
   * with specified labels for the gray and white matters.
   */
  void
  SetSegmentationImage(const InputImageType * seg)
  {
    this->SetNthInput(0, const_cast<InputImageType *>(seg));
  }

  /**
   * Get the segmentation image.
   */
  const InputImageType *
  GetSegmentationImage() const
  {
    return this->GetInput(0);
  }


  /**
   * Get the label image.
   */
  const RealImageType *
  GetThicknessPriorImage() const
  {
    return this->m_ThicknessPriorImage;
  }
  /**
   * Set the label image.
   */
  void
  SetThicknessPriorImage(RealImagePointer seg)
  {
    this->m_ThicknessPriorImage = seg;
  }


  /**
   * Set the grey matter probability image.
   */
  void
  SetGrayMatterProbabilityImage(const RealImageType * gm)
  {
    this->SetNthInput(1, const_cast<RealImageType *>(gm));
    this->Modified();
  }

  /**
   * Get the grey matter probability image.
   */
  const RealImageType *
  GetGrayMatterProbabilityImage() const
  {
    return static_cast<const RealImageType *>(this->ProcessObject::GetInput(1));
  }

  /**
   * Set the white matter probability image.
   */
  void
  SetWhiteMatterProbabilityImage(const RealImageType * wm)
  {
    this->SetNthInput(2, const_cast<RealImageType *>(wm));
    this->Modified();
  }

  /**
   * Get the grey matter probability image.
   */
  const RealImageType *
  GetWhiteMatterProbabilityImage() const
  {
    return static_cast<const RealImageType *>(this->ProcessObject::GetInput(2));
  }

  /**
   * Get the warped white matter probability map.
   */
  const RealImageType *
  GetWarpedWhiteMatterProbabilityImage() const
  {
    return static_cast<const RealImageType *>(this->ProcessObject::GetOutput(1));
  }

  /**
   * Set/Get the maximum number of registration iterations.  Default = 50.
   */
  itkNonVirtualSetMacro(MaximumNumberOfIterations, unsigned int);
  itkNonVirtualGetConstMacro(MaximumNumberOfIterations, unsigned int);

  /**
   * Set/Get the variance for time regularization.
   */
  itkNonVirtualSetMacro(TimeSmoothingVariance, RealType);
  itkNonVirtualGetConstMacro(TimeSmoothingVariance, RealType);

  /**
   * Set/Get the maximum number of inversion iterations.  Default = 20.
   */
  itkNonVirtualSetMacro(MaximumNumberOfInvertDisplacementFieldIterations, unsigned int);
  itkNonVirtualGetConstMacro(MaximumNumberOfInvertDisplacementFieldIterations, unsigned int);

  /**
   * Set/Get the gray matter label in the segmentation image.  Default = 2.
   */
  itkNonVirtualSetMacro(GrayMatterLabel, LabelType);
  itkNonVirtualGetConstMacro(GrayMatterLabel, LabelType);

  /**
   * Set/Get the white matter label in the segmentation image.  Default = 3.
   */
  itkNonVirtualSetMacro(WhiteMatterLabel, LabelType);
  itkNonVirtualGetConstMacro(WhiteMatterLabel, LabelType);

  /**
   * Set/Get the convergence threshold.  Default = 0.001.
   */
  itkNonVirtualSetMacro(ConvergenceThreshold, RealType);
  itkNonVirtualGetConstMacro(ConvergenceThreshold, RealType);

  /**
   * Set/Get the convergence window size. Convergence is determined by fitting a
   * line to the normalized energy profile of the last N iterations (where
   * N is specified by the window size) and determining the slope which is
   * then compared with the convergence threshold.  Default = 10.
   */
  itkNonVirtualSetMacro(ConvergenceWindowSize, unsigned int);
  itkNonVirtualGetConstMacro(ConvergenceWindowSize, unsigned int);

  /**
   * Set/Get the thickness prior estimate---provides a constraint on the
   * final thickness measurement.  References in the literature
   * give a normal thickness of typically 3 mm with normal range
   * from ~2 mm in the calcarine cortex to ~4 mm in the precentral
   * gyrus. Default = 10 mm.
   */
  itkNonVirtualSetMacro(ThicknessPriorEstimate, RealType);
  itkNonVirtualGetConstMacro(ThicknessPriorEstimate, RealType);

  /**
   * Set/Get the gradient step size.  Default = 0.025.
   */
  itkNonVirtualSetClampMacro(InitialGradientStep, RealType, 0, NumericTraits<RealType>::max());
  itkNonVirtualGetConstMacro(InitialGradientStep, RealType);

  /**
   * Get the current gradient step.
   */
  itkNonVirtualGetConstMacro(CurrentGradientStep, RealType);

  /**
   * Set/Get the gradient smoothing sigma (in mm).  Default = 1.0.
   */
  itkNonVirtualSetClampMacro(GradientSmoothingSigma, RealType, 0, NumericTraits<RealType>::max());
  itkNonVirtualGetConstMacro(GradientSmoothingSigma, RealType);

  /**
   * Set/Get the smoothing variance for the total and hit images (in voxel units squared).  Default = 1.0.
   */
  itkNonVirtualSetClampMacro(SmoothingVariance, RealType, 0, NumericTraits<RealType>::max());
  itkNonVirtualGetConstMacro(SmoothingVariance, RealType);

  /**
   * Set/Get the isotropic mesh spacing for smoothing the velocity field (in mm).  Default = 5.75.
   */
  itkNonVirtualSetClampMacro(BSplineSmoothingIsotropicMeshSpacing, RealType, 0, NumericTraits<RealType>::max());
  itkNonVirtualGetConstMacro(BSplineSmoothingIsotropicMeshSpacing, RealType);

  /**
   * Set/Get the B-spline smoothing variance for the velocity field (in voxels).  Default = 1.5.
   */
  itkNonVirtualSetClampMacro(SmoothingVelocityFieldVariance, RealType, 0, NumericTraits<RealType>::max());
  itkNonVirtualGetConstMacro(SmoothingVelocityFieldVariance, RealType);

  /**
   * Set/Get the number of integration points.  Default = 10.
   */
  itkNonVirtualSetMacro(NumberOfIntegrationPoints, unsigned int);
  itkNonVirtualGetConstMacro(NumberOfIntegrationPoints, unsigned int);

  /**
   * Set/Get the sparse image neighborhood radius.  Default = 2.
   */
  itkNonVirtualSetMacro(SparseImageNeighborhoodRadius, unsigned int);
  itkNonVirtualGetConstMacro(SparseImageNeighborhoodRadius, unsigned int);

  /**
   * Set/Get the option to restrict deformation along the last dimension.  Default = false.
   */
  itkNonVirtualSetMacro(RestrictDeformation, bool);
  itkNonVirtualGetConstMacro(RestrictDeformation, bool);
  itkNonVirtualBooleanMacro(RestrictDeformation);

  /**
   * Set/Get the temporal points values.  Default = no special value.
   */
  void
  SetTimePoints(std::vector<RealType> timePoints)
  {
    this->m_TimePoints = timePoints;
    this->Modified();
  }
  itkNonVirtualGetConstMacro(TimePoints, std::vector<RealType>);

  /**
   * Set/Get the option to use masked smoothing.  Default = false.
   */
  itkNonVirtualSetMacro(UseMaskedSmoothing, bool);
  itkNonVirtualGetConstMacro(UseMaskedSmoothing, bool);
  itkNonVirtualBooleanMacro(UseMaskedSmoothing);

  /**
   * Set/Get the option to use B-spline smoothing.  Default = false.
   */
  itkNonVirtualSetMacro(UseBSplineSmoothing, bool);
  itkNonVirtualGetConstMacro(UseBSplineSmoothing, bool);
  itkNonVirtualBooleanMacro(UseBSplineSmoothing);

  /**
   * Set/Get the option to include the cumulative velocity fields as output.  Default = false.
   */
  itkNonVirtualSetMacro(IncludeCumulativeVelocityFields, bool);
  itkNonVirtualGetConstMacro(IncludeCumulativeVelocityFields, bool);
  itkNonVirtualBooleanMacro(IncludeCumulativeVelocityFields);

  /**
   * Get the number of elapsed iterations.  This is a helper function for
   * reporting observations.
   */
  itkNonVirtualGetConstMacro(ElapsedIterations, unsigned int);

  /**
   * Get the current energy.  This is a helper function for reporting
   * observations.
   */
  itkNonVirtualGetConstMacro(CurrentEnergy, RealType);

  /**
   * Get the current convergence measurement.  This is a helper function for
   * reporting observations.
   */
  itkNonVirtualGetConstMacro(CurrentConvergenceMeasurement, RealType);

  /**
   * Get thickness image.
   */
  OutputImageType * GetThicknessImage();

  /**
   * Get forward cumulative velocity
   */
  CumulativeVelocityFieldType * GetForwardCumulativeVelocityField();

  /**
   * Get inverse cumulative velocity
   */
  CumulativeVelocityFieldType * GetInverseCumulativeVelocityField();

  /** Standard itk::ProcessObject subclass method. */
  using DataObjectPointerArraySizeType = ProcessObject::DataObjectPointerArraySizeType;
  using Superclass::MakeOutput;
  DataObjectPointer
  MakeOutput(DataObjectPointerArraySizeType idx) override;

protected:
  DiReCTImageFilter();
  ~DiReCTImageFilter() override;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  GenerateOutputInformation() override;

  void
  ThreadedGenerateData(const RegionType &, ThreadIdType) override{};

  void
  GenerateData() override;

private:
  /**
   * Private function to determine if a voxel is in the mask.
   */
  bool
  IsInsideMask(IndexType index)
  {
    return (this->GetSegmentationImage()->GetPixel(index) == this->m_GrayMatterLabel |
            this->GetSegmentationImage()->GetPixel(index) == this->m_WhiteMatterLabel);
  }

  /**
   * Private function for extracting regions (e.g. gray or white).
   */
  InputImagePointer
  ExtractRegion(const InputImageType *, LabelType);

  /**
   * Private function for extracting regional contours (e.g. gray or white).
   */
  InputImagePointer
  ExtractRegionalContours(const InputImageType *, LabelType);

  /**
   * Private function for making thickness image.
   */
  void MakeThicknessImage(RealImagePointer, RealImagePointer, InputImagePointer, RealImagePointer);

  /**
   * Private function for warping an image.
   */
  RealImagePointer
  WarpImage(const RealImageType *, const DisplacementFieldType *);

  /**
   * Private function for smoothing the deformation field.
   */
  DisplacementFieldPointer
  GaussianSmoothDisplacementField(const DisplacementFieldType *, const RealType);

  /**
   * Private function for smoothing the deformation field.
   */
  DisplacementFieldPointer
  MaskedGaussianSmoothDisplacementField(const DisplacementFieldType *);

  /**
   * Private function for smoothing the deformation field.
   */
  DisplacementFieldPointer
  BSplineSmoothDisplacementField(const DisplacementFieldType *, const RealType);

  /**
   * Private function for smoothing the image.
   */
  RealImagePointer
  SmoothImage(const RealImageType *, const RealType);

  RealType     m_ThicknessPriorEstimate;
  RealType     m_GradientSmoothingSigma;
  RealType     m_SmoothingVariance;
  RealType     m_SmoothingVelocityFieldVariance;
  RealType     m_BSplineSmoothingIsotropicMeshSpacing;
  RealType     m_InitialGradientStep;
  RealType     m_CurrentGradientStep;
  unsigned int m_NumberOfIntegrationPoints;

  unsigned int m_SparseImageNeighborhoodRadius;

  LabelType m_GrayMatterLabel;
  LabelType m_WhiteMatterLabel;

  unsigned int m_ElapsedIterations;
  unsigned int m_MaximumNumberOfIterations;
  unsigned int m_MaximumNumberOfInvertDisplacementFieldIterations;
  RealType     m_CurrentEnergy;
  RealType     m_CurrentConvergenceMeasurement;
  RealType     m_ConvergenceThreshold;
  unsigned int m_ConvergenceWindowSize;

  RealImagePointer m_ThicknessPriorImage;

  bool                  m_UseBSplineSmoothing;
  bool                  m_UseMaskedSmoothing;
  bool                  m_RestrictDeformation;
  bool                  m_IncludeCumulativeVelocityFields;
  SparseMatrixType      m_SparseMatrix;
  RealImagePointer      m_SparseMatrixIndexImage;
  std::vector<RealType> m_TimePoints;
  RealType              m_TimeSmoothingVariance;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkDiReCTImageFilter.hxx"
#endif

#endif
