/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: itkDiReCTImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/11/15 23:46:06 $
  Version:   $Revision: 1.16 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDiReCTImageFilter_h
#define __itkDiReCTImageFilter_h

#include "itkImageToImageFilter.h"

#include "itkVector.h"

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
 * \author Nicholas J. Tustison
 *
 * \par REFERENCE
 * S. R. Das, B. B. Avants, M. Grossman, and J. C. Gee, "Registration based
 * cortical thickness measurement," Neuroimage 2009, 45:867--879.
 *
 */

template <class TInputImage, class TOutputImage>
class DiReCTImageFilter :
  public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef DiReCTImageFilter                             Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Extract dimension from input and output image. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TInputImage::ImageDimension );

  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage                         InputImageType;
  typedef typename InputImageType::Pointer    InputImagePointer;
  typedef typename InputImageType::PixelType  InputPixelType;
  typedef InputPixelType                      LabelType;
  typedef typename InputImageType::RegionType RegionType;
  typedef typename InputImageType::IndexType  IndexType;
  typedef typename IndexType::IndexValueType  IndexValueType;
  typedef TOutputImage                        OutputImageType;
  typedef TOutputImage                        RealImageType;

  typedef typename OutputImageType::PixelType       RealType;
  typedef typename RealImageType::Pointer           RealImagePointer;
  typedef Vector<RealType, ImageDimension>          VectorType;
  typedef Image<VectorType, ImageDimension>         DisplacementFieldType;
  typedef typename DisplacementFieldType::Pointer   DisplacementFieldPointer;
  typedef typename VectorType::ValueType            VectorValueType;
  typedef typename DisplacementFieldType::PointType PointType;

  /**
   * Set the segmentation image.  The segmentation image is a labeled image
   * with specified labels for the gray and white matters.
   */
  void SetSegmentationImage( const InputImageType *seg )
  {
    this->SetNthInput( 0, const_cast<InputImageType *>( seg ) );
  }

  /**
   * Get the segmentation image.
   */
  const InputImageType * GetSegmentationImage() const
  {
    return this->GetInput( 0 );
  }


  /**
   * Get the label image.
   */
  const RealImageType * GetThicknessPriorImage() const
  {
    return this->m_ThicknessPriorImage;
  }
  /**
   * Set the label image.
   */
  void SetThicknessPriorImage( RealImagePointer seg )
  {
    this->m_ThicknessPriorImage = seg;
  }


  /**
   * Set the grey matter probability image.
   */
  void SetGrayMatterProbabilityImage( const RealImageType *gm )
  {
    this->SetNthInput( 1, const_cast<RealImageType *>( gm ) );
    this->Modified();
  }

  /**
   * Get the grey matter probability image.
   */
  const RealImageType * GetGrayMatterProbabilityImage() const
  {
    return static_cast<const RealImageType *>( this->ProcessObject::GetInput( 1 ) );
  }

  /**
   * Set the white matter probability image.
   */
  void SetWhiteMatterProbabilityImage( const RealImageType *wm )
  {
    this->SetNthInput( 2, const_cast<RealImageType *>( wm ) );
    this->Modified();
  }

  /**
   * Get the grey matter probability image.
   */
  const RealImageType * GetWhiteMatterProbabilityImage() const
  {
    return static_cast<const RealImageType *>( this->ProcessObject::GetInput( 2 ) );
  }

  /**
   * Get the warped white matter probability map.
   */
  const RealImageType * GetWarpedWhiteMatterProbabilityImage() const
    {
    return static_cast<const RealImageType *>( this->ProcessObject::GetOutput( 1 ) );
    }

  /**
   * Set the maximum number of registration iterations.  Default = 50.
   */
  itkSetMacro( MaximumNumberOfIterations, unsigned int );

  /**
   * Get the maximum number of registration iterations.  Default = 50.
   */
  itkGetConstMacro( MaximumNumberOfIterations, unsigned int );

  /**
   * Set the maximum number of registration iterations.  Default = 20.
   */
  itkSetMacro( MaximumNumberOfInvertDisplacementFieldIterations, unsigned int );

  /**
   * Get the maximum number of registration iterations.  Default = 20.
   */
  itkGetConstMacro( MaximumNumberOfInvertDisplacementFieldIterations, unsigned int );

  /**
   * Set the gray matter label in the segmentation image.  Default = 2.
   */
  itkSetMacro( GrayMatterLabel, LabelType );

  /**
   * Get the gray matter label in the segmentation image.  Default = 2.
   */
  itkGetConstMacro( GrayMatterLabel, LabelType );

  /**
   * Set the white matter label in the segmentation image.  Default = 3.
   */
  itkSetMacro( WhiteMatterLabel, LabelType );

  /**
   * Get the white matter label in the segmentation image.  Default = 3.
   */
  itkGetConstMacro( WhiteMatterLabel, LabelType );

  /**
   * Set the convergence threshold.  Default = 0.001.
   */
  itkSetMacro( ConvergenceThreshold, RealType );

  /**
   * Get the convergence threshold.  Default = 0.001.
   */
  itkGetConstMacro( ConvergenceThreshold, RealType );

  /**
   * Set the convergence window size. Convergence is determined by fitting a
   * line to the normalized energy profile of the last N iterations (where
   * N is specified by the window size) and determining the slope which is
   * then compared with the convergence threshold.  Default = 10.
   */
  itkSetMacro( ConvergenceWindowSize, unsigned int );

  /**
   * Get the convergence window size. Convergence is determined by fitting a
   * line to the normalized energy profile of the last N iterations (where
   * N is specified by the window size) and determining the slope which is
   * then compared with the convergence threshold. Default = 10.
   */
  itkGetConstMacro( ConvergenceWindowSize, unsigned int );

  /**
   * Set the thickness prior estimate---provides a constraint on the
   * final thickness measurement.  References in the literature
   * give a normal thickness of typically 3 mm with normal range
   * from ~2 mm in the calcarine cortex to ~4 mm in the precentral
   * gyrus. Default = 10 mm.
   */
  itkSetMacro( ThicknessPriorEstimate, RealType );

  /**
   * Get the thickness prior estimate---provides a constraint on the
   * final thickness measurement.  References in the literature
   * give a normal thickness of typically 3 mm with normal range
   * from ~2 mm in the calcarine cortex to ~4 mm in the precentral
   * gyrus.  Default = 10 mm.
   */
  itkGetConstMacro( ThicknessPriorEstimate, RealType );

  /**
   * Set the gradient step size.  Default = 0.025.
   */
  itkSetClampMacro( InitialGradientStep, RealType, 0, NumericTraits<RealType>::max() );

  /**
   * Get the gradient step.  Default = 0.025.
   */
  itkGetConstMacro( InitialGradientStep, RealType );

  /**
   * Get the current gradient step.
   */
  itkGetConstMacro( CurrentGradientStep, RealType );

  /**
   * Set the smoothing sigma for the velocity field (in voxels).  Default = 1.5.
   */
  itkSetClampMacro( SmoothingVelocityFieldVariance, RealType, 0, NumericTraits<RealType>::max() );

  /**
   * Get the smoothing sigma for the velocity field (in voxels).  Default = 1.5.
   */
  itkGetConstMacro( SmoothingVelocityFieldVariance, RealType );

  /**
   * Set the smoothing sigma for the total and hit images (in voxels).  Default = 1.0.
   */
  itkSetClampMacro( SmoothingVariance, RealType, 0, NumericTraits<RealType>::max() );

  /**
   * Set the smoothing sigma for the total and hit images (in voxels).  Default = 1.5.
   */
  itkGetConstMacro( SmoothingVariance, RealType );

  /**
   * Set the number of integration points.  Default = 10.
   */
  itkSetMacro( NumberOfIntegrationPoints, unsigned int  );

  /**
   * Get the number of integration points.  Default = 10.
   */
  itkGetConstMacro( NumberOfIntegrationPoints, unsigned int );

  /**
   * Get the number of elapsed iterations.  This is a helper function for
   * reporting observations.
   */
  itkGetConstMacro( ElapsedIterations, unsigned int );

  /**
   * Get the current energy.  This is a helper function for reporting
   * observations.
   */
  itkGetConstMacro( CurrentEnergy, RealType );

  /**
   * Get the current convergence measurement.  This is a helper function for
   * reporting observations.
   */
  itkGetConstMacro( CurrentConvergenceMeasurement, RealType );

protected:

  DiReCTImageFilter();
  virtual ~DiReCTImageFilter();

  void PrintSelf( std::ostream& os, Indent indent ) const;

  void ThreadedGenerateData( const RegionType &, ThreadIdType )
  {
  };

  void GenerateData();

private:

  /**
   * Private function for extracting regions (e.g. gray or white).
   */
  InputImagePointer ExtractRegion( const InputImageType *, LabelType );

  /**
   * Private function for extracting regional contours (e.g. gray or white).
   */
  InputImagePointer ExtractRegionalContours( const InputImageType *, LabelType );

  /**
   * Private function for making thickness image.
   */
  void MakeThicknessImage( RealImagePointer, RealImagePointer, InputImagePointer, RealImagePointer );

  /**
   * Private function for warping an image.
   */
  RealImagePointer WarpImage( const RealImageType *, const DisplacementFieldType * );

  /**
   * Private function for inverting the deformation field.
   */
  void InvertDisplacementField( const DisplacementFieldType *, DisplacementFieldType * );

  /**
   * Private function for smoothing the deformation field.
   */
  DisplacementFieldPointer SmoothDisplacementField( const DisplacementFieldType *, const RealType );

  /**
   * Private function for smoothing the image.
   */
  RealImagePointer SmoothImage( const RealImageType *, const RealType );

  RealType     m_ThicknessPriorEstimate;
  RealType     m_SmoothingVariance;
  RealType     m_SmoothingVelocityFieldVariance;
  RealType     m_InitialGradientStep;
  RealType     m_CurrentGradientStep;
  unsigned int m_NumberOfIntegrationPoints;

  LabelType m_GrayMatterLabel;
  LabelType m_WhiteMatterLabel;

  unsigned int m_ElapsedIterations;
  unsigned int m_MaximumNumberOfIterations;
  unsigned int m_MaximumNumberOfInvertDisplacementFieldIterations;
  RealType     m_CurrentEnergy;
  RealType     m_CurrentConvergenceMeasurement;
  RealType     m_ConvergenceThreshold;
  unsigned int m_ConvergenceWindowSize;

  RealImagePointer  m_ThicknessPriorImage;

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiReCTImageFilter.hxx"
#endif

#endif
