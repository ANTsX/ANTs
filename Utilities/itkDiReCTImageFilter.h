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
 * In addition to specifying the input point set, one must specify the number
 * of control points.  The specified number of control points must be
 * > m_SplineOrder.  If one wishes to use the multilevel component of
 * this algorithm, one must also specify the number of levels in the
 * hierarchy.  If this is desired, the number of control points becomes
 * the number of control points for the coarsest level.  The algorithm
 * then increases the number of control points at each level so that
 * the B-spline n-D grid is refined to twice the previous level.
 *
 * \author Nicholas J. Tustison
 *
 * Contributed by Nicholas J. Tustison, Brian B. Avants
 * in the Insight Journal paper:
 *
 * \par REFERENCE
 * S. R. Das, B. B. Avants, M. Grossman, and J. C. Gee, "Registration based
 * cortical thickness measurement," Neuroimage 2009, 45:867--879.
 *
 * \par REFERENCE
 * S. E. Jones, B. R. Buchbinder, and I Aharon, "Three-dimensional mapping
 * of Cortical thickness using Laplace's Equation." Human Brian Mapping 2000,
 * 11:12-32.
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
  typedef TInputImage                        InputImageType;
  typedef typename InputImageType::Pointer   InputImagePointer;
  typedef typename InputImageType::PixelType InputPixelType;
  typedef TOutputImage                       OutputImageType;

  typedef double RealType;
  typedef Image<RealType,
                itkGetStaticConstMacro( ImageDimension )>   RealImageType;
  typedef typename RealImageType::Pointer RealImagePointer;
  typedef Vector<RealType,
                 itkGetStaticConstMacro( ImageDimension )>   VectorType;
  typedef Image<VectorType,
                itkGetStaticConstMacro( ImageDimension )>   VectorImageType;
  typedef typename VectorImageType::Pointer   VectorImagePointer;
  typedef typename VectorType::ValueType      VectorValueType;
  typedef typename VectorImageType::PointType PointType;

  /**
   * Set the segmentation image.  The segmentation image is a labeled image
   * with voxel values of '0' for background, '1' for csf, '2' for gm, and
   * '3' for wm.
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
    return static_cast<const RealImageType *>(
      this->ProcessObject::GetInput( 1 ) );
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
    return static_cast<const RealImageType *>(
      this->ProcessObject::GetInput( 2 ) );
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
   * Set the gray matter label in the segmentation image.
   */
  itkSetMacro( GrayMatterLabel, unsigned int );

  /**
   * Get the gray matter label in the segmentation image.
   */
  itkGetConstMacro( GrayMatterLabel, unsigned int );

  /**
   * Set the white matter label in the segmentation image.
   */
  itkSetMacro( WhiteMatterLabel, unsigned int );

  /**
   * Get the white matter label in the segmentation image.
   */
  itkGetConstMacro( WhiteMatterLabel, unsigned int );

  /**
   * Set the convergence threshold.  Convergence is determined by the
   */
  itkSetMacro( ConvergenceThreshold, RealType );

  /**
   * Get the convergence threshold.  Convergence is determined by the
   */
  itkGetConstMacro( ConvergenceThreshold, RealType );

  /**
   * Set the convergence threshold.  Convergence is determined by the
   */
  itkSetMacro( ConvergenceWindowSize, unsigned int );

  /**
   * Get the convergence threshold.  Convergence is determined by the
   */
  itkGetConstMacro( ConvergenceWindowSize, unsigned int );

  /**
   * Set the thickness prior estimate---provides a constraint on the
   * final thickness measurement.  References in the literature
   * give a normal thickness of typically 3 mm with normal range
   * from ~2 mm in the calcarine cortex to ~4 mm in the precentral
   * gyrus. Default = 6 mm.
   */
  itkSetMacro( ThicknessPriorEstimate, RealType );

  /**
   * Get the thickness prior estimate---provides a constraint on the
   * final thickness measurement.  References in the literature
   * give a normal thickness of typically 3 mm with normal range
   * from ~2 mm in the calcarine cortex to ~4 mm in the precentral
   * gyrus.  Default = 6 mm.
   */
  itkGetConstMacro( ThicknessPriorEstimate, RealType );

  /**
   * Set the gradient step size.  Default = 0.5.
   */
  itkSetClampMacro( GradientStep, RealType, 0, NumericTraits<RealType>::max() );

  /**
   * Get the gradient step.  Default = 0.5.
   */
  itkGetConstMacro( GradientStep, RealType );

  /**
   * Set the smoothing sigma.  Default = 1.0.
   */
  itkSetClampMacro( SmoothingSigma, RealType, 0, NumericTraits<RealType>::max() );

  /**
   * Get the smoothing sigma.  Default = 1.0.
   */
  itkGetConstMacro( SmoothingSigma, RealType );

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

  void GenerateData();

private:

  /**
   */
  InputImagePointer ExtractRegion( const InputImageType *, unsigned int );

  /**
   */
  InputImagePointer ExtractRegionalContours( const InputImageType *, unsigned int );

  /**
   */
  RealImagePointer WarpImage( const RealImageType *, const VectorImageType * );

  /**
   */
  void InvertDeformationField( const VectorImageType *, VectorImageType * );

  /**
   */
  VectorImagePointer SmoothDeformationField( const VectorImageType *, const RealType );

  RealType     m_ThicknessPriorEstimate;
  RealType     m_SmoothingSigma;
  RealType     m_GradientStep;
  unsigned int m_NumberOfIntegrationPoints;

  unsigned int m_GrayMatterLabel;
  unsigned int m_WhiteMatterLabel;

  unsigned int m_ElapsedIterations;
  unsigned int m_MaximumNumberOfIterations;
  RealType     m_CurrentEnergy;
  RealType     m_CurrentConvergenceMeasurement;
  RealType     m_ConvergenceThreshold;
  unsigned int m_ConvergenceWindowSize;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiReCTImageFilter.txx"
#endif

#endif
