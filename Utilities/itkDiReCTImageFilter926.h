/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: itkDiReCTImageFilter926.h,v $
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
#ifndef __itkDiReCTImageFilter926_h
#define __itkDiReCTImageFilter926_h

#include "itkImageToImageFilter.h"

#include "itkPointSet.h"
#include "itkVector.h"

namespace itk
{
/** \class DiReCTImageFilter926
 * \brief Diffeomorphic Registration-based Cortical Thickness measurement.
 *
 * To estimate the cortical thickness, the following inputs are required:
 *   - Segmentation image in which the csf, grey matter, and white matter voxels
 *     are all labeled with values of 1, 2, and 3, respectively.
 *   - Corresponding grey matter and white matter probability maps.
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
class DiReCTImageFilter926 :
  public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef DiReCTImageFilter926                          Self;
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
  typedef typename InputImageType::RegionType RegionType;
  typedef typename InputImageType::IndexType  IndexType;
  typedef typename IndexType::IndexValueType  IndexValueType;
  typedef TOutputImage                        OutputImageType;
  typedef TOutputImage                        RealImageType;

  typedef typename OutputImageType::PixelType RealType;
  typedef typename RealImageType::Pointer     RealImagePointer;
  typedef Vector<RealType,
                 itkGetStaticConstMacro( ImageDimension )>   VectorType;
  typedef Image<VectorType,
                itkGetStaticConstMacro( ImageDimension )>   VectorImageType;
  typedef typename VectorImageType::Pointer   VectorImagePointer;
  typedef typename VectorType::ValueType      VectorValueType;
  typedef typename VectorImageType::PointType PointType;

  typedef PointSet<RealType,
                   itkGetStaticConstMacro( ImageDimension )>      SparseImageType;
  typedef typename SparseImageType::Pointer SparseImagePointer;
  typedef PointSet<VectorType,
                   itkGetStaticConstMacro( ImageDimension )>      SparseVectorImageType;
  typedef typename SparseVectorImageType::Pointer SparseVectorImagePointer;
  typedef typename SparseImageType::PointType     SparseIndexType;

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
   * Set the gray matter label in the segmentation image.  Default = 2.
   */
  itkSetMacro( GrayMatterLabel, unsigned int );

  /**
   * Get the gray matter label in the segmentation image.  Default = 2.
   */
  itkGetConstMacro( GrayMatterLabel, unsigned int );

  /**
   * Set the white matter label in the segmentation image.  Default = 3.
   */
  itkSetMacro( WhiteMatterLabel, unsigned int );

  /**
   * Get the white matter label in the segmentation image.  Default = 3.
   */
  itkGetConstMacro( WhiteMatterLabel, unsigned int );

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
   * Set the gradient step size.  Default = 0.5.
   */
  itkSetClampMacro( GradientStep, RealType, 0, NumericTraits<RealType>::max() );

  /**
   * Get the gradient step.  Default = 0.5.
   */
  itkGetConstMacro( GradientStep, RealType );

  /**
   * Set the smoothing sigma.  Default = 1.5.
   */
  itkSetClampMacro( SmoothingSigma, RealType, 0, NumericTraits<RealType>::max() );

  /**
   * Get the smoothing sigma.  Default = 1.5.
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

  DiReCTImageFilter926();
  virtual ~DiReCTImageFilter926();

  void PrintSelf( std::ostream& os, Indent indent ) const;

  void ThreadedGenerateData( const RegionType &, ThreadIdType )
  {
  };

  void GenerateData();

private:

  /**
   * Private function for extracting regions (e.g. gray or white).
   */
  InputImagePointer ExtractRegion( const InputImageType *, unsigned int );

  /**
   * Private function for extracting regional contours (e.g. gray or white).
   */
  InputImagePointer ExtractRegionalContours( const InputImageType *, unsigned int );

  /**
   * Private function for warping an image.
   */
  RealImagePointer WarpImage( const RealImageType *, const VectorImageType * );

  /**
   * Private function for inverting the deformation field.
   */
  void InvertDisplacementField( const VectorImageType *, VectorImageType * );

  /**
   * Private function for smoothing the deformation field.
   */
  VectorImagePointer SmoothDisplacementField( const VectorImageType *, const RealType );

  /**
   * Private function for converting a scalar image to a sparse image
   */
  SparseImagePointer ConvertRealImageToSparseImage( const RealImageType *, const InputImageType * );

  /**
   * Private function for converting a sparse image to a scalar image
   */
  RealImagePointer ConvertSparseImageToRealImage( const SparseImageType *, const InputImageType * );

  /**
   * Private function for converting a vector image to a sparse vector image
   */
  SparseVectorImagePointer ConvertVectorImageToSparseVectorImage(const VectorImageType *, const InputImageType * );

  /**
   * Private function for converting a sparse vector image to a vector image
   */
  VectorImagePointer ConvertSparseVectorImageToVectorImage(const SparseVectorImageType *, const InputImageType * );

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
#include "itkDiReCTImageFilter926.hxx"
#endif

#endif
