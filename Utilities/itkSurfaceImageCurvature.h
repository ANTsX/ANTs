/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _SurfaceImageCurvature_h
#define _SurfaceImageCurvature_h

#include "itkNeighborhoodIterator.h"

#include "itkSurfaceCurvatureBase.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"

namespace itk
{
/** \class SurfaceImageCurvature
 *
 * This class takes a surface as input and creates a local
 * geometric frame for each surface point.
 *
 *
 */
template <typename TSurface>
class SurfaceImageCurvature final : public SurfaceCurvatureBase<TSurface, 3>
{
public:
  /** Standard class typedefs. */
  typedef SurfaceImageCurvature          Self;
  typedef SurfaceCurvatureBase<TSurface> Superclass;
  typedef SmartPointer<Self>             Pointer;
  typedef SmartPointer<const Self>       ConstPointer;

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(SurfaceImageCurvature);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Image related types. */
  typedef typename TSurface::PixelType PixelType;
  enum
  {
    ImageDimension = TSurface::ImageDimension
  };
  typedef Image<PixelType, Self::ImageDimension> ImageType;
  typedef typename ImageType::IndexType                            IndexType;
  typedef typename ImageType::SpacingType                          SpacingType;
  typedef typename ImageType::SizeType                             SizeType;
  typedef ImageRegionIteratorWithIndex<ImageType>                  ImageIteratorType;
  /** Image dimension. */
  static constexpr unsigned int SurfaceDimension = TSurface::ImageDimension;

  typedef typename Superclass::RealType   RealType;
  typedef typename Superclass::PointType  VectorType;
  typedef typename Superclass::PointType  FixedVectorType;
  typedef typename Superclass::PointType  PointType;
  typedef typename Superclass::MatrixType MatrixType;
  typedef typename ImageType::PointType   ImagePointType;

  typedef Image<PixelType, Self::ImageDimension> OutputImageType;
  typedef ImageRegionIteratorWithIndex<OutputImageType>            OutputImageIteratorType;

  typedef typename OutputImageType::Pointer OutputImagePointer;

  typedef Image<MatrixType, Self::ImageDimension> FrameImageType;

  /** Gradient filtering */
  typedef CovariantVector<RealType, Self::ImageDimension>        GradientPixelType;
  typedef Image<GradientPixelType, Self::ImageDimension>         GradientImageType;
  typedef itk::VectorLinearInterpolateImageFunction<GradientImageType, RealType>   VectorInterpolatorType;
  typedef SmartPointer<GradientImageType>                                          GradientImagePointer;
  typedef GradientRecursiveGaussianImageFilter<OutputImageType, GradientImageType> GradientImageFilterType;
  typedef GradientImageFilter<OutputImageType>                                     GradientImageFilterType2;
  typedef typename GradientImageFilterType::Pointer                                GradientImageFilterPointer;

  typedef NeighborhoodIterator<ImageType> NeighborhoodIteratorType;

  /** Find all points within some distance of the origin.
   * The argument gives the number of times to apply the
   * mean shift algorithm to find the best neighborhood.
   */
  void
  FindNeighborhood(unsigned int numMeanShifts = 0) override;

  void
  FindEuclideanNeighborhood(PointType p);

  void
  FindGeodesicNeighborhood();

  /** This applies one of the algorithms for finding the local curvature
      and frame.  The default is joshi. */
  void
  ComputeFrameOverDomain(unsigned int which = 0) override;

  ImageType *
  GetInput();

  virtual void
  SetInputImage(typename ImageType::Pointer & input);
  OutputImageType *
  GetOutput();

  /** Apply the level set curvature equation over the whole image */
  void
  LevelSetMeanCurvature();

  /** Use the gradient of the image to estimate the normals everywhere.
   *  Also compute area, if boolean is set.
   */
  void
  EstimateNormalsFromGradient();

  /** Use the Weingarten map to estimate the curvature.*/
  void
  WeingartenMap();

  /** Use the Weingarten map to estimate the curvature.*/
  void
  WeingartenMapGradients();

  /** Computes a neighborhood surface area function everywhere*/
  void
  ComputeSurfaceArea();

  /** Use the gradient estimated normal to get the local frame.
      Requires call to SetNormal to find tangents. */
  void EstimateFrameFromGradient(IndexType);
  void EstimateFrameFromGradient(ImagePointType);

  /** Implemented version of virtual function from parent class.
      Here, we just sum the computed function, held in CurvatureImage,
      over a neighborhood */
  RealType
  IntegrateFunctionOverNeighborhood(bool norm = false) override;

  /** Get the neighborhood integral for every surface point.*/
  RealType
  IntegrateFunctionOverSurface(bool norm = false);

  /** Postprocess the curvature function by, e.g., gaussian
      smoothing of the curvature (and perhaps frame)
      in the local neighbhorhood. */
  void
  PostProcessGeometry();

  itkSetMacro(NeighborhoodRadius, RealType);
  itkGetMacro(NeighborhoodRadius, RealType);

  itkSetMacro(UseLabel, bool);
  itkGetMacro(UseLabel, bool);

  itkSetMacro(kSign, float);

  itkSetMacro(Threshold, float);

  itkGetMacro(FunctionImage, OutputImagePointer);
  itkSetMacro(FunctionImage, OutputImagePointer);

  void
  ProcessLabelImage();

  float
  CurvatureAtIndex(IndexType index)
  {
    PointType p;
    this->m_FunctionImage->TransformIndexToPhysicalPoint(index, p);
//    for (unsigned int k = 0; k < ImageDimension; k++) {
//      p[k] = (RealType)index[k];
//    }
    this->SetOrigin(p);
    this->EstimateFrameFromGradient(index);
    this->FindNeighborhood();
    this->WeingartenMap();
    float fval = this->m_GaussianKappa;
    float kpix;
    fval = this->m_MeanKappa;
    if (fabs(fval) > 1)
    {
      fval /= fval;
    }
    kpix = kpix + m_kSign * fval;

    return kpix;
  }

  inline RealType
  innerProduct(PointType v1, PointType v2)
  {
    return (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]);
  }

  inline bool
  IsValidIndex(IndexType ind)
  {
    for (unsigned int i = 0; i < ImageDimension; i++)
    {
      float shifted = ind[i];
      if (shifted < 1 || shifted > (this->m_ImageSize[i] - 1))
      {
        return false;
      }
    }
    return true;
  }

  inline bool
  IsValidSurface(PixelType pix, IndexType ind)
  {
    for (unsigned int i = 0; i < ImageDimension; i++)
    {
      float shifted = ind[i];
      if (shifted < 1 || shifted > (this->m_ImageSize[i] - 1))
      {
        return false;
      }
    }

    if (this->m_UseLabel)
    {
      if (itk::Math::FloatAlmostEqual(pix, this->m_SurfaceLabel))
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      if (pix > static_cast<PixelType>(this->m_Threshold))
      {
        return true;
      }
      else
      {
        return false;
      }
    }
  }

protected:
  SurfaceImageCurvature();
  ~SurfaceImageCurvature() override = default;

  void CopyImageToFunctionImage(OutputImagePointer, OutputImagePointer);

  /** This function changes the values of the label image for use with
      the fast marching image filter. */
private:
  PixelType                                m_SurfaceLabel;
  OutputImagePointer                       m_FunctionImage;
  RealType                                 m_NeighborhoodRadius;
  SizeType                                 m_ImageSize;
  GradientImagePointer                     m_GradientImage;
  NeighborhoodIteratorType                 m_ti;
  NeighborhoodIteratorType                 m_ti2;
  bool                                     m_UseLabel;
  float                                    m_kSign;
  float                                    m_Threshold;
  float                                    m_Area;
  RealType                                 m_MinSpacing;
  typename VectorInterpolatorType::Pointer m_Vinterp;
  SpacingType                              m_Spacing;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkSurfaceImageCurvature.hxx"
#endif

#endif
