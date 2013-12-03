/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: itkSurfaceImageCurvature.h,v $
  Language:  C++
  Date:      $Date: 2008/11/15 23:46:06 $
  Version:   $Revision: 1.12 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

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
class SurfaceImageCurvature :
  public SurfaceCurvatureBase<TSurface,
                              3>
{
public:

  /** Standard class typedefs. */
  typedef SurfaceImageCurvature          Self;
  typedef SurfaceCurvatureBase<TSurface> Superclass;
  typedef SmartPointer<Self>             Pointer;
  typedef SmartPointer<const Self>       ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(SurfaceImageCurvature, SurfaceCurvatureBase);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Image related types. */
  typedef typename TSurface::PixelType PixelType;
  enum { ImageDimension = TSurface::ImageDimension };
  typedef Image<PixelType, itkGetStaticConstMacro(ImageDimension)>
    ImageType;
  typedef typename ImageType::IndexType           IndexType;
  typedef typename ImageType::SizeType            SizeType;
  typedef ImageRegionIteratorWithIndex<ImageType> ImageIteratorType;
  /** Image dimension. */
  itkStaticConstMacro(SurfaceDimension, unsigned int, TSurface::ImageDimension);

  typedef typename Superclass::RealType   RealType;
  typedef typename Superclass::PointType  VectorType;
  typedef typename Superclass::PointType  FixedVectorType;
  typedef typename Superclass::PointType  PointType;
  typedef typename Superclass::MatrixType MatrixType;

  typedef  Image<PixelType, itkGetStaticConstMacro(ImageDimension)>
    OutputImageType;
  typedef ImageRegionIteratorWithIndex<OutputImageType> OutputImageIteratorType;

  typedef typename OutputImageType::Pointer OutputImagePointer;

  typedef Image<MatrixType, itkGetStaticConstMacro(ImageDimension)>
    FrameImageType;

  /** Gradient filtering */
  typedef CovariantVector<RealType,
                          itkGetStaticConstMacro(ImageDimension)> GradientPixelType;
  typedef Image<GradientPixelType,
                itkGetStaticConstMacro(ImageDimension)> GradientImageType;
  typedef SmartPointer<GradientImageType> GradientImagePointer;
  typedef GradientRecursiveGaussianImageFilter<OutputImageType, GradientImageType>
    GradientImageFilterType;
  typedef GradientImageFilter<OutputImageType>
    GradientImageFilterType2;
  typedef typename GradientImageFilterType::Pointer GradientImageFilterPointer;

  typedef NeighborhoodIterator<ImageType> NeighborhoodIteratorType;

  /** Find all points within some distance of the origin.
    * The argument gives the number of times to apply the
    * mean shift algorithm to find the best neighborhood.
    */
  void FindNeighborhood(unsigned int numMeanShifts = 0);

  void  FindEuclideanNeighborhood(PointType p);

  void  FindGeodesicNeighborhood();

  /** This applies one of the algorithms for finding the local curvature
      and frame.  The default is joshi. */
  void ComputeFrameOverDomain(unsigned int which = 0);

  ImageType * GetInput();

  virtual void SetInputImage(typename ImageType::Pointer & input);
  OutputImageType * GetOutput();

  /** Apply the level set curvature equation over the whole image */
  void LevelSetMeanCurvature();

  /** Use the gradient of the image to estimate the normals everywhere.
   *  Also compute area, if boolean is set.
  */
  void EstimateNormalsFromGradient();

  /** Use the Weingarten map to estimate the curvature.*/
  void WeingartenMap();

  /** Computes a neighborhood surface area function everywhere*/
  void ComputeSurfaceArea();

  /** Use the gradient estimated normal to get the local frame.
      Requires call to SetNormal to find tangents. */
  void EstimateFrameFromGradient(IndexType);

  /** Implemented version of virtual function from parent class.
      Here, we just sum the computed function, held in CurvatureImage,
      over a neighborhood */
  RealType IntegrateFunctionOverNeighborhood(bool norm = false);

  /** Get the neighborhood integral for every surface point.*/
  RealType IntegrateFunctionOverSurface(bool norm = false);

  /** Postprocess the curvature function by, e.g., gaussian
      smoothing of the curvature (and perhaps frame)
      in the local neighbhorhood. */
  void PostProcessGeometry();

  itkSetMacro(NeighborhoodRadius, RealType);
  itkGetMacro(NeighborhoodRadius, RealType);

  itkSetMacro(UseLabel, bool);
  itkGetMacro(UseLabel, bool);

  itkSetMacro(kSign, float);
  itkSetMacro(Sigma, float);

  itkSetMacro(Threshold, float);

  itkGetMacro(FunctionImage, OutputImagePointer);
  itkSetMacro(FunctionImage, OutputImagePointer);

  void ProcessLabelImage();

  float CurvatureAtIndex(IndexType index)
  {
    PointType p;

    for( unsigned int k = 0; k < ImageDimension; k++ )
      {
      p[k] = (RealType) index[k];
      }
    this->SetOrigin(p);
    this->EstimateFrameFromGradient(index);
    this->FindNeighborhood();
    this->WeingartenMap();
    float fval = this->m_GaussianKappa;
    float kpix;
    fval = this->m_MeanKappa;
    if( fabs(fval) > 1 )
      {
      fval /= fval;
      }
    kpix = kpix + m_kSign * fval;

    return kpix;
  }

  inline bool IsValidSurface(PixelType pix, IndexType /* ind */)
  {
    //   std::cout << "m_UseLabel "<< m_UseLabel << " pix " << pix << " Thresh " << m_Threshold << " ind " << ind
    // <<
    // std::endl;
    if( this->m_UseLabel )
      {
      if( pix == this->m_SurfaceLabel )
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
      if( pix > this->m_Threshold )
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
  ~SurfaceImageCurvature()
  {
  };

  void CopyImageToFunctionImage( OutputImagePointer, OutputImagePointer);

  /** This function changes the values of the label image for use with
      the fast marching image filter. */
private:

  PixelType                m_SurfaceLabel;
  OutputImagePointer       m_FunctionImage;
  RealType                 m_NeighborhoodRadius;
  SizeType                 m_ImageSize;
  GradientImagePointer     m_GradientImage;
  NeighborhoodIteratorType m_ti;
  NeighborhoodIteratorType m_ti2;
  bool                     m_UseLabel;
  float                    m_kSign;
  float                    m_Sigma;
  float                    m_Threshold;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSurfaceImageCurvature.hxx"
#endif

#endif
