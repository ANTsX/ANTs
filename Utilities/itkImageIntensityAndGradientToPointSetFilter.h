/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkImageIntensityAndGradientToPointSetFilter_h
#define itkImageIntensityAndGradientToPointSetFilter_h

#include "itkCovariantVector.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkImage.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkPointSet.h"

namespace itk
{
/** \class ImageIntensityAndGradientToPointSetFilter
 * \brief
 * Reads a file and creates an ikt point set.
 *
 */
template <typename TInputImage, typename TMaskImage, typename TOutputMesh>
class ImageIntensityAndGradientToPointSetFilter final : public MeshSource<TOutputMesh>
{
public:
  /** Standard "Self" typedef. */
  typedef ImageIntensityAndGradientToPointSetFilter Self;
  typedef MeshSource<TOutputMesh>                   Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Extract dimension from the input image. */
  static constexpr unsigned int Dimension = TInputImage::ImageDimension;

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(ImageIntensityAndGradientToPointSetFilter);

  /** Hold on to the type information specified by the template parameters. */
  typedef TInputImage                        InputImageType;
  typedef typename InputImageType::PixelType InputImagePixelType;

  typedef TMaskImage MaskImageType;

  typedef TOutputMesh                         OutputMeshType;
  typedef typename OutputMeshType::MeshTraits MeshTraits;
  typedef typename OutputMeshType::Superclass PointSetType;
  typedef typename OutputMeshType::PointType  PointType;
  typedef typename MeshTraits::PixelType      PointSetPixelType;

  typedef CovariantVector<InputImagePixelType, Dimension> GradientPixelType;
  typedef Image<GradientPixelType, Dimension>             GradientImageType;

  typedef ConstNeighborhoodIterator<GradientImageType>       ConstNeighborhoodIteratorType;
  typedef typename ConstNeighborhoodIteratorType::RadiusType NeighborhoodRadiusType;

  typedef GradientRecursiveGaussianImageFilter<InputImageType, GradientImageType> GradientFilterType;

  /**
   * Set/Get the input image.
   */
  void
  SetInputImage(const InputImageType * image)
  {
    this->SetNthInput(0, const_cast<InputImageType *>(image));
  }
  void
  SetInput1(const InputImageType * image)
  {
    this->SetInputImage(image);
  }

  /**
   * Get mask image function.
   */
  const InputImageType *
  GetInputImage() const
  {
    return static_cast<const InputImageType *>(this->ProcessObject::GetInput(0));
  }

  /**
   * Set mask image function.
   */
  void
  SetMaskImage(const MaskImageType * mask)
  {
    this->SetNthInput(1, const_cast<MaskImageType *>(mask));
  }
  void
  SetInput2(const MaskImageType * mask)
  {
    this->SetMaskImage(mask);
  }

  /**
   * Get mask image function.
   */
  const MaskImageType *
  GetMaskImage() const
  {
    return static_cast<const MaskImageType *>(this->ProcessObject::GetInput(1));
  }

  void
  Update() final;

  /**
   * Set/Get sigma for the gradient recursive gaussian image filter
   */
  itkSetMacro(Sigma, double);
  itkGetConstMacro(Sigma, double);

  /**
   * Set/Get boolean for gradient calculation.
   */
  itkSetMacro(UseCentralDifferenceFunction, bool);
  itkGetConstMacro(UseCentralDifferenceFunction, bool);

  /**
   * Set/Get neighborhood radius
   */
  itkSetMacro(NeighborhoodRadius, NeighborhoodRadiusType);
  itkGetConstMacro(NeighborhoodRadius, NeighborhoodRadiusType);

protected:
  ImageIntensityAndGradientToPointSetFilter();
  ~ImageIntensityAndGradientToPointSetFilter() override = default;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  GenerateData() final;

  double m_Sigma;

  NeighborhoodRadiusType m_NeighborhoodRadius;

  bool m_UseCentralDifferenceFunction;

private:
  ImageIntensityAndGradientToPointSetFilter(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented

  void
  ReadPoints();
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkImageIntensityAndGradientToPointSetFilter.hxx"
#endif

#endif // _itkImageIntensityAndGradientToPointSetFilter_h
