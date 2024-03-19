#ifndef ANTSAllocImage_h
#define ANTSAllocImage_h
#include "itkImageBase.h"
#include "itkImage.h"

/** Allocate an image based only on region */
template <typename ImageType>
typename ImageType::Pointer
AllocImage(const typename ImageType::RegionType & region)
{
  typename ImageType::Pointer rval = ImageType::New();
  rval->SetRegions(region);
  rval->AllocateInitialized();
  return rval;
}

/** Allocate an image based only on size */
template <typename ImageType>
typename ImageType::Pointer
AllocImage(const typename ImageType::SizeType & size)
{
  typename ImageType::RegionType region;
  region.SetSize(size);
  return AllocImage<ImageType>(region);
}

template <typename ImageType>
typename ImageType::Pointer
AllocImage(const typename ImageType::RegionType & region, const typename ImageType::PixelType & init)
{
  typename ImageType::Pointer rval = AllocImage<ImageType>(region);
  rval->FillBuffer(init);
  return rval;
}

/** Allocate image based on region, spaciing, origin & directions */
template <typename ImageType>
typename ImageType::Pointer
AllocImage(const typename ImageType::RegionType &    region,
           const typename ImageType::SpacingType &   spacing,
           const typename ImageType::PointType &     origin,
           const typename ImageType::DirectionType & directions)
{
  typename ImageType::Pointer rval = AllocImage<ImageType>(region);
  rval->SetSpacing(spacing);
  rval->SetOrigin(origin);
  rval->SetDirection(directions);
  return rval;
}

/** Allocate image based on region, spacing, origin & directions
 *   then initialize with initValue
 */
template <typename ImageType>
typename ImageType::Pointer
AllocImage(const typename ImageType::RegionType &    region,
           const typename ImageType::SpacingType &   spacing,
           const typename ImageType::PointType &     origin,
           const typename ImageType::DirectionType & directions,
           const typename ImageType::PixelType       initValue)
{
  typename ImageType::Pointer rval = AllocImage<ImageType>(region, spacing, origin, directions);
  rval->FillBuffer(initValue);
  return rval;
}

/** Allocate an image based on size, spacing, origin, directions */
template <typename ImageType>
typename ImageType::Pointer
AllocImage(const typename ImageType::SizeType &      size,
           const typename ImageType::SpacingType &   spacing,
           const typename ImageType::PointType &     origin,
           const typename ImageType::DirectionType & directions)
{
  typename ImageType::RegionType region;
  region.SetSize(size);
  return AllocImage<ImageType>(region, spacing, origin, directions);
}

/** Allocate image based on size,spacing,origin&directions
 *   then initialize with initial Pixel.
 */
template <typename ImageType>
typename ImageType::Pointer
AllocImage(const typename ImageType::SizeType &      size,
           const typename ImageType::SpacingType &   spacing,
           const typename ImageType::PointType &     origin,
           const typename ImageType::DirectionType & directions,
           const typename ImageType::PixelType       initValue)
{
  typename ImageType::RegionType region;
  region.SetSize(size);
  return AllocImage<ImageType>(region, spacing, origin, directions, initValue);
}

/** Allocate an image based on an exemplar image. */
template <typename ImageType>
typename ImageType::Pointer
AllocImage(const typename itk::ImageBase<ImageType::ImageDimension> * exemplar)
{
  typename ImageType::Pointer rval = ImageType::New();
  // it may be the case that the output image might have a different
  // number of PixelComponents than the exemplar, so only copy this
  // information.
  // rval->CopyInformation(exemplar);
  rval->SetLargestPossibleRegion(exemplar->GetLargestPossibleRegion());
  rval->SetBufferedRegion(exemplar->GetBufferedRegion());
  rval->SetRequestedRegion(exemplar->GetRequestedRegion());
  rval->SetSpacing(exemplar->GetSpacing());
  rval->SetOrigin(exemplar->GetOrigin());
  rval->SetDirection(exemplar->GetDirection());
  rval->AllocateInitialized();
  return rval;
}

/** Allocate an image based on an exemplar image, and then init all
 *  voxels to the suppled fillValue
 */
template <typename ImageType>
typename ImageType::Pointer
AllocImage(const typename itk::ImageBase<ImageType::ImageDimension> * exemplar,
           const typename ImageType::PixelType &                      fillValue)
{
  typename ImageType::Pointer rval = AllocImage<ImageType>(exemplar);
  rval->FillBuffer(fillValue);
  return rval;
}

#endif // ANTSAllocImage_h
