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
#ifndef itkNonLocalPatchBasedImageFilter_h
#define itkNonLocalPatchBasedImageFilter_h

#include "itkImageToImageFilter.h"

#include "itkConstNeighborhoodIterator.h"

namespace itk {

/**
 * \class NonLocalPatchBasedImageFilter
 * \brief Implementation of a non-local upsampling (i.e., superresolution) image filter.
 *
 * \ingroup ITKFiltering
 */

template<typename TInputImage, class TOutputImage = TInputImage>
class NonLocalPatchBasedImageFilter :
  public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef NonLocalPatchBasedImageFilter                 Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Runtime information support. */
  itkTypeMacro( NonLocalPatchBasedImageFilter, ImageToImageFilter );

  /** Standard New method. */
  itkNewMacro( Self );

  /** ImageDimension constants */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TInputImage::ImageDimension );

  /** Some convenient typedefs. */
  typedef TInputImage                                    InputImageType;
  typedef typename InputImageType::PixelType             InputPixelType;
  typedef typename InputImageType::Pointer               InputImagePointer;
  typedef std::vector<InputImagePointer>                 InputImageList;
  typedef std::vector<InputImageList>                    InputImageSetList;
  typedef typename InputImageType::RegionType            RegionType;

  typedef TOutputImage                                   OutputImageType;
  typedef typename OutputImageType::PixelType            OutputPixelType;

  typedef std::vector<InputPixelType>                    InputImagePixelVectorType;

  typedef float                                          RealType;
  typedef Image<RealType, ImageDimension>                RealImageType;
  typedef typename RealImageType::Pointer                RealImagePointer;
  typedef typename RealImageType::IndexType              IndexType;

  typedef Neighborhood<InputPixelType, ImageDimension>         NeighborhoodType;
  typedef SizeValueType                                        NeighborhoodSizeType;

  typedef ConstNeighborhoodIterator<InputImageType>            ConstNeighborhoodIteratorType;
  typedef typename ConstNeighborhoodIteratorType::RadiusType   NeighborhoodRadiusType;
  typedef typename ConstNeighborhoodIteratorType::OffsetType   NeighborhoodOffsetType;

  typedef std::vector<NeighborhoodOffsetType>                  NeighborhoodOffsetListType;
  /**
   * Neighborhood patch similarity metric enumerated type
   */
  enum SimilarityMetricType {
    PEARSON_CORRELATION,
    MEAN_SQUARES
  };

  /**
   * Get/set neighborhood search radius.
   * Default = 3x3x...
   */
  itkSetMacro( NeighborhoodSearchRadius, NeighborhoodRadiusType );
  itkGetConstMacro( NeighborhoodSearchRadius, NeighborhoodRadiusType );

  /**
   * Get/set neighborhood search size.
   */
  itkSetMacro( NeighborhoodSearchSize, NeighborhoodSizeType );
  itkGetConstMacro( NeighborhoodSearchSize, NeighborhoodSizeType );

  /**
   * Get/set neighborhood search offset list.
   */
  virtual void SetNeighborhoodSearchOffsetList( const NeighborhoodOffsetListType list )
    {
    this->m_NeighborhoodSearchOffsetList = list;
    this->Modified();
    }
  itkGetConstMacro( NeighborhoodSearchOffsetList, NeighborhoodOffsetListType );

  /**
   * Get/set neighborhood patch radius.
   * Default = 1x1x...
   */
  itkSetMacro( NeighborhoodPatchRadius, NeighborhoodRadiusType );
  itkGetConstMacro( NeighborhoodPatchRadius, NeighborhoodRadiusType );

  /**
   * Get/set neighborhood patch size.
   */
  itkSetMacro( NeighborhoodPatchSize, NeighborhoodSizeType );
  itkGetConstMacro( NeighborhoodPatchSize, NeighborhoodSizeType );

  /**
   * Get/set neighborhood patch offset list.
   */
  virtual void SetNeighborhoodPatchOffsetList( const NeighborhoodOffsetListType list )
    {
    this->m_NeighborhoodPatchOffsetList = list;
    this->Modified();
    }
  itkGetConstMacro( NeighborhoodPatchOffsetList, NeighborhoodOffsetListType );

  /**
   * Enumerated type for neighborhood similarity.  Default = MEAN_SQUARES
   */
  itkSetMacro( SimilarityMetric, SimilarityMetricType );
  itkGetConstMacro( SimilarityMetric, SimilarityMetricType );

protected:

  NonLocalPatchBasedImageFilter();
  ~NonLocalPatchBasedImageFilter() {}

  void BeforeThreadedGenerateData() ITK_OVERRIDE;

  void PrintSelf( std::ostream & os, Indent indent ) const ITK_OVERRIDE;

  RealType ComputeNeighborhoodPatchSimilarity( const InputImageList &, const IndexType, const InputImagePixelVectorType &, const bool );

  InputImagePixelVectorType VectorizeImageListPatch( const InputImageList &, const IndexType, const bool );

  InputImagePixelVectorType VectorizeImagePatch( const InputImagePointer, const IndexType, const bool );

  void GetMeanAndStandardDeviationOfVectorizedImagePatch( const InputImagePixelVectorType &, RealType &, RealType & );

  itkSetMacro( TargetImageRegion, RegionType );
  itkGetConstMacro( TargetImageRegion, RegionType );

  SimilarityMetricType                                 m_SimilarityMetric;

  SizeValueType                                        m_NeighborhoodSearchSize;
  NeighborhoodRadiusType                               m_NeighborhoodSearchRadius;
  NeighborhoodOffsetListType                           m_NeighborhoodSearchOffsetList;

  SizeValueType                                        m_NeighborhoodPatchSize;
  NeighborhoodRadiusType                               m_NeighborhoodPatchRadius;
  NeighborhoodOffsetListType                           m_NeighborhoodPatchOffsetList;

  RegionType                                           m_TargetImageRegion;


private:

  NonLocalPatchBasedImageFilter( const Self& ) ITK_DELETED_FUNCTION;
  void operator=( const Self& ) ITK_DELETED_FUNCTION;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkNonLocalPatchBasedImageFilter.hxx"
#endif

#endif
