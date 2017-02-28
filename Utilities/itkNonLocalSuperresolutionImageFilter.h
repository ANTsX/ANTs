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
#ifndef itkNonLocalSuperresolutionImageFilter_h
#define itkNonLocalSuperresolutionImageFilter_h

#include "itkImageToImageFilter.h"

#include "itkConstNeighborhoodIterator.h"

namespace itk {

/**
 * \class NonLocalSuperresolutionImageFilter
 * \brief Implementation of a non-local upsampling (i.e., superresolution) image filter.

 *
 * \author Jose V. Manjon with ITK porting by Nick Tustison
 *
 * Contributed by
 *
 * \par REFERENCE
 *
 * José V. Manjón, Pierrick Coupe, Antonio Buades, Vladimir Fonov, Louis Collins and
 * Montserrat Robles.
 * "Non-local MRI Upsampling",
 * Medical Image Analysis, 14:784-792, 2010.
 *
 * José V. Manjón, Pierrick Coupe, Antonio Buades, Louis Collins and Montserrat Robles.
 * "MRI Superresolution Using Self-Similarity and Image Priors",
 * International Journal of Biomedical Imaging, 2010.
 *
 * \ingroup ITKFiltering
 */

template<typename TInputImage, class TOutputImage = TInputImage>
class NonLocalSuperresolutionImageFilter :
  public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef NonLocalSuperresolutionImageFilter            Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Runtime information support. */
  itkTypeMacro( NonLocalSuperresolutionImageFilter, ImageToImageFilter );

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

  typedef std::vector<RealType>                          ScaleLevelsArrayType;

  typedef ConstNeighborhoodIterator<InputImageType>            ConstNeighborhoodIteratorType;
  typedef typename ConstNeighborhoodIteratorType::RadiusType   NeighborhoodRadiusType;
  typedef typename ConstNeighborhoodIteratorType::OffsetType   NeighborhoodOffsetType;

  /**
   * Neighborhood patch similarity metric enumerated type
   */
  enum SimilarityMetricType {
    PEARSON_CORRELATION,
    MEAN_SQUARES
  };

  /**
   * Set the low resolution input image to be refined.
   */
  void SetLowResolutionInputImage( const InputImageType *image )
    {
    this->SetInput( const_cast<InputImageType *>( image ) );
    }
  void SetInput1( const InputImageType *image )
    {
    this->SetLowResolutionInputImage( image );
    }

  /**
   * Get the low resolution input image to be upsampled.
   */
  const InputImageType* GetLowResolutionInputImage() const
    {
    return static_cast<const InputImageType*>( this->ProcessObject::GetInput( 0 ) );
    }

  /**
   * Set the high resolution input image to be used as the reference input.
   */
  void SetHighResolutionReferenceImage( const InputImageType *image )
    {
    this->SetNthInput( 1, const_cast<InputImageType *>( image ) );
    }
  void SetInput2( const InputImageType *image )
    {
    this->SetHighResolutionReferenceImage( image );
    }

  /**
   * Get the high resolution input image to be used as the reference input.
   */
  const InputImageType* GetHighResolutionReferenceImage() const
    {
    return static_cast<const InputImageType*>( this->ProcessObject::GetInput( 1 ) );
    }

  /**
   * Type of the Interpolator Base class
   */
  typedef InterpolateImageFunction<InputImageType, RealType> InterpolatorType;

  /**
   * Set the interpolator.
   */
  itkSetObjectMacro( Interpolator, InterpolatorType );

  /**
   * Get the interpolator.
   */
  itkGetModifiableObjectMacro( Interpolator, InterpolatorType );

  /**
   * Set/get the interpolator.
   */
  itkSetMacro( IntensityDifferenceSigma, RealType );
  itkGetConstMacro( IntensityDifferenceSigma, RealType );

  /**
   * Get the interpolator.
   */
  itkSetMacro( PatchSimilaritySigma, RealType );
  itkGetConstMacro( PatchSimilaritySigma, RealType );

  /**
   * Neighborhood search radius.
   * Default = 3x3x...
   */
  itkSetMacro( NeighborhoodSearchRadius, NeighborhoodRadiusType );
  itkGetConstMacro( NeighborhoodSearchRadius, NeighborhoodRadiusType );

  /**
   * Neighborhood block radius.
   * Default = 1x1x...
   */
  itkSetMacro( NeighborhoodPatchRadius, NeighborhoodRadiusType );
  itkGetConstMacro( NeighborhoodPatchRadius, NeighborhoodRadiusType );

  /**
   * Enumerated type for neighborhood similarity.  Default = MEAN_SQUARES
   */
  itkSetMacro( SimilarityMetric, SimilarityMetricType );
  itkGetConstMacro( SimilarityMetric, SimilarityMetricType );

  /**
   * Scale levels for .
   */
  itkSetMacro( ScaleLevels, ScaleLevelsArrayType );
  itkGetConstMacro( ScaleLevels, ScaleLevelsArrayType );

  /**
   * Get the number of elapsed iterations.  This is a helper function for
   * reporting observations.
   */
  itkGetConstMacro( CurrentIteration, SizeValueType );

protected:
  NonLocalSuperresolutionImageFilter();
  ~NonLocalSuperresolutionImageFilter() {}

  void PrintSelf( std::ostream & os, Indent indent ) const ITK_OVERRIDE;

  void GenerateData() ITK_OVERRIDE;

  void ThreadedGenerateData( const RegionType &, ThreadIdType ) ITK_OVERRIDE;

  void BeforeThreadedGenerateData() ITK_OVERRIDE;

  void AfterThreadedGenerateData() ITK_OVERRIDE;

  void AllocateOutputs() ITK_OVERRIDE;

  void VerifyInputInformation() ITK_OVERRIDE;

  void GenerateOutputInformation() ITK_OVERRIDE;

  void GenerateInputRequestedRegion() ITK_OVERRIDE;

private:

  NonLocalSuperresolutionImageFilter( const Self& ) ITK_DELETE_FUNCTION;
  void operator=( const Self& ) ITK_DELETE_FUNCTION;

  void PerformMeanCorrection();

  RealType ComputeNeighborhoodPatchSimilarity( const InputImageList &, const IndexType, const InputImagePixelVectorType &, const bool );

  InputImagePixelVectorType VectorizeImageListPatch( const InputImageList &, const IndexType, const bool );

  InputImagePixelVectorType VectorizeImagePatch( const InputImagePointer, const IndexType, const bool );

  void GetMeanAndStandardDeviationOfVectorizedImagePatch( const InputImagePixelVectorType &, RealType &, RealType & );

  typedef std::pair<unsigned int, RealType>           DistanceIndexType;
  typedef std::vector<DistanceIndexType>              DistanceIndexVectorType;

  struct DistanceIndexComparator
    {
    bool operator () ( const DistanceIndexType& left, const DistanceIndexType& right )
      {
      return left.second < right.second;
      }
    };

  RealType                                             m_Epsilon;

  SimilarityMetricType                                 m_SimilarityMetric;

  RealImagePointer                                     m_WeightSumImage;

  RegionType                                           m_TargetImageRequestedRegion;

  typename InputImageType::Pointer                     m_InterpolatedLowResolutionInputImage;

  typename InterpolatorType::Pointer                   m_Interpolator;

  SizeValueType                                        m_NeighborhoodSearchSize;
  NeighborhoodRadiusType                               m_NeighborhoodSearchRadius;
  std::vector<NeighborhoodOffsetType>                  m_NeighborhoodSearchOffsetList;

  NeighborhoodRadiusType                               m_NeighborhoodPatchRadius;
  SizeValueType                                        m_NeighborhoodPatchSize;
  std::vector<NeighborhoodOffsetType>                  m_NeighborhoodPatchOffsetList;

  RealType                                             m_PatchSimilaritySigma;
  RealType                                             m_IntensityDifferenceSigma;

  ScaleLevelsArrayType                                 m_ScaleLevels;
  SizeValueType                                        m_CurrentIteration;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkNonLocalSuperresolutionImageFilter.hxx"
#endif

#endif
