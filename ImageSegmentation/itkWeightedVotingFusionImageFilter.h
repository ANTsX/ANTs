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
#ifndef __itkWeightedVotingFusionImageFilter_h
#define __itkWeightedVotingFusionImageFilter_h

#include "itkNonLocalPatchBasedImageFilter.h"

#include "itkConstNeighborhoodIterator.h"

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include <vector>
#include <map>
#include <set>

namespace itk
{

/** \class WeightedVotingFusionImageFilter
 * \brief Implementation of the joint label fusion and joint intensity fusion algorithm
 *
 * \author Paul Yushkevich with modifications by Brian Avants and Nick Tustison
 *
 * \par REFERENCE
 *
 * H. Wang, J. W. Suh, S. Das, J. Pluta, C. Craige, P. Yushkevich,
 * "Multi-atlas segmentation with joint label fusion," IEEE Trans.
 * on Pattern Analysis and Machine Intelligence, 35(3), 611-623, 2013.
 *
 * H. Wang and P. A. Yushkevich, "Multi-atlas segmentation with joint
 * label fusion and corrective learning--an open source implementation,"
 * Front. Neuroinform., 2013.
 *
 * \ingroup ImageSegmentation
 */

template <typename TInputImage, typename TOutputImage>
class WeightedVotingFusionImageFilter final : public NonLocalPatchBasedImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef WeightedVotingFusionImageFilter                          Self;
  typedef NonLocalPatchBasedImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                                       Pointer;
  typedef SmartPointer<const Self>                                 ConstPointer;

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(WeightedVotingFusionImageFilter);

  itkNewMacro(Self);

  /** ImageDimension constants */
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;


  /** Some convenient typedefs. */
  typedef TInputImage                           InputImageType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename InputImageType::PixelType    InputImagePixelType;

  typedef typename Superclass::InputImageList    InputImageList;
  typedef typename Superclass::InputImageSetList InputImageSetList;

  typedef typename Superclass::InputImagePixelVectorType InputImagePixelVectorType;

  typedef TOutputImage                        OutputImageType;
  typedef typename OutputImageType::PixelType LabelType;
  typedef std::set<LabelType>                 LabelSetType;
  typedef Image<LabelType, ImageDimension>    LabelImageType;
  typedef typename LabelImageType::Pointer    LabelImagePointer;
  typedef std::vector<LabelImagePointer>      LabelImageList;

  typedef Image<unsigned int, ImageDimension> CountImageType;

  typedef LabelImageType                  MaskImageType;
  typedef typename MaskImageType::Pointer MaskImagePointer;

  typedef typename InputImageType::RegionType RegionType;
  typedef typename InputImageType::SizeType   SizeType;
  typedef typename InputImageType::IndexType  IndexType;

  typedef Image<float, ImageDimension>           ProbabilityImageType;
  typedef typename ProbabilityImageType::Pointer ProbabilityImagePointer;

  typedef double RealType;

  typedef vnl_matrix<RealType> MatrixType;
  typedef vnl_vector<RealType> VectorType;

  typedef std::map<LabelType, ProbabilityImagePointer> LabelPosteriorProbabilityMap;
  typedef std::map<LabelType, LabelImagePointer>       LabelExclusionMap;
  typedef std::vector<ProbabilityImagePointer>         VotingWeightImageList;

  typedef typename Superclass::ConstNeighborhoodIteratorType ConstNeighborhoodIteratorType;
  typedef typename Superclass::NeighborhoodRadiusType        NeighborhoodRadiusType;
  typedef typename Superclass::NeighborhoodOffsetType        NeighborhoodOffsetType;
  typedef typename Superclass::NeighborhoodOffsetListType    NeighborhoodOffsetListType;

  typedef typename SizeType::SizeValueType       RadiusValueType;
  typedef Image<RadiusValueType, ImageDimension> RadiusImageType;
  typedef typename RadiusImageType::Pointer      RadiusImagePointer;

  /**
   * Set the multimodal target image
   */
  void
  SetTargetImage(InputImageList imageList)
  {
    this->m_TargetImage = imageList;
    this->UpdateInputs();
  }

  /**
   * Add an atlas (multi-modal image + segmentation)
   * imageList is a vector of image pointers
   * *segmentation is the pointer for the single segmentation corresponding to imageList
   */
  void
  AddAtlas(InputImageList imageList, LabelImageType * segmentation = nullptr)
  {
    this->m_AtlasImages.push_back(imageList);

    if (this->m_NumberOfAtlasModalities == 0)
    {
      itkDebugMacro("Setting the number of modalities to " << this->m_NumberOfAtlasModalities);
      this->m_NumberOfAtlasModalities = imageList.size();
    }
    else if (this->m_NumberOfAtlasModalities != imageList.size())
    {
      itkExceptionMacro("The number of atlas multimodal images is not equal to " << this->m_NumberOfAtlasModalities);
    }
    this->m_NumberOfAtlases++;

    if (segmentation != nullptr)
    {
      this->m_AtlasSegmentations.push_back(segmentation);
      this->m_NumberOfAtlasSegmentations++;
    }

    this->UpdateInputs();
  }

  /**
   * Set mask image function.  If a binary mask image is specified, only
   * those input image voxels corresponding with mask image values equal
   * to one are used.
   */
  void
  SetMaskImage(MaskImageType * mask)
  {
    this->m_MaskImage = mask;
    this->UpdateInputs();
  }

  /**
   * Add a label exclusion map
   */
  void
  AddLabelExclusionImage(LabelType label, LabelImageType * exclusionImage)
  {
    this->m_LabelExclusionImages[label] = exclusionImage;
    this->UpdateInputs();
  }

  /**
   * Get the number of modalities used in determining the optimal label fusion
   * or optimal fused image.
   */
  itkGetConstMacro(NumberOfAtlasModalities, unsigned int);

  /**
   * Get the label set.
   */
  itkGetConstMacro(LabelSet, LabelSetType);

  /**
   * Set/Get the local search neighborhood radius image.
   */
  void
  SetNeighborhoodSearchRadiusImage(RadiusImageType * image)
  {
    this->m_NeighborhoodSearchRadiusImage = image;
  }

  /**
   * Set/Get the Alpha parameter---the regularization weight added to the matrix Mx for
   * the inverse.  Default = 0.1.
   */
  itkSetMacro(Alpha, RealType);
  itkGetConstMacro(Alpha, RealType);

  /**
   * Set/Get the Beta parameter---exponent for mapping intensity difference to joint error.
   * Default = 2.0.
   */
  itkSetMacro(Beta, RealType);
  itkGetConstMacro(Beta, RealType);

  /** Set the requested region */
  void
  GenerateInputRequestedRegion() override;

  /**
   * Boolean for retaining the posterior images. This can have a negative effect
   * on memory use, so it should only be done if one wishes to save the posterior
   * maps. The posterior maps (n = number of labels) give the probability of each
   * voxel in the target image belonging to each label.  Default = false.
   */
  itkSetMacro(RetainLabelPosteriorProbabilityImages, bool);
  itkGetConstMacro(RetainLabelPosteriorProbabilityImages, bool);
  itkBooleanMacro(RetainLabelPosteriorProbabilityImages);

  /**
   * Boolean for retaining the voting weights images.  This can have a negative effect
   * on memory use, so it should only be done if one wishes to save the voting weight
   * maps.  The voting weight maps (n = number of atlases) gives the contribution of
   * a particular atlas to the final label/intensity fusion.
   */
  itkSetMacro(RetainAtlasVotingWeightImages, bool);
  itkGetConstMacro(RetainAtlasVotingWeightImages, bool);
  itkBooleanMacro(RetainAtlasVotingWeightImages);

  /**
   * Boolean for constraining the weights to be positive and sum to 1.  We use
   * an implementation of the algorithm based on the algorithm by Lawson, Charles L.;
   * Hanson, Richard J. (1995). Solving Least Squares Problems. SIAM.
   */
  itkSetMacro(ConstrainSolutionToNonnegativeWeights, bool);
  itkGetConstMacro(ConstrainSolutionToNonnegativeWeights, bool);
  itkBooleanMacro(ConstrainSolutionToNonnegativeWeights);

  /**
   * Get the current state for progress reporting.
   */
  itkGetConstMacro(IsWeightedAveragingComplete, bool);

  /**
   * Get the posterior probability image corresponding to a label.
   */
  ProbabilityImagePointer
  GetLabelPosteriorProbabilityImage(LabelType label)
  {
    if (this->m_RetainLabelPosteriorProbabilityImages)
    {
      if (std::find(this->m_LabelSet.begin(), this->m_LabelSet.end(), label) != this->m_LabelSet.end())
      {
        return this->m_LabelPosteriorProbabilityImages[label];
      }
      else
      {
        itkDebugMacro("Not returning a label posterior probability image.  Requested label not found.");
        return nullptr;
      }
    }
    else
    {
      itkDebugMacro("Not returning a label posterior probability image.  These images were not saved.");
      return nullptr;
    }
  }

  /**
   * Get the voting weight image corresponding to an atlas.
   */
  ProbabilityImagePointer
  GetAtlasVotingWeightImage(unsigned int n)
  {
    if (this->m_RetainAtlasVotingWeightImages)
    {
      if (n < this->m_NumberOfAtlases)
      {
        return this->m_AtlasVotingWeightImages[n];
      }
      else
      {
        itkDebugMacro("Not returning a voting weight image.  Requested index is greater than the number of atlases.");
        return nullptr;
      }
    }
    else
    {
      itkDebugMacro("Not returning a voting weight image.  These images were not saved.");
      return nullptr;
    }
  }

  /**
   * Get the joint intensity fusion output image
   */
  ProbabilityImagePointer
  GetJointIntensityFusionImage(unsigned int n)
  {
    if (n < this->m_NumberOfAtlasModalities)
    {
      return this->m_JointIntensityFusionImage[n];
    }
    else
    {
      itkDebugMacro(
        "Not returning a joint intensity fusion image.  Requested index is greater than the number of modalities.");
      return nullptr;
    }
  }

protected:
  WeightedVotingFusionImageFilter();
  ~WeightedVotingFusionImageFilter() override = default;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  ThreadedGenerateData(const RegionType &, ThreadIdType) override;

  void
  BeforeThreadedGenerateData() override;

  void
  AfterThreadedGenerateData() override;

  void
  GenerateData() override;

private:
  void
  ThreadedGenerateDataForWeightedAveraging(const RegionType &, ThreadIdType);

  void
  ThreadedGenerateDataForReconstruction(const RegionType &, ThreadIdType);

  VectorType
  NonNegativeLeastSquares(const MatrixType &, const VectorType &, const RealType);

  void
  UpdateInputs();

  typedef std::pair<unsigned int, RealType> DistanceIndexType;
  typedef std::vector<DistanceIndexType>    DistanceIndexVectorType;

  struct DistanceIndexComparator
  {
    bool
    operator()(const DistanceIndexType & left, const DistanceIndexType & right)
    {
      return left.second < right.second;
    }
  };

  bool m_IsWeightedAveragingComplete;

  /** Input variables   */
  InputImageList    m_TargetImage;
  InputImageSetList m_AtlasImages;
  LabelImageList    m_AtlasSegmentations;
  LabelExclusionMap m_LabelExclusionImages;
  MaskImagePointer  m_MaskImage;

  typename CountImageType::Pointer m_CountImage;

  LabelSetType  m_LabelSet;
  SizeValueType m_NumberOfAtlases;
  SizeValueType m_NumberOfAtlasSegmentations;
  SizeValueType m_NumberOfAtlasModalities;

  std::map<RadiusValueType, NeighborhoodOffsetListType> m_NeighborhoodSearchOffsetSetsMap;

  RealType m_Alpha;
  RealType m_Beta;

  bool m_RetainLabelPosteriorProbabilityImages;
  bool m_RetainAtlasVotingWeightImages;
  bool m_ConstrainSolutionToNonnegativeWeights;

  ProbabilityImagePointer m_WeightSumImage;

  RadiusImagePointer m_NeighborhoodSearchRadiusImage;

  /** Output variables     */
  LabelPosteriorProbabilityMap m_LabelPosteriorProbabilityImages;
  VotingWeightImageList        m_AtlasVotingWeightImages;

  InputImageList m_JointIntensityFusionImage;
};

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkWeightedVotingFusionImageFilter.hxx"
#endif

#endif
