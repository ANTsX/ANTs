/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
  https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt
  for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __antsAtroposSegmentationImageFilter_h
#define __antsAtroposSegmentationImageFilter_h

#include "itkImageToImageFilter.h"

#include "antsListSampleFunction.h"
#include "antsListSampleToListSampleFilter.h"

#include "itkArray.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkFixedArray.h"
#include "itkListSample.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkNeighborhoodIterator.h"
#include "itkPointSet.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkVector.h"

#include <algorithm>
#include <vector>
#include <map>
#include <utility>

namespace itk
{
namespace ants
{
/** \class AtroposSegmentationImageFilter
 * \brief Atropos:  A Priori Classification with Registration Initialized
 *  Template Assistance
 *
 * This filter provides an Expectation-Maximization framework for statistical
 * segmentation where the intensity profile of each class is modeled as a
 * mixture model and spatial smoothness is enforced by an MRF prior.
 *
 * Initial labeling can be performed by otsu thresholding, kmeans clustering,
 * a set of user-specified prior probability images, or a prior label image.
 * If specified, the latter two initialization options are also used as
 * priors in the MRF update step.
 *
 * The assumed labeling is such that classes are assigned consecutive
 * indices 1, 2, 3, etc.  Label 0 is reserved for the background when a
 * mask is specified.
 *
 */

template <typename TInputImage,
          typename TMaskImage = Image<unsigned char, TInputImage::ImageDimension>,
          class TClassifiedImage = TMaskImage>
class AtroposSegmentationImageFilter final : public ImageToImageFilter<TInputImage, TClassifiedImage>
{
public:
  /** Standard class typdedefs. */
  typedef AtroposSegmentationImageFilter                    Self;
  typedef ImageToImageFilter<TInputImage, TClassifiedImage> Superclass;
  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(AtroposSegmentationImageFilter);

  /** Dimension of the images. */
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int ClassifiedImageDimension = TClassifiedImage::ImageDimension;
  static constexpr unsigned int MaskImageDimension = TMaskImage::ImageDimension;

  /** Typedef support of input types. */
  typedef TInputImage                   ImageType;
  typedef typename ImageType::PixelType PixelType;
  typedef typename ImageType::IndexType IndexType;
  typedef typename ImageType::SizeType  SizeType;

  typedef TMaskImage                              MaskImageType;
  typedef typename MaskImageType::PixelType       MaskLabelType;
  typedef TClassifiedImage                        ClassifiedImageType;
  typedef typename ClassifiedImageType::Pointer   ClassifiedImagePointer;
  typedef typename ClassifiedImageType::PixelType LabelType;

  /** Some convenient typedefs. */
  typedef float                           RealType;
  typedef Image<RealType, ImageDimension> RealImageType;
  typedef typename RealImageType::Pointer RealImagePointer;

  typedef FixedArray<unsigned, ImageDimension> ArrayType;
  typedef PointSet<RealType, 1>                SparseImageType;
  typedef typename SparseImageType::Pointer    SparseImagePointer;

  /** Mixture model component typedefs */
  typedef Array<RealType>                                                      MeasurementVectorType;
  typedef typename itk::Statistics::ListSample<MeasurementVectorType>          SampleType;
  typedef SmartPointer<SampleType>                                             SamplePointer;
  typedef ants::Statistics::ListSampleFunction<SampleType, RealType, RealType> LikelihoodFunctionType;
  typedef typename LikelihoodFunctionType::Pointer                             LikelihoodFunctionPointer;
  typedef typename LikelihoodFunctionType::ListSampleWeightArrayType           WeightArrayType;

  typedef std::vector<LabelType>                 PartialVolumeLabelSetType;
  typedef std::vector<PartialVolumeLabelSetType> PartialVolumeClassesType;

  /** Outlier handling typedefs */
  typedef ants::Statistics::ListSampleToListSampleFilter<SampleType, SampleType> OutlierHandlingFilterType;

  /** Randomizer typedefs */
  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator RandomizerType;
  typedef RandomizerType::IntegerType                            RandomizerSeedType;

  /** B-spline fitting typedefs */
  typedef Vector<RealType, 1>                                                      ScalarType;
  typedef Image<ScalarType, ImageDimension>                                        ScalarImageType;
  typedef PointSet<ScalarType, ImageDimension>                                     PointSetType;
  typedef BSplineScatteredDataPointSetToImageFilter<PointSetType, ScalarImageType> BSplineFilterType;
  typedef typename BSplineFilterType::PointDataImageType                           ControlPointLatticeType;
  typedef typename ControlPointLatticeType::Pointer                                ControlPointLatticePointer;
  typedef std::vector<ControlPointLatticePointer>                                  ControlPointLatticeContainerType;

  /** Initialization typedefs */
  enum InitializationStrategyType
  {
    Random,
    KMeans,
    Otsu,
    PriorProbabilityImages,
    PriorLabelImage
  };

  typedef std::pair<RealType, RealType>            LabelParametersType;
  typedef std::map<LabelType, LabelParametersType> LabelParameterMapType;
  typedef Array<RealType>                          ParametersType;

  /** Posterior probability formulation typedefs */
  enum PosteriorProbabilityFormulationType
  {
    Socrates,
    Plato,
    Aristotle,
    Sigmoid
  };

  // ivars Set/Get functionality

  /**
   * Set the number of tissue classes which is clamped from below at 2.
   * Default = 3.
   */
  itkSetClampMacro(NumberOfTissueClasses, LabelType, 2, NumericTraits<unsigned int>::max());

  /**
   * Get the number of segmentation classes.
   */
  itkGetConstMacro(NumberOfTissueClasses, unsigned int);

  /**
   * Set the partial-volume-label set one at a time.
   */
  void AddPartialVolumeLabelSet(PartialVolumeLabelSetType);

  /**
   * Get the number of partial volume classes
   */
  itkGetConstMacro(NumberOfPartialVolumeClasses, unsigned int);

  /**
   * The user can specify whether or not to use the partial volume likelihoods,
   * in which case the partial volume class is considered separate from the
   * tissue classes.  Alternatively, one can use the MRF only to handle
   * partial volume in which case, partial volume voxels are not considered
   * as separate classes.
   */
  itkSetMacro(UsePartialVolumeLikelihoods, bool);

  /**
   * The user can specify whether or not to use the partial volume likelihoods,
   * in which case the partial volume class is considered separate from the
   * tissue classes.  Alternatively, one can use the MRF only to handle
   * partial volume in which case, partial volume voxels are not considered
   * as separate classes.
   */
  itkGetConstMacro(UsePartialVolumeLikelihoods, bool);

  /**
   * The user can specify whether or not to use the partial volume likelihoods,
   * in which case the partial volume class is considered separate from the
   * tissue classes.  Alternatively, one can use the MRF only to handle
   * partial volume in which case, partial volume voxels are not considered
   * as separate classes.
   */
  itkBooleanMacro(UsePartialVolumeLikelihoods);

  /**
   * Set the maximum number of iterations.  The algorithm terminates at either
   * the maximum number of iterations or when the convergence threshold has
   * been exceeded.  Default = 5.
   */
  itkSetMacro(MaximumNumberOfIterations, unsigned int);

  /**
   * Get the maximum number of iterations.
   */
  itkGetConstMacro(MaximumNumberOfIterations, unsigned int);

  /**
   * Set the convergence threshold.  The algorithm terminates at either
   * the maximum number of iterations or when the convergence threshold has
   * been exceeded.  Default = 0.001.
   */
  itkSetMacro(ConvergenceThreshold, RealType);

  /**
   * Get the convergence threshold
   */
  itkGetConstMacro(ConvergenceThreshold, RealType);

  /**
   * Get the current convergence posterior probability for comparison with the
   * convergence threshold at each iteration.  Also is used to report progress.
   */
  itkGetConstMacro(CurrentPosteriorProbability, RealType);

  /**
   * Get the current number of iterations for comparison with the maximum
   * number of iterations.  Also is used to report progress.
   */
  itkGetConstMacro(ElapsedIterations, unsigned int);

  /**
   * Set the MRF smoothing parameter (sometimes designated as \beta in the
   * literature) which is clamped at 0.0 from below.  Greater values increase
   * cause a greater spatial homogeneity in the final labeling.  Default value
   * = 0.3.
   */
  itkSetClampMacro(MRFSmoothingFactor, RealType, NumericTraits<RealType>::ZeroValue(), NumericTraits<RealType>::max());

  /**
   * Get the MRF smoothing parameter.
   */
  itkGetConstMacro(MRFSmoothingFactor, RealType);

  /**
   * Set the MRF smoothing radius.  The radius can be set independently in each
   * dimension.  Also note that the each neighbor's contribution in calculating
   * the MRF-based prior is weighted by its distance to the center voxel.
   * Default value = 1^(ImageDimension).
   */
  itkSetMacro(MRFRadius, ArrayType);

  /**
   * Get the MRF smoothing radius.
   */
  itkGetConstMacro(MRFRadius, ArrayType);

  /**
   * Set the MRF neighborhood-defining image.
   */
  itkSetObjectMacro(MRFCoefficientImage, RealImageType);

  /**
   * Get the MRF neighborhood-defining image.
   */
  itkGetConstObjectMacro(MRFCoefficientImage, RealImageType);

  /**
   * Set the annealing temperature for ICM asynchronous updating.  For values
   * different from unity, the posterior probability is exponentiated by the
   * inverse of the annealing temperature, i.e. posterior probability \prop
   * (likelihoods * priors)^(1/T) where T is specified annealing temperature
   * raised to the number of elapsed iterations.  Default value = 1.0.
   */
  itkSetClampMacro(InitialAnnealingTemperature,
                   RealType,
                   NumericTraits<RealType>::ZeroValue(),
                   NumericTraits<RealType>::max());

  /**
   * Get the initial annealing temperature.  For values
   * different from unity, the posterior probability is exponentiated by the
   * inverse of the annealing temperature, i.e. posterior probability \prop
   * (likelihoods * priors)^(1/T) where T is specified annealing temperature
   * raised to the number of elapsed iterations.  Default value = 1.0.
   */
  itkGetConstMacro(InitialAnnealingTemperature, RealType);

  /**
   * Set the minimum annealing temperature for ICM asynchronous updating.
   * Typically, the algorithm becomes unstable for values < 0.1.  Default value
   * = 0.1.
   */
  itkSetClampMacro(MinimumAnnealingTemperature,
                   RealType,
                   NumericTraits<RealType>::ZeroValue(),
                   NumericTraits<RealType>::max());

  /**
   * Get the minimum annealing temperature for ICM asynchronous updating.
   * Typically, the algorithm becomes unstable for values < 0.1.  Default value
   * = 0.1.
   */
  itkGetConstMacro(MinimumAnnealingTemperature, RealType);

  /**
   * Set the annealing rate for ICM asynchronous updating.  For values
   * different from unity, the posterior probability is exponentiated by the
   * inverse of the annealing temperature, i.e. posterior probability \prop
   * (likelihoods * priors)^(1/T) where T is specified annealing temperature
   * raised to the number of elapsed iterations. Default value = 1.0.
   */
  itkSetClampMacro(AnnealingRate, RealType, NumericTraits<RealType>::ZeroValue(), NumericTraits<RealType>::OneValue());

  /**
   * Set the annealing rate for ICM asynchronous updating.  For values
   * different from unity, the posterior probability is exponentiated by the
   * inverse of the annealing temperature, i.e. posterior probability \prop
   * (likelihoods * priors)^(1/T) where T is specified annealing temperature
   * raised to the number of elapsed iterations. Default value = 1.0.
   */
  itkGetConstMacro(AnnealingRate, RealType);

  /**
   * Set initialization unsigned integer for random number generator.  Default is to
   * initialize randomly using the system clock.
   */
  void
  SetRandomizerInitializationSeed(const RandomizerSeedType);

  /**
   * Set initialization unsigned integer for random number generator.  Default is to
   * initialize randomly using the system clock.
   */
  itkGetConstMacro(RandomizerInitializationSeed, RandomizerSeedType);

  /**
   * Set the initialization strategy.  Initialization can occur without prior
   * information using kmeans or otsu thresholding or with prior information
   * using prior label images or prior probability images.  Default is Kmeans.
   */
  itkSetMacro(InitializationStrategy, InitializationStrategyType);

  /**
   * Get the initialization strategy.
   */
  itkGetConstMacro(InitializationStrategy, InitializationStrategyType);

  /**
   * Set the initial kmeans parameters.  For kmeans initialization, one can
   * set the initial cluster centers (for the first intensity image only).
   */
  itkSetMacro(InitialKMeansParameters, ParametersType);

  /**
   * Get the initial kmeans parameters.
   */
  itkGetConstMacro(InitialKMeansParameters, ParametersType);

  /**
   * Set the posterior probability formulation type.  This flexibility is more
   * for developmental experimentation.  Most applications will use the
   * default.  Default = Socrates.
   */
  itkSetMacro(PosteriorProbabilityFormulation, PosteriorProbabilityFormulationType);

  /**
   * Get the posterior probability formulation.
   */
  itkGetConstMacro(PosteriorProbabilityFormulation, PosteriorProbabilityFormulationType);

  /**
   * Set whether or not to use the mixture model proportions.  These proportion
   * parameters form part of the finite mixture formulation.  Default = true.
   */
  itkSetMacro(UseMixtureModelProportions, bool);

  /**
   * Get the value of boolean describing whether or not the mixture model
   * proportions should be used.  Default = true.
   */
  itkGetConstMacro(UseMixtureModelProportions, bool);

  /**
   * Set the value of the boolean parameter dictating whether or not memory
   * usage should be minimized.  Memory minimization takes more time per
   * iteration but the resulting memory footprint allows one to perform
   * problems with a large number of classes.  Default value = false.
   */
  itkSetMacro(MinimizeMemoryUsage, bool);

  /**
   * Get the value of the boolean parameter dictating whether or not memory
   * usage should be minimized.
   */
  itkGetConstMacro(MinimizeMemoryUsage, bool);

  /**
   * Set the value of the boolean parameter dictating whether or not memory
   * usage should be minimized.  Memory minimization takes more time per
   * iteration but the resulting memory footprint allows one to perform
   * problems with a large number of classes.  Default value = false.
   */
  itkBooleanMacro(MinimizeMemoryUsage);

  /**
   * Set the prior probability threshold value.  This determines what pixel
   * values are included in the sparse representation of the prior probability
   * images.  Default value =
   */
  itkSetClampMacro(ProbabilityThreshold,
                   RealType,
                   NumericTraits<RealType>::ZeroValue(),
                   NumericTraits<RealType>::OneValue());

  /**
   * Get the prior probability threshold value.
   */
  itkGetConstMacro(ProbabilityThreshold, RealType);

  // The following parameters are used for adaptive smoothing of one or more of
  // the intensity images.

  /**
   * Set the spline order of the adaptive smoothing.  Default = 3.
   */
  itkSetMacro(SplineOrder, unsigned int);

  /**
   * Get the spline order of the adaptive smoothing.
   */
  itkGetConstMacro(SplineOrder, unsigned int);

  /**
   * Set the number of fitting levels for the adaptive smoothing.  Default = 6.
   */
  itkSetMacro(NumberOfLevels, ArrayType);

  /**
   * Get the number of fitting levels for the adaptive smoothing.
   */
  itkGetConstMacro(NumberOfLevels, ArrayType);

  /**
   * Set the control point grid size for the adaptive smoothing.  Default =
   * 4^(ImageDimension) for a mesh size of 1^(ImageDimension).
   */
  itkSetMacro(NumberOfControlPoints, ArrayType);

  /**
   * Get the control point grid size.
   */
  itkGetConstMacro(NumberOfControlPoints, ArrayType);

  /**
   * Set the adaptive smoothing weight clamped between 0 and 1 which weights
   * between using just the intensity image (weight = 0) and using only the
   * full smoothed image (weight = 1).  Each intensity input image uses a
   * seperate smoothing weight value.
   */
  void
  SetAdaptiveSmoothingWeight(unsigned int idx, RealType weight)
  {
    RealType clampedWeight =
      std::min(NumericTraits<RealType>::OneValue(), std::max(NumericTraits<RealType>::ZeroValue(), weight));

    if (idx >= this->m_AdaptiveSmoothingWeights.size())
    {
      this->m_AdaptiveSmoothingWeights.resize(idx + 1);
      this->m_AdaptiveSmoothingWeights[idx] = clampedWeight;
      this->Modified();
    }
    if (!itk::Math::FloatAlmostEqual(this->m_AdaptiveSmoothingWeights[idx], weight))
    {
      this->m_AdaptiveSmoothingWeights[idx] = clampedWeight;
      this->Modified();
    }
  }

  /**
   * Get the adaptive smoothing weight for a specific intensity image.
   */
  RealType
  GetAdaptiveSmoothingWeight(unsigned int idx)
  {
    if (idx < this->m_AdaptiveSmoothingWeights.size())
    {
      return this->m_AdaptiveSmoothingWeights[idx];
    }
    else
    {
      return 0;
    }
  }

  /**
   * Set the prior label parameters.  For each class/label for label propagation
   * the boundary probabilty value is set and the exponential decay parameter.
   * The prior labeled regions are weighted linearly from the boundary to the
   * center of the region (as defined by either a Euclidean or Geodesic
   * distance) whereas outside the region, the probability value is modulated
   * by an exponential characterized by the decay parameter.
   */
  void
  SetPriorLabelParameterMap(LabelParameterMapType m)
  {
    this->m_PriorLabelParameterMap = m;
    this->Modified();
  }

  /**
   * Get the prior label parameters.
   */
  LabelParameterMapType
  GetPriorLabelParameterMap()
  {
    return this->m_PriorLabelParameterMap;
  }

  /**
   * Set the prior probability weight.  Determines what percentage of the
   * prior probability information should be included in the posterior
   * probability information.
   */
  itkSetClampMacro(PriorProbabilityWeight, RealType, NumericTraits<RealType>::ZeroValue(), static_cast<RealType>(1.e9));

  /**
   * Get the prior probability weight.
   */
  itkGetConstMacro(PriorProbabilityWeight, RealType);

  /**
   * Set a prior probability image (numbered between 1,...,numberOfClasses).
   */
  void
  SetPriorProbabilityImage(unsigned int whichClass, RealImageType * prior);

  /**
   * Get a prior probability image (numbered between 1,...,numberOfClasses).
   */
  RealImagePointer
  GetPriorProbabilityImage(unsigned int whichClass) const;

  /**
   * Set the prior label image which is assumed to have intensity values \in
   * {1,...,numberOfClasses}
   */
  void
  SetPriorLabelImage(const ClassifiedImageType * prior);

  /**
   * Get the prior label image.
   */
  const ClassifiedImageType *
  GetPriorLabelImage() const;

  /**
   * Get the number of intensity images used during the segmentation process.
   */
  itkGetConstMacro(NumberOfIntensityImages, unsigned int);

  /**
   * Set the input intensity image (numbered between 1,...,numberOfClasses)
   */
  void
  SetIntensityImage(unsigned int which, const ImageType * image);

  /**
   * Get the input intensity image (numbered between 1,...,numberOfClasses)
   */
  const ImageType *
  GetIntensityImage(unsigned int which) const;

  /**
   * Set the mask image.  The regional mask defines the domain of the
   * segmentation.
   */
  void
  SetMaskImage(const MaskImageType * mask);

  /**
   * Get the mask image.
   */
  const MaskImageType *
  GetMaskImage() const;

  /**
   * Set the label propagation type.  Euclidean distance uses the Maurer distance
   * transform to calculate the distance transform image. Otherwise the fast
   * marching filter is used to produce the geodesic distance.  The former option
   * is faster but for non-Euclidean shapes (such as the cortex), it might be
   * more accurate to use the latter option.  Default = false.
   */
  itkSetMacro(UseEuclideanDistanceForPriorLabels, bool);

  /**
   * Get the label propagation type.
   */
  itkGetConstMacro(UseEuclideanDistanceForPriorLabels, bool);

  /**
   * Set the label propagation type.  Euclidean distance uses the Maurer distance
   * transform to calculate the distance transform image. Otherwise the fast
   * marching filter is used to produce the geodesic distance.  The former option
   * is faster but for non-Euclidean shapes (such as the cortex), it might be
   * more accurate to use the latter option.  Default = false.
   */
  itkBooleanMacro(UseEuclideanDistanceForPriorLabels);

  /**
   * Set the outlier handling filter.  This takes the intensity samples from the
   * input images and modifies the sample such that the outlier effects of the
   * sample points are removed.  Default = nullptr.
   */
  itkSetObjectMacro(OutlierHandlingFilter, OutlierHandlingFilterType);

  /**
   * Get the outlier handling filter.
   */
  itkGetModifiableObjectMacro(OutlierHandlingFilter, OutlierHandlingFilterType);

  /**
   * Set the likelihood function for a specified class.  These functions are
   * traditionally characterized as parametric, i.e. Gaussian, or nonparametric.
   * A likelihood function must be assigned for each class.
   */
  void
  SetLikelihoodFunction(unsigned int n, LikelihoodFunctionType * prob)
  {
    if (n < this->m_MixtureModelComponents.size() && this->m_MixtureModelComponents[n] != prob)
    {
      this->m_MixtureModelComponents[n] = prob;
      this->Modified();
    }
    else if (n >= this->m_MixtureModelComponents.size())
    {
      this->m_MixtureModelComponents.resize(n + 1);
      this->m_MixtureModelComponents[n] = prob;
      this->Modified();
    }
  }

  /**
   * Get the likelihood function for a specified class.
   */
  LikelihoodFunctionType *
  GetLikelihoodFunction(unsigned int n)
  {
    if (n < this->m_MixtureModelComponents.size())
    {
      return this->m_MixtureModelComponents[n].GetPointer();
    }
    else
    {
      return nullptr;
    }
  }

  /**
   * Get the likelihood image for a specified class.  Note that this function
   * facilitates looking at the likelihood image for the user but is not used
   * internally during the optimization of the segmentation solution.
   */
  RealImagePointer
  GetLikelihoodImage(unsigned int);

  /**
   * Get the posterior probability image.  This function provides the soft
   * classification results.
   */
  RealImagePointer
  GetPosteriorProbabilityImage(unsigned int);

  /**
   * Get the smooth intensity image.  Available when adaptive smoothing is
   * enabled.
   */
  RealImagePointer
  GetSmoothIntensityImageFromPriorImage(unsigned int, unsigned int);

  /**
   * Get the distance prior probability image.  Available when prior images are
   * used.
   */
  RealImagePointer
  GetDistancePriorProbabilityImage(unsigned int);

  /**
   * Boolean variable governing the update scheme.  If set to true (default), an
   * asynchronous approach to updating the class labels is performed which
   * has theoretical convergence properties.
   */
  itkBooleanMacro(UseAsynchronousUpdating);

  /**
   * Boolean variable governing the update scheme.  If set to true (default), an
   * asynchronous approach to updating the class labels is performed which
   * has theoretical convergence properties.
   */
  itkSetMacro(UseAsynchronousUpdating, bool);

  /**
   * Boolean variable governing the update scheme.  If set to true (default), an
   * asynchronous approach to updating the class labels is performed which
   * has theoretical convergence properties.
   */
  itkGetConstMacro(UseAsynchronousUpdating, bool);

  /**
   * Set the number of maximum allowed ICM iterations.  When asynchronous
   * updating is used, at each iteration an ICM optimization step occurs in
   * which the posterior is maximized.  During the ICM optimization, monotonic
   * increase of the posterior is guaranteed during the iterations forming
   * the optimization which usually maximizes out between 5-10 iterations.
   */
  itkSetMacro(MaximumNumberOfICMIterations, unsigned int);

  /**
   * Get the number of maximum allowed ICM iterations.  When asynchronous
   * updating is used, at each iteration an ICM optimization step occurs in
   * which the posterior is maximized.  During the ICM optimization, monotonic
   * increase of the posterior is guaranteed during the iterations forming
   * the optimization which usually maximizes out between 5-10 iterations.
   */
  itkGetConstMacro(MaximumNumberOfICMIterations, unsigned int);

  /**
   * Get the ICM code image.
   */
  ClassifiedImagePointer
  GetICMCodeImage()
  {
    return this->m_ICMCodeImage;
  };

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(SameDimensionCheck1, (Concept::SameDimension<ImageDimension, ClassifiedImageDimension>));
  itkConceptMacro(SameDimensionCheck2, (Concept::SameDimension<ImageDimension, MaskImageDimension>));
  /** End concept checking */
#endif
protected:
  AtroposSegmentationImageFilter();
  ~AtroposSegmentationImageFilter() override;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  GenerateData() override;

private:
  AtroposSegmentationImageFilter(const Self &) = delete;
  void
  operator=(const Self &) = delete;

  /**
   * Initialize the segmentation labeling.
   */
  void
  GenerateInitialClassLabeling();

  /**
   * Initialize labeling using otsu thresholding on the first input image.
   */
  void
  GenerateInitialClassLabelingWithOtsuThresholding();

  /**
   * Initialize labeling using kmeans classification.
   */
  void
  GenerateInitialClassLabelingWithKMeansClustering();

  /**
   * Initialize labeling using prior probability images.
   */
  void
  GenerateInitialClassLabelingWithPriorProbabilityImages();

  /**
   * Update the class labeling at each iteration using asynchronous ICM updating.
   * and return the max posterior probability.
   */
  RealType
  UpdateClassLabeling();

  /**
   * Compute the ICM code image for asynchronous updating to ensure that the
   * the local MRF neighborhoods are updated independently.  For more information
   * see notes for the bool variable m_UseAsynchronousUpdating.
   */
  void
  ComputeICMCodeImage();

  /**
   * This function returns a set of samples for each class such that each
   * measurement vector of the returned SampleType corresponds to a single
   * voxel across the set of auxiliary and input images.
   */
  typename SampleType::Pointer
  GetScalarSamples();

  /**
   * Calculate the local posterior probability.
   */
  RealType
  CalculateLocalPosteriorProbability(RealType, RealType, RealType, RealType, RealType, IndexType, unsigned int);

  void
  EvaluateMRFNeighborhoodWeights(ConstNeighborhoodIterator<ClassifiedImageType>, Array<RealType> &);

  RealType PerformLocalLabelingUpdate(NeighborhoodIterator<ClassifiedImageType>);

  // ivars

  unsigned int             m_NumberOfTissueClasses;
  unsigned int             m_NumberOfPartialVolumeClasses;
  PartialVolumeClassesType m_PartialVolumeClasses;
  bool                     m_UsePartialVolumeLikelihoods;

  unsigned int m_NumberOfIntensityImages;

  unsigned int m_ElapsedIterations;
  unsigned int m_MaximumNumberOfIterations;
  RealType     m_CurrentPosteriorProbability;
  RealType     m_ConvergenceThreshold;

  std::vector<LikelihoodFunctionPointer> m_MixtureModelComponents;
  Array<RealType>                        m_MixtureModelProportions;

  InitializationStrategyType m_InitializationStrategy;
  ParametersType             m_InitialKMeansParameters;

  PosteriorProbabilityFormulationType m_PosteriorProbabilityFormulation;
  bool                                m_UseMixtureModelProportions;
  RealType                            m_InitialAnnealingTemperature;
  RealType                            m_MinimumAnnealingTemperature;
  RealType                            m_AnnealingRate;

  typename OutlierHandlingFilterType::Pointer m_OutlierHandlingFilter;

  typename RandomizerType::Pointer m_Randomizer;

  ArrayType        m_MRFRadius;
  RealType         m_MRFSmoothingFactor;
  RealImagePointer m_MRFCoefficientImage;

  typename ClassifiedImageType::SpacingType m_ImageSpacing;

  unsigned int           m_MaximumICMCode;
  ClassifiedImagePointer m_ICMCodeImage;
  bool                   m_UseAsynchronousUpdating;
  unsigned int           m_MaximumNumberOfICMIterations;
  RandomizerSeedType     m_RandomizerInitializationSeed;

  std::vector<RealType>           m_AdaptiveSmoothingWeights;
  RealType                        m_PriorProbabilityWeight;
  LabelParameterMapType           m_PriorLabelParameterMap;
  RealType                        m_ProbabilityThreshold;
  std::vector<RealImagePointer>   m_PriorProbabilityImages;
  std::vector<SparseImagePointer> m_PriorProbabilitySparseImages;

  unsigned int                                  m_SplineOrder;
  ArrayType                                     m_NumberOfLevels;
  ArrayType                                     m_NumberOfControlPoints;
  std::vector<ControlPointLatticeContainerType> m_ControlPointLattices;

  RealImagePointer m_SumDistancePriorProbabilityImage;
  RealImagePointer m_SumPosteriorProbabilityImage;
  bool             m_MinimizeMemoryUsage;

  bool                          m_UseEuclideanDistanceForPriorLabels;
  std::vector<RealImagePointer> m_DistancePriorProbabilityImages;
  std::vector<RealImagePointer> m_PosteriorProbabilityImages;

  itk::Array<unsigned long> m_LabelVolumes;

  std::vector<const ImageType *> m_IntensityImages;

  typename ClassifiedImageType::ConstPointer m_PriorLabelImage;
  typename MaskImageType::ConstPointer       m_MaskImage;

  // inline functions to help with the sparse image creation

  inline typename RealImageType::IndexType
  NumberToIndex(unsigned long number, const SizeType size) const
  {
    IndexType k;

    k[0] = 1;
    for (unsigned int i = 1; i < ImageDimension; i++)
    {
      k[i] = size[i - 1] * k[i - 1];
    }
    IndexType index;
    for (unsigned int i = 0; i < ImageDimension; i++)
    {
      index[ImageDimension - i - 1] = static_cast<unsigned long>(number / k[ImageDimension - i - 1]);
      number %= k[ImageDimension - i - 1];
    }
    return index;
  }

  inline unsigned long
  IndexToNumber(const IndexType k, const SizeType size) const
  {
    unsigned long number = k[0];

    for (unsigned int i = 1; i < ImageDimension; i++)
    {
      unsigned long s = 1;
      for (unsigned int j = 0; j < i; j++)
      {
        s *= size[j];
      }
      number += s * k[i];
    }
    return number;
  }
};
} // namespace ants
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "antsAtroposSegmentationImageFilter.hxx"
#endif

#endif
