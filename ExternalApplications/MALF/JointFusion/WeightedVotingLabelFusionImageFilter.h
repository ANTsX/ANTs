#ifndef __WeightedVotingLabelFusionImageFilter_h_
#define __WeightedVotingLabelFusionImageFilter_h_

#include "itkImageToImageFilter.h"
#include "itkConstNeighborhoodIterator.h"

template <class TInputImage, class TOutputImage>
class WeightedVotingLabelFusionImageFilter : public itk::ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef WeightedVotingLabelFusionImageFilter               Self;
  typedef itk::ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef itk::SmartPointer<Self>                            Pointer;
  typedef itk::SmartPointer<const Self>                      ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageToImageFilter, ImageSource);

  itkNewMacro(Self);

  /** Superclass typedefs. */
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
  typedef typename Superclass::OutputImagePixelType  OutputImagePixelType;

  /** Some convenient typedefs. */
  typedef TInputImage                           InputImageType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename InputImageType::PixelType    InputImagePixelType;

  typedef typename InputImageType::RegionType RegionType;
  typedef typename InputImageType::SizeType   SizeType;
  typedef typename InputImageType::IndexType  IndexType;

  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Set target image */
  void SetTargetImage(InputImageType *image)
  {
    m_Target = image; UpdateInputs();
  }

  /** Add an atlas */
  void AddAtlas(InputImageType *grey, InputImageType *seg)
  {
    m_Atlases.push_back(grey);
    m_AtlasSegs.push_back(seg);
    UpdateInputs();
  }

  void AddExclusionMap(InputImagePixelType label, InputImageType *excl)
  {
    m_Exclusions[label] = excl;
    UpdateInputs();
  }

  /** Set the parameters */
  itkSetMacro(SearchRadius, SizeType);
  itkGetMacro(SearchRadius, SizeType);

  itkSetMacro(PatchRadius, SizeType);
  itkGetMacro(PatchRadius, SizeType);

  itkSetMacro(Alpha, double);
  itkGetMacro(Alpha, double);

  itkSetMacro(Beta, double);
  itkGetMacro(Beta, double);

  /** Set the requested region */
  void GenerateInputRequestedRegion();

  /**
   * Whether the posterior maps should be retained. This can have a negative effect
   * on memory use, so it should only be done if one wishes to save the posterior
   * maps. The posterior maps given the probability of each voxel in the target image
   * belonging to each label.
   */
  itkSetMacro(RetainPosteriorMaps, bool)
  itkGetMacro(RetainPosteriorMaps, bool)

  typedef itk::Image<float, InputImageDimension>                    PosteriorImage;
  typedef typename PosteriorImage::Pointer                          PosteriorImagePtr;
  typedef typename std::map<InputImagePixelType, PosteriorImagePtr> PosteriorMap;

  /**
   * Get the posterior maps (if they have been retained)
   */
  const PosteriorMap & GetPosteriorMaps()
  {
    return m_PosteriorMap;
  }

  void GenerateData();

protected:

  WeightedVotingLabelFusionImageFilter()
  {
    m_Alpha = 0.01;
    m_Beta = 2;
    m_RetainPosteriorMaps = false;
  }

  ~WeightedVotingLabelFusionImageFilter()
  {
  }

private:

  typedef itk::Neighborhood<InputImagePixelType, InputImageDimension> HoodType;
  typedef itk::ConstNeighborhoodIterator<InputImageType>              NIter;

  double PatchSimilarity(const InputImagePixelType * const psearch, const InputImagePixelType * const normtrg,
                         const size_t n, int * offsets, InputImagePixelType & sum_psearch,
                         InputImagePixelType & ssq_psearch);

  void ComputeOffsetTable(const InputImageType * const image, const SizeType & radius, int * *offset, size_t & nPatch,
                          int * *manhattan = NULL);

  void UpdateInputs();

  void PatchStats(const InputImagePixelType * p, size_t n, int * offsets,
                  InputImagePixelType & mean, InputImagePixelType & std);

  double JointErrorEstimate(const InputImagePixelType * t, const InputImagePixelType * a1,
                            const InputImagePixelType * a2, size_t n, int *offsets);

  SizeType m_SearchRadius;
  SizeType m_PatchRadius;

  double m_Alpha;
  double m_Beta;

  typedef std::vector<InputImagePointer>                   InputImageList;
  typedef std::map<InputImagePixelType, InputImagePointer> ExclusionMap;

  // Posterior maps
  PosteriorMap m_PosteriorMap;

  // Whether they are retained
  bool m_RetainPosteriorMaps;

  // Organized lists of inputs
  InputImagePointer m_Target;
  InputImageList    m_AtlasSegs;
  InputImageList    m_Atlases;
  ExclusionMap      m_Exclusions;
};

#endif
