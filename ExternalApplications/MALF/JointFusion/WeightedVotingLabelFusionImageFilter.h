#ifndef __WeightedVotingLabelFusionImageFilter_h_
#define __WeightedVotingLabelFusionImageFilter_h_

#include "itkImageToImageFilter.h"
#include "itkConstNeighborhoodIterator.h"

template <class TInputImage, class TOutputImage>
class WeightedVotingLabelFusionImageFilter : public itk::ImageToImageFilter <TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef WeightedVotingLabelFusionImageFilter Self;
  typedef itk::ImageToImageFilter<TInputImage, TOutputImage>  Superclass;
  typedef itk::SmartPointer<Self>         Pointer;
  typedef itk::SmartPointer<const Self>   ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageToImageFilter,ImageSource);

  itkNewMacro(Self);

  /** Superclass typedefs. */
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
  typedef typename Superclass::OutputImagePixelType  OutputImagePixelType;

  /** Some convenient typedefs. */
  typedef TInputImage                             InputImageType;
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename InputImageType::PixelType      InputImagePixelType;

  typedef typename InputImageType::RegionType     RegionType;
  typedef typename InputImageType::SizeType       SizeType;
  typedef typename InputImageType::IndexType      IndexType;

  typedef std::vector<int>  LabelSetType;
  typedef std::vector<InputImagePointer> InputImageList;
  typedef std::vector<InputImageList> InputMultiModList;


  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
                                                                               
  /** Set target image */
  void SetTargetImage(InputImageList image)
    { m_Target = image; UpdateInputs(); }

//  void SetTargetImage(InputImageType *image)
//    { m_Target = image; UpdateInputs(); }

  /** Add an atlas */
  void AddAtlas(InputImageList grey, InputImageType *seg)
    {
    for (size_t i=0;i<grey.size();i++){
      m_Atlases.push_back(grey[i]);
    }
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

  itkSetMacro(Modality, int);
  itkGetMacro(Modality, int);

  itkSetMacro(PatchRadius, SizeType);
  itkGetMacro(PatchRadius, SizeType);

  itkSetMacro(Alpha, double);
  itkGetMacro(Alpha, double);

  itkSetMacro(Beta, double);
  itkGetMacro(Beta, double);

  //itkSetMacro(GroupID, std::vector<int>);
  virtual void SetGroupID(const std::vector<int> _arg)
    {
      this->m_GroupID = _arg;
      this->Modified();
    }
  itkGetMacro(GroupID, std::vector<int>);

  // itkSetMacro(GroupWeight, std::vector<double>);
  virtual void SetGroupWeight(const std::vector<double> _arg)
    {
      this->m_GroupWeight = _arg;
      this->Modified();
    }
  itkGetMacro(GroupWeight, std::vector<double>);


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

  itkSetMacro(RetainVotingWeight, bool)
  itkGetMacro(RetainVotingWeight, bool)


  typedef itk::Image<float, InputImageDimension> PosteriorImage;
  typedef typename PosteriorImage::Pointer PosteriorImagePtr;
  typedef typename std::map<InputImagePixelType, PosteriorImagePtr> PosteriorMap;
                                                                    
  /**
   * Get the posterior maps (if they have been retained)
   */
  const PosteriorMap &GetPosteriorMaps()
    { return m_PosteriorMap; }


  /**
   * Get the voting weights (if they have been retained)
   */
  const PosteriorMap &GetVotingWeight()
    { return m_VotingWeight; }



  void GenerateData();
 
protected:

  WeightedVotingLabelFusionImageFilter() 
    { 
    m_Alpha=0.1; 
    m_Beta=2; 
    m_RetainPosteriorMaps = false;
    m_RetainVotingWeight = false;
    }
  ~WeightedVotingLabelFusionImageFilter() {}

private:

  typedef itk::Neighborhood<InputImagePixelType, InputImageDimension> HoodType;
  typedef itk::ConstNeighborhoodIterator<InputImageType> NIter;

  double PatchSimilarity(
    const InputImagePixelType *psearch, const InputImagePixelType *pnormtrg, int offset,
    size_t n, int *offsets, InputImagePixelType &psearchSum, InputImagePixelType &psearchSSQ);

  void ComputeOffsetTable(
    const InputImageType *image, const SizeType &radius, 
    int **offset, size_t &nPatch, int **manhattan = NULL);

  void UpdateInputs();

  void PatchStats(const InputImagePixelType *p, size_t n, int *offsets, InputImagePixelType &mean, InputImagePixelType &sd);

  double JointErrorEstimate(const InputImagePixelType *t, const InputImagePixelType *a1, const InputImagePixelType *a2, size_t n, int *offsets);

  SizeType m_SearchRadius, m_PatchRadius;
  int m_Modality;
  double m_Alpha, m_Beta;

  typedef std::map<InputImagePixelType, InputImagePointer> ExclusionMap;
  
  typedef std::vector<int> ArrayInt;

  ArrayInt m_GroupID;
  std::vector<double> m_GroupWeight;


  // Posterior maps
  PosteriorMap m_PosteriorMap;
  PosteriorMap m_VotingWeight;

  // Whether they are retained
  bool m_RetainPosteriorMaps;
  bool m_RetainVotingWeight;

  // Organized lists of inputs
  InputImageList m_Target;
  InputImageList m_Atlases;
  InputImageList m_AtlasSegs;
  ExclusionMap m_Exclusions;

};


#endif
