// ===========================================================================
//
//                            PUBLIC DOMAIN NOTICE
//            Office of High Performance Computing and Communications
//
//  This software/database is a "United States Government Work" under the
//  terms of the United States Copyright Act.  It was written as part of
//  the author's official duties as a United States Government employee and
//  thus cannot be copyrighted.  This software is freely available
//  to the public for use. The National Library of Medicine and the U.S.
//  Government have not placed any restriction on its use or reproduction.
//
//  Although all reasonable efforts have been taken to ensure the accuracy
//  and reliability of the software and data, the NLM and the U.S.
//  Government do not and cannot warrant the performance or results that
//  may be obtained by using this software or data. The NLM and the U.S.
//  Government disclaim all warranties, express or implied, including
//  warranties of performance, merchantability or fitness for any particular
//  purpose.
//
//  Please cite the author in any work or product based on this material.
//
// ===========================================================================
//

#ifndef __itkMultiScaleLaplacianBlobDetectorImageFilter_h
#define __itkMultiScaleLaplacianBlobDetectorImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkGaussianSpatialObject.h"

namespace itk
{
/** \class ScaleSpaceBlobSpatialObject
 * \brief a spatial object to represent a gaussian blob with size
 *
 * The superclass parameter "Sigma" is the size of the blob if it was
 * a gaussian.
 * \ingroup ITKImageScaleSpace
 **/
template <unsigned int TDimension = 3>
class ScaleSpaceBlobSpatialObject final : public GaussianSpatialObject<TDimension>
{
public:
  typedef ScaleSpaceBlobSpatialObject       Self;
  typedef double                            ScalarType;
  typedef SmartPointer<Self>                Pointer;
  typedef SmartPointer<const Self>          ConstPointer;
  typedef GaussianSpatialObject<TDimension> Superclass;
  typedef SmartPointer<Superclass>          SuperclassPointer;
  typedef typename Superclass::PointType    PointType;
  typedef float                             RealPixelType;
  typedef Image<RealPixelType, TDimension>  RealImageType;
  typedef typename RealImageType::IndexType CenterType;

  static constexpr unsigned int NumberOfDimensions = TDimension;

  itkNewMacro(Self);
  itkOverrideGetNameOfClassMacro(ScaleSpaceBlobSpatialObject);

  /** Set/Get the normalized laplacian value of the extrema */
  itkGetMacro(ScaleSpaceValue, double);
  itkSetMacro(ScaleSpaceValue, double);

  /** The radius of the object if it is a solid hyper-sphere */
  double
  GetObjectRadius(void) const
  {
    return this->GetSigmaInObjectSpace() * itk::Math::sqrt2;
  }

  /** The sigma of the laplacian where the extrema occoured */
  double
  GetScaleSpaceSigma(void) const
  {
    return this->GetSigmaInObjectSpace() / (std::sqrt(TDimension / 2.0));
  }

  /** The location where the extrema occoured */
  itkGetMacro(Center, CenterType);
  itkSetMacro(Center, CenterType);

private:
  double     m_ScaleSpaceValue;
  double     m_ObjectRadius;
  CenterType m_Center;
};

/** \class MultiScaleLaplacianBlobDetectorImageFilter
 * \brief Performs detection of "blob" by searching for local extrema
 * of the Laplacian of the Gassian (LoG) in pixel-space and
 * scale-space.
 *
 * Define a blob.
 *
 * Explain notation of scale space used
 *
 *
 *
 * \sa LaplacianRecursiveGaussianImageFilter
 *
 * \ingroup ITKImageScaleSpace
 *
 * \author Bradley Lowekamp
 */
template <typename TInputImage>
class MultiScaleLaplacianBlobDetectorImageFilter final : public ImageToImageFilter<TInputImage, TInputImage>
{
public:
  // todo figure out a better base class

  /** Standard class typedefs. */
  typedef MultiScaleLaplacianBlobDetectorImageFilter   Self;
  typedef ImageToImageFilter<TInputImage, TInputImage> Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** label image typedef */
  typedef float                                                   BlobRadiusType;
  typedef itk::Image<BlobRadiusType, TInputImage::ImageDimension> BlobRadiusImageType;
  typedef typename BlobRadiusImageType::Pointer                   BlobRadiusImagePointer;

  /** Blob typedef */
  typedef ScaleSpaceBlobSpatialObject<TInputImage::ImageDimension> BlobType;
  typedef typename BlobType::Pointer                               BlobPointer;
  typedef std::vector<BlobPointer>                                 BlobsListType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(MultiScaleLaplacianBlobDetectorImageFilter);

  /** Typedef to images */
  typedef TInputImage                           InputImageType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;

  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

  /** Set/Get the starting value to search scale-space.
   *
   * T is equivalent to varinace or sigma squared.
   **/
  itkSetMacro(StartT, double);
  itkGetMacro(StartT, double);

  /** Set/Get the ending value to search scale-space.
   */
  itkSetMacro(EndT, double);
  itkGetMacro(EndT, double);

  /** Set/Get the number of steps per doubling of sigma sampled.
   */
  itkSetMacro(StepsPerOctave, double);
  itkGetMacro(StepsPerOctave, double);

  /** Set/Get the number of blobs to find
   */
  itkSetMacro(NumberOfBlobs, size_t);
  itkGetMacro(NumberOfBlobs, size_t);

  /** Set/Get the label image
   */
  itkSetMacro(BlobRadiusImage, BlobRadiusImagePointer);
  itkGetMacro(BlobRadiusImage, BlobRadiusImagePointer);

  /** Pseudo-output
   * Get the list of circles. This recomputes the circles
   */
  BlobsListType &
  GetBlobs(void)
  {
    return this->m_BlobList;
  }

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  // / \todo need to concepts checked
  /** End concept checking */
#endif
protected:
  /** Internal Real ImageType */
  typedef float                                             RealPixelType;
  typedef Image<RealPixelType, TInputImage::ImageDimension> RealImageType;
  typedef typename RealImageType::Pointer                   RealImagePointer;
  typedef typename RealImageType::ConstPointer              RealImageConstPointer;

  MultiScaleLaplacianBlobDetectorImageFilter();

  // not defined or implemented as default works
  // virtual ~MultiScaleLaplacianBlobDetectorImageFilter( void ) {}

  void
  GenerateData() override;

  void
  ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId) override;

private:
  MultiScaleLaplacianBlobDetectorImageFilter(const Self &) = delete;
  void
  operator=(const Self &) = delete;

  class Blob
  {
  public:
    Blob(void) = default;

    Blob(typename RealImageType::IndexType center, double sigma, RealPixelType value)
      : m_Center(center)
      , m_Sigma(sigma)
      , m_Value(value)
    {}

    Blob(const Blob & b)
      : m_Center(b.m_Center)
      , m_Sigma(b.m_Sigma)
      , m_Value(b.m_Value)
    {}

    Blob &
    operator=(const Blob & b)
    {
      this->m_Center = b.m_Center;
      this->m_Sigma = b.m_Sigma;
      this->m_Value = b.m_Value;
      return *this;
    }

    typename RealImageType::IndexType m_Center;
    double                            m_Sigma;
    RealPixelType                     m_Value;
  };

  static bool
  BlobValueGreaterCompare(Blob & a, Blob & b)
  {
    return a.m_Value > b.m_Value;
  }

  static bool
  BlobValueLesserCompare(Blob & a, Blob & b)
  {
    return a.m_Value < b.m_Value;
  }

  typedef std::vector<Blob> BlobHeapType;

  // private member variable go here
  RealImageConstPointer     m_LaplacianImage[3];
  std::vector<BlobHeapType> m_BlobHeapPerThread;
  RealPixelType             m_GlobalMinimalBestBlobValue;

  double m_CurrentSigma;
  float  m_CurrentProgress;

  // private IVAR
  size_t m_NumberOfBlobs;

  BlobsListType m_BlobList;

  unsigned int m_StepsPerOctave;
  double       m_StartT;
  double       m_EndT;

  BlobRadiusImagePointer m_BlobRadiusImage;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkMultiScaleLaplacianBlobDetectorImageFilter.hxx"
#endif

#endif // __itkMultiScaleLaplacianBlobDetectorImageFilter_h
