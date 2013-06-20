#ifndef ITKDEFORMATIONFIELDFROMMULTITRANSFORMFILTER_H_
#define ITKDEFORMATIONFIELDFROMMULTITRANSFORMFILTER_H_

#include "itkWarpImageMultiTransformFilter.h"

namespace itk
{
template <
  class TOutputImage,
  class TDisplacementField,
  class TTransform
  >
class DisplacementFieldFromMultiTransformFilter :
  public         WarpImageMultiTransformFilter<TOutputImage, TOutputImage, TDisplacementField, TTransform>
{
public:
  /** Standard class typedefs. */
  typedef TOutputImage                                                                             TInputImage;
  typedef DisplacementFieldFromMultiTransformFilter                                                Self;
  typedef WarpImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform> Superclass;
  typedef SmartPointer<Self>                                                                       Pointer;
  typedef SmartPointer<const Self>                                                                 ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( DisplacementFieldFromMultiTransformFilter, WarpImageMultiTransformFilter );

  /** Typedef to describe the output image region type. */
  typedef typename TOutputImage::RegionType OutputImageRegionType;

  /** Inherit some types from the superclass. */
  typedef typename Superclass::InputImageType         InputImageType;
  typedef typename Superclass::InputImagePointer      InputImagePointer;
  typedef typename Superclass::OutputImageType        OutputImageType;
  typedef typename Superclass::OutputImagePointer     OutputImagePointer;
  typedef typename Superclass::InputImageConstPointer InputImageConstPointer;
  typedef typename OutputImageType::IndexType         IndexType;
  typedef typename OutputImageType::SizeType          SizeType;
  typedef typename OutputImageType::PixelType         PixelType;
  typedef typename OutputImageType::SpacingType       SpacingType;

  /** Determine the image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension );
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension );
  itkStaticConstMacro(DisplacementFieldDimension, unsigned int,
                      TDisplacementField::ImageDimension );

  /** Displacement field typedef support. */
  typedef TDisplacementField                        DisplacementFieldType;
  typedef typename DisplacementFieldType::Pointer   DisplacementFieldPointer;
  typedef typename DisplacementFieldType::PixelType DisplacementType;
  typedef typename DisplacementType::ValueType      DisplacementScalarValueType;

  typedef typename Superclass::PointType PointType;

  virtual void GenerateInputRequestedRegion()
  {
    Superclass::GenerateInputRequestedRegion();
  }

  virtual void BeforeThreadedGenerateData()
  {
  };

  virtual void AfterThreadedGenerateData()
  {
  };

  virtual void GenerateOutputInformation()
  {
    // call the superclass's implementation
    Superclass::GenerateOutputInformation();
  };
protected:
  DisplacementFieldFromMultiTransformFilter() : Superclass()
  {
    this->SetNumberOfRequiredInputs( 0 );
    const DisplacementScalarValueType kMaxDisp = itk::NumericTraits<DisplacementScalarValueType>::max();
    Superclass::m_EdgePaddingValue.Fill(kMaxDisp);
  }

  ~DisplacementFieldFromMultiTransformFilter()
  {
  };
  void PrintSelf(std::ostream& os, Indent indent) const
  {
    Superclass::PrintSelf(os, indent);
  };

  /** WarpImageMultiTransformFilter is implemented as a multi-threaded filter.
   * As such, it needs to provide and implementation for
   * ThreadedGenerateData(). */
  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                            ThreadIdType threadId )
  {
    OutputImagePointer outputPtr = this->GetOutput();

    // support progress methods/callbacks
    ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels() );

    // iterator for the output image
    ImageRegionIteratorWithIndex<OutputImageType> outputIt(outputPtr, outputRegionForThread);

//        int cnt = 0;
    while( !outputIt.IsAtEnd() )
      {
      PointType point1, point2;

      // get the output image index
      IndexType index = outputIt.GetIndex();
      outputPtr->TransformIndexToPhysicalPoint( index, point1 );

      const bool isinside = this->MultiTransformPoint(point1, point2, Superclass::m_bFirstDeformNoInterp, index);

      if( isinside )
        {
        PixelType value;
        for( int ii = 0; ii < OutputImageType::ImageDimension; ii++ )
          {
          value[ii] = point2[ii] - point1[ii];
          }

        outputIt.Set( value );
        }
      else
        {
        PixelType                         value;
        const DisplacementScalarValueType kMaxDisp = itk::NumericTraits<DisplacementScalarValueType>::max();
        for( int ii = 0; ii < OutputImageType::ImageDimension; ii++ )
          {
          value[ii] = kMaxDisp;
          }
        outputIt.Set( value );
        }

      ++outputIt;
      }

    progress.CompletedPixel();
  };
};
} // end namespace itk
#endif /*ITKDEFORMATIONFIELDFROMMULTITRANSFORMFILTER_H_*/
