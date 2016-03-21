#ifndef __itkLabelSelectionImageAdaptor_h
#define __itkLabelSelectionImageAdaptor_h

#include "itkImageAdaptor.h"

namespace itk
{
namespace Accessor
{
/** \class LabelSelectionPixelAccessor
 * \brief Return a binary mask of the selected label
 *
 * LabelSelectionPixelAccessor is templated over an internal type and an
 * external type representation. This class cast the input
 * applies the function to it and cast the result according
 * to the types defined as template parameters
 *
 * \ingroup ImageAdaptors
 * \ingroup ITKImageAdaptors
 */
template< class TInternalType, class TExternalType >
class ITK_EXPORT LabelSelectionPixelAccessor
{
public:
  /** External typedef. It defines the external aspect
   * that this class will exhibit. */
  typedef TExternalType ExternalType;

  /** Internal typedef. It defines the internal real
   * representation of data. */
  typedef TInternalType InternalType;

	void SetAcceptedValue(TInternalType value) { m_AcceptedValue = value; }

  inline TExternalType Get(const TInternalType & input) const
  {
    return (TExternalType)(
             ( input == m_AcceptedValue ) ? 1 : 0 );
  }
protected:
	TInternalType m_AcceptedValue;
};
} // end namespace Accessor

/** \class LabelSelectionImageAdaptor
 * \brief Presents a label image as a binary image of one label
 *
 * Additional casting is performed according to the input and output image
 * types following C++ default casting rules.
 *
 * \ingroup ImageAdaptors
 * \ingroup ITKImageAdaptors
 */
template< class TImage, class TOutputPixelType >
class ITK_EXPORT LabelSelectionImageAdaptor:public
  ImageAdaptor< TImage,
                Accessor::LabelSelectionPixelAccessor<
                  typename TImage::PixelType,
                  TOutputPixelType >   >
{
public:
  /** Standard class typedefs. */
  typedef LabelSelectionImageAdaptor Self;
  typedef ImageAdaptor< TImage, Accessor::LabelSelectionPixelAccessor<
                          typename TImage::PixelType,
                          TOutputPixelType > >  Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(LabelSelectionImageAdaptor, ImageAdaptor);

	void SetAcceptedValue(typename TImage::PixelType value) {
		this->GetPixelAccessor().SetAcceptedValue(value);
	}

protected:
  LabelSelectionImageAdaptor() {}
  virtual ~LabelSelectionImageAdaptor() {}

private:
  LabelSelectionImageAdaptor(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented
};
} // end namespace itk

#endif
