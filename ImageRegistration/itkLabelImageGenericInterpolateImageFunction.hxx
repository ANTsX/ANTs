#ifndef __itkLabelImageGenericInterpolateImageFunction_hxx
#define __itkLabelImageGenericInterpolateImageFunction_hxx

#include "itkLabelImageGenericInterpolateImageFunction.h"
#include <itkImageRegionConstIterator.h>

namespace itk
{

template<typename TInputImage, template<class, typename> class TInterpolator , typename TCoordRep>
LabelImageGenericInterpolateImageFunction<TInputImage, TInterpolator, TCoordRep>
::LabelImageGenericInterpolateImageFunction()
{
}

template<typename TInputImage, template<class, typename> class TInterpolator , typename TCoordRep>
void LabelImageGenericInterpolateImageFunction<TInputImage,TInterpolator, TCoordRep>
::SetInputImage( const TInputImage *image ) {
	/** We have one adaptor and one interpolator per label to keep the class thread-safe:
	 *  changing the adaptor's accepted value wouldn't work when called from a multi-threaded filter */
	if (image) {
		m_Labels.clear();
		typedef itk::ImageRegionConstIterator<TInputImage> IteratorType;
		IteratorType it(image,image->GetLargestPossibleRegion());
		it.GoToBegin();
		for (; !it.IsAtEnd(); ++it) {
			m_Labels.insert(it.Get());
		}
		m_InternalInterpolators.clear();
		m_LabelSelectionAdaptors.clear();
		for(typename LabelSetType::const_iterator i=m_Labels.begin(); i != m_Labels.end(); ++i) {
			typename LabelSelectionAdaptorType::Pointer adapt = LabelSelectionAdaptorType::New();
			// This adaptor doesn't implement Set() so this should be safe
			adapt->SetImage(const_cast<TInputImage*>(image));
			adapt->SetAcceptedValue(*i);
			m_LabelSelectionAdaptors.push_back(adapt);
			typename InternalInterpolatorType::Pointer interp = InternalInterpolatorType::New();
			interp->SetInputImage(adapt);
			m_InternalInterpolators.push_back(interp);
		}
	}
	Superclass::SetInputImage(image);
}

template<typename TInputImage, template<class, typename> class TInterpolator , typename TCoordRep>
typename LabelImageGenericInterpolateImageFunction<TInputImage, TInterpolator, TCoordRep>
::OutputType
LabelImageGenericInterpolateImageFunction<TInputImage, TInterpolator, TCoordRep>
::EvaluateAtContinuousIndex( const ContinuousIndexType & cindex, OutputType * grad ) const
{
	/** Interpolate the binary mask corresponding to each label and return the label
	 * with the highest value */
	double value=0;
	typename TInputImage::PixelType best_label=0;
	typename LabelSetType::const_iterator it;
	int i;
	for(it=m_Labels.begin(), i = 0;
		it != m_Labels.end(); ++it, ++i) {
		double tmp = m_InternalInterpolators[i]->EvaluateAtContinuousIndex(cindex);
		if( tmp > value) {
			value = tmp;
			best_label = (*it);
		}
	}
		return best_label;
}

} // namespace itk

#endif
