/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGaussianListSampleFunction.txx,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGaussianListSampleFunction_txx
#define __itkGaussianListSampleFunction_txx

#include "itkGaussianListSampleFunction.h"
#include "itkCovarianceCalculator.h"
#include "itkMeanCalculator.h"
#include "itkWeightedCovarianceCalculator.h"
#include "itkWeightedMeanCalculator.h"

namespace itk
{
namespace Statistics
{
template <class TListSample, class TOutput, class TCoordRep>
GaussianListSampleFunction<TListSample, TOutput, TCoordRep>
::GaussianListSampleFunction()
{
  this->m_Gaussian = GaussianType::New();
}

template <class TListSample, class TOutput, class TCoordRep>
GaussianListSampleFunction<TListSample, TOutput, TCoordRep>
::~GaussianListSampleFunction()
{
}

template <class TListSample, class TOutput, class TCoordRep>
void
GaussianListSampleFunction<TListSample, TOutput, TCoordRep>
::SetInputListSample( const InputListSampleType * ptr )
{
  this->m_ListSample = ptr;

  if( this->m_Weights.Size() == this->m_ListSample->Size() )
    {
    typedef typename Statistics::
      WeightedMeanCalculator<InputListSampleType> MeanCalculatorType;
    typename MeanCalculatorType::Pointer meanCalculator =
      MeanCalculatorType::New();
    meanCalculator->SetWeights( &this->m_Weights );
    meanCalculator->SetInputSample( this->GetInputListSample() );
    meanCalculator->Update();

    typedef typename Statistics::
      WeightedCovarianceCalculator<InputListSampleType> CovarianceCalculatorType;
    typename CovarianceCalculatorType::Pointer covarianceCalculator =
      CovarianceCalculatorType::New();
    covarianceCalculator->SetWeights( &this->m_Weights );
    covarianceCalculator->SetMean( meanCalculator->GetOutput() );
    covarianceCalculator->SetInputSample( this->GetInputListSample() );
    covarianceCalculator->Update();

    this->m_Gaussian->SetMean( *meanCalculator->GetOutput() );
    this->m_Gaussian->SetCovariance( *covarianceCalculator->GetOutput() );
    }
  else
    {
    typedef Statistics::MeanCalculator<InputListSampleType> MeanCalculatorType;
    typename MeanCalculatorType::Pointer meanCalculator =
      MeanCalculatorType::New();
    meanCalculator->SetInputSample( this->GetInputListSample() );
    meanCalculator->Update();

    typedef Statistics::CovarianceCalculator<InputListSampleType>
      CovarianceCalculatorType;
    typename CovarianceCalculatorType::Pointer covarianceCalculator =
      CovarianceCalculatorType::New();
    covarianceCalculator->SetMean( meanCalculator->GetOutput() );
    covarianceCalculator->SetInputSample( this->GetInputListSample() );
    covarianceCalculator->Update();

    this->m_Gaussian->SetMean( *meanCalculator->GetOutput() );
    this->m_Gaussian->SetCovariance( *covarianceCalculator->GetOutput() );
    }
}

template <class TListSample, class TOutput, class TCoordRep>
TOutput
GaussianListSampleFunction<TListSample, TOutput, TCoordRep>
::Evaluate( const InputMeasurementVectorType & measurement ) const
{
  return this->m_Gaussian->Evaluate( measurement );
}

/**
 * Standard "PrintSelf" method
 */
template <class TListSample, class TOutput, class TCoordRep>
void
GaussianListSampleFunction<TListSample, TOutput, TCoordRep>
::PrintSelf(
  std::ostream& os,
  Indent indent) const
{
  os << indent << "mean = " << this->m_Gaussian->GetMean() << ", ";

  typename GaussianType::CovarianceType covariance =
    this->m_Gaussian->GetCovariance();
  os << "covariance = [";
  for( unsigned int r = 0; r < covariance.Rows(); r++ )
    {
    for( unsigned int c = 0; c < covariance.Cols() - 1; c++ )
      {
      os << covariance( r, c ) << ", ";
      }
    if( r == covariance.Rows() - 1 )
      {
      os << covariance( r, covariance.Cols() - 1 ) << "]" << std::endl;
      }
    else
      {
      os << covariance( r, covariance.Cols() - 1 ) << "; ";
      }
    }
}
} // end of namespace Statistics
} // end of namespace itk

#endif
