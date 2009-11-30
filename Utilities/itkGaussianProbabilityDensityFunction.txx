/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGaussianProbabilityDensityFunction.txx,v $
  Language:  C++
  Date:      $Date: 2008/11/15 23:46:06 $
  Version:   $Revision: 1.17 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGaussianProbabilityDensityFunction_txx
#define __itkGaussianProbabilityDensityFunction_txx

#include "itkGaussianProbabilityDensityFunction.h"

#include "itkDecomposeTensorFunction.h"

namespace itk
{
namespace Statistics
{
template <class TMeasurementVector>
GaussianProbabilityDensityFunction<TMeasurementVector>
::GaussianProbabilityDensityFunction()
{
  this->m_Mean.SetSize( this->GetMeasurementVectorSize() );
  this->m_Mean.Fill( 0.0 );

  this->SetSigma( 1.0 );
  this->m_UseAnisotropicCovariance = false;

  this->m_GenerateRandomSamples = false;

  this->m_Randomizer = GeneratorType::New();
  this->m_Randomizer->SetSeed( static_cast<ITK_UINT32>( 0 ) );
}

template <class TMeasurementVector>
void
GaussianProbabilityDensityFunction<TMeasurementVector>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Mean: ";
  os << this->m_Mean << std::endl;

  if( this->m_UseAnisotropicCovariance )
    {
    os << indent << "Covariance: " << std::endl;
    os << this->m_Covariance.GetVnlMatrix();
    os << indent << "InverseCovariance: " << std::endl;
    os << indent << this->m_InverseCovariance.GetVnlMatrix();
    }
  else
    {
    os << indent << "Sigma: " << this->m_Sigma << std::endl;
    }

  os << indent << "Prefactor: " << this->m_PreFactor << std::endl;
}

template <class TMeasurementVector>
void
GaussianProbabilityDensityFunction<TMeasurementVector>
::SetCovariance( MatrixType cov )
{
  if( this->m_Covariance != cov )
    {
    // Sanity check
    if( cov.GetVnlMatrix().rows() != cov.GetVnlMatrix().cols() )
      {
      itkExceptionMacro( << "Covariance matrix must be square" );
      }
    if( this->GetMeasurementVectorSize() )
      {
      if( cov.GetVnlMatrix().rows() != this->GetMeasurementVectorSize() )
        {
        itkExceptionMacro( "Length of measurement vectors in the sample must be"
                           << " the same as the size of the covariance." );
        }
      }
    else
      {
      this->SetMeasurementVectorSize( cov.GetVnlMatrix().rows() );
      }

    this->m_Covariance = cov;

    this->m_IsCovarianceZero = this->m_Covariance.GetVnlMatrix().is_zero();

    if( !this->m_IsCovarianceZero )
      {
      // allocate the memory for this->m_InverseCovariance matrix
      this->m_InverseCovariance.GetVnlMatrix() =
        vnl_matrix_inverse<RealType>( this->m_Covariance.GetVnlMatrix() );

      // the determinant of the covaraince matrix
      typedef DecomposeTensorFunction<MatrixType, RealType> DecomposerType;
      typename DecomposerType::Pointer decomposer = DecomposerType::New();
      RealType det = decomposer->EvaluateDeterminant( this->m_Covariance );

      // calculate coefficient C of multivariate gaussian
      this->m_PreFactor = 1.0 / ( vcl_sqrt( det )
                                  * vcl_pow( sqrt( 2.0 * vnl_math::pi ),
                                             RealType( this->GetMeasurementVectorSize() ) ) );

      // calculate eigendecomposition

      if( this->m_GenerateRandomSamples )
        {
        typename DecomposerType::OutputMatrixType eigValues;
        decomposer->EvaluateSymmetricEigenDecomposition( this->m_Covariance,
                                                         eigValues, this->m_EigenVectors );

        this->m_EigenValues.SetSize( this->GetMeasurementVectorSize() );
        for( unsigned int i = 0; i < this->m_EigenValues.GetSize(); i++ )
          {
          this->m_EigenValues[i] = eigValues[i][i];
          }
        }
      }

    this->m_UseAnisotropicCovariance = true;
    this->Modified();
    }
}

template <class TMeasurementVector>
typename GaussianProbabilityDensityFunction<TMeasurementVector>::MatrixType
GaussianProbabilityDensityFunction<TMeasurementVector>
::GetCovariance()
{
  if( this->m_UseAnisotropicCovariance )
    {
    return this->m_Covariance;
    }
  else
    {
    MatrixType covariance( this->GetMeasurementVectorSize(),
                           this->GetMeasurementVectorSize() );
    covariance.SetIdentity();
    covariance *= ( this->m_Sigma * this->m_Sigma );
    return covariance;
    }
}

template <class TMeasurementVector>
typename GaussianProbabilityDensityFunction<TMeasurementVector>::MatrixType
GaussianProbabilityDensityFunction<TMeasurementVector>
::GetInverseCovariance()
{
  if( this->m_UseAnisotropicCovariance )
    {
    return this->m_InverseCovariance;
    }
  else
    {
    MatrixType invCovariance( this->GetMeasurementVectorSize(),
                              this->GetMeasurementVectorSize() );
    invCovariance.SetIdentity();
    invCovariance *= ( 1.0 / ( this->m_Sigma * this->m_Sigma ) );
    return invCovariance;
    }
}

template <class TMeasurementVector>
inline typename GaussianProbabilityDensityFunction<TMeasurementVector>::RealType
GaussianProbabilityDensityFunction<TMeasurementVector>
::Evaluate( const MeasurementVectorType & measurement ) const
{
  RealType temp;

  const MeasurementVectorSizeType measurementVectorSize
    = this->GetMeasurementVectorSize();
  MeanType tempVector;
  MeasurementVectorTraits::SetLength( tempVector, measurementVectorSize );

  // Compute |y - mean |
  for( unsigned int i = 0; i < measurementVectorSize; i++ )
    {
    tempVector[i] = measurement[i] - this->m_Mean[i];
    }

  if( this->m_UseAnisotropicCovariance )
    {
    MeanType tempVector2;
    MeasurementVectorTraits::SetLength( tempVector2, measurementVectorSize );
    // Compute |y - mean | * inverse(cov)
    for( unsigned int i = 0; i < measurementVectorSize; i++ )
      {
      temp = 0;
      for( unsigned int j = 0; j < measurementVectorSize; j++ )
        {
        temp += tempVector[j]
          * this->m_InverseCovariance.GetVnlMatrix().get(j, i);
        }
      tempVector2[i] = temp;
      }

    // Compute |y - mean | * inverse(cov) * |y - mean|^T
    temp = 0;
    for( unsigned int i = 0; i < measurementVectorSize; i++ )
      {
      temp += tempVector2[i] * tempVector[i];
      }

    return this->m_PreFactor * vcl_exp( -0.5 * temp );
    }
  else
    {
    temp = 0;
    for( unsigned int i = 0; i < measurementVectorSize; i++ )
      {
      temp += ( tempVector[i] * tempVector[i] );
      }
    return this->m_PreFactor
           * vcl_exp( -0.5 * temp / ( vnl_math_sqr( this->m_Sigma ) ) );
    }
}

template <class TMeasurementVector>
typename GaussianProbabilityDensityFunction<TMeasurementVector>
::MeasurementVectorType
GaussianProbabilityDensityFunction<TMeasurementVector>
::GenerateRandomSample()
{
  MeasurementVectorType measurement;

  measurement.Fill( 0 );

  if( this->m_GenerateRandomSamples )
    {
    if( this->m_UseAnisotropicCovariance )
      {
      MeanType sample( this->GetMeasurementVectorSize() );
      for( unsigned int d = 0; d < this->GetMeasurementVectorSize(); d++ )
        {
        sample[d] = this->m_Randomizer->GetNormalVariate(
            0.0, this->m_EigenValues[d] );
        }

      sample = this->m_EigenVectors * sample + this->m_Mean;
      for( unsigned int d = 0; d < this->GetMeasurementVectorSize(); d++ )
        {
        measurement[d] = sample[d];
        }
      }
    else
      {
      MeanType sample( this->GetMeasurementVectorSize() );
      for( unsigned int d = 0; d < this->GetMeasurementVectorSize(); d++ )
        {
        sample[d] = this->m_Randomizer->GetNormalVariate(
            0.0, vnl_math_sqr( this->m_Sigma )  );
        }
      sample += this->m_Mean;
      for( unsigned int d = 0; d < this->GetMeasurementVectorSize(); d++ )
        {
        measurement[d] = sample[d];
        }
      }
    }
  return measurement;
}
} // end namespace Statistics
} // end of namespace itk

#endif
