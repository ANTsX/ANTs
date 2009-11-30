/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGaussianProbabilityDensityFunction.h,v $
  Language:  C++
  Date:      $Date: 2008/11/15 23:46:06 $
  Version:   $Revision: 1.16 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGaussianProbabilityDensityFunction_h
#define __itkGaussianProbabilityDensityFunction_h

#include "itkArray.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkVariableSizeMatrix.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include "vnl/algo/vnl_determinant.h"
#include "vnl/vnl_math.h"

#include "itkMatrix.h"
#include "itkDensityFunction.h"

namespace itk
{
namespace Statistics
{
/** \class GaussianProbabilityDensityFunction
 * \brief GaussianProbabilityDensityFunction class represents Gaussian
 *  Density Function.
 *
 * This class keeps parameter to define Gaussian Density Function  and has
 * method to return the probability density of an instance (pattern).
 * If the all element of the covariance matrix is zero the "usual" density
 * calculations ignored. if the measurement vector to be evaluated is equal to
 * the mean, then the Evaluate method will return maximum value of
 * RealType and return 0 for others
 *
 */

template <class TMeasurementVector>
class ITK_EXPORT GaussianProbabilityDensityFunction :
  public         DensityFunction<TMeasurementVector>
{
public:
  /** Standard class typedefs */
  typedef GaussianProbabilityDensityFunction  Self;
  typedef DensityFunction<TMeasurementVector> Superclass;
  typedef SmartPointer<Self>                  Pointer;
  typedef SmartPointer<const Self>            ConstPointer;

  /** Strandard macros */
  itkTypeMacro( GaussianProbabilityDensityFunction, DensityFunction );
  itkNewMacro( Self );

  /** Typedef alias for the measurement vectors */
  typedef TMeasurementVector MeasurementVectorType;

  /** Length of each measurement vector */
  typedef typename
    Superclass::MeasurementVectorSizeType            MeasurementVectorSizeType;

  typedef double RealType;

  /** Type of the mean vector */
  typedef Array<RealType> MeanType;

  /** Type of the covariance matrix */
  typedef VariableSizeMatrix<RealType> MatrixType;
  typedef Array<RealType>              EigenValuesType;
  typedef VariableSizeMatrix<RealType> EigenVectorsType;

  typedef typename Statistics
    ::MersenneTwisterRandomVariateGenerator         GeneratorType;

  /** Sets the mean */
  void SetMean( const MeanType mean )
  {
    if( this->GetMeasurementVectorSize() )
      {
      if( mean.Size() != this->GetMeasurementVectorSize() )
        {
        itkExceptionMacro( "GaussianProbabilityDensityFunction::SetMean() "
                           << "Size of measurement vectors in the sample must the same as "
                           << "the size of the mean." );
        }
      }
    else
      {
      this->SetMeasurementVectorSize( mean.Size() );
      }

    if( this->m_Mean != mean )
      {
      this->m_Mean = mean;
      this->Modified();
      }
  }

  /** Gets the mean */
  const MeanType GetMean() const
  {
    return this->m_Mean;
  }

  /** Sets the covariance matrix.
   * Also, this function calculates inverse covariance and pre factor of
   * Gaussian Distribution to speed up GetProbability */
  void SetCovariance( MatrixType cov );

  /** Gets the covariance matrix */
  MatrixType GetCovariance();

  MatrixType GetInverseCovariance();

  void SetSigma( RealType sigma )
  {
    if( this->m_Sigma != sigma )
      {
      if( sigma <= 0 )
        {
        itkExceptionMacro( "Sigma must be greater than 0." );
        }
      this->m_Sigma = sigma;
      this->m_UseAnisotropicCovariance = false;

      this->m_PreFactor = 1.0 / ( this->m_Sigma
                                  * vcl_pow( sqrt( 2.0 * vnl_math::pi ),
                                             RealType( this->GetMeasurementVectorSize() ) ) );

      this->Modified();
      }
  }

  const RealType GetSigma() const
  {
    return this->m_Sigma;
  }

  itkSetMacro( GenerateRandomSamples, bool );
  itkGetConstMacro( GenerateRandomSamples, bool );
  itkBooleanMacro( GenerateRandomSamples );

  /** Gets the probability density of a measurement vector. */
  RealType Evaluate( const MeasurementVectorType & measurement ) const;

  /** Gets the probability density of a measurement vector. */
  MeasurementVectorType GenerateRandomSample();

protected:
  GaussianProbabilityDensityFunction( void );
  virtual ~GaussianProbabilityDensityFunction( void )
  {
  }

  void PrintSelf( std::ostream& os, Indent indent ) const;

private:
  MeanType   m_Mean;               // mean
  MatrixType m_Covariance;         // covariance matrix
  RealType   m_Sigma;              // sigma - use for isotropic gaussian

  bool             m_GenerateRandomSamples;
  EigenValuesType  m_EigenValues;
  EigenVectorsType m_EigenVectors;

  // if isotropic covariance is used, we don't want to make
  // unnecessary calculations.
  bool m_UseAnisotropicCovariance;

  // inverse covariance matrix which is automatically calculated
  // when covariace matrix is set.  This speed up the GetProbability()
  MatrixType m_InverseCovariance;

  // pre_factor which is automatically calculated
  // when covariace matirx is set.  This speeds up the GetProbability()
  RealType m_PreFactor;

  /** if the all element of the given covarinace is zero, then this
   * value set to true */
  bool m_IsCovarianceZero;

  typename GeneratorType::Pointer m_Randomizer;
};
} // end of namespace Statistics
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGaussianProbabilityDensityFunction.txx"
#endif

#endif
