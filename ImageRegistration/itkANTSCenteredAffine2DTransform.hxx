/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkANTSCenteredAffine2DTransform_hxx
#define _itkANTSCenteredAffine2DTransform_hxx

// #include "itkRigid2DTransform.h"
#include "itkANTSCenteredAffine2DTransform.h"
#include "vnl/algo/vnl_qr.h"

namespace itk
{
// Constructor with default arguments
// template<typename TScalarType>
// ANTSCenteredAffine2DTransform<TScalarType>::
// ANTSCenteredAffine2DTransform():
//   Superclass(OutputSpaceDimension, ParametersDimension)
// {
//   m_Angle = NumericTraits< TScalarType >::ZeroValue();
// }
template <typename TScalarType>
ANTSCenteredAffine2DTransform<TScalarType>::ANTSCenteredAffine2DTransform() :
  Superclass(ParametersDimension)
{
  m_Angle = NumericTraits<TScalarType>::ZeroValue();
  m_S1 = NumericTraits<TScalarType>::OneValue();
  m_S2 = NumericTraits<TScalarType>::OneValue();
  m_K = NumericTraits<TScalarType>::ZeroValue();
}

// Constructor with arguments
template <typename TScalarType>
ANTSCenteredAffine2DTransform<TScalarType>::ANTSCenteredAffine2DTransform( unsigned int spaceDimension,
                                                                           unsigned int parametersDimension) :
  Superclass(spaceDimension, parametersDimension)
{
  m_Angle = NumericTraits<TScalarType>::ZeroValue();
  m_S1 = NumericTraits<TScalarType>::OneValue();
  m_S2 = NumericTraits<TScalarType>::OneValue();
  m_K = NumericTraits<TScalarType>::ZeroValue();
}

// Destructor
template <typename TScalarType>
ANTSCenteredAffine2DTransform<TScalarType>::
~ANTSCenteredAffine2DTransform()
= default;

// Print self
template <typename TScalarType>
void
ANTSCenteredAffine2DTransform<TScalarType>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Angle       = " << m_Angle        << std::endl;
  os << indent << "S1          = " << m_S1 << std::endl;
  os << indent << "S2          = " << m_S2 << std::endl;
  os << indent << "K           = " << m_K << std::endl;
  //  os << indent << "S2          = " << m_ << std::endl;
}

// Set the rotation matrix
// template<typename TScalarType>
// void
// ANTSCenteredAffine2DTransform<TScalarType>::
// SetMatrix(const MatrixType & matrix )
// {
//   itkDebugMacro("setting  m_Matrix  to " << matrix );
//   // The matrix must be orthogonal otherwise it is not
//   // representing a valid rotaion in 2D space
//   typename MatrixType::InternalMatrixType test =
//     matrix.GetVnlMatrix() * matrix.GetTranspose();

//   const double tolerance = 1e-10;
//   if( !test.is_identity( tolerance ) )
//     {
//     itk::ExceptionObject ex(__FILE__,__LINE__,"Attempt to set a Non-Orthogonal matrix",ITK_LOCATION);
//     throw ex;
//     }

//   this->SetVarMatrix( matrix );
//   this->ComputeOffset();
//   this->ComputeMatrixParameters();
//   this->Modified();

// }

/** Compute the Angle from the Rotation Matrix */
template <typename TScalarType>
void
ANTSCenteredAffine2DTransform<TScalarType>
::ComputeMatrixParameters( void )
{
  typedef vnl_matrix_fixed<TScalarType, 2U, 2U> TMatrix;

  TMatrix A, Q, R;

  A = this->GetMatrix().GetVnlMatrix();
  vnl_qr<ScalarType> myqr(A);

  R = myqr.Q();   // Q() is the rotation
  Q = myqr.R();   // R() is the upper triangluar

  TMatrix dq;
  for( unsigned i = 0; i < 2; i++ )
    {
    dq(i, i) = (Q(i, i) >= 0) ? 1 : -1;
    }

  R = R * dq;
  Q = dq * Q;

//    std::cout << "A=" << A << std::endl;
//    std::cout << "Q=" << Q << std::endl;
//    std::cout << "R=" << R << std::endl;
//    std::cout << "dq=" << dq << std::endl;

  m_Angle = std::acos(R[0][0]);

  if( this->GetMatrix()[1][0] < itk::NumericTraits<TScalarType>::ZeroValue() )
    {
    m_Angle = -m_Angle;
    }

  m_S1 = Q[0][0];
  m_S2 = Q[1][1];
  m_K = Q[0][1] / Q[0][0];

  this->ComputeMatrix();

  if( static_cast<double>( this->GetMatrix()[1][0] ) -
      static_cast<double>( std::sin(m_Angle) ) > 0.000001 )
    {
    itkWarningMacro("Bad Rotation Matrix " << this->GetMatrix() );
    }
}

// // Compose with a translation
// template<typename TScalarType>
// void
// ANTSCenteredAffine2DTransform<TScalarType>::
// Translate(const OffsetType &offset, bool)
// {
//   OutputVectorType newOffset = this->GetOffset();
//   newOffset += offset;
//   this->SetOffset(newOffset);
// }

// Create and return an inverse transformation
// template<typename TScalarType>
// void
// ANTSCenteredAffine2DTransform<TScalarType>::
// CloneInverseTo( Pointer & result ) const
// {
//   result = New();
//   result->SetCenter( this->GetCenter() );  // inverse have the same center
//   result->SetAngle( -this->GetAngle() );
//   result->SetTranslation( -( this->GetInverseMatrix() * this->GetTranslation() ) );
// }

// Create and return a clone of the transformation
template <typename TScalarType>
void
ANTSCenteredAffine2DTransform<TScalarType>::CloneTo( Pointer & result ) const
{
  result = New();
  result->SetCenter( this->GetCenter() );
  result->SetAngle( this->GetAngle() );
  result->SetS1( this->GetS1() );
  result->SetS2( this->GetS2() );
  result->SetK( this->GetK() );
  result->SetTranslation( this->GetTranslation() );
}

// Reset the transform to an identity transform
template <typename TScalarType>
void
ANTSCenteredAffine2DTransform<TScalarType>::SetIdentity( void )
{
  this->Superclass::SetIdentity();
  m_Angle = NumericTraits<TScalarType>::ZeroValue();
  m_S1 = NumericTraits<TScalarType>::OneValue();
  m_S2 = NumericTraits<TScalarType>::OneValue();
  m_K = NumericTraits<TScalarType>::ZeroValue();
}

// Set the angle of rotation
template <typename TScalarType>
void
ANTSCenteredAffine2DTransform<TScalarType>
::SetAngle(TScalarType angle)
{
  m_Angle = angle;
  this->ComputeMatrix();
  this->ComputeOffset();
  this->Modified();
}

template <typename TScalarType>
void
ANTSCenteredAffine2DTransform<TScalarType>
::SetS1(TScalarType S1)
{
  m_S1 = S1;
  this->ComputeMatrix();
  this->ComputeOffset();
  this->Modified();
}

template <typename TScalarType>
void
ANTSCenteredAffine2DTransform<TScalarType>
::SetS2(TScalarType S2)
{
  m_S2 = S2;
  this->ComputeMatrix();
  this->ComputeOffset();
  this->Modified();
}

template <typename TScalarType>
void
ANTSCenteredAffine2DTransform<TScalarType>
::SetK(TScalarType K)
{
  m_K = K;
  this->ComputeMatrix();
  this->ComputeOffset();
  this->Modified();
}

// Set the angle of rotation
template <typename TScalarType>
void
ANTSCenteredAffine2DTransform<TScalarType>
::SetAngleInDegrees(TScalarType angle)
{
  const TScalarType angleInRadians = angle * std::atan(1.0) / 45.0;

  this->SetAngle( angleInRadians );
}

// Compute the matrix from the angle
template <typename TScalarType>
void
ANTSCenteredAffine2DTransform<TScalarType>
::ComputeMatrix( void )
{
  const double ca = std::cos(m_Angle );
  const double sa = std::sin(m_Angle );

  const double s1 = m_S1;
  const double s2 = m_S2;
  const double k = m_K;

  MatrixType rotationMatrix;

  rotationMatrix[0][0] = ca; rotationMatrix[0][1] = -sa;
  rotationMatrix[1][0] = sa; rotationMatrix[1][1] = ca;

  MatrixType scaleMatrix;
  scaleMatrix[0][0] = s1; scaleMatrix[0][1] = 0;
  scaleMatrix[1][0] = 0; scaleMatrix[1][1] = s2;

  MatrixType shearMatrix;
  shearMatrix[0][0] = 1; shearMatrix[0][1] = k;
  shearMatrix[1][0] = 0; shearMatrix[1][1] = 1;

  MatrixType varMatrix;

  varMatrix = rotationMatrix * scaleMatrix * shearMatrix;

  this->SetVarMatrix( varMatrix );
}

// Set Parameters
template <typename TScalarType>
void
ANTSCenteredAffine2DTransform<TScalarType>::SetParameters( const ParametersType & parameters )
{
  itkDebugMacro( << "Setting parameters " << parameters );

  // Set angle/s1/s2/k
  this->SetVarAngle( parameters[0] );
  this->SetVarS1( parameters[1] );
  this->SetVarS2( parameters[2] );
  this->SetVarK( parameters[3] );

  // Set the center
  InputPointType center;
  for( unsigned int i = 0; i < OutputSpaceDimension; i++ )
    {
    center[i] = parameters[i + 4];
    }
  this->SetVarCenter( center );

  // Set translation
  OutputVectorType translation;
  for( unsigned int i = 0; i < OutputSpaceDimension; i++ )
    {
    translation[i] = parameters[i + 6];
    }
  this->SetVarTranslation( translation );

  // Update matrix and offset
  this->ComputeMatrix();
  this->ComputeOffset();

  // Modified is always called since we just have a pointer to the
  // parameters and cannot know if the parameters have changed.
  this->Modified();

  itkDebugMacro(<< "After setting parameters ");
}

// Get Parameters
template <typename TScalarType>
const typename ANTSCenteredAffine2DTransform<TScalarType>::ParametersType
& ANTSCenteredAffine2DTransform<TScalarType>::
GetParameters( void ) const
  {
  itkDebugMacro( << "Getting parameters ");

  // Get the angle/s1/s2/k
  this->m_Parameters[0] = this->GetAngle();
  this->m_Parameters[1] = this->GetS1();
  this->m_Parameters[2] = this->GetS2();
  this->m_Parameters[3] = this->GetK();
  // Get the center
  for( unsigned int i = 0; i < OutputSpaceDimension; i++ )
    {
    this->m_Parameters[i + 4] = this->GetCenter()[i];
    }
  // Get the translation
  for( unsigned int i = 0; i < OutputSpaceDimension; i++ )
    {
    this->m_Parameters[i + 6] = this->GetTranslation()[i];
    }

  itkDebugMacro(<< "After getting parameters " << this->m_Parameters );

  return this->m_Parameters;
  }

// // Compute transformation Jacobian
// template<typename TScalarType>
// const typename ANTSCenteredAffine2DTransform<TScalarType>::JacobianType &
// ANTSCenteredAffine2DTransform<TScalarType>::
// GetJacobian( const InputPointType & p ) const
// {
//
//    const double ca = std::cos(this->GetAngle() );
//    const double sa = std::sin(this->GetAngle() );
//    const double s1 = m_S1;
//    const double s2 = m_S2;
//    const double k = m_K;
//
//    this->m_Jacobian.Fill(0.0);
//
//    const double cx = this->GetCenter()[0];
//    const double cy = this->GetCenter()[1];
//
//    // derivatives with respect to the angle
//    // this->m_Jacobian[0][0] = -sa * ( p[0] - cx ) - ca * ( p[1] - cy );
//    // this->m_Jacobian[1][0] =  ca * ( p[0] - cx ) - sa * ( p[1] - cy );
//
//    double pxoff = (p[0]-cx)+k*(p[1]-cy);
//    double pyoff = p[1]-cy;
//
//    // wrt. theta
//    this->m_Jacobian[0][0] = s1*( pxoff )*(-sa) + s2*( pyoff )*(-ca);
//    this->m_Jacobian[1][0] = s1*( pxoff )*(ca) + s2*( pyoff )*(-sa);
//
//    // wrt. s1/s2
//    this->m_Jacobian[0][1] = ca * pxoff;
//    this->m_Jacobian[0][2] = -sa * pyoff;
//
//    this->m_Jacobian[1][1] = sa * pxoff;
//    this->m_Jacobian[1][2] = ca * pyoff;
//
//    // wrt. k
//    this->m_Jacobian[0][3] = ca * s1 * pyoff;
//    this->m_Jacobian[1][3] = sa * s1 * pyoff;
//
//    // wrt. cx/cy
//    this->m_Jacobian[0][4] = - s1*ca + 1.0;
//    this->m_Jacobian[0][5] = -(k*s1*ca - s2*sa);
//    this->m_Jacobian[1][4] = - s1*sa;
//    this->m_Jacobian[1][5] = -(k*s1*sa + s2*ca) + 1.0;
//
//
//    // wrt. t1/t2
//    this->m_Jacobian[0][6] = 1.0;
//    this->m_Jacobian[1][7] = 1.0;
//
//    // compute derivatives for the translation part
//    // unsigned int blockOffset = 1;
//    // for(unsigned int dim=0; dim < OutputSpaceDimension; dim++ )
//    //  {
//    //  this->m_Jacobian[ dim ][ blockOffset + dim ] = 1.0;
//    //  }
//
//    return this->m_Jacobian;
//
// }

template <typename TScalarType>
void
ANTSCenteredAffine2DTransform<TScalarType>::ComputeJacobianWithRespectToParameters(const InputPointType  & p,
                                                                                   JacobianType & j) const
{
  const double ca = std::cos(this->GetAngle() );
  const double sa = std::sin(this->GetAngle() );
  const double s1 = m_S1;
  const double s2 = m_S2;
  const double k = m_K;

  j.SetSize( this->GetOutputSpaceDimension(), this->GetNumberOfLocalParameters() );
  j.Fill(0.0);

  const double cx = this->GetCenter()[0];
  const double cy = this->GetCenter()[1];

  // derivatives with respect to the angle
  // j[0][0] = -sa * ( p[0] - cx ) - ca * ( p[1] - cy );
  // j[1][0] =  ca * ( p[0] - cx ) - sa * ( p[1] - cy );

  double pxoff = ( static_cast<double>( p[0] ) - cx) + k * ( static_cast<double>( p[1] ) - cy );
  double pyoff = static_cast<double>( p[1] ) - cy;

  // wrt. theta
  j[0][0] = s1 * ( pxoff ) * (-sa) + s2 * ( pyoff ) * (-ca);
  j[1][0] = s1 * ( pxoff ) * (ca) + s2 * ( pyoff ) * (-sa);

  // wrt. s1/s2
  j[0][1] = ca * pxoff;
  j[0][2] = -sa * pyoff;

  j[1][1] = sa * pxoff;
  j[1][2] = ca * pyoff;

  // wrt. k
  j[0][3] = ca * s1 * pyoff;
  j[1][3] = sa * s1 * pyoff;

  // wrt. cx/cy
  j[0][4] = -s1 * ca + 1.0;
  j[0][5] = -(k * s1 * ca - s2 * sa);
  j[1][4] = -s1 * sa;
  j[1][5] = -(k * s1 * sa + s2 * ca) + 1.0;

  // wrt. t1/t2
  j[0][6] = 1.0;
  j[1][7] = 1.0;
}

// Back transform a point
template <typename TScalarType>
typename ANTSCenteredAffine2DTransform<TScalarType>::InputPointType
ANTSCenteredAffine2DTransform<TScalarType>::BackTransform(const OutputPointType & point) const
{
  itkWarningMacro(
    <<
    "BackTransform(): This method is slated to be removed from ITK.  Instead, please use GetInverse() to generate an inverse transform and then perform the transform using that inverted transform.");
  return this->GetInverseMatrix() * (point - this->GetOffset() );
}

// Back transform a vector
template <typename TScalarType>
typename ANTSCenteredAffine2DTransform<TScalarType>::InputVectorType
ANTSCenteredAffine2DTransform<TScalarType>::BackTransform(const OutputVectorType & vect ) const
{
  itkWarningMacro(
    <<
    "BackTransform(): This method is slated to be removed from ITK.  Instead, please use GetInverse() to generate an inverse transform and then perform the transform using that inverted transform.");
  return this->GetInverseMatrix() * vect;
}

// Back transform a vnl_vector
template <typename TScalarType>
typename ANTSCenteredAffine2DTransform<TScalarType>::InputVnlVectorType
ANTSCenteredAffine2DTransform<TScalarType>::BackTransform(const OutputVnlVectorType & vect ) const
{
  itkWarningMacro(
    <<
    "BackTransform(): This method is slated to be removed from ITK.  Instead, please use GetInverse() to generate an inverse transform and then perform the transform using that inverted transform.");
  return this->GetInverseMatrix() * vect;
}

// Back Transform a CovariantVector
template <typename TScalarType>
typename ANTSCenteredAffine2DTransform<TScalarType>::InputCovariantVectorType
ANTSCenteredAffine2DTransform<TScalarType>::BackTransform(const OutputCovariantVectorType & vect) const
{
  itkWarningMacro(
    <<
    "BackTransform(): This method is slated to be removed from ITK.  Instead, please use GetInverse() to generate an inverse transform and then perform the transform using that inverted transform.");
  return this->GetMatrix() * vect;
}
} // namespace

#endif
