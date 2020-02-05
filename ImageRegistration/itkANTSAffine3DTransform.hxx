#ifndef _itkANTSAffine3DTransform_hxx
#define _itkANTSAffine3DTransform_hxx

#include "itkANTSAffine3DTransform.h"
#include "vnl/algo/vnl_qr.h"

namespace itk
{
// Constructor with default arguments
template <typename TScalarType>
ANTSAffine3DTransform<TScalarType>::ANTSAffine3DTransform() :
  Superclass(ParametersDimension)
{
  m_Rotation = VnlQuaternionType(0, 0, 0, 1);   // axis * std::sin(t/2), std::cos(t/2)
  m_S1 = NumericTraits<TScalarType>::OneValue();
  m_S2 = NumericTraits<TScalarType>::OneValue();
  m_S3 = NumericTraits<TScalarType>::OneValue();
  m_K1 = NumericTraits<TScalarType>::ZeroValue();
  m_K2 = NumericTraits<TScalarType>::ZeroValue();
  m_K3 = NumericTraits<TScalarType>::ZeroValue();
}

// Constructor with default arguments
template <typename TScalarType>
ANTSAffine3DTransform<TScalarType>::ANTSAffine3DTransform(unsigned int outputSpaceDimension,
                                                          unsigned int parametersDimension) :
  Superclass(outputSpaceDimension, parametersDimension)
{
  m_Rotation = VnlQuaternionType(0, 0, 0, 1);   // axis * std::sin(t/2), std::cos(t/2)
  m_S1 = NumericTraits<TScalarType>::OneValue();
  m_S2 = NumericTraits<TScalarType>::OneValue();
  m_S3 = NumericTraits<TScalarType>::OneValue();
  m_K1 = NumericTraits<TScalarType>::ZeroValue();
  m_K2 = NumericTraits<TScalarType>::ZeroValue();
  m_K3 = NumericTraits<TScalarType>::ZeroValue();
}

// // Constructor with explicit arguments
// template<typename TScalarType>
// ANTSAffine3DTransform<TScalarType>::
// ANTSAffine3DTransform( const MatrixType & matrix,
//                           const OutputVectorType & offset ) :
//   Superclass(matrix, offset)
// {
//   this->ComputeMatrixParameters();
// }

// Print self
template <typename TScalarType>
void ANTSAffine3DTransform<TScalarType>::PrintSelf(std::ostream & os,
                                                   Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Rotation:    " << m_Rotation << std::endl;
  os << indent << "S1, S2, S3: " << m_S1 << ", " << m_S2 << ", " << m_S3
     << std::endl;
  os << indent << "K1, K2, K3: " << m_K1 << ", " << m_K2 << ", " << m_K3
     << std::endl;
}

// Set rotation
template <typename TScalarType>
void ANTSAffine3DTransform<TScalarType>::SetRotation(
  const VnlQuaternionType & rotation)
{
  m_Rotation = rotation;

  this->ComputeMatrix();
  this->ComputeOffset();
  this->Modified();

  return;
}

template <typename TScalarType>
void ANTSAffine3DTransform<TScalarType>::SetS1(const TScalarType S1)
{
  m_S1 = S1;

  this->ComputeMatrix();
  this->ComputeOffset();
  this->Modified();

  return;
}

template <typename TScalarType>
void ANTSAffine3DTransform<TScalarType>::SetS2(const TScalarType S2)
{
  m_S2 = S2;

  this->ComputeMatrix();
  this->ComputeOffset();
  this->Modified();

  return;
}

template <typename TScalarType>
void ANTSAffine3DTransform<TScalarType>::SetS3(const TScalarType S3)
{
  m_S3 = S3;

  this->ComputeMatrix();
  this->ComputeOffset();
  this->Modified();

  return;
}

template <typename TScalarType>
void ANTSAffine3DTransform<TScalarType>::SetK1(const TScalarType K1)
{
  m_K1 = K1;

  this->ComputeMatrix();
  this->ComputeOffset();
  this->Modified();

  return;
}

template <typename TScalarType>
void ANTSAffine3DTransform<TScalarType>::SetK2(const TScalarType K2)
{
  m_K2 = K2;

  this->ComputeMatrix();
  this->ComputeOffset();
  this->Modified();

  return;
}

template <typename TScalarType>
void ANTSAffine3DTransform<TScalarType>::SetK3(const TScalarType K3)
{
  m_K3 = K3;

  this->ComputeMatrix();
  this->ComputeOffset();
  this->Modified();

  return;
}

// Set the parameters in order to fit an Identity transform
template <typename TScalarType>
void ANTSAffine3DTransform<TScalarType>::SetIdentity(void)
{
  m_Rotation = VnlQuaternionType(0, 0, 0, 1);
  m_S1 = NumericTraits<TScalarType>::OneValue();
  m_S2 = NumericTraits<TScalarType>::OneValue();
  m_S3 = NumericTraits<TScalarType>::OneValue();
  m_K1 = NumericTraits<TScalarType>::ZeroValue();
  m_K2 = NumericTraits<TScalarType>::ZeroValue();
  m_K3 = NumericTraits<TScalarType>::ZeroValue();
  this->Superclass::SetIdentity();
}

// Set Parameters
template <typename TScalarType>
void ANTSAffine3DTransform<TScalarType>::SetParameters(
  const ParametersType & parameters)
{
  OutputVectorType translation;

  // Transfer the quaternion part
  unsigned int par = 0;

  for( unsigned int j = 0; j < 4; j++ )
    {
    m_Rotation[j] = parameters[par];
    ++par;
    }

  m_S1 = parameters[par++];
  m_S2 = parameters[par++];
  m_S3 = parameters[par++];
  m_K1 = parameters[par++];
  m_K2 = parameters[par++];
  m_K3 = parameters[par++];

  this->ComputeMatrix();
  // Transfer the constant part
  for( unsigned int i = 0; i < SpaceDimension; i++ )
    {
    translation[i] = parameters[par];
    ++par;
    }
  this->SetVarTranslation(translation);

  this->ComputeOffset();

  // Modified is always called since we just have a pointer to the
  // parameters and cannot know if the parameters have changed.
  this->Modified();
}

// Set Parameters
template <typename TScalarType>
const typename ANTSAffine3DTransform<TScalarType>::ParametersType
& ANTSAffine3DTransform<TScalarType>::GetParameters() const {
  VnlQuaternionType quaternion = this->GetRotation();
  OutputVectorType  translation = this->GetTranslation();

  // Transfer the quaternion part
  unsigned int par = 0;
  for( unsigned int j = 0; j < 4; j++ )
    {
    this->m_Parameters[par] = quaternion[j];
    ++par;
    }

  this->m_Parameters[par++] = m_S1;
  this->m_Parameters[par++] = m_S2;
  this->m_Parameters[par++] = m_S3;
  this->m_Parameters[par++] = m_K1;
  this->m_Parameters[par++] = m_K2;
  this->m_Parameters[par++] = m_K3;
  // Transfer the constant part
  for( unsigned int i = 0; i < SpaceDimension; i++ )
    {
    this->m_Parameters[par] = translation[i];
    ++par;
    }

  return this->m_Parameters;
  }

// // Get parameters
// template<typename TScalarType>
// const typename ANTSAffine3DTransform<TScalarType>::JacobianType &
// ANTSAffine3DTransform<TScalarType>::GetJacobian(
//        const InputPointType & p) const {
//
//    this->m_Jacobian.Fill(0.0);
//
//    TScalarType c1 = this->GetCenter()[0];
//    TScalarType c2 = this->GetCenter()[1];
//    TScalarType c3 = this->GetCenter()[2];
//    TScalarType s1 = this->m_S1;
//    TScalarType s2 = this->m_S2;
//    TScalarType s3 = this->m_S3;
//    TScalarType k1 = this->m_K1;
//    TScalarType k2 = this->m_K2;
//    TScalarType k3 = this->m_K3;
//    TScalarType x1 = p[0];
//    TScalarType x2 = p[1];
//    TScalarType x3 = p[2];
//
//    // z1,z2,z3 is the S*K*point p
//    TScalarType w1 = (x1 - c1) + k1 * (x2 - c2) + k2 * (x3 - c3);
//    TScalarType w2 = (x2 - c2) + k3 * (x3 - c2);
//    TScalarType w3 = (x3 - c3);
//
//    TScalarType z1 = s1 * w1;
//    TScalarType z2 = s2 * w2;
//    TScalarType z3 = s3 * w3;
//
//    // compute Jacobian with respect to quaternion parameters
//    this->m_Jacobian[0][0] = 2.0
//            * (m_Rotation.x() * z1 + m_Rotation.y() * z2 + m_Rotation.z() * z3);
//    this->m_Jacobian[0][1] =
//            2.0
//                    * (-m_Rotation.y() * z1 + m_Rotation.x() * z2
//                            + m_Rotation.r() * z3);
//    this->m_Jacobian[0][2] =
//            2.0
//                    * (-m_Rotation.z() * z1 - m_Rotation.r() * z2
//                            + m_Rotation.x() * z3);
//    this->m_Jacobian[0][3] =
//            -2.0
//                    * (-m_Rotation.r() * z1 + m_Rotation.z() * z2
//                            - m_Rotation.y() * z3);
//
//    this->m_Jacobian[1][0] = -this->m_Jacobian[0][1];
//    this->m_Jacobian[1][1] = this->m_Jacobian[0][0];
//    this->m_Jacobian[1][2] = this->m_Jacobian[0][3];
//    this->m_Jacobian[1][3] = -this->m_Jacobian[0][2];
//
//    this->m_Jacobian[2][0] = -this->m_Jacobian[0][2];
//    this->m_Jacobian[2][1] = -this->m_Jacobian[0][3];
//    this->m_Jacobian[2][2] = this->m_Jacobian[0][0];
//    this->m_Jacobian[2][3] = this->m_Jacobian[0][1];
//
//    // get rotation matrix first
//    // this is done to compensate for the transposed representation
//    // between VNL and ITK
//    VnlQuaternionType conjugateRotation = m_Rotation.conjugate();
//    MatrixType newMatrix;
//    newMatrix = conjugateRotation.rotation_matrix_transpose();
//
//    TScalarType r11 = newMatrix[0][0];
//    TScalarType r12 = newMatrix[0][1];
//    TScalarType r13 = newMatrix[0][2];
//    TScalarType r21 = newMatrix[1][0];
//    TScalarType r22 = newMatrix[1][1];
//    TScalarType r23 = newMatrix[1][2];
//    TScalarType r31 = newMatrix[2][0];
//    TScalarType r32 = newMatrix[2][1];
//    TScalarType r33 = newMatrix[2][2];
//
//    // compute Jacobian wrt S1/S2/S3
//    this->m_Jacobian[0][4] = r11 * w1;
//    this->m_Jacobian[0][5] = r12 * w2;
//    this->m_Jacobian[0][6] = r13 * w3;
//    this->m_Jacobian[1][4] = r21 * w1;
//    this->m_Jacobian[1][5] = r22 * w2;
//    this->m_Jacobian[1][6] = r23 * w3;
//    this->m_Jacobian[2][4] = r31 * w1;
//    this->m_Jacobian[2][5] = r32 * w2;
//    this->m_Jacobian[2][6] = r33 * w3;
//
//    // compute Jacobian wrt K1/K2/K3
//    this->m_Jacobian[0][7] = r11 * s1 * (x2 - c2);
//    this->m_Jacobian[0][8] = r11 * s1 * (x3 - c3);
//    this->m_Jacobian[0][9] = r12 * s2 * (x3 - c3);
//    this->m_Jacobian[1][7] = r21 * s1 * (x2 - c2);
//    this->m_Jacobian[1][8] = r21 * s1 * (x3 - c3);
//    this->m_Jacobian[1][9] = r22 * s2 * (x3 - c3);
//    this->m_Jacobian[2][7] = r31 * s1 * (x2 - c2);
//    this->m_Jacobian[2][8] = r31 * s1 * (x3 - c3);
//    this->m_Jacobian[2][9] = r32 * s2 * (x3 - c3);
//
//    // for test purpose
//    // FORCE ONLY DO DERIVATIVE ON S
//    //   for(int ii=0; ii <3; ii++){
//    //     for(int jj=0; jj<=9;jj++){
//    //       this->m_Jacobian[ii][jj] = 0.0;
//    //     }
//    //   }
//    //   this->m_Jacobian[0][4] = r11 * w1;
//    //   this->m_Jacobian[1][4] = r21 * w1;
//    //   this->m_Jacobian[2][4] = r31 * w1;
//
//    // compute derivatives for the translation part
//    unsigned int blockOffset = 10;
//    for (unsigned int dim = 0; dim < SpaceDimension; dim++) {
//        this->m_Jacobian[dim][blockOffset + dim] = 1.0;
//    }
//
//    return this->m_Jacobian;
//
//    //   // compute derivatives with respect to rotation
//    //   this->m_Jacobian.Fill(0.0);
//
//    //   const TScalarType x = p[0] - this->GetCenter()[0];
//    //   const TScalarType y = p[1] - this->GetCenter()[1];
//    //   const TScalarType z = p[2] - this->GetCenter()[2];
//
//    //   // compute Jacobian with respect to quaternion parameters
//    //   this->m_Jacobian[0][0] =   2.0 * (  m_Rotation.x() * x + m_Rotation.y() * y
//    //                               + m_Rotation.z() * z );
//    //   this->m_Jacobian[0][1] =   2.0 * (- m_Rotation.y() * x + m_Rotation.x() * y
//    //                               + m_Rotation.r() * z );
//    //   this->m_Jacobian[0][2] =   2.0 * (- m_Rotation.z() * x - m_Rotation.r() * y
//    //                               + m_Rotation.x() * z );
//    //   this->m_Jacobian[0][3] = - 2.0 * (- m_Rotation.r() * x + m_Rotation.z() * y
//    //                               - m_Rotation.y() * z );
//
//    //   this->m_Jacobian[1][0] = - this->m_Jacobian[0][1];
//    //   this->m_Jacobian[1][1] =   this->m_Jacobian[0][0];
//    //   this->m_Jacobian[1][2] =   this->m_Jacobian[0][3];
//    //   this->m_Jacobian[1][3] = - this->m_Jacobian[0][2];
//
//    //   this->m_Jacobian[2][0] = - this->m_Jacobian[0][2];
//    //   this->m_Jacobian[2][1] = - this->m_Jacobian[0][3];
//    //   this->m_Jacobian[2][2] =   this->m_Jacobian[0][0];
//    //   this->m_Jacobian[2][3] =   this->m_Jacobian[0][1];
//
//    //   // compute derivatives for the translation part
//    //   unsigned int blockOffset = 4;
//    //   for(unsigned int dim=0; dim < SpaceDimension; dim++ )
//    //     {
//    //     this->m_Jacobian[ dim ][ blockOffset + dim ] = 1.0;
//    //     }
//
//    //   return this->m_Jacobian;
//
// }

template <typename TScalarType>
void ANTSAffine3DTransform<TScalarType>
::ComputeJacobianWithRespectToParameters(
  const InputPointType  & p,
  JacobianType & j) const
{
  j.SetSize( this->GetOutputSpaceDimension(), this->GetNumberOfLocalParameters() );
  j.Fill(0.0);

  TScalarType c1 = this->GetCenter()[0];
  TScalarType c2 = this->GetCenter()[1];
  TScalarType c3 = this->GetCenter()[2];
  TScalarType s1 = this->m_S1;
  TScalarType s2 = this->m_S2;
  TScalarType s3 = this->m_S3;
  TScalarType k1 = this->m_K1;
  TScalarType k2 = this->m_K2;
  TScalarType k3 = this->m_K3;
  TScalarType x1 = p[0];
  TScalarType x2 = p[1];
  TScalarType x3 = p[2];

  // z1,z2,z3 is the S*K*point p
  TScalarType w1 = (x1 - c1) + k1 * (x2 - c2) + k2 * (x3 - c3);
  TScalarType w2 = (x2 - c2) + k3 * (x3 - c3);
  TScalarType w3 = (x3 - c3);

  TScalarType z1 = s1 * w1;
  TScalarType z2 = s2 * w2;
  TScalarType z3 = s3 * w3;

  // compute Jacobian with respect to quaternion parameters
  j[0][0] = static_cast<TScalarType>( 2.0 )
    * (m_Rotation.x() * z1 + m_Rotation.y() * z2 + m_Rotation.z() * z3);
  j[0][1] =
   static_cast<TScalarType>( 2.0 )
    * (-m_Rotation.y() * z1 + m_Rotation.x() * z2
       + m_Rotation.r() * z3);
  j[0][2] =
    static_cast<TScalarType>( 2.0 )
    * (-m_Rotation.z() * z1 - m_Rotation.r() * z2
       + m_Rotation.x() * z3);
  j[0][3] =
    static_cast<TScalarType>( -2.0 )
    * (-m_Rotation.r() * z1 + m_Rotation.z() * z2
       - m_Rotation.y() * z3);

  j[1][0] = -j[0][1];
  j[1][1] = j[0][0];
  j[1][2] = j[0][3];
  j[1][3] = -j[0][2];

  j[2][0] = -j[0][2];
  j[2][1] = -j[0][3];
  j[2][2] = j[0][0];
  j[2][3] = j[0][1];

  // get rotation matrix first
  // this is done to compensate for the transposed representation
  // between VNL and ITK
  VnlQuaternionType conjugateRotation = m_Rotation.conjugate();
  MatrixType        newMatrix;
  newMatrix = conjugateRotation.rotation_matrix_transpose();

  TScalarType r11 = newMatrix[0][0];
  TScalarType r12 = newMatrix[0][1];
  TScalarType r13 = newMatrix[0][2];
  TScalarType r21 = newMatrix[1][0];
  TScalarType r22 = newMatrix[1][1];
  TScalarType r23 = newMatrix[1][2];
  TScalarType r31 = newMatrix[2][0];
  TScalarType r32 = newMatrix[2][1];
  TScalarType r33 = newMatrix[2][2];

  // compute Jacobian wrt S1/S2/S3
  j[0][4] = r11 * w1;
  j[0][5] = r12 * w2;
  j[0][6] = r13 * w3;
  j[1][4] = r21 * w1;
  j[1][5] = r22 * w2;
  j[1][6] = r23 * w3;
  j[2][4] = r31 * w1;
  j[2][5] = r32 * w2;
  j[2][6] = r33 * w3;

  // compute Jacobian wrt K1/K2/K3
  j[0][7] = r11 * s1 * (x2 - c2);
  j[0][8] = r11 * s1 * (x3 - c3);
  j[0][9] = r12 * s2 * (x3 - c3);
  j[1][7] = r21 * s1 * (x2 - c2);
  j[1][8] = r21 * s1 * (x3 - c3);
  j[1][9] = r22 * s2 * (x3 - c3);
  j[2][7] = r31 * s1 * (x2 - c2);
  j[2][8] = r31 * s1 * (x3 - c3);
  j[2][9] = r32 * s2 * (x3 - c3);

  // for test purpose
  // FORCE ONLY DO DERIVATIVE ON S
  //   for(int ii=0; ii <3; ii++){
  //     for(int jj=0; jj<=9;jj++){
  //       j[ii][jj] = 0.0;
  //     }
  //   }
  //   j[0][4] = r11 * w1;
  //   j[1][4] = r21 * w1;
  //   j[2][4] = r31 * w1;

  // compute derivatives for the translation part
  unsigned int blockOffset = 10;
  for( unsigned int dim = 0; dim < SpaceDimension; dim++ )
    {
    j[dim][blockOffset + dim] = 1.0;
    }

  //   // compute derivatives with respect to rotation
  //   j.Fill(0.0);

  //   const TScalarType x = p[0] - this->GetCenter()[0];
  //   const TScalarType y = p[1] - this->GetCenter()[1];
  //   const TScalarType z = p[2] - this->GetCenter()[2];

  //   // compute Jacobian with respect to quaternion parameters
  //   j[0][0] =   2.0 * (  m_Rotation.x() * x + m_Rotation.y() * y
  //                               + m_Rotation.z() * z );
  //   j[0][1] =   2.0 * (- m_Rotation.y() * x + m_Rotation.x() * y
  //                               + m_Rotation.r() * z );
  //   j[0][2] =   2.0 * (- m_Rotation.z() * x - m_Rotation.r() * y
  //                               + m_Rotation.x() * z );
  //   j[0][3] = - 2.0 * (- m_Rotation.r() * x + m_Rotation.z() * y
  //                               - m_Rotation.y() * z );

  //   j[1][0] = - j[0][1];
  //   j[1][1] =   j[0][0];
  //   j[1][2] =   j[0][3];
  //   j[1][3] = - j[0][2];

  //   j[2][0] = - j[0][2];
  //   j[2][1] = - j[0][3];
  //   j[2][2] =   j[0][0];
  //   j[2][3] =   j[0][1];

  //   // compute derivatives for the translation part
  //   unsigned int blockOffset = 4;
  //   for(unsigned int dim=0; dim < SpaceDimension; dim++ )
  //     {
  //     j[ dim ][ blockOffset + dim ] = 1.0;
  //     }

  //   return j;
}

// template<typename TScalarType>
// const typename ANTSAffine3DTransform< TScalarType >::InverseMatrixType &
// ANTSAffine3DTransform<TScalarType>::
// GetInverseMatrix() const
// {
//   // If the transform has been modified we recompute the inverse
//   if(this->InverseMatrixIsOld())
//     {
//     InverseMatrixType newMatrix;
//     VnlQuaternionType conjugateRotation = m_Rotation.conjugate();
//     VnlQuaternionType inverseRotation = conjugateRotation.inverse();
//     newMatrix = inverseRotation.rotation_matrix_transpose();
//     this->SetVarInverseMatrix(newMatrix);
//     }
//   return this->GetVarInverseMatrix();
// }

template <typename TScalarType>
typename ANTSAffine3DTransform<TScalarType>::MatrixType ANTSAffine3DTransform<
  TScalarType>::ComputeMyRotationMatrix()
{
  VnlQuaternionType conjugateRotation = m_Rotation.conjugate();
  // this is done to compensate for the transposed representation
  // between VNL and ITK
  MatrixType R;

  R = conjugateRotation.rotation_matrix_transpose();
  return R;
}

template <typename TScalarType>
void ANTSAffine3DTransform<TScalarType>::ComputeMatrix()
{
  VnlQuaternionType conjugateRotation = m_Rotation.conjugate();
  // this is done to compensate for the transposed representation
  // between VNL and ITK
  MatrixType R;

  R = conjugateRotation.rotation_matrix_transpose();

  MatrixType S;
  S.Fill(NumericTraits<TScalarType>::ZeroValue());
  S[0][0] = this->m_S1;
  S[1][1] = this->m_S2;
  S[2][2] = this->m_S3;

  MatrixType K;
  K.Fill(NumericTraits<TScalarType>::ZeroValue());
  K[0][0] = NumericTraits<TScalarType>::OneValue();
  K[0][1] = this->m_K1;
  K[0][2] = this->m_K2;
  K[1][1] = NumericTraits<TScalarType>::OneValue();
  K[1][2] = this->m_K3;
  K[2][2] = NumericTraits<TScalarType>::OneValue();

  MatrixType newMatrix;

  newMatrix = R * S * K;

  this->SetVarMatrix(newMatrix);
}

template <typename TScalarType>
void ANTSAffine3DTransform<TScalarType>::ComputeMatrixParameters()
{
//   VnlQuaternionType quat(this->GetMatrix().GetVnlMatrix());
//   m_Rotation = quat;

//     std::cout << "compute para: to be done!" << std::endl;

//     InternalMatrixType A, Q, R;

  typedef vnl_matrix_fixed<TScalarType, 3U, 3U> TMatrix;

  TMatrix A, R, Q;

  A = this->GetMatrix().GetVnlMatrix();
  vnl_qr<ScalarType> myqr(A);

  Q = myqr.Q();   // Q() is the rotation
  R = myqr.R();   // R() is the upper triangluar

  // songgang: anyone added this???
  //      this is not necessary, for the mirror case
  //      the scale factor could not negative
  //      normalize R
  // need to run this, otherwise, identity matrix have negative scale!!!

//     // force diagnoal of rotation matrix to be positive
//     TMatrix dq(3,3,0);
//     for(unsigned i=0;i<3;i++){
//         dq(i,i) = (R(i,i)>=0)? 1 : -1;
//     }
//     Q = Q * dq;
//     R = dq * R;

  // force trace of rotation matrix be maximum possible by multiplying
  // a diagonal (+/-1 +/-1 +/-1) whose determinant is always positive +1.

  double trA = Q(0, 0) + Q(1, 1) + Q(2, 2);  // 1, 1, 1
  double trB = Q(0, 0) - Q(1, 1) - Q(2, 2);  // 1, -1, -1
  double trC = -Q(0, 0) + Q(1, 1) - Q(2, 2); // -1, 1, -1
  double trD = -Q(0, 0) - Q(1, 1) + Q(2, 2); // -1, -1, 1
  double maxTr = trA;                        // find the maximum of all terms;
  if( trB > maxTr )
    {
    maxTr = trB;     // dividing by the maximum makes
    }
  if( trC > maxTr )
    {
    maxTr = trC;     // the computations more stable
    }
  if( trD > maxTr )
    {
    maxTr = trD;     // and avoid division by zero
    }
  if( itk::Math::FloatAlmostEqual( maxTr, trB ) )
    {
    TMatrix dq;
    dq(0, 0) = 1;
    dq(1, 1) = -1;
    dq(2, 2) = -1;
    Q = Q * dq;
    R = dq * R;
    }

  if( itk::Math::FloatAlmostEqual( maxTr, trC ) )
    {
    TMatrix dq;
    dq(0, 0) = -1;
    dq(1, 1) = 1;
    dq(2, 2) = -1;
    Q = Q * dq;
    R = dq * R;
    }

  if( itk::Math::FloatAlmostEqual( maxTr, trD ) )
    {
    TMatrix dq;
    dq(0, 0) = -1;
    dq(1, 1) = -1;
    dq(2, 2) = 1;
    Q = Q * dq;
    R = dq * R;
    }

  double tr = 1 + Q(0, 0) + Q(1, 1) + Q(2, 2);
  double s, r, u, v, w;

  if( tr > 0 )
    {
    s = 0.5 / sqrt(tr);
    r = 0.25 / s;
    u = static_cast<double>(Q(2, 1) - Q(1, 2) ) * s;
    v = static_cast<double>(Q(0, 2) - Q(2, 0) ) * s;
    w = static_cast<double>(Q(1, 0) - Q(0, 1) ) * s;
    }
  else if( Q(0, 0) > Q(1, 1) && Q(0, 0) > Q(2, 2) )
    {
    s = 2 * sqrt(1 + Q(0, 0) - Q(1, 1) - Q(2, 2) );
    u = 0.25 * s;
    v = static_cast<double>(Q(0, 1) + Q(1, 0) ) / s;
    w = static_cast<double>(Q(0, 2) + Q(2, 0) ) / s;
    r = static_cast<double>(Q(1, 2) - Q(2, 1) ) / s;
    }
  else if( Q(0, 0) > Q(1, 1) )
    {
    s = 2 * sqrt(1 + Q(1, 1) - Q(0, 0) - Q(2, 2) );
    u = static_cast<double>( Q(0, 1) + Q(1, 0) ) / s;
    v = 0.25 * s;
    w = static_cast<double>(Q(1, 2) + Q(2, 1) ) / s;
    r = static_cast<double>(Q(0, 2) - Q(2, 0) ) / s;
    }
  else
    {
    s = 2 * sqrt(1 + Q(2, 2) - Q(0, 0) - Q(1, 1) );
    u = static_cast<double>(Q(0, 2) + Q(2, 0) ) / s;
    v = static_cast<double>(Q(1, 2) + Q(2, 1) ) / s;
    w = 0.25 * s;
    r = static_cast<double>(Q(0, 1) - Q(1, 0) ) / s;
    }

  std::cout << "A=" << A << std::endl;
  std::cout << "rotation R" << Q << std::endl;
  std::cout << "upper R" << R << std::endl;
  std::cout << "s=" << s << " u=" << u << " v=" << v << " w" << w << " r="
                   << r << std::endl;

  m_Rotation = VnlQuaternionType(u, v, w, r);

  std::cout << "m_Rotation from vnl" << VnlQuaternionType(u, v, w, r)
                   << std::endl;

  m_S1 = R(0, 0);
  m_S2 = R(1, 1);
  m_S3 = R(2, 2);

  m_K1 = R(0, 1) / R(0, 0);
  m_K2 = R(0, 2) / R(0, 0);
  m_K3 = R(1, 2) / R(1, 1);

  // std::cout << "before: this->GetMatrix(): " << this->GetMatrix();

  this->ComputeMatrix();

  // std::cout << "after: this->GetMatrix(): " << this->GetMatrix();

//     std::cout << "A=" << A << std::endl;
//     std::cout << "R=" << R << std::endl;
//     std::cout << "R=" << R << std::endl;
//     std::cout << "dq=" << dq << std::endl;
}
} // namespace

#endif
