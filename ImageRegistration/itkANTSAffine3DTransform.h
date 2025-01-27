#ifndef __itkANTSAffine3DTransform_h
#define __itkANTSAffine3DTransform_h

#include <iostream>
#include "itkRigid3DTransform.h"
#include "vnl/vnl_quaternion.h"

namespace itk
{
/** \brief ANTSAffine3DTransform of a vector space (e.g. space coordinates).
 *
 * This transform applies a rotation and translation to the space
 *
 * \ingroup Transforms
 */
template <typename TScalarType = double>
// Data type for scalars (float or double)
class ANTSAffine3DTransform final : public MatrixOffsetTransformBase<TScalarType, 3, 3>
//        public Rigid3DTransform< TScalarType >
{
public:
  /** Standard class typedefs.   */
  typedef ANTSAffine3DTransform Self;
  //  typedef Rigid3DTransform< TScalarType >     Superclass;
  typedef MatrixOffsetTransformBase<TScalarType, 3, 3> Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

  /** Run-time type information (and related methods).   */
  // itkOverrideGetNameOfClassMacro( ANTSAffine3DTransform);
  itkOverrideGetNameOfClassMacro(ANTSAffine3DTransform);

  /** Dimension of parameters   */
  static constexpr unsigned int InputSpaceDimension = 3;
  static constexpr unsigned int OutputSpaceDimension = 3;
  static constexpr unsigned int SpaceDimension = 3;
  static constexpr unsigned int ParametersDimension = 13;

  /** Parameters Type   */
  typedef typename Superclass::ParametersType            ParametersType;
  typedef typename Superclass::JacobianType              JacobianType;
  typedef typename Superclass::ScalarType                ScalarType;
  typedef typename Superclass::InputPointType            InputPointType;
  typedef typename Superclass::OutputPointType           OutputPointType;
  typedef typename Superclass::InputVectorType           InputVectorType;
  typedef typename Superclass::OutputVectorType          OutputVectorType;
  typedef typename Superclass::InputVnlVectorType        InputVnlVectorType;
  typedef typename Superclass::OutputVnlVectorType       OutputVnlVectorType;
  typedef typename Superclass::InputCovariantVectorType  InputCovariantVectorType;
  typedef typename Superclass::OutputCovariantVectorType OutputCovariantVectorType;
  typedef typename Superclass::MatrixType                MatrixType;
  typedef typename Superclass::InverseMatrixType         InverseMatrixType;
  typedef typename Superclass::CenterType                CenterType;
  typedef typename Superclass::OffsetType                OffsetType;
  typedef typename Superclass::TranslationType           TranslationType;

  /** VnlQuaternion type.  */
  typedef vnl_quaternion<TScalarType> VnlQuaternionType;

  /** Compute the Jacobian Matrix of the transformation at one point */
  /** Set the rotation of the rigid transform.
   * This method sets the rotation of a ANTSAffine3DTransform to a
   * value specified by the user. */
  void
  SetRotation(const VnlQuaternionType & rotation);

  void
  SetS1(const TScalarType S1);

  void
  SetS2(const TScalarType S2);

  void
  SetS3(const TScalarType S3);

  void
  SetK1(const TScalarType K1);

  void
  SetK2(const TScalarType K2);

  void
  SetK3(const TScalarType K3);

  /** Get the rotation from an ANTSAffine3DTransform.
   * This method returns the value of the rotation of the
   * ANTSAffine3DTransform.   **/
  const VnlQuaternionType &
  GetRotation() const
  {
    return m_Rotation;
  }

  MatrixType
  ComputeMyRotationMatrix();

  itkGetConstReferenceMacro(S1, TScalarType);
  itkGetConstReferenceMacro(S2, TScalarType);
  itkGetConstReferenceMacro(S3, TScalarType);
  itkGetConstReferenceMacro(K1, TScalarType);
  itkGetConstReferenceMacro(K2, TScalarType);
  itkGetConstReferenceMacro(K3, TScalarType);

  /** Set the parameters to the IdentityTransform */
  void
  SetIdentity() override;

  /** Set the transformation from a container of parameters.
   * This is typically used by optimizers.
   * There are 7 parameters. The first four represents the
   * quaternion and the last three represents the
   * offset. */
  void
  SetParameters(const ParametersType & parameters) override;

  const ParametersType &
  GetParameters() const override;

  //  /** Compute the Jacobian of the transformation.
  //   * This method computes the Jacobian matrix of the transformation.
  //   * given point or vector, returning the transformed point or
  //   * vector. The rank of the Jacobian will also indicate if the transform
  //   * is invertible at this point. */
  //  const JacobianType & GetJacobian(const InputPointType  &point ) const;

  /** Compute the Jacobian of the transformation
   *
   * This method computes the Jacobian matrix of the transformation.
   * given point or vector, returning the transformed point or
   * vector. The rank of the Jacobian will also indicate if the transform
   * is invertible at this point.
   * Get local Jacobian for the given point
   * \c j will sized properly as needed.
   */
  void
  ComputeJacobianWithRespectToParameters(const InputPointType & p, JacobianType & j) const override;

protected:
  /*   ANTSAffine3DTransform(const MatrixType &matrix, */
  /*                            const OutputVectorType &offset); */
  ANTSAffine3DTransform(unsigned int outputDims, unsigned int paramDims);
  ANTSAffine3DTransform();
  ~ANTSAffine3DTransform() override = default;

  void
  ComputeMatrix() override;

  void
  ComputeMatrixParameters() override;

  void
  SetVarRotation(const VnlQuaternionType & rotation)
  {
    m_Rotation = rotation;
  };

  //  const InverseMatrixType & GetInverseMatrix() const;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  ANTSAffine3DTransform(const Self &) = delete;
  void
  operator=(const Self &) = delete;

  /** Rotation of the transformation. */
  VnlQuaternionType m_Rotation;

  /** added affine parameters **/
  TScalarType m_S1;
  TScalarType m_S2;
  TScalarType m_S3;
  TScalarType m_K1;
  TScalarType m_K2;
  TScalarType m_K3;
}; // class ANTSAffine3DTransform
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkANTSAffine3DTransform.hxx"
#endif

#endif /* __itkANTSAffine3DTransform_h */
