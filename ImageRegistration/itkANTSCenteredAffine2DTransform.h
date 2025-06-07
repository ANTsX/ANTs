#ifndef __itkANTSCenteredAffine2DTransform_h
#define __itkANTSCenteredAffine2DTransform_h

#include <iostream>
#include "itkMatrixOffsetTransformBase.h"
#include "itkMacro.h"

namespace itk
{
template <typename TScalarType = double>
// Data type for scalars (float or double)
// class Rigid2DTransform :
class ANTSCenteredAffine2DTransform final
  : public MatrixOffsetTransformBase<TScalarType, 2, 2> // Dimensions of input and output spaces
{
public:
  /** Standard class typedefs. */
  //  typedef Rigid2DTransform Self;
  typedef ANTSCenteredAffine2DTransform                Self;
  typedef MatrixOffsetTransformBase<TScalarType, 2, 2> Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** Run-time type information (and related methods). */
  // itkOverrideGetNameOfClassMacro( Rigid2DTransform);
  itkOverrideGetNameOfClassMacro(ANTSCenteredAffine2DTransform);

  /** New macro for creation of through a Smart Pointer */
  itkNewMacro(Self);

  /** Dimension of the space. */
  static constexpr unsigned int InputSpaceDimension = 2;
  static constexpr unsigned int OutputSpaceDimension = 2;
  // static constexpr unsigned int ParametersDimension = 3;
  static constexpr unsigned int ParametersDimension = 8;

  /** Scalar type. */
  typedef typename Superclass::ScalarType ScalarType;

  /** Parameters type. */
  typedef typename Superclass::ParametersType ParametersType;

  /** Jacobian type. */
  typedef typename Superclass::JacobianType JacobianType;

  // / Standard matrix type for this class
  typedef typename Superclass::MatrixType MatrixType;

  // / Standard vector type for this class
  typedef typename Superclass::OffsetType OffsetType;

  // / Standard vector type for this class
  typedef typename Superclass::InputVectorType  InputVectorType;
  typedef typename Superclass::OutputVectorType OutputVectorType;

  // / Standard covariant vector type for this class
  typedef typename Superclass::InputCovariantVectorType  InputCovariantVectorType;
  typedef typename Superclass::OutputCovariantVectorType OutputCovariantVectorType;

  // / Standard vnl_vector type for this class
  typedef typename Superclass::InputVnlVectorType  InputVnlVectorType;
  typedef typename Superclass::OutputVnlVectorType OutputVnlVectorType;

  // / Standard coordinate point type for this class
  typedef typename Superclass::InputPointType  InputPointType;
  typedef typename Superclass::OutputPointType OutputPointType;

  /**
   * Set the rotation Matrix of a Rigid2D Transform
   *
   * This method sets the 2x2 matrix representing the rotation
   * in the transform.  The Matrix is expected to be orthogonal
   * with a certain tolerance.
   *
   * \warning This method will throw an exception is the matrix
   * provided as argument is not orthogonal.
   *
   * \sa MatrixOffsetTransformBase::SetMatrix()
   **/
  // virtual void SetMatrix( const MatrixType & matrix );

  /**
   * Set/Get the rotation matrix. These methods are old and are
   * retained for backward compatibility. Instead, use SetMatrix()
   * GetMatrix().
   **/
  // virtual void SetRotationMatrix(const MatrixType &matrix)
  //    { this->SetMatrix( matrix ); }
  //   const MatrixType & GetRotationMatrix() const
  //    { return this->GetMatrix(); }

  /**
   * Compose the transformation with a translation
   *
   * This method modifies self to include a translation of the
   * origin.  The translation is precomposed with self if pre is
   * true, and postcomposed otherwise.
   **/
  // void Translate(const OffsetType &offset, bool pre=false);

  /**
   * Back transform by an rigid transformation.
   *
   * The BackTransform() methods are slated to be removed from ITK.
   * Instead, please use GetInverse() or CloneInverseTo() to generate
   * an inverse transform and  then perform the transform using that
   * inverted transform.
   **/
  inline InputPointType
  BackTransform(const OutputPointType & point) const;

  inline InputVectorType
  BackTransform(const OutputVectorType & vector) const;

  inline InputVnlVectorType
  BackTransform(const OutputVnlVectorType & vector) const;

  inline InputCovariantVectorType
  BackTransform(const OutputCovariantVectorType & vector) const;

  /** Set/Get the angle of rotation in radians */
  void
  SetAngle(TScalarType angle);

  void
  SetS1(TScalarType S1);

  void
  SetS2(TScalarType S2);

  void
  SetK(TScalarType K);

  itkGetConstReferenceMacro(Angle, TScalarType);
  itkGetConstReferenceMacro(S1, TScalarType);
  itkGetConstReferenceMacro(S2, TScalarType);
  itkGetConstReferenceMacro(K, TScalarType);

  /** Set the angle of rotation in degrees. */
  void
  SetAngleInDegrees(TScalarType angle);

  /** Set/Get the angle of rotation in radians. These methods
   * are old and are retained for backward compatibility.
   * Instead, use SetAngle() and GetAngle(). */
  // void SetRotation(TScalarType angle)
  //  { this->SetAngle(angle); }
  // virtual const TScalarType & GetRotation() const
  //  { return m_Angle; }

  /** Set the transformation from a container of parameters
   * This is typically used by optimizers.
   * There are 3 parameters. The first one represents the
   * angle of rotation in radians and the last two represents the translation.
   * The center of rotation is fixed.
   *
   * \sa Transform::SetParameters()
   * \sa Transform::SetFixedParameters() */
  void
  SetParameters(const ParametersType & parameters) override;

  /** Get the parameters that uniquely define the transform
   * This is typically used by optimizers.
   * There are 3 parameters. The first one represents the
   * angle or rotation in radians and the last two represents the translation.
   * The center of rotation is fixed.
   *
   * \sa Transform::GetParameters()
   * \sa Transform::GetFixedParameters() */
  const ParametersType &
  GetParameters() const override;

  /** This method computes the Jacobian matrix of the transformation
   * at a given input point.
   *
   * \sa Transform::GetJacobian() */
  // const JacobianType & GetJacobian(const InputPointType  &point ) const;

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

  /**
   * This method creates and returns a new ANTSCenteredAffine2DTransform object
   * which is the inverse of self.
   **/
  // void CloneInverseTo( Pointer & newinverse ) const;

  /**
   * This method creates and returns a new ANTSCenteredAffine2DTransform object
   * which has the same parameters.
   **/
  void
  CloneTo(Pointer & clone) const;

  /** Reset the parameters to create and identity transform. */
  void
  SetIdentity() override;

protected:
  // Rigid2DTransform();
  ANTSCenteredAffine2DTransform();

  //  Rigid2DTransform( unsigned int outputSpaceDimension,
  //                    unsigned int parametersDimension);
  ANTSCenteredAffine2DTransform(unsigned int outputSpaceDimension, unsigned int parametersDimension);

  //  ~Rigid2DTransform();
  ~ANTSCenteredAffine2DTransform() override;

  /**
   * Print contents of an ANTSCenteredAffine2DTransform
   **/
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /** Compute the matrix from angle. This is used in Set methods
   * to update the underlying matrix whenever a transform parameter
   * is changed. */
  void
  ComputeMatrix() final;

  /** Compute the angle from the matrix. This is used to compute
   * transform parameters from a given matrix. This is used in
   * MatrixOffsetTransformBase::Compose() and
   * MatrixOffsetTransformBase::GetInverse(). */
  void
  ComputeMatrixParameters() override;

  /** Update angle without recomputation of other internal variables. */
  void
  SetVarAngle(TScalarType angle)
  {
    m_Angle = angle;
  }

  void
  SetVarS1(TScalarType S1)
  {
    m_S1 = S1;
  }

  void
  SetVarS2(TScalarType S2)
  {
    m_S2 = S2;
  }

  void
  SetVarK(TScalarType K)
  {
    m_K = K;
  }

private:
  ANTSCenteredAffine2DTransform(const Self &) = delete;
  void
  operator=(const Self &) = delete;

  TScalarType m_Angle;
  TScalarType m_S1;
  TScalarType m_S2;
  TScalarType m_K;
}; // class ANTSCenteredAffine2DTransform
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkANTSCenteredAffine2DTransform.hxx"
#endif

#endif /* __itkANTSCenteredAffine2DTransform_h */
