#ifndef antsNearestRotation_h
#define antsNearestRotation_h

// Nearest-rotation (special-orthogonal polar factor) helper used by the
// principal-axis initializers (antsAffineInitializer, antsAI, ImageMath) to
// solve Wahba's problem. Prefers the Eigen-backed itk::Math::SVD when the ITK in
// use provides it, otherwise vnl_svd; the result is backend-invariant either way.
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/algo/vnl_svd.h"
#include "vnl/algo/vnl_svd_fixed.h"
#include "vnl/algo/vnl_determinant.h"

#if __has_include(<itkMathSVD.h>)
#  include "itkMathSVD.h"
#  ifndef ANTS_HAS_ITK_MATH_SVD
#    define ANTS_HAS_ITK_MATH_SVD 1
#  endif
#endif

namespace ants
{
// Nearest rotation (special-orthogonal polar factor) of A: the R in SO(n)
// minimizing ||R - A||_F. With A = U diag(W) V^T it is R = U diag(1,..,1,d) V^T
// where d = det(U V^T) is +/-1, forcing det(R) = +1. Unlike a raw U V^T (which may
// be a reflection) or an ad-hoc diagonal-sign fix, this is the true nearest
// rotation and is invariant to the singular-vector sign convention, so it is
// identical across SVD backends.
template <typename TMatrix, typename T>
TMatrix
NearestRotationImpl(const TMatrix & U, const TMatrix & V)
{
  TMatrix            rotation = U * V.transpose();
  const T            d = static_cast<T>(vnl_determinant(rotation));
  const unsigned int last = U.cols() - 1;
  TMatrix            scaledU = U;
  for (unsigned int i = 0; i < U.rows(); ++i)
  {
    scaledU(i, last) *= d;
  }
  return scaledU * V.transpose();
}

template <typename T>
vnl_matrix<T>
NearestRotation(const vnl_matrix<T> & A)
{
#if defined(ANTS_HAS_ITK_MATH_SVD)
  const auto s = itk::Math::SVD(A);
  return NearestRotationImpl<vnl_matrix<T>, T>(s.U, s.V);
#else
  vnl_svd<T> s(A);
  return NearestRotationImpl<vnl_matrix<T>, T>(s.U(), s.V());
#endif
}

template <typename T, unsigned int VDim>
vnl_matrix_fixed<T, VDim, VDim>
NearestRotation(const vnl_matrix_fixed<T, VDim, VDim> & A)
{
  using MatrixType = vnl_matrix_fixed<T, VDim, VDim>;
#if defined(ANTS_HAS_ITK_MATH_SVD)
  const auto s = itk::Math::SVD(A);
  return NearestRotationImpl<MatrixType, T>(s.U, s.V);
#else
  vnl_svd_fixed<T, VDim, VDim> s(A);
  return NearestRotationImpl<MatrixType, T>(s.U(), s.V());
#endif
}
} // namespace ants

#endif // antsNearestRotation_h
