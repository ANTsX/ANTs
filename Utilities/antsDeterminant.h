#ifndef antsDeterminant_h
#define antsDeterminant_h

// ITK 6 deprecates vnl/algo/vnl_determinant.h in favor of the Eigen-backed
// itk::bridge::Math::Determinant; ITK 5 has only the former.
#if __has_include("itkBridgeMathDeterminant.h")
#  include "itkBridgeMathDeterminant.h"
#else
#  include "vnl/algo/vnl_determinant.h"
#endif

namespace ants
{
template <typename TMatrix>
auto
Determinant(const TMatrix & m)
{
#if __has_include("itkBridgeMathDeterminant.h")
  return itk::bridge::Math::Determinant(m);
#else
  return vnl_determinant(m);
#endif
}
} // namespace ants

#endif // antsDeterminant_h
