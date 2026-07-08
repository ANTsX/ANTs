/*=========================================================================
 *
 *  Regression test for the Wahba principal-axis rotation used by
 *  antsAffineInitializer / antsAI. Those tools build an attitude matrix
 *  B = sum_k outer(moving_axis_k, fixed_axis_k) and recover the aligning
 *  rotation from its SVD. B is rank-deficient by construction, so the proper
 *  special-orthogonal polar factor (ants::NearestRotation) must be used; an
 *  ad-hoc "flip negative diagonal" correction fails ~48% of the time. This
 *  test constructs consistent problems with a known ground-truth rotation and
 *  asserts exact recovery, orthogonality and det = +1. It is fixture-free and
 *  deterministic, and exercises the vnl_svd path ANTs builds against.
 *
 *=========================================================================*/
#include "antsNearestRotation.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "vnl/algo/vnl_determinant.h"
#include <cmath>
#include <iostream>
#include <cstdlib>

namespace
{
using T = double;

vnl_matrix<T>
identity(unsigned int n)
{
  vnl_matrix<T> I(n, n, 0.0);
  I.set_identity();
  return I;
}

// 3-D rotation about a unit axis by angle (Rodrigues' formula).
vnl_matrix<T>
rotation3D(T ax, T ay, T az, T angle)
{
  vnl_vector<T> u(3);
  u[0] = ax;
  u[1] = ay;
  u[2] = az;
  u.normalize();
  vnl_matrix<T> K(3, 3, 0.0);
  K(0, 1) = -u[2];
  K(0, 2) = u[1];
  K(1, 0) = u[2];
  K(1, 2) = -u[0];
  K(2, 0) = -u[1];
  K(2, 1) = u[0];
  return identity(3) + K * std::sin(angle) + K * K * (1.0 - std::cos(angle));
}

vnl_matrix<T>
rotation2D(T angle)
{
  vnl_matrix<T> R(2, 2);
  R(0, 0) = std::cos(angle);
  R(0, 1) = -std::sin(angle);
  R(1, 0) = std::sin(angle);
  R(1, 1) = std::cos(angle);
  return R;
}

// Build the ANTs Wahba attitude matrix for fixed axes = the first (n-1) standard
// basis vectors and moving axes = Rtrue * fixed, then assert NearestRotation
// recovers Rtrue and is a proper rotation.
bool
checkRecovers(const vnl_matrix<T> & Rtrue, const std::string & name)
{
  const unsigned int n = Rtrue.rows();
  vnl_matrix<T>      B(n, n, 0.0);
  for (unsigned int k = 0; k < n - 1; ++k) // (n-1) axis pairs, as antsAffineInitializer builds
  {
    vnl_vector<T> fixedAxis(n, 0.0);
    fixedAxis[k] = 1.0;
    const vnl_vector<T> movingAxis = Rtrue * fixedAxis;
    B += outer_product(movingAxis, fixedAxis);
  }

  const vnl_matrix<T> R = ants::NearestRotation(B);

  const T recoverError = (R - Rtrue).frobenius_norm();
  const T orthoError = (R.transpose() * R - identity(n)).frobenius_norm();
  const T detError = std::abs(vnl_determinant(R) - 1.0);

  const bool ok = (recoverError < 1e-9) && (orthoError < 1e-9) && (detError < 1e-9);
  std::cout << "  " << (ok ? "PASS" : "FAIL") << "  " << name << "  |R-Rtrue|=" << recoverError
            << "  |R^T R - I|=" << orthoError << "  |det-1|=" << detError << "\n";
  return ok;
}
} // namespace

int
main()
{
  unsigned int failures = 0;
  std::cout << std::scientific;

  // 2-D: rank-1 attitude matrix (single axis pair).
  for (const T deg : { 15.0, 30.0, 45.0, 80.0, 170.0 })
  {
    const T angle = deg * 3.14159265358979323846 / 180.0;
    if (!checkRecovers(rotation2D(angle), "2D rot " + std::to_string(static_cast<int>(deg)) + "deg"))
    {
      ++failures;
    }
  }

  // 3-D: rank-2 attitude matrix (two axis pairs) about several axes/angles.
  const double cases[][4] = { { 0, 0, 1, 30 }, { 1, 0, 0, 45 },     { 0, 1, 0, 60 },
                              { 1, 1, 1, 35 }, { 1, -2, 0.5, 100 }, { -1, 0.3, 2, 155 } };
  for (const auto & c : cases)
  {
    const T angle = c[3] * 3.14159265358979323846 / 180.0;
    if (!checkRecovers(rotation3D(c[0], c[1], c[2], angle),
                       "3D axis(" + std::to_string(c[0]) + "," + std::to_string(c[1]) + "," + std::to_string(c[2]) +
                         ") " + std::to_string(static_cast<int>(c[3])) + "deg"))
    {
      ++failures;
    }
  }

  std::cout << (failures ? "\nFAILED (" + std::to_string(failures) + ")\n" : "\nAll Wahba rotations recovered.\n");
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
