/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _SurfaceCurvatureBase_h
#define _SurfaceCurvatureBase_h

#include <vnl/vnl_matops.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_cross.h>
// #include <vnl/algo/vnl_rpoly_roots.h>
#include "itkObject.h"
#include "itkProcessObject.h"

#include "itkVectorContainer.h"
#include "itkCastImageFilter.h"

namespace itk
{
/** \class SurfaceCurvatureBase
 *
 * This class takes a surface as input and creates a local
 * geometric frame for each surface point.
 *
 * \note The origin of a neighborhood is always taken to be
 *       the first point entered into and the
 *       last point stored in the list.
 */
template <typename TSurface, unsigned int TDimension = 3>
class SurfaceCurvatureBase : public ProcessObject
{
public:
  /** Standard class typedefs. */
  using Self = SurfaceCurvatureBase<TSurface, TDimension>;
  using Superclass = ProcessObject;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(SurfaceCurvatureBase);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Image types. */
  using SurfaceType = TSurface;

  /** Image dimension. */
  //  static constexpr unsigned int ImageDimension = TImage::ImageDimension;
  static constexpr unsigned int ImageDimension = TDimension;
  static constexpr unsigned int SurfaceDimension = TDimension;

  using RealType = float;
  using VectorType = vnl_vector<RealType>;
  using FixedVectorType = vnl_vector_fixed<RealType, (Self::ImageDimension)>;
  using PointType = vnl_vector_fixed<RealType, (Self::ImageDimension)>;
  using MatrixType = vnl_matrix<double>;
  using PointContainerType = std::vector<PointType>;
  using FunctionContainerType = std::vector<float>;

  /** Set input parameter file */
  itkSetStringMacro(ParameterFileName);

  /** Set input parameter file */
  itkGetStringMacro(ParameterFileName);

  /** Fill the point list with the points neighboring the current origin */
  virtual void
  FindNeighborhood(unsigned int temp = 0);

  /** A Euclidean relative of L.D. Griffin's compactness.*/
  RealType
  ComputeMeanEuclideanDistance();

  /** */
  void
  ComputeAveragePoint();

  void ProjectToTangentPlane(PointType);

  void
  EigenDecomposition(MatrixType D);

  /** Estimate the plane tangent to a point using the neighborhood
   *   of the point.  This is, in general, an over-constrained least
   *   squares fitting problem.
   */
  void EstimateTangentPlane(PointType);
  void WeightedEstimateTangentPlane(PointType);
  void TestEstimateTangentPlane(PointType);

  /** This function sets the reference tangent arbitrarily.
   *   It can be overridden in case there is a better practical approach.
   */
  void
  ChooseReferenceTangent();

  /** This function fills the weight and angle vectors for
   * a given neighborhood.
   */
  virtual void
  ComputeWeightsAndDirectionalKappaAndAngles(PointType origin);

  /** */
  void
  ComputeFrame(PointType origin);

  /** */
  void
  ComputeFrameAndKappa(PointType origin);

  void
  ShimshoniFrame(PointType origin);

  /** */
  void
  ComputeJoshiFrame(PointType origin);

  /** */
  void
  EstimateCurvature(RealType w1 = 3. / 8., RealType w2 = 1. / 8., RealType w3 = 1. / 8., RealType w4 = 3. / 8.);

  /** Use the Besl and Jain analytical curvature computation.
   * We use a least square polynomial fit to the local neighborhood
   * to estimate the mean and gaussian curvature.
   */
  void JainMeanAndGaussianCurvature(PointType);

  /** This function returns a weight given a distance
   *  It may be the identity function, a normalization
   *  or a gaussianization of the input distance. */
  virtual RealType GetWeight(PointType, PointType);

  /** This function returns the angle between the reference tangent
      and the projection onto the tangent plane of the vector between
      the neighborhood focus and its neighbor. */
  virtual RealType
  GetTheta(PointType Q, PointType origin);

  /** Estimate the directional curvature using Shimshoni's method (eq 6).*/
  virtual void EstimateDirectionalCurvature(PointType, PointType);

  void
  PrintFrame();

  virtual void
  ComputeFrameOverDomain(unsigned int /* which */ = 3){};

  RealType
  ErrorEstimate(PointType, RealType sign = 1.0);

  unsigned int
  CharacterizeSurface();

  itkSetMacro(Origin, PointType);
  itkGetMacro(Origin, PointType);

  itkSetMacro(AveragePoint, PointType);
  itkGetMacro(AveragePoint, PointType);

  itkGetMacro(Normal, FixedVectorType);

  itkGetMacro(MeanKappa, RealType);
  itkGetMacro(Sigma, RealType);
  itkSetMacro(Sigma, RealType);

  itkGetMacro(UseGeodesicNeighborhood, bool);
  itkSetMacro(UseGeodesicNeighborhood, bool);

  /** Set normal estimates a 3D frame from a given normal */
  void SetFrameFromNormal(FixedVectorType);

  /** We use the cross-product of the tangents times the image spacing
      to get the local area. */
  RealType ComputeLocalArea(FixedVectorType);

  /** We estimate the integral as a sum, assuming the local
      area (from compute local area) scales the value of the
      function at the pixel.  See http://mathworld.wolfram.com/SurfaceIntegral.html*/
  virtual RealType
  IntegrateFunctionOverNeighborhood(bool /* norm */ = false)
  {
    return 0;
  }

  void
  SwitchNormalSign()
  {
    m_Normal *= (-1.0);
  }

  // for conjugate harmonic function
  float
  dstarUestimate();

  // get this from the local frame - 1st order vs 2nd order shape operator
  void
  EstimateMetricTensor();

protected:
  SurfaceCurvatureBase();
  ~SurfaceCurvatureBase() override = default;

  /** Holds the value of Pi. */
  double m_Pi;

  bool m_Debug;

  /** Data structures to contain points. */
  PointContainerType m_PointList;

  /** Data structures to contain function
      values associated with points. */
  FunctionContainerType m_FunctionValueList;

  /** This list contains the projection of the vectors onto
      the tangent plane (T_i Shimshoni). */
  PointContainerType m_TangentProjectionList;

  /** The point that is the origin of the current neighborhood. */
  PointType m_Origin;
  PointType m_AveragePoint;
  PointType m_PlaneVector;

  /** Data that represents single vectors */
  FixedVectorType m_ArbitraryTangent;
  FixedVectorType m_Normal;
  FixedVectorType m_Tangent1;
  FixedVectorType m_Tangent2;
  RealType        m_dX; // size in local x dir
  RealType        m_dY; // size in local y dir
  FixedVectorType m_MetricTensor;

  VectorType m_ThetaVector;
  VectorType m_WeightVector;
  VectorType m_DirectionalKappaVector;

  /** Data for representing neighborhoods and derived from the vector frames*/

  /** Approximate directional curvature */
  RealType m_DirectionalKappa;

  RealType m_MetricTensorDeterminant;
  /** Approximate principal curvature 1*/
  RealType m_Kappa1;
  /** Approximate principal curvature 2*/
  RealType m_Kappa2;
  RealType m_GaussianKappa;
  RealType m_MeanKappa;

  /** Solution weights eq. (7) (8) Shimshoni */
  RealType m_A;
  RealType m_B;
  RealType m_C;

  /** True Eigenvector weights */
  RealType m_W1;
  RealType m_W2;

  /** Eigenvalues */
  RealType m_Eval0;
  RealType m_Eval1;
  RealType m_Eval2;

  unsigned int m_CurrentNeighborhoodPointIndex;

  std::string m_ParameterFileName;

  /** We use this to avoid computing the frame in degenerate cases. */
  RealType m_TotalDKap;

  RealType m_TotalArea;

  RealType m_Sigma;

  bool m_UseGeodesicNeighborhood;

private:
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkSurfaceCurvatureBase.hxx"
#endif

#endif
