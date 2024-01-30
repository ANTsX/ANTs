/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _SurfaceCurvatureBase_hxx
#define _SurfaceCurvatureBase_hxx

#include <vcl_compiler.h>
#include <iostream>
#include <cmath>
#include <iostream>
#include <vnl/vnl_real_polynomial.h>
// #include <vnl/algo/vnl_rpoly_roots.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_vector_fixed.h>
#include <itkMath.h>


namespace itk
{
template <typename TSurface, unsigned int TDimension>
SurfaceCurvatureBase<TSurface, TDimension>::SurfaceCurvatureBase()
{
  m_Origin.fill(0.0);
  m_ArbitraryTangent.fill(0.);
  m_Normal.fill(0.);
  m_Tangent1.fill(0.);
  m_Tangent2.fill(0.);

  m_DirectionalKappa = 0.0;
  m_Kappa1 = 0.0;
  m_Kappa2 = 0.0;
  m_GaussianKappa = 0.0;
  m_A = 0.0;
  m_B = 0.0;
  m_C = 0.0;
  m_W1 = 3. / 8.;
  m_W2 = 1. / 8.;
  m_Eval0 = 0.0;
  m_Eval1 = 0.0;
  m_Eval2 = 0.0;
  m_CurrentNeighborhoodPointIndex = 0;
  m_ParameterFileName = "";

  m_Pi = 3.14159265358979323846;
  m_Debug = true;
  m_Debug = false;

  m_UseGeodesicNeighborhood = false;

  m_TotalArea = 0.0;
  m_Sigma = 1.0F;
}

template <typename TSurface, unsigned int TDimension>
void
SurfaceCurvatureBase<TSurface, TDimension>::ComputeAveragePoint()
{
  m_AveragePoint.fill(0.0);
  unsigned int npts = m_PointList.size();
  unsigned int dim = SurfaceDimension;
  for (unsigned int i = 0; i < npts; i++)
  {
    for (unsigned int j = 0; j < dim; j++)
    {
      m_AveragePoint[j] += m_PointList[i][j];
    }
  }

  m_AveragePoint /= (RealType)npts;
}

template <typename TSurface, unsigned int TDimension>
void
SurfaceCurvatureBase<TSurface, TDimension>::ProjectToTangentPlane(PointType dif)
{
  MatrixType id(ImageDimension, ImageDimension);

  id.fill(0.0);
  unsigned int i = 0, j = 0;

  MatrixType NN(ImageDimension, ImageDimension);
  for (i = 0; i < ImageDimension; i++)
  {
    id(i, i) = 1.0;
    m_PlaneVector[i] = 0.0;
    for (j = 0; j < ImageDimension; j++)
    {
      NN(i, j) = id(i, j) - m_Normal(i) * m_Normal(j);
      m_PlaneVector[i] += NN(i, j) * dif(j);
    }
  }
}

template <typename TSurface, unsigned int TDimension>
void
SurfaceCurvatureBase<TSurface, TDimension>::SetFrameFromNormal(FixedVectorType N)
{
  // WORKS ONLY FOR 3D!

  RealType mag = N.magnitude();

  if (!itk::Math::FloatAlmostEqual(mag, itk::NumericTraits<RealType>::ZeroValue()))
  {
    m_Normal = N / mag;
  }
  else
  {
    m_Tangent1.fill(0.0);
    m_Tangent2.fill(0.0);
    return;
  }

  // now estimate the tangent plane
  if (fabs(m_Normal[0]) > 1.e-1)
  {
    float norm = itk::NumericTraits<float>::OneValue() / m_Normal[0];
    m_Tangent1[1] = norm * N[0];
    m_Tangent1[0] = -1.0f * norm * (N[1] + N[2]);
    m_Tangent1[2] = norm * N[0];
  }
  else if (fabs(m_Normal[1]) > 1.e-1)
  {
    float norm = itk::NumericTraits<float>::OneValue() / m_Normal[1];
    m_Tangent1[0] = norm * N[1];
    m_Tangent1[1] = -1.0f * norm * (N[0] + N[2]);
    m_Tangent1[2] = norm * N[1];
  }
  else if (fabs(m_Normal[2]) > 1.e-1)
  {
    float norm = itk::NumericTraits<float>::OneValue() / m_Normal[2];
    m_Tangent1[0] = norm * N[2];
    m_Tangent1[2] = -1.0f * norm * (N[0] + N[1]);
    m_Tangent1[1] = norm * N[2];
  }
  m_Tangent1 /= m_Tangent1.magnitude();

  m_Tangent2 = vnl_cross_3d(m_Normal, m_Tangent1);
  m_Tangent2 /= m_Tangent2.magnitude();

  this->ChooseReferenceTangent();
}

template <typename TSurface, unsigned int TDimension>
void
SurfaceCurvatureBase<TSurface, TDimension>::EigenDecomposition(MatrixType D)
{
  // Compute estimated frame using eigensystem of D'*D
  {
    vnl_symmetric_eigensystem<double> eig(D.transpose() * D);
    for (unsigned int j = 0; j < SurfaceDimension; j++)
    {
      m_Normal[j] = eig.get_eigenvector(0)[j];
      m_Tangent1[j] = eig.get_eigenvector(1)[j];
      m_Tangent2[j] = eig.get_eigenvector(2)[j];
    }

    //    this->SetFrameFromNormal(m_Normal);
    m_Eval0 = eig.get_eigenvalue(0);
    m_Eval1 = eig.get_eigenvalue(1);
    m_Eval2 = eig.get_eigenvalue(2);
    if (m_Debug)
    {
      vnl_vector<double> a = eig.get_eigenvector(0);
      std::cout << "Eig residual = " << (D * a).magnitude() << std::endl;
      a = eig.get_eigenvector(1);
      std::cout << "Eig residual = " << (D * a).magnitude() << std::endl;
      a = eig.get_eigenvector(2);
      std::cout << "Eig residual = " << (D * a).magnitude() << std::endl;
    }
  }
}

template <typename TSurface, unsigned int TDimension>
void
SurfaceCurvatureBase<TSurface, TDimension>::EstimateTangentPlane(PointType origin)
{
  // Build cov matrix D
  unsigned int npts = m_PointList.size() - 1;
  unsigned int dim = SurfaceDimension;
  MatrixType   D(npts, dim);

  for (unsigned int i = 0; i < npts; i++)
  {
    for (unsigned int j = 0; j < dim; j++)
    {
      D(i, j) = (m_PointList[i][j] - origin(j)) * (m_PointList[i][j] - origin(j));
    }
  }

  this->EigenDecomposition(D);
  m_dX = m_Tangent1.magnitude();
  m_dY = m_Tangent2.magnitude();
  m_Tangent1 /= m_dX;
  m_Tangent2 /= m_dY;
  this->ChooseReferenceTangent();
}

template <typename TSurface, unsigned int TDimension>
void
SurfaceCurvatureBase<TSurface, TDimension>::WeightedEstimateTangentPlane(PointType origin)
{
  // Build cov matrix D
  RealType     twi = 0.0;
  unsigned int npts = m_PointList.size() - 1;
  unsigned int dim = SurfaceDimension;
  MatrixType   D(npts, dim);

  for (unsigned int i = 0; i < npts; i++)
  {
    PointType t = m_PointList[i] - origin;
    RealType  wi = t.magnitude();
    if (!itk::Math::FloatAlmostEqual(wi, itk::NumericTraits<RealType>::ZeroValue()))
    {
      wi = itk::NumericTraits<RealType>::OneValue() / wi;
    }
    twi += wi;
    for (unsigned int j = 0; j < dim; j++)
    {
      D(i, j) = wi * (m_PointList[i][j] - origin(j)) * (m_PointList[i][j] - origin(j));
    }
  }

  this->EigenDecomposition(D / twi);
  this->ChooseReferenceTangent();
}

template <typename TSurface, unsigned int TDimension>
void SurfaceCurvatureBase<TSurface, TDimension>::ComputeFrame(PointType /* origin */)
{
  // Build cov matrix D
  unsigned int npts = m_PointList.size() - 1;
  unsigned int dim = SurfaceDimension;

  RealType twi = 1., wi = 0.;

  MatrixType D;
  bool       method1 = true;

  if (method1)
  {
    twi = 0.0;
    D.set_size(dim, dim);
    D.fill(0.0);
    for (unsigned int i = 0; i < npts; i++)
    {
      PointType p = m_PointList[i] - m_Origin;
      wi = 1. / p.magnitude();
      twi += wi;
      for (unsigned int j = 0; j < dim; j++)
      {
        for (unsigned int k = 0; k < dim; k++)
        {
          D(j, k) += m_TangentProjectionList[i][j] * m_TangentProjectionList[i][k]; // *wi;
        }
      }
    }
  }
  else
  {
    D.set_size(npts, dim);
    for (unsigned int i = 0; i < npts; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
      {
        D(i, j) = m_TangentProjectionList[i][j] * m_TangentProjectionList[i][j] * m_DirectionalKappaVector[i] *
                  m_WeightVector[i];
      }
    }
  }
  this->EigenDecomposition(D / twi);
  this->ChooseReferenceTangent();
}

template <typename TSurface, unsigned int TDimension>
void SurfaceCurvatureBase<TSurface, TDimension>::ComputeFrameAndKappa(PointType /* origin */)
{
  // Build cov matrix D
  unsigned int npts = m_PointList.size() - 1;
  unsigned int dim = SurfaceDimension;
  /*
    if (m_TotalDKap/(float)m_PointList.size() < 0.05)
    {
      m_Kappa1=0.0;
      m_Kappa2=0.0;
      m_A=0.;
      m_B=0.;
      m_C=0.0;
      m_GaussianKappa=0.0;
      m_TotalDKap=0.0;
      return;
    }
  */
  MatrixType D;

  //  dim=2;
  D.set_size(dim, dim);
  D.fill(0.0);
  for (unsigned int i = 0; i < npts; i++)
  {
    for (unsigned int j = 0; j < dim; j++)
    {
      for (unsigned int k = 0; k < dim; k++)
      {
        D(j, k) += static_cast<double>(m_TangentProjectionList[i][j] * m_TangentProjectionList[i][k] *
                                       m_WeightVector[i] * m_DirectionalKappaVector[i]);
      }
    }
  }

  vnl_symmetric_eigensystem<double> eig(D.transpose() * D);
  m_Eval0 = eig.get_eigenvalue(0);
  m_Eval1 = eig.get_eigenvalue(1);
  m_Eval2 = eig.get_eigenvalue(2);
}

template <typename TSurface, unsigned int TDimension>
void
SurfaceCurvatureBase<TSurface, TDimension>::ShimshoniFrame(PointType origin)
{
  this->ComputeWeightsAndDirectionalKappaAndAngles(origin);
  this->ComputeFrameAndKappa(origin);
  this->EstimateCurvature(m_A, m_B, m_B, m_C);
  //      this->EstimateCurvature();
}

template <typename TSurface, unsigned int TDimension>
void
SurfaceCurvatureBase<TSurface, TDimension>::JainMeanAndGaussianCurvature(PointType origin)
{
  // Build cov matrix D
  unsigned int npts = m_PointList.size() - 1;

  if (npts < 7)
  {
    m_MeanKappa = 0;
    m_GaussianKappa = 0;
    m_Kappa1 = 0;
    m_Kappa2 = 0;
    return;
  }

  vnl_vector<double> dists(npts);

  vnl_vector<double> wts(npts);

  MatrixType D;

  D.set_size(npts, 6); // each row contains [u^2 , uv, v^2, u, v, 1] for point p
  D.fill(0.0);

  float totwt = 0.0;
  float wt = 0.0;

  bool robust = true;
  robust = false;
  unsigned int i = 0;
  for (i = 0; i < npts; i++)
  {
    PointType Q = m_PointList[i];

    PointType Dif = Q - origin;

    float u1 = 0.0;
    float u2 = 0.0;
    wt = itk::NumericTraits<float>::OneValue() / static_cast<float>(Dif.magnitude());
    wts[i] = wt;
    //   totwt+=wt;
    for (unsigned int pp = 0; pp < SurfaceDimension; pp++)
    {
      u1 += Dif[pp] * m_Tangent1[pp];
      u2 += Dif[pp] * m_Tangent2[pp];
    }

    // co-ordinates in tangent plane
    //   RealType sign=1.0;
    //   if (u1 < 0.) sign=-1.0;
    //   u1=sqrt(fabs(u1))*sign;
    //   if (u2 < 0.) sign=-1.0; else sign=1.0;
    //   u2=sqrt(fabs(u2))*sign;

    dists[i] = 0.0;
    PointType tanproj;
    for (unsigned int jj = 0; jj < SurfaceDimension; jj++)
    {
      m_PlaneVector[jj] = u1 * m_Tangent1[jj] + u2 * m_Tangent2[jj]; // tangent projection of Q
      tanproj[jj] = origin[jj] + m_PlaneVector[jj];
    }

    PointType temp = Q - tanproj;
    dists[i] = temp.magnitude(); // sqrt(dists[i]);

    //   if (robust) dists[i]/=wt;

    D(i, 0) = u1 * u1;
    D(i, 1) = 2.0f * u1 * u2;
    D(i, 2) = u2 * u2;
    D(i, 3) = u1;
    D(i, 4) = u2;
    D(i, 5) = 1.0f;
  }

  if (robust)
  {
    totwt = 0.0;
    float sigma = m_Sigma;
    for (unsigned int tt = 0; tt < npts; tt++)
    {
      wt = exp(-1.0f * static_cast<float>(wts[tt]) / sigma);
      dists[tt] *= static_cast<double>(wt);
      totwt += wt;
      //      std::cout << " wt " << wt << std::endl;
    }
    if (totwt > 0)
    {
      dists /= totwt;
    }
  }

  vnl_svd<double>    svd(D);
  vnl_vector<double> a = svd.solve(dists);

  //  std::cout << " a vec " << a << std::endl;

  //  for (int tt=0; tt<6; tt++)
  //    if (a[tt] < 1.e-6) a[tt]=0.0;
  /*
    double fu = a(0);  //a(3);
    double fv = a(2);  //a(4);
    double fuu = a(3); //2.*a(0);
    double fuv = a(1); //a(1);
    double fvv = a(4); //2.*a(2);
  */
  double fu = a(3);
  double fv = a(4);
  double fuu = 2. * a(0);
  double fuv = a(1);
  double fvv = 2. * a(2);

  m_MeanKappa = static_cast<RealType>((1.0 + fv * fv) * fuu - 2.0 * fu * fv * fuv + (1.0 + fu * fu) * fvv);
  m_MeanKappa /= static_cast<RealType>(2.0 * std::pow(1.0 + fu * fu + fv * fv, 1.5));

  m_GaussianKappa = (fuu * fvv - fuv * fuv) / ((1.0 + fu * fu + fv * fv) * (1.0 + fu * fu + fv * fv));

  m_Kappa1 = m_MeanKappa + static_cast<float>(std::sqrt(m_MeanKappa * m_MeanKappa - m_GaussianKappa));
  m_Kappa2 = m_MeanKappa - static_cast<float>(std::sqrt(m_MeanKappa * m_MeanKappa - m_GaussianKappa));

  if (itk::Math::FloatAlmostEqual(origin[0], 20.0f) && itk::Math::FloatAlmostEqual(origin[1], 65.0f) &&
      itk::Math::FloatAlmostEqual(origin[2], 39.0f))
  {
    std::cout << " surf params " << a << std::endl;
    std::cout << " k1 " << m_Kappa1 << " k2 " << m_Kappa2 << std::endl;
    this->PrintFrame();
    for (unsigned int jj = 0; jj < npts; jj++)
    {
      std::cout << " point " << m_PointList[jj];
      std::cout << " dist " << dists[jj] << std::endl;
    }
  }

  // NOTE: THIS IS THE ERROR ESTIMATE
  // std::cout << "SVD residual = " << (D * a).magnitude() << std::endl;
  // std::cout << " a " << a << std::endl;
  //  MatrixType C(2,2);
  //  C(0,0)=a(0);   C(1,0)=a(1);   C(0,1)=a(1);   C(1,1)=a(2);
  //  vnl_symmetric_eigensystem<double> eig(C);
  //  m_Kappa1=eig.get_eigenvalue(0);
  //  m_Kappa2=eig.get_eigenvalue(1);
  //
  //  m_MeanKappa=0.5*(m_Kappa1+m_Kappa2);
  //  m_GaussianKappa = m_Kappa1*m_Kappa2;
  //
  //  std::cout << " eval 0 " << m_Kappa1 << " eval 1 " << m_Kappa2 << std::endl;
}

// end jain frame

template <typename TSurface, unsigned int TDimension>
void
SurfaceCurvatureBase<TSurface, TDimension>::ComputeJoshiFrame(PointType origin)
{
  // Build cov matrix D
  unsigned int npts = m_PointList.size() - 1;

  if (npts < 5)
  {
    m_MeanKappa = 0;
    m_GaussianKappa = 0;
    m_Kappa1 = 0;
    m_Kappa2 = 0;
    return;
  }

  vnl_vector<double> dists(npts);

  MatrixType D;

  D.set_size(npts, 3);
  D.fill(0.0);
  // float totw = 0;
  for (unsigned int i = 0; i < npts; i++)
  {
    PointType Q = m_PointList[i];

    float u1 = 0.0;
    float u2 = 0.0;
    float sqdif = 0;
    for (unsigned int ijj = 0; ijj < SurfaceDimension; ijj++)
    {
      u1 += (Q[ijj] - origin[ijj]) * m_Tangent1[ijj];
      u2 += (Q[ijj] - origin[ijj]) * m_Tangent2[ijj];
      sqdif += ((Q[ijj] - origin[ijj]) * (Q[ijj] - origin[ijj]));
    }
    // totw += 1.0f / sqdif;
    //
    // co-ordinates in tangent plane
    //   RealType sign=1.0;
    //   if (u1 < 0.) sign=-1.0;
    //   u1=sqrt(fabs(u1))*sign;
    //   if (u2 < 0.) sign=-1.0; else sign=1.0;
    //   u2=sqrt(fabs(u2))*sign;
    //
    //
    //   RealType pvMag=0;
    //   for (unsigned int i=0; i<SurfaceDimension; i++)
    //   {
    //      m_PlaneVector[i]=u1*m_Tangent1[i]+u2*m_Tangent2[i];
    //      pvMag+=m_PlaneVector[i]*m_PlaneVector[i];
    //   }
    //
    PointType tanproj;
    for (unsigned int ijj = 0; ijj < SurfaceDimension; ijj++)
    {
      m_PlaneVector[ijj] = u1 * m_Tangent1[ijj] + u2 * m_Tangent2[ijj]; // tangent projection of Q
      tanproj[ijj] = origin[ijj] + m_PlaneVector[ijj];
    }

    PointType temp = Q - tanproj;
    dists[i] = temp.magnitude(); // /sqdif;//*exp(-1.*sqdif/2);//sqrt(dists[i]);

    D(i, 0) = u1 * u1;
    D(i, 1) = 2.0f * u1 * u2;
    D(i, 2) = u2 * u2;
  }

  vnl_svd<double>    svd(D);
  vnl_vector<double> c = svd.solve(dists);

  // NOTE: THIS IS THE ERROR ESTIMATE
  // std::cout << "SVD residual = " << (D * c).magnitude() << std::endl;
  // std::cout << " C " << c << std::endl;

  MatrixType C(2, 2);
  C(0, 0) = c(0);
  C(1, 0) = c(1);
  C(0, 1) = c(1);
  C(1, 1) = c(2);

  vnl_symmetric_eigensystem<double> eig(C);

  m_Kappa1 = eig.get_eigenvalue(0);
  m_Kappa2 = eig.get_eigenvalue(1);

  m_MeanKappa = 0.5f * (m_Kappa1 + m_Kappa2);
  m_GaussianKappa = m_Kappa1 * m_Kappa2;
  //  std::cout << " eval 0 " << m_Kappa1 << " eval 1 " << m_Kappa2 << std::endl;
}

template <typename TSurface, unsigned int TDimension>
typename SurfaceCurvatureBase<TSurface, TDimension>::RealType
SurfaceCurvatureBase<TSurface, TDimension>::ComputeLocalArea(FixedVectorType spacing)
{
  PointType t1;
  PointType t2;

  for (unsigned int i = 0; i < SurfaceDimension; i++)
  {
    t1[i] = m_Tangent1[i] * spacing[i];
    t2[i] = m_Tangent2[i] * spacing[i];
  }

  PointType temp = vnl_cross_3d(t1, t2);

  return temp.magnitude();
}

template <typename TSurface, unsigned int TDimension>
unsigned int
SurfaceCurvatureBase<TSurface, TDimension>::CharacterizeSurface()
{
  float th = 1.e-6;

  if (static_cast<float>(fabs(m_GaussianKappa)) < th * th)
  {
    m_GaussianKappa = 0.0;
  }
  if (static_cast<float>(fabs(m_MeanKappa)) < th)
  {
    m_MeanKappa = 0.0;
  }

  if (m_MeanKappa > 0 && m_GaussianKappa > 0)
  {
    return 1; // bowl
  }
  else if (m_MeanKappa < 0 && m_GaussianKappa > 0)
  {
    return 2; // inv bowl
  }
  else if (m_MeanKappa > 0 && m_GaussianKappa < 0)
  {
    return 3; // pos saddle
  }
  else if (m_MeanKappa < 0 && m_GaussianKappa < 0)
  {
    return 4; // neg saddle
  }
  else if (m_MeanKappa > 0 && itk::Math::FloatAlmostEqual(m_GaussianKappa, 0.0f))
  {
    return 5; // pos U
  }
  else if (m_MeanKappa < 0 && itk::Math::FloatAlmostEqual(m_GaussianKappa, 0.0f))
  {
    return 6; // neg U
  }
  else if (itk::Math::FloatAlmostEqual(m_MeanKappa, 0.0f) && itk::Math::FloatAlmostEqual(m_GaussianKappa, 0.0f))
  {
    return 7; // flat
  }
  else if (itk::Math::FloatAlmostEqual(m_MeanKappa, 0.0f) && m_GaussianKappa < 0)
  {
    return 8; // perfect saddle
  }
  else
  {
    return 0;
  }

  // old below
  if (m_MeanKappa > 0 && m_GaussianKappa > 0)
  {
    return 8;
  }
  else if (m_MeanKappa > 0 && itk::Math::FloatAlmostEqual(m_GaussianKappa, 0.0f))
  {
    return 7;
  }
  else if (m_MeanKappa > 0 && m_GaussianKappa < 0)
  {
    return 6;
  }
  else if (itk::Math::FloatAlmostEqual(m_MeanKappa, 0.0f) && itk::Math::FloatAlmostEqual(m_GaussianKappa, 0.0f))
  {
    return 5;
  }
  else if (itk::Math::FloatAlmostEqual(m_MeanKappa, 0.0f) && m_GaussianKappa < 0)
  {
    return 4;
  }
  else if (m_MeanKappa < 0 && m_GaussianKappa > 0)
  {
    return 3;
  }
  else if (m_MeanKappa < 0 && itk::Math::FloatAlmostEqual(m_GaussianKappa, 0.0f))
  {
    return 2;
  }
  else if (m_MeanKappa < 0 && m_GaussianKappa < 0)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

template <typename TSurface, unsigned int TDimension>
typename SurfaceCurvatureBase<TSurface, TDimension>::RealType
SurfaceCurvatureBase<TSurface, TDimension>::ComputeMeanEuclideanDistance()
{
  double       totdist = 0.0;
  unsigned int ct = 0;
  unsigned int pls = m_PointList.size();

  for (unsigned int i = 0; i < pls; i++)
  {
    for (unsigned int j = 0; j < pls; j++)
    {
      PointType p = m_PointList[j] - m_PointList[i];
      totdist += static_cast<double>(p.magnitude());
      ct++;
    }
  }

  RealType meandist = static_cast<RealType>(totdist) / static_cast<RealType>(ct);
  // 4 pi r^2 = a,  r^2 = a/(4pi)
  RealType spheremeandist = 1.0; // sqrt((float)ct/(4.0*3.1418));
  RealType compactness = meandist / spheremeandist;

  return 10.0f - compactness;
}

template <typename TSurface, unsigned int TDimension>
void
SurfaceCurvatureBase<TSurface, TDimension>::TestEstimateTangentPlane(PointType origin)
{
  // Read points from stdin
  MatrixType pts;

  std::cout << " input points " << std::endl;

  std::cin >> pts;

  // Build cov matrix D
  int        npts = pts.rows();
  int        dim = pts.columns();
  MatrixType D(npts, dim);

  // get average of points
  RealType  count = 0.0;
  PointType avgpt;
  avgpt.fill(0);
  for (int i = 0; i < npts; i++)
  {
    for (int j = 0; j < dim; j++)
    {
      avgpt(j) += pts(i, j);
      count += static_cast<RealType>(1.0);
    }
  }

  avgpt /= count;
  if (m_Debug)
  {
    std::cout << " avg " << avgpt << std::endl;
  }

  origin = avgpt;
  for (int i = 0; i < npts; i++)
  {
    for (int j = 0; j < dim; j++)
    {
      D(i, j) = (pts(i, j) - origin(j)) * (pts(i, j) - origin(j));
    }
    // D(i,dim) = 1;
  }

  // 1. Compute using SVD
  {
    vnl_svd<double>    svd(D);
    vnl_vector<double> a = svd.nullvector();
    std::cout << "SVD residual = " << (D * a).magnitude() << std::endl;
    std::cout << "SVD normal " << a << std::endl;
  }

  // 2. Compute using eigensystem of D'*D
  {
    vnl_symmetric_eigensystem<double> eig(D.transpose() * D);
    vnl_vector<double>                a = eig.get_eigenvector(0);
    std::cout << "Eig residual = " << (D * a).magnitude() << std::endl;
    std::cout << " normal  " << eig.get_eigenvector(0) << std::endl;
    std::cout << "Eigvec 1  " << eig.get_eigenvector(1) << std::endl;
    std::cout << "Eigvec 2  " << eig.get_eigenvector(2) << std::endl;
    std::cout << "Eigval normal  " << eig.get_eigenvalue(0) << std::endl;
    std::cout << "Eigval 1  " << eig.get_eigenvalue(1) << std::endl;
    std::cout << "Eigval 2  " << eig.get_eigenvalue(2) << std::endl;
  }
}

template <typename TSurface, unsigned int TDimension>
void
SurfaceCurvatureBase<TSurface, TDimension>::FindNeighborhood(unsigned int /* temp */)
{
  // Read points from stdin
  MatrixType pts;

  std::cout << " input points " << std::endl;

  std::cin >> pts;

  // Build cov matrix D
  unsigned int npts = pts.rows();
  unsigned int dim = pts.columns();

  m_AveragePoint.fill(0.0);
  for (unsigned int i = 0; i < npts; i++)
  {
    PointType pt;
    for (unsigned int j = 0; j < dim; j++)
    {
      pt(j) = pts(i, j);
      m_AveragePoint[j] += static_cast<RealType>(pts(i, j));
    }
    if (i == 0)
    {
      m_Origin = pt;
    }
    m_PointList.insert(m_PointList.begin(), pt);
  }

  m_AveragePoint /= (RealType)npts;
  if (m_Debug)
  {
    std::cout << " point list size " << m_PointList.size() << std::endl;
    for (unsigned int i = 0; i < m_PointList.size(); i++)
    {
      std::cout << " point  " << m_PointList[i];
    }
    std::cout << std::endl;
  }
}

/** This function sets the reference tangent arbitrarily.
 *   It can be overridden in case there is a better practical approach.
 */
template <typename TSurface, unsigned int TDimension>
void
SurfaceCurvatureBase<TSurface, TDimension>::ChooseReferenceTangent()
{
  float w = 1.;
  float w2 = (1.0f - w);

  // m_ArbitraryTangent[0]=w;  m_ArbitraryTangent[1]=1.-w;  m_ArbitraryTangent[2]=0.;
  m_ArbitraryTangent = w * m_Tangent2 + (w2 * m_Tangent1);
  m_ArbitraryTangent /= m_ArbitraryTangent.magnitude();
  if (m_Debug)
  {
    std::cout << " arb tan " << m_ArbitraryTangent << std::endl;
  }
}

template <typename TSurface, unsigned int TDimension>
void
SurfaceCurvatureBase<TSurface, TDimension>::ComputeWeightsAndDirectionalKappaAndAngles(PointType origin)
{
  unsigned int sz = m_PointList.size();

  if (sz > 1)
  {
    sz = sz - 1;
  }
  m_TangentProjectionList.clear();
  m_DirectionalKappaVector.set_size(sz);
  m_DirectionalKappaVector.fill(0.0);
  m_WeightVector.set_size(sz);
  m_WeightVector.fill(0.0);
  m_ThetaVector.set_size(sz);
  m_ThetaVector.fill(0.0);

  RealType w1 = 0, w2 = 0, totalWeight = 0.0;
  RealType difMag = 0, sqdifMag = 0;

  m_TotalDKap = 0.0;

  m_A = 0.0;
  m_B = 0.0;
  m_C = 0.0;
  for (unsigned int ii = 0; ii < sz; ii++)
  {
    PointType Q = m_PointList[ii];

    PointType Dif = Q - origin;

    w1 = 0.0;
    w2 = 0.0;
    sqdifMag = 0.0;
    m_DirectionalKappa = 0.0;
    for (unsigned int i = 0; i < SurfaceDimension; i++)
    {
      m_DirectionalKappa += (m_Normal[i] * Dif[i]);
      w1 += Dif[i] * m_Tangent1[i];
      w2 += Dif[i] * m_Tangent2[i];
      sqdifMag += Dif[i] * Dif[i];
    }

    if (!itk::Math::FloatAlmostEqual(sqdifMag, 0.0f) && Q != origin)
    {
      //   m_PlaneVector[0]=w1;  m_PlaneVector[1]=w2;

      difMag = sqrt(sqdifMag);
      m_ThetaVector[ii] = GetTheta(Q, origin);
      m_TangentProjectionList.insert(m_TangentProjectionList.begin(), m_PlaneVector);

      totalWeight += 1.0f / difMag;
      m_WeightVector[ii] = 1.0f / difMag;
      m_TotalDKap += m_DirectionalKappa;
      m_DirectionalKappaVector[ii] = m_DirectionalKappa;
    }
    else
    {
      m_WeightVector[ii] = 0.0;
      m_DirectionalKappaVector[ii] = 0.0;
      PointType pp;
      pp.fill(0.0);
      m_TangentProjectionList.insert(m_TangentProjectionList.begin(), pp);
    }

    // if( m_Origin[0] == 95 && m_Origin[1] == 94 && m_Origin[2] == 63 )
    //   {
    //   std::cout << " tdk " << m_TotalDKap << " nor " << m_Normal << " dk "
    //                    << m_DirectionalKappa << " dif " << Dif << " mpv " << m_PlaneVector << std::endl;
    //   // m_Debug=true;
    //   }
    double costheta = cos(m_ThetaVector[ii]);
    double sintheta = sin(m_ThetaVector[ii]);
    double weight = m_WeightVector[ii];
    m_A += static_cast<float>(weight * costheta * costheta * costheta * costheta);
    m_B += static_cast<float>(weight * sintheta * sintheta * costheta * costheta);
    m_C += static_cast<float>(weight * sintheta * sintheta * sintheta * sintheta);
  }

  // if( m_Origin[0] == 54 && m_Origin[1] == 54 && m_Origin[2] == 63 )
  //   {
  //   std::cout << m_Origin << " tdk " << m_TotalDKap << " nor " << m_Normal << std::endl;
  //   // m_Debug=true;
  //   }
  // Now compute A, B, C
  m_A /= totalWeight;
  m_B /= totalWeight;
  m_C /= totalWeight;
  m_WeightVector /= totalWeight;

  if (m_Debug)
  {
    std::cout << " weight " << m_WeightVector << std::endl;
    std::cout << " theta " << m_ThetaVector << std::endl;
    std::cout << " dkap  " << m_DirectionalKappaVector << std::endl;
    std::cout << " A " << m_A << " B " << m_B << " C " << m_C << std::endl;
  }
}

/** This function returns a weight given a distance
 *  It may be the identity function, a normalization
 *  or a gaussianization of the input distance.
 */
template <typename TSurface, unsigned int TDimension>
typename SurfaceCurvatureBase<TSurface, TDimension>::RealType
SurfaceCurvatureBase<TSurface, TDimension>::GetWeight(PointType p1, PointType p2)
{
  PointType d = p1 - p2;

  return d.magnitude();
}

/** This function returns the angle between the reference tangent
    and the projection onto the tangent plane of the vector between
    the neighborhood focus and its neighbor.
*/
template <typename TSurface, unsigned int TDimension>
typename SurfaceCurvatureBase<TSurface, TDimension>::RealType
SurfaceCurvatureBase<TSurface, TDimension>::GetTheta(PointType Q, PointType origin)
{
  RealType  theta = 0.0;
  PointType Dif = Q - origin;

  //   float dm=Dif.magnitude();
  //   if (dm==0) dm=1.0;
  //   Dif/=dm;

  float u1 = 0.0;
  float u2 = 0.0;

  for (unsigned int i = 0; i < SurfaceDimension; i++)
  {
    u1 += Dif[i] * m_Tangent1[i];
    u2 += Dif[i] * m_Tangent2[i];
  }

  // co-ordinates in tangent plane
  //   RealType sign=1.0;
  //   if (u1 < 0.) sign=-1.0;
  //   u1=sqrt(fabs(u1))*sign;
  //   if (u2 < 0.) sign=-1.0; else sign=1.0;
  //   u2=sqrt(fabs(u2))*sign;

  RealType tot = 1; // fabs(u1)+fabs(u2);
                    //   if (tot ==0)
  {
    //    std::cout << " tan1 " << m_Tangent1 << " tan2 " << m_Tangent2 << std::endl;
    //    std::cout << " norm " << m_Normal << " u1 " << u1 << " u2 " << u2 << std::endl;
    //    std::cout << " Origin " << m_Origin << " Q " << Q << std::endl;
    //    tot=1.;
  }
  RealType pvMag = 0;
  RealType ip = 0.0;
  RealType arbMag = 0.0;
  for (unsigned int i = 0; i < SurfaceDimension; i++)
  {
    m_PlaneVector[i] = u1 / (tot)*m_Tangent1[i] + u2 / (tot)*m_Tangent2[i];
    ip += m_PlaneVector[i] * m_ArbitraryTangent[i];
    pvMag += m_PlaneVector[i] * m_PlaneVector[i];
    arbMag += m_ArbitraryTangent[i] * m_ArbitraryTangent[i];
  }

  //   ProjectToTangentPlane(Dif);
  //   m_PlaneVector/=m_PlaneVector.magnitude();
  //   if (pvMag!=0)m_PlaneVector/=sqrt(pvMag);

  if (m_Debug)
  {
    std::cout << " m_PlaneVector " << m_PlaneVector << " dif " << Dif << std::endl;
  }

  float temp = (sqrt(pvMag) * sqrt(arbMag));
  if (!itk::Math::FloatAlmostEqual(temp, 0.0f))
  {
    theta = acos(ip / temp);
  }
  else
  {
    theta = 0.0; // FIXME?
  }
  return theta;
}

/** */
template <typename TSurface, unsigned int TDimension>
void
SurfaceCurvatureBase<TSurface, TDimension>::EstimateDirectionalCurvature(PointType Q, PointType P)
{
  PointType Dif = Q - P;

  m_DirectionalKappa = 0.0;
  for (unsigned int i = 0; i < SurfaceDimension; i++)
  {
    m_DirectionalKappa += (m_Normal[i] * Dif[i]);
  }

  m_DirectionalKappa = 2.0f * m_DirectionalKappa / Dif.magnitude();
}

/** */
template <typename TSurface, unsigned int TDimension>
void
SurfaceCurvatureBase<TSurface, TDimension>::EstimateCurvature(RealType w1, RealType w2, RealType w3, RealType w4)
{
  /*
    if (m_TotalDKap/(float)m_PointList.size() < 0.05)
    {
      m_Kappa1=0.0;
      m_Kappa2=0.0;
      m_A=0.;
      m_B=0.;
      m_C=0.0;
      m_GaussianKappa=0.0;
      m_TotalDKap=0.0;
      return;
    }*/

  RealType denom = (w4 - w2 * w3 / w1);

  if (std::fabs(denom) > static_cast<RealType>(1.e-12))
  {
    m_Kappa2 = (m_Eval2 - m_Eval1 * w3 / w1) / denom;
  }
  else
  {
    w1 = m_W1;
    w2 = m_W2;
    w3 = m_W2;
    w4 = m_W1;
    denom = (w4 - w2 * w3 / w1);
    m_Kappa2 = (m_Eval2 - m_Eval1 * w3 / w1) / denom;
  }
  m_Kappa1 = (m_Eval1 - w2 * m_Kappa1) / w1;

  m_MeanKappa = 0.5f * (m_Kappa1 + m_Kappa2);
  m_GaussianKappa = m_Kappa1 * m_Kappa2;

  //  std::cout << " eval test " << w1*m_Kappa1 + w2*m_Kappa2 << " e1 " << m_Eval1 << std::endl;
  //  std::cout << " eval test " << w3*m_Kappa1 + w4*m_Kappa2 << " e2 " << m_Eval2 << std::endl;
}

/** */
template <typename TSurface, unsigned int TDimension>
void
SurfaceCurvatureBase<TSurface, TDimension>::PrintFrame()
{
  std::cout << " normal " << m_Normal << std::endl;
  std::cout << " t1 " << m_Tangent1 << std::endl;
  std::cout << " t2 " << m_Tangent2 << std::endl;
  std::cout << " k1 " << m_Kappa1 << " k2 " << m_Kappa2 << std::endl;
  std::cout << " e0 " << m_Eval0 << " e1 " << m_Eval1 << " e2 " << m_Eval2 << std::endl;
  std::cout << " A " << m_A << " B " << m_B << " C " << m_C << std::endl;
  std::cout << std::endl;
}

/** */
template <typename TSurface, unsigned int TDimension>
typename SurfaceCurvatureBase<TSurface, TDimension>::RealType
SurfaceCurvatureBase<TSurface, TDimension>::ErrorEstimate(PointType origin, RealType sign)
{
  unsigned int npts = m_PointList.size() - 1;
  unsigned int dim = SurfaceDimension;

  //  float error=0.0;
  float toterror = 0.0;

  double    qpt = 0.0;
  double    qpnkpt = 0.0;
  double    costheta;
  double    sintheta;
  double    kp;
  RealType  theta;
  PointType Q, est;

  for (unsigned int pp = 0; pp < npts; pp++)
  {
    Q = m_PointList[pp];

    if (Q != origin)
    {
      theta = this->GetTheta(Q, origin);
      costheta = cos(theta);
      sintheta = sin(theta);
      kp = static_cast<double>(m_Kappa1) * costheta * costheta + static_cast<double>(m_Kappa2) * sintheta * sintheta;

      if (!itk::Math::FloatAlmostEqual(kp, 0.0))
      {
        PointType qp = Q - origin;

        FixedVectorType normal = m_Normal * sign;

        float w1 = 0, w2 = 0;
        for (unsigned int j = 0; j < dim; j++)
        {
          qpt += qp[j] * normal[j] * kp;
          w1 += qp[j] * m_Tangent1[j];
          w2 += qp[j] * m_Tangent2[j];
          qpnkpt += qp[j] * m_PlaneVector[j];
        }

        vnl_vector<double> f2(4);
        f2(0) = 0.5 * kp * kp;
        f2(1) = 0.0;
        f2(2) = 1. - qpnkpt;
        f2(3) = -1. * qpt;

        /** commenting out until we can get vnl_rpoly_roots working
            vnl_rpoly_roots roots(f2);
            // Evaluate results
            //vnl_real_polynomial p(f2);
            //for(unsigned int i = 0; i < p.degree(); i++)
            //  vnl_test_assert("Root residual", std::abs(p.evaluate(roots[i])) < 1e-12);

            float minrel=9.e9;
            float mins=0.0;
            float s=0.0;
            vnl_vector<double> so=roots.real();
        //    for (s=0.1; s<=2.0; s=s+.01)
            for (unsigned int ind=0; ind<so.size(); ind++)
            {
              s=fabs(so[ind]);
              est= s*m_PlaneVector+ (float)(s*s*0.5*kp*kp)*normal;
              PointType dif = est - qp;
              error=dif.magnitude();
              if (error < minrel)
              {
                minrel=error;
                mins=s;
              }
            }
          if (m_Debug){
              std::cout << " single point error " << minrel << " mins " << s<< std::endl;
              std::cout << " est pt " << est << " pt " << qp << std::endl << std::endl;
            }
            toterror+=minrel;
        */
      }
    }
  }

  if (m_Debug)
  {
    std::cout << " tot error " << 1.0f / (static_cast<float>(npts + 1)) * toterror << std::endl;
  }

  return 1.0f / (static_cast<float>(npts + 1)) * toterror;
}

template <typename TSurface, unsigned int TDimension>
void
SurfaceCurvatureBase<TSurface, TDimension>::EstimateMetricTensor()
{
  float a = 0;

  for (unsigned int i = 0; i < 3; i++)
  {
    a += m_Tangent1[i] * m_Tangent1[i];
  }
  m_MetricTensor[0] = a;

  float b = 0;
  for (unsigned int i = 0; i < 3; i++)
  {
    b += m_Tangent2[i] * m_Tangent1[i];
  }
  m_MetricTensor[1] = b;

  float c = 0;
  for (unsigned int i = 0; i < 3; i++)
  {
    c += m_Tangent2[i] * m_Tangent2[i];
  }
  m_MetricTensor[2] = c;

  m_MetricTensorDeterminant = sqrt(a * c - b * b);

  if (m_MetricTensorDeterminant < 0.0f)
  {
    std::cout << "bad metric tensor ";
  }
}

template <typename TSurface, unsigned int TDimension>
float
SurfaceCurvatureBase<TSurface, TDimension>::dstarUestimate()
{
  // this->FindNeighborhood();  make sure this also fills the function values functionvaluelist
  this->EstimateTangentPlane(m_Origin);

  // below, basically replace the heights with the function values

  // Build cov matrix D
  int npts = m_PointList.size();

  if (npts < 4)
  {
    return 0.0;
  }

  vnl_vector<double> dists(npts);
  vnl_vector<double> surfx(npts);
  vnl_vector<double> surfy(npts);
  vnl_vector<double> surfz(npts);
  vnl_vector<double> funcvals(npts);
  vnl_vector<double> wts(npts);
  MatrixType         D;

  D.set_size(npts, 3); // each row contains [ u, v, 1] for point p
  D.fill(0.0);

  // float totwt = 0.0;
  float totfunc = 0.0;
  float wt = 0.0;
  float meanu1 = 0;
  float meanu2 = 0;
  for (int i = 0; i < npts; i++)
  {
    PointType Q = m_PointList[i];
    PointType Dif = Q - m_Origin;

    surfx[i] = Dif[0];
    surfy[i] = Dif[1];
    surfz[i] = Dif[2];

    float u1 = 0.0;
    float u2 = 0.0;
    wt = Dif.magnitude();
    wts[i] = wt;
    // totwt += wt;
    for (unsigned int pp = 0; pp < SurfaceDimension; pp++)
    {
      u1 += Dif[pp] * m_Tangent1[pp];
      u2 += Dif[pp] * m_Tangent2[pp];
    }

    dists[i] = 0.0;
    PointType tanproj;
    for (unsigned int jj = 0; jj < SurfaceDimension; jj++)
    {
      m_PlaneVector[jj] = u1 * m_Tangent1[jj] + u2 * m_Tangent2[jj]; // tangent projection of Q
      tanproj[jj] = m_Origin[jj] + m_PlaneVector[jj];
    }

    PointType temp = Q - tanproj;
    dists[i] = temp.magnitude();
    funcvals[i] = m_FunctionValueList[i];
    totfunc += m_FunctionValueList[i];

    D(i, 0) = u1;
    D(i, 1) = u2;
    D(i, 2) = 1.0;

    meanu1 += static_cast<float>(fabs(u1));
    meanu2 += static_cast<float>(fabs(u2));
  }

  meanu1 /= (float)npts;
  meanu2 /= (float)npts;

  float mx;
  if (meanu1 > meanu2)
  {
    mx = meanu1;
  }
  else
  {
    mx = meanu2;
  }
  meanu1 /= mx;
  meanu2 /= mx;

  vnl_svd<double> svd(D);

  float              aa, bb, cc;
  vnl_vector<double> a = svd.solve(surfx);
  vnl_vector<double> b = svd.solve(surfy);
  vnl_vector<double> c = svd.solve(surfz);

  aa = a(0) * a(0) + b(0) * b(0) + c(0) * c(0);
  bb = a(1) * a(0) + b(1) * b(0) + c(1) * c(0);
  cc = a(1) * a(1) + b(1) * b(1) + c(1) * c(1);

  vnl_vector<double> func = svd.solve(funcvals);
  auto               Ux = static_cast<float>(func(0));
  auto               Uy = static_cast<float>(func(1));
  float              dx = meanu1; // m_dX;//m_Eval1;//m_dX;
  float              dy = meanu2; // m_dY;//m_Eval2;//m_dY;
  m_MetricTensor[0] = aa;
  m_MetricTensor[1] = bb;
  m_MetricTensor[2] = cc;
  float denom = sqrt(aa * cc - bb * bb);
  //  denom=1.0;
  //  std::cout << " denom " << denom << " dx " << dx << " dy " << dy << std::endl;
  // std::cout << " a " << aa << " b " << bb << " c " << cc << std::endl;

  float dstarU = 1.0f / denom * ((bb * Ux - aa * Uy) * dx + (cc * Ux - bb * Uy) * dy);

  // std::cout << " dstarU " << dstarU << std::endl;

  return dstarU;
}
} // namespace itk

#endif
