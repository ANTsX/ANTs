#ifndef TensorFunctions_cxx
#define TensorFunctions_cxx

// #include "itkSmallStrainDiffusionTensorReorientationImageFilter.h"
// #include "itkVectorIndexSelectionCastImageFilter.h"

#include "antsUtilities.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkVector.h"
#include "itkVersor.h"
#include "itkVariableSizeMatrix.h"
#include "itkDecomposeTensorFunction2.h"
#include "itkRotationMatrixFromVectors.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "itkMatrix.h"
#include "itkVariableSizeMatrix.h"
namespace matHelper
{
template <typename TFloat, unsigned int dim>
unsigned int
Rows(const vnl_matrix_fixed<TFloat, dim, dim> &)
{
  return dim;
}

template <typename TFloat, unsigned int dim>
unsigned int
Columns(const vnl_matrix_fixed<TFloat, dim, dim> &)
{
  return dim;
}

template <typename TFloat, unsigned int rows, unsigned int columns>
unsigned int
Rows(const itk::Matrix<TFloat, rows, columns> &)
{
  return rows;
}

template <typename TFloat, unsigned int rows, unsigned int columns>
unsigned int
Columns(const itk::Matrix<TFloat, rows, columns> &)
{
  return columns;
}

template <typename TFloat>
unsigned int
Rows(const itk::VariableSizeMatrix<TFloat> & mat)
{
  return mat.Rows();
}

template <typename TFloat>
unsigned int
Columns(const itk::VariableSizeMatrix<TFloat> & mat)
{
  return mat.Cols();
}

} // namespace matHelper

template <typename TensorType, typename MatrixType>
void
Vector2Matrix(TensorType & dtv, MatrixType & dtm)
{
  // dtm(0, 0) = dtv[0];
  // dtm(0, 1) = dtm(1, 0) = dtv[1];
  // dtm(0, 2) = dtm(2, 0) = dtv[2];
  // dtm(1, 1) = dtv[3];
  // dtm(1, 2) = dtm(2, 1) = dtv[4];
  // dtm(2, 2) = dtv[5];
  unsigned int tensorIndex = 0;

  for (unsigned i = 0; i < matHelper::Rows(dtm); ++i)
  {
    for (unsigned j = i; j < matHelper::Columns(dtm); ++j, ++tensorIndex)
    {
      dtm(i, j) = dtm(j, i) = dtv[tensorIndex];
    }
  }
}

template <typename TensorType, typename MatrixType>
MatrixType
Vector2Matrix(TensorType dtv)
{
  MatrixType dtm(3, 3);

  Vector2Matrix<TensorType, MatrixType>(dtv, dtm);

  return dtm;
}

template <typename TensorType, typename MatrixType>
TensorType
Matrix2Vector(MatrixType dtm)
{
  TensorType dtv;

  // dtv[0] = dtm(0, 0);
  // dtv[1] = dtm(0, 1);
  // dtv[2] = dtm(0, 2);
  // dtv[3] = dtm(1, 1);
  // dtv[4] = dtm(1, 2);
  // dtv[5] = dtm(2, 2);
  unsigned int tensorIndex = 0;

  for (unsigned i = 0; i < dtm.rows(); ++i)
  {
    for (unsigned j = i; j < dtm.cols(); ++j, ++tensorIndex)
    {
      dtv[tensorIndex] = dtm(i, j);
    }
  }
  return dtv;
}

template <typename TensorType, typename MatrixType>
void
EigenAnalysis(TensorType dtv, MatrixType & evals, MatrixType & evecs)
{
  MatrixType dtm = Vector2Matrix<TensorType, MatrixType>(dtv);

  itk::DecomposeTensorFunction2<MatrixType, typename MatrixType::ValueType, MatrixType> decomposer;

  decomposer.EvaluateSymmetricEigenDecomposition(dtm, evals, evecs);
}

template <typename TensorType, typename VectorType>
float
DiffusionCoefficient(TensorType dtv, VectorType direction, bool normalized = false)
{
  vnl_matrix<double> tensor(3, 3);
  tensor[0][0] = dtv[0];
  tensor[1][0] = tensor[0][1] = dtv[1];
  tensor[2][0] = tensor[0][2] = dtv[2];
  tensor[1][1] = dtv[3];
  tensor[1][2] = tensor[2][1] = dtv[4];
  tensor[2][2] = dtv[5];

  vnl_matrix<double> tangentmat(3, 1);
  tangentmat(0, 0) = direction[0];
  tangentmat(1, 0) = direction[1];
  tangentmat(2, 0) = direction[2];

  vnl_matrix<double> fddmat = tangentmat.transpose() * tensor * tangentmat;
  double             fdd = fddmat(0, 0);

  if (normalized > 0)
  {
    fdd = fdd / (tensor[0][0] + tensor[1][1] + tensor[2][2]);
  }

  return static_cast<float>(fdd);
}

namespace tensorHelper
{

template <typename TFloat, unsigned int NDim>
unsigned int
size(const itk::SymmetricSecondRankTensor<TFloat, NDim> &)
{
  return NDim;
}

template <typename TFloat, unsigned int NDim>
unsigned int
size(const itk::Vector<TFloat, NDim> &)
{
  return NDim;
}

} // namespace tensorHelper
template <typename TensorType>
TensorType
TensorLogAndExp(TensorType dtv, bool takelog, bool & success)
{
  float eps = 1.e-12, mag = 0;

  for (unsigned int jj = 0; jj < tensorHelper::size(dtv); jj++)
  {
    float ff = dtv[jj];
    mag += ff * ff;
    if (std::isnan(ff) || std::isinf(ff))
    {
      dtv.Fill(0); // dtv[0]=eps;   dtv[3]=eps;  dtv[5]=eps;
      success = false;
      return dtv;
    }
  }
  mag = sqrt(mag);

  // if(  dtv[1] == 0 && dtv[2] == 0 && dtv[4] == 0 )
  //   {
  //   success = false;
  //   return dtv;
  //   }
  // rewrite above to avoid going beyond tensor array bounds.
  {
    unsigned int indices[3] = { 1, 2, 4 };
    unsigned int jj;
    for (jj = 0; indices[jj] < tensorHelper::size(dtv); ++jj)
    {
      if (!itk::Math::FloatAlmostEqual(static_cast<double>(dtv[indices[jj]]), 0.0))
      {
        break;
      }
    }
    if (jj == 3)
    {
      success = false;
      return dtv;
    }
  }

  if (mag < eps)
  {
    success = false;
    return dtv;
  }
  // typedef vnl_matrix<double>        MatrixType;
  typedef itk::VariableSizeMatrix<typename TensorType::ValueType> MatrixType;

  // MatrixType DT = Vector2Matrix<TensorType,MatrixType>(dtv);
  MatrixType D;
  MatrixType V;
  EigenAnalysis<TensorType, MatrixType>(dtv, D, V);
  double e1 = D(0, 0);
  double e2 = D(1, 1);
  double e3 = D(2, 2);
  // float peigeps=1.e-12;

  // if( fabs(e3) < eps )
  //   {
  //   // success = false;
  //   // std::cout << "-4" << std::flush;
  //   // return dtv;
  //   }

  MatrixType eigmat(3, 3);
  eigmat.Fill(0);
  if (takelog)
  {
    if (e1 < 0)
    {
      e1 = e2;
    }
    if (e3 < 0)
    {
      e3 = e2;
    }
    eigmat(0, 0) = log(fabs(e1));
    eigmat(1, 1) = log(fabs(e2));
    eigmat(2, 2) = log(fabs(e3));
  }
  else // take exp
  {
    eigmat(0, 0) = exp(e1);
    eigmat(1, 1) = exp(e2);
    eigmat(2, 2) = exp(e3);
  }

  if (std::isnan(eigmat(0, 0)) || std::isnan(eigmat(1, 1)) || std::isnan(eigmat(2, 2)))
  {
    dtv.Fill(0);
    success = false;
    return dtv;
  }

  typedef typename MatrixType::InternalMatrixType VnlMatrixType;
  VnlMatrixType                                   DTrec = V.GetVnlMatrix() * eigmat.GetVnlMatrix() * V.GetTranspose();
  TensorType                                      dtv2 = Matrix2Vector<TensorType, VnlMatrixType>(DTrec);

  return dtv2;
}

template <typename TensorType>
TensorType
TensorLog(TensorType dtv, bool success = true)
{
  return TensorLogAndExp<TensorType>(dtv, true, success);
}

template <typename TensorType>
TensorType
TensorExp(TensorType dtv, bool /* takelog */, bool success = true)
{
  return TensorLogAndExp<TensorType>(dtv, false, success);
}

template <typename TensorType>
bool
IsRealTensor(TensorType dtv)
{
  bool isreal = true;

  for (unsigned int i = 0; i < 6; i++)
  {
    if (!std::isfinite(dtv[i]))
    {
      isreal = false;
    }
  }
  return isreal;
}

template <typename TensorType>
float
GetTensorFA(TensorType dtv)
{

  // Check for zero diffusion (probably background) and return zero FA
  // if that's the case
  typename TensorType::ValueType sum = dtv[0] + dtv[3] + dtv[5];
  if (itk::Math::FloatAlmostEqual(sum, itk::NumericTraits<typename TensorType::ValueType>::ZeroValue()))
  {
    return itk::NumericTraits<float>::ZeroValue();
  }

  typedef vnl_matrix<double> MatrixType;
  MatrixType                 DT(3, 3);
  DT.fill(0);
  DT(0, 0) = dtv[0];
  DT(1, 1) = dtv[3];
  DT(2, 2) = dtv[5];
  DT(1, 0) = DT(0, 1) = dtv[1];
  DT(2, 0) = DT(0, 2) = dtv[2];
  DT(2, 1) = DT(1, 2) = dtv[4];
  vnl_symmetric_eigensystem<double> eig(DT);
  double                            e1 = (eig.D(0, 0));
  double                            e2 = (eig.D(1, 1));
  double                            e3 = (eig.D(2, 2));
  if (e1 < 0)
  {
    e1 = e2;
  }
  if (e3 < 0)
  {
    e3 = e2;
  }
  // compute variance of e's
  double emean = (e1 + e2 + e3) / 3.0;
  double numer = sqrt((e1 - emean) * (e1 - emean) + (e2 - emean) * (e2 - emean) + (e3 - emean) * (e3 - emean));
  double denom = sqrt(e1 * e1 + e2 * e2 + e3 * e3);
  double fa = sqrt(3.0 / 2.0) * numer / denom;
  return fa;
}

template <typename TensorType>
float
GetTensorFANumerator(TensorType dtv)
{
  typedef vnl_matrix<double> MatrixType;
  MatrixType                 DT(3, 3);
  DT.fill(0);
  DT(0, 0) = dtv[0];
  DT(1, 1) = dtv[3];
  DT(2, 2) = dtv[5];
  DT(1, 0) = DT(0, 1) = dtv[1];
  DT(2, 0) = DT(0, 2) = dtv[2];
  DT(2, 1) = DT(1, 2) = dtv[4];
  vnl_symmetric_eigensystem<double> eig(DT);
  double                            e1 = (eig.D(0, 0));
  double                            e2 = (eig.D(1, 1));
  double                            e3 = (eig.D(2, 2));
  if (e1 < 0)
  {
    e1 = e2;
  }
  if (e3 < 0)
  {
    e3 = e2;
  }
  // compute variance of e's
  double emean = (e1 + e2 + e3) / 3.0;
  float  numer = sqrt((e1 - emean) * (e1 - emean) + (e2 - emean) * (e2 - emean) + (e3 - emean) * (e3 - emean));
  return numer;
}

template <typename TensorType>
float
GetTensorFADenominator(TensorType dtv)
{
  typedef vnl_matrix<double> MatrixType;
  MatrixType                 DT(3, 3);
  DT.fill(0);
  DT(0, 0) = dtv[0];
  DT(1, 1) = dtv[3];
  DT(2, 2) = dtv[5];
  DT(1, 0) = DT(0, 1) = dtv[1];
  DT(2, 0) = DT(0, 2) = dtv[2];
  DT(2, 1) = DT(1, 2) = dtv[4];
  vnl_symmetric_eigensystem<double> eig(DT);
  double                            e1 = (eig.D(0, 0));
  double                            e2 = (eig.D(1, 1));
  double                            e3 = (eig.D(2, 2));
  if (e1 < 0)
  {
    e1 = e2;
  }
  if (e3 < 0)
  {
    e3 = e2;
  }
  // compute variance of e's
  float denom = sqrt(e1 * e1 + e2 * e2 + e3 * e3);
  return denom;
}

template <typename TVectorType, typename TTensorType>
float
GetMetricTensorCost(TVectorType dpath, TTensorType dtv, unsigned int matrixpower)
{
  typedef vnl_matrix<double> MatrixType;
  MatrixType                 DT(3, 3);
  DT.fill(0);
  DT(0, 0) = dtv[0];
  DT(1, 1) = dtv[3];
  DT(2, 2) = dtv[5];
  DT(1, 0) = DT(0, 1) = dtv[1];
  DT(2, 0) = DT(0, 2) = dtv[2];
  DT(2, 1) = DT(1, 2) = dtv[4];

  vnl_symmetric_eigensystem<double> eig(DT);
  double                            e1 = (eig.D(0, 0));
  double                            e2 = (eig.D(1, 1));
  double                            e3 = (eig.D(2, 2));
  double                            etot = e1 + e2 + e3;
  if (itk::Math::FloatAlmostEqual(etot, 0.0))
  {
    etot = 1.0;
  }

  MatrixType vec(3, 1);
  vec(0, 0) = dpath[0];
  vec(1, 0) = dpath[1];
  vec(2, 0) = dpath[2];
  MatrixType inv = vnl_matrix_inverse<double>(DT).inverse();
  for (unsigned int lo = 1; lo < matrixpower; lo++)
  {
    inv = inv * inv;
  }

  MatrixType sol = vec.transpose() * inv * vec;
  float      cost = sol(0, 0); // /etot;

  return sqrt(cost);
}

template <typename TVectorType, typename TTensorType>
TVectorType
ChangeTensorByVector(TVectorType dpath, TTensorType dtv, float epsilon)
{
  typedef vnl_matrix<double> MatrixType;
  MatrixType                 DT(3, 3);
  DT.fill(0);
  DT(0, 0) = dtv[0];
  DT(1, 1) = dtv[3];
  DT(2, 2) = dtv[5];
  DT(1, 0) = DT(0, 1) = dtv[1];
  DT(2, 0) = DT(0, 2) = dtv[2];
  DT(2, 1) = DT(1, 2) = dtv[4];
  vnl_symmetric_eigensystem<double> eig(DT);
  double                            e3 = (eig.D(0, 0));
  double                            e2 = (eig.D(1, 1));
  double                            e1 = (eig.D(2, 2));
  double                            etot = e1 + e2 + e3;
  if (itk::Math::FloatAlmostEqual(etot, 0.0))
  {
    etot = 1.0;
  }

  MatrixType vec(3, 1);
  vec(0, 0) = dpath[0];
  vec(1, 0) = dpath[1];
  vec(2, 0) = dpath[2];

  MatrixType evec1(3, 1); // biggest
  evec1(0, 0) = eig.V(0, 2);
  evec1(1, 0) = eig.V(1, 2);
  evec1(2, 0) = eig.V(2, 2);
  MatrixType evec2(3, 1); // middle
  evec2(0, 0) = eig.V(0, 1);
  evec2(1, 0) = eig.V(1, 1);
  evec2(2, 0) = eig.V(2, 1);
  MatrixType evec3(3, 1); // smallest
  evec3(0, 0) = eig.V(0, 0);
  evec3(1, 0) = eig.V(1, 0);
  evec3(2, 0) = eig.V(2, 0);

  float temp;
  temp = (vec.transpose() * evec1)(0, 0);
  temp = sqrt(temp * temp);
  e1 *= (1.0 - static_cast<double>(epsilon * temp));
  if (e1 < 1.e-11)
  {
    e1 = 1.e-11;
  }

  temp = (vec.transpose() * evec2)(0, 0);
  temp = sqrt(temp * temp);
  e2 *= (1.0 - static_cast<double>(epsilon * temp));
  if (e2 < 1.e-11)
  {
    e2 = 1.e-11;
  }

  temp = (vec.transpose() * evec3)(0, 0);
  temp = sqrt(temp * temp);
  e3 *= (1.0 - static_cast<double>(epsilon * temp));
  if (e3 < 1.e-11)
  {
    e3 = 1.e-11;
  }

  DT = (evec3 * evec3.transpose()) * e3 + (evec2 * evec2.transpose()) * e2 + (evec1 * evec1.transpose()) * e1;

  itk::Vector<float, 6> newtens;
  newtens[0] = DT(0, 0);
  newtens[3] = DT(1, 1);
  newtens[5] = DT(2, 2);
  newtens[1] = DT(0, 1);
  newtens[2] = DT(0, 2);
  newtens[4] = DT(2, 1);

  return newtens;
}

template <typename TTensorType>
float
GetTensorADC(TTensorType dtv, unsigned int opt = 0)
{
  typename TTensorType::ValueType sum = dtv[0] + dtv[3] + dtv[5];

  if (opt <= 1)
  {
    return (sum / static_cast<typename TTensorType::ValueType>(3.0));
  }
  float eps = 1.e-9, mag = 0;
  for (unsigned int jj = 0; jj < 6; jj++)
  {
    float ff = dtv[jj];
    mag += ff * ff;
    if (std::isnan(ff) || std::isinf(ff))
    {
      return itk::NumericTraits<float>::ZeroValue();
    }
  }
  mag = sqrt(mag);

  if (itk::Math::FloatAlmostEqual(dtv[1], itk::NumericTraits<typename TTensorType::ValueType>::ZeroValue()) &&
      itk::Math::FloatAlmostEqual(dtv[2], itk::NumericTraits<typename TTensorType::ValueType>::ZeroValue()) &&
      itk::Math::FloatAlmostEqual(dtv[4], itk::NumericTraits<typename TTensorType::ValueType>::ZeroValue()))
  {
    return itk::NumericTraits<float>::ZeroValue();
  }
  if (mag < eps)
  {
    return itk::NumericTraits<float>::ZeroValue();
  }

  typedef vnl_matrix<double> MatrixType;
  MatrixType                 DT(3, 3);
  DT.fill(0);
  DT(0, 0) = dtv[0];
  DT(1, 1) = dtv[3];
  DT(2, 2) = dtv[5];
  DT(1, 0) = DT(0, 1) = dtv[1];
  DT(2, 0) = DT(0, 2) = dtv[2];
  DT(2, 1) = DT(1, 2) = dtv[4];
  //  if (takelog )std::cout << " TAKING LOG " << std::endl;  elsestd::cout << "TAKING EXP " << std::endl;
  // std::cout << " dtv " << dtv << std::endl;
  vnl_symmetric_eigensystem<double> eig(DT);
  double                            e1 = (eig.D(0, 0));
  double                            e2 = (eig.D(1, 1));
  double                            e3 = (eig.D(2, 2));

  /*
  opt  return
  1    ADC
  2    radial diffusion
  3    e1
  4    e2
  5    e3
  */

  if (opt <= 1)
  {
    return (e1 + e2 + e3) / 3.0;
  }
  else if (opt == 2)
  {
    return (e2 + e1) / 2.0; // radial diffusion
  }
  else if (opt == 3) // e1
  {
    return e1;
  }
  else if (opt == 4) // e2
  {
    return e2;
  }
  else if (opt == 5) // e3
  {
    return e3;
  }
  else
  {
    return (e1 + e2 + e3) / 3.0;
  }
}

template <typename TTensorType>
itk::RGBPixel<unsigned char>
GetTensorRGB(TTensorType dtv)
{
  typedef TTensorType TensorType;

  itk::RGBPixel<unsigned char> zero;
  zero.Fill(0);
  float eps = 1.e-9, mag = 0;
  for (unsigned int jj = 0; jj < 6; jj++)
  {
    float ff = dtv[jj];
    mag += ff * ff;
    if (std::isnan(ff) || std::isinf(ff))
    {
      return zero;
    }
  }
  mag = sqrt(mag);

  itk::RGBPixel<unsigned char> rgb;

  if (itk::Math::FloatAlmostEqual(dtv[1], itk::NumericTraits<typename TTensorType::ValueType>::ZeroValue()) &&
      itk::Math::FloatAlmostEqual(dtv[2], itk::NumericTraits<typename TTensorType::ValueType>::ZeroValue()) &&
      itk::Math::FloatAlmostEqual(dtv[4], itk::NumericTraits<typename TTensorType::ValueType>::ZeroValue()))
  {
    return zero;
  }
  if (mag < eps)
  {
    return zero;
  }

  typedef itk::VariableSizeMatrix<typename TTensorType::ValueType> EigenMatrixType;
  EigenMatrixType                                                  evals(3, 3);
  EigenMatrixType                                                  evecs(3, 3);

  EigenAnalysis<TensorType, EigenMatrixType>(dtv, evals, evecs);
  float fa = GetTensorFA<TensorType>(dtv);

  rgb[0] = (unsigned char)(std::fabs(evecs(0, 2)) * fa * 255);
  rgb[1] = (unsigned char)(std::fabs(evecs(1, 2)) * fa * 255);
  rgb[2] = (unsigned char)(std::fabs(evecs(2, 2)) * fa * 255);

  return rgb;
}

template <typename TTensorType>
itk::RGBPixel<float>
GetTensorPrincipalEigenvector(TTensorType dtv)
{
  itk::RGBPixel<float> zero;

  zero.Fill(0);
  float eps = 1.e-9, mag = 0;
  for (unsigned int jj = 0; jj < 6; jj++)
  {
    float ff = dtv[jj];
    mag += ff * ff;
    if (std::isnan(ff) || std::isinf(ff))
    {
      return zero;
    }
  }
  mag = sqrt(mag);

  if (itk::Math::FloatAlmostEqual(dtv[1], itk::NumericTraits<typename TTensorType::ValueType>::ZeroValue()) &&
      itk::Math::FloatAlmostEqual(dtv[2], itk::NumericTraits<typename TTensorType::ValueType>::ZeroValue()) &&
      itk::Math::FloatAlmostEqual(dtv[4], itk::NumericTraits<typename TTensorType::ValueType>::ZeroValue()))
  {
    return zero;
  }
  if (mag < eps)
  {
    return zero;
  }

  typedef vnl_matrix<double> MatrixType;
  MatrixType                 DT(3, 3);
  DT.fill(0);
  DT(0, 0) = dtv[0];
  DT(1, 1) = dtv[3];
  DT(2, 2) = dtv[5];
  DT(1, 0) = DT(0, 1) = dtv[1];
  DT(2, 0) = DT(0, 2) = dtv[2];
  DT(2, 1) = DT(1, 2) = dtv[4];
  //  if (takelog )std::cout << " TAKING LOG " << std::endl;  elsestd::cout << "TAKING EXP " << std::endl;
  // std::cout << " dtv " << dtv << std::endl;
  vnl_symmetric_eigensystem<double> eig(DT);

  itk::RGBPixel<float> rgb;

  // float xx = dtv[0];
  // float xy = dtv[1];
  // float xz = dtv[2];
  // float yy = dtv[3];
  // float yz = dtv[4];
  // float zz = dtv[5];
  // float isp = (xx * xx + yy * yy + zz * zz + 2.0f * (xy * xy + xz * xz + yz * yz));

  // float fa = 0.0;
  // if (isp > 0.0f)
  // {
  //   float trace = dtv[0];
  //   trace += dtv[3];
  //   trace += dtv[5];
  //   float anisotropy = 3.0f * isp - trace * trace;
  //   fa = (std::sqrt(anisotropy / (2.0f * isp)));
  // }

  // rgb[0]=eig.V(2,0)*fa*255;//+eig.V(1,0)*e2;
  // rgb[1]=eig.V(2,1)*fa*255;//+eig.V(1,1)*e2;
  //  rgb[2]=eig.V(2,2)*fa*255;//+eig.V(1,2)*e2;

  // biggest evec
  rgb[0] = eig.V(0, 2); // +eig.V(1,0)*e2;
  rgb[1] = eig.V(1, 2); // +eig.V(1,1)*e2;
  rgb[2] = eig.V(2, 2); // +eig.V(1,2)*e2;

  return rgb;
  mag = rgb[0] * rgb[0] + rgb[1] * rgb[1] + rgb[2] * rgb[2];

  mag = sqrt(mag);
  rgb[0] = rgb[0] / mag;
  rgb[1] = rgb[1] / mag;
  rgb[2] = rgb[2] / mag;

  return rgb;
}

template <typename TTensorType>
itk::Vector<float>
GetTensorPrincipalEigenvector(TTensorType dtv, unsigned int whichvec)
{
  itk::Vector<float, 3> zero;

  zero.Fill(0);
  float eps = 1.e-9, mag = 0;
  for (unsigned int jj = 0; jj < 6; jj++)
  {
    float ff = dtv[jj];
    mag += ff * ff;
    if (std::isnan(ff) || std::isinf(ff))
    {
      return zero;
    }
  }
  mag = sqrt(mag);

  if (itk::Math::FloatAlmostEqual(dtv[1], itk::NumericTraits<typename TTensorType::ValueType>::ZeroValue()) &&
      itk::Math::FloatAlmostEqual(dtv[2], itk::NumericTraits<typename TTensorType::ValueType>::ZeroValue()) &&
      itk::Math::FloatAlmostEqual(dtv[4], itk::NumericTraits<typename TTensorType::ValueType>::ZeroValue()))
  {
    return zero;
  }
  if (mag < eps)
  {
    return zero;
  }

  typedef vnl_matrix<double> MatrixType;
  MatrixType                 DT(3, 3);
  DT.fill(0);
  DT(0, 0) = dtv[0];
  DT(1, 1) = dtv[3];
  DT(2, 2) = dtv[5];
  DT(1, 0) = DT(0, 1) = dtv[1];
  DT(2, 0) = DT(0, 2) = dtv[2];
  DT(2, 1) = DT(1, 2) = dtv[4];
  //  if (takelog )std::cout << " TAKING LOG " << std::endl;  else std::cout << "TAKING EXP " << std::endl;
  // std::cout << " dtv " << dtv << std::endl;
  vnl_symmetric_eigensystem<double> eig(DT);

  itk::Vector<float, 3> rgb;

  float xx = dtv[0];
  float xy = dtv[1];
  float xz = dtv[2];
  float yy = dtv[3];
  float yz = dtv[4];
  float zz = dtv[5];
  float isp = (xx * xx + yy * yy + zz * zz + 2.0f * (xy * xy + xz * xz + yz * yz));

  if (isp > 0.0f)
  {
    float trace = dtv[0];
    trace += dtv[3];
    trace += dtv[5];
  }

  rgb[0] = eig.V(0, whichvec); // +eig.V(1,0)*e2;
  rgb[1] = eig.V(1, whichvec); // +eig.V(1,1)*e2;
  rgb[2] = eig.V(2, whichvec); // +eig.V(1,2)*e2;

  return rgb;
}

template <typename TTensorType>
static float
GetMetricTensorCost(itk::Vector<float, 3> dpath, TTensorType dtv)
{
  typedef vnl_matrix<double> MatrixType;
  MatrixType                 DT(3, 3);
  DT.fill(0);
  DT(0, 0) = dtv[0];
  DT(1, 1) = dtv[3];
  DT(2, 2) = dtv[5];
  DT(1, 0) = DT(0, 1) = dtv[1];
  DT(2, 0) = DT(0, 2) = dtv[2];
  DT(2, 1) = DT(1, 2) = dtv[4];
  vnl_symmetric_eigensystem<double> eig(DT);
  double                            e1 = (eig.D(0, 0));
  double                            e2 = (eig.D(1, 1));
  double                            e3 = (eig.D(2, 2));
  double                            etot = e1 + e2 + e3;
  if (itk::Math::AlmostEquals(etot, 0.0))
  {
    etot = 1;
  }

  MatrixType vec(3, 1);
  vec(0, 0) = dpath[0];
  vec(1, 0) = dpath[1];
  vec(2, 0) = dpath[2];
  MatrixType inv = vnl_matrix_inverse<double>(DT).inverse();

  MatrixType sol = vec.transpose() * inv * vec;
  float      cost = sol(0, 0) / etot;

  return cost;
}

template <typename TensorType, typename VersorTensorType, typename MatrixType>
VersorTensorType
VersorTensor(TensorType dtv)
{
  typedef itk::VariableSizeMatrix<typename MatrixType::ValueType> EigenMatrixType;
  typedef itk::Vector<typename MatrixType::ValueType, 3>          VectorType;

  EigenMatrixType D;
  EigenMatrixType V;
  EigenAnalysis<TensorType, EigenMatrixType>(dtv, D, V);

  VectorType e3;
  e3[0] = V(0, 2);
  e3[1] = V(1, 2);
  e3[2] = V(2, 2);
  VectorType e2;
  e2[0] = V(0, 1);
  e2[1] = V(1, 1);
  e2[2] = V(2, 1);

  VectorType xAxis;
  xAxis[0] = 1;
  xAxis[1] = 0;
  xAxis[2] = 0;
  VectorType yAxis;
  yAxis[0] = 0;
  yAxis[1] = 1;
  yAxis[2] = 0;

  MatrixType R1 = RotationMatrixFromVectors<VectorType, MatrixType>(e3, xAxis);
  e2 = R1 * e2;
  MatrixType R2 = RotationMatrixFromVectors<VectorType, MatrixType>(e2, yAxis);

  MatrixType R = R2 * R1;

  itk::Versor<typename MatrixType::ValueType> versor;
  versor.SetMatrix(R);

  VersorTensorType dtq;
  dtq[0] = D(0, 0);
  dtq[1] = D(1, 1);
  dtq[2] = D(2, 2);
  dtq[3] = versor[0];
  dtq[4] = versor[1];
  dtq[5] = versor[2];

  return dtq;
}

template <typename TensorType, typename VersorTensorType, typename MatrixType>
VersorTensorType
VersorTensor(TensorType dtv, MatrixType frame)
{
  typedef itk::VariableSizeMatrix<typename MatrixType::ValueType> EigenMatrixType;
  typedef itk::Vector<typename MatrixType::ValueType, 3>          VectorType;

  EigenMatrixType D;
  EigenMatrixType V;
  EigenAnalysis<TensorType, EigenMatrixType>(dtv, D, V);

  VectorType e3;
  e3[0] = V(0, 2);
  e3[1] = V(1, 2);
  e3[2] = V(2, 2);
  VectorType e2;
  e2[0] = V(0, 1);
  e2[1] = V(1, 1);
  e2[2] = V(2, 1);

  VectorType xAxis;
  xAxis[0] = frame(0, 0);
  xAxis[1] = frame(1, 0);
  xAxis[2] = frame(2, 0);
  VectorType yAxis;
  yAxis[0] = frame(0, 1);
  yAxis[1] = frame(1, 1);
  yAxis[2] = frame(2, 1);

  MatrixType R1 = RotationMatrixFromVectors<VectorType, MatrixType>(e3, xAxis);
  e2 = R1 * e2;
  MatrixType R2 = RotationMatrixFromVectors<VectorType, MatrixType>(e2, yAxis);

  MatrixType R = R2 * R1;

  itk::Versor<typename MatrixType::ValueType> versor;
  versor.SetMatrix(R);

  VersorTensorType dtq;
  dtq[0] = D(0, 0);
  dtq[1] = D(1, 1);
  dtq[2] = D(2, 2);
  dtq[3] = versor[0];
  dtq[4] = versor[1];
  dtq[5] = versor[2];

  return dtq;
}

#endif
