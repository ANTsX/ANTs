/*=========================================================================


  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
  https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt
  for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include <vnl/vnl_random.h>
#include <vnl/vnl_trace.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_generalized_eigensystem.h>

namespace itk
{
namespace ants
{
template <typename TInputImage, typename TRealType>
antsMatrixUtilities<TInputImage, TRealType>::antsMatrixUtilities()
{
  this->m_Debug = false;
  this->m_PinvTolerance = 1.e-2;
  this->m_PercentVarianceForPseudoInverse = 0.9;
  this->m_MaskImageP = nullptr;
  this->m_MaskImageQ = nullptr;
  this->m_MaskImageR = nullptr;
}

template <typename TInputImage, typename TRealType>
typename antsMatrixUtilities<TInputImage, TRealType>::MatrixType
antsMatrixUtilities<TInputImage, TRealType>::NormalizeMatrix(
  typename antsMatrixUtilities<TInputImage, TRealType>::MatrixType p)
{
  vnl_random randgen(time(nullptr));
  MatrixType np(p.rows(), p.columns());

  for (unsigned long i = 0; i < p.columns(); i++)
  {
    VectorType wpcol = p.get_column(i);
    VectorType wpcol2 = wpcol - wpcol.mean();
    double     sd = wpcol2.squared_magnitude();
    sd = sqrt(sd / (p.rows() - 1));
    if (sd <= 0)
    {
      if (this->m_Debug)
      {
        std::cout << " bad-row " << i << wpcol << std::endl;
      }
      for (unsigned long j = 0; j < wpcol.size(); j++)
      {
        wpcol2(j) = randgen.drand32();
      }
      wpcol2 = wpcol2 - wpcol2.mean();
      sd = wpcol2.squared_magnitude();
      sd = sqrt(sd / (p.rows() - 1));
    }
    wpcol = wpcol2 / sd;
    np.set_column(i, wpcol);
  }
  return np;
}

template <typename TInputImage, typename TRealType>
typename antsMatrixUtilities<TInputImage, TRealType>::MatrixType
antsMatrixUtilities<TInputImage, TRealType>::VNLPseudoInverse(
  typename antsMatrixUtilities<TInputImage, TRealType>::MatrixType rin,
  bool                                                             take_sqrt)
{
  double       pinvTolerance = this->m_PinvTolerance;
  MatrixType   dd = rin;
  unsigned int ss = dd.rows();

  if (dd.rows() > dd.columns())
  {
    ss = dd.columns();
  }
  vnl_svd<RealType> eig(dd, pinvTolerance);
  for (unsigned int j = 0; j < ss; j++)
  {
    RealType eval = eig.W(j, j);
    if (eval > pinvTolerance) // FIXME -- check tolerances against matlab pinv
    {
      eig.W(j, j) = 1 / (eval); // need eval for inv cov
      if (take_sqrt)
      {
        eig.W(j, j) = 1 / sqrt(eval);
      }
    }
    else
    {
      eig.W(j, j) = 0;
    }
  }
  return (eig.recompose()).transpose();
}

template <typename TInputImage, typename TRealType>
typename antsMatrixUtilities<TInputImage, TRealType>::VectorType
antsMatrixUtilities<TInputImage, TRealType>::GetCovMatEigenvector(
  typename antsMatrixUtilities<TInputImage, TRealType>::MatrixType rin,
  unsigned int                                                     which)
{
  double     pinvTolerance = this->m_PinvTolerance;
  MatrixType dd = this->NormalizeMatrix(rin);
  MatrixType cov = dd * dd.transpose();

  cov.set_identity();
  TRealType regularization = 1.e-3;
  cov = cov * regularization + rin * rin.transpose();
  vnl_svd<RealType> eig(cov, pinvTolerance);
  VectorType        vec1 = eig.U().get_column(which);
  VectorType        vec2 = eig.V().get_column(which);
  //  std::cout <<" W 1 " << eig.W(0,0) << " W 2 " << eig.W(1,1) << std::endl;
  if (vec2.size() == rin.rows())
  {
    return vec2;
  }
  else
  {
    return vec1;
  }
}

template <typename TInputImage, typename TRealType>
typename antsMatrixUtilities<TInputImage, TRealType>::MatrixType
antsMatrixUtilities<TInputImage, TRealType>::GetCovMatEigenvectors(
  typename antsMatrixUtilities<TInputImage, TRealType>::MatrixType rin)
{
  double     pinvTolerance = this->m_PinvTolerance;
  MatrixType dd = this->NormalizeMatrix(rin);
  MatrixType cov = dd * dd.transpose();

  cov.set_identity();
  TRealType regularization = 1.e-3;
  cov = cov * regularization + rin * rin.transpose();
  vnl_svd<RealType> eig(cov, pinvTolerance);
  VectorType        vec1 = eig.U().get_column(0);
  VectorType        vec2 = eig.V().get_column(0);
  double            trace = vnl_trace<double>(cov);
  double            evalsum = 0;
  for (unsigned int i = 0; i < cov.rows(); i++)
  {
    evalsum += eig.W(i, i);
    std::cout << " variance-explained-eval " << i << " = " << evalsum / trace * 100 << std::endl;
  }
  //  std::cout <<" W 1 " << eig.W(0,0) << " W 2 " << eig.W(1,1) << std::endl;
  if (vec2.size() == rin.rows())
  {
    return eig.V();
  }
  else
  {
    return eig.U();
  }
}
} // namespace ants
} // namespace itk
