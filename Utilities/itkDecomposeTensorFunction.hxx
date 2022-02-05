/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDecomposeTensorFunction_hxx
#define __itkDecomposeTensorFunction_hxx


#include "vnl/algo/vnl_cholesky.h"
#include "vnl/algo/vnl_qr.h"
#include "vnl/algo/vnl_real_eigensystem.h"
#include "vnl/algo/vnl_svd.h"
#include "vnl/algo/vnl_svd_economy.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/vnl_det.h"

namespace itk
{
template <typename TInput, typename TRealType, typename TOutput>
DecomposeTensorFunction<TInput, TRealType, TOutput>::DecomposeTensorFunction() = default;

template <typename TInput, typename TRealType, typename TOutput>
void
DecomposeTensorFunction<TInput, TRealType, TOutput>::EvaluateEigenDecomposition(InputMatrixType &  M,
                                                                                OutputMatrixType & D,
                                                                                OutputMatrixType & V)
{
  unsigned int RowDimensions = M.Rows();
  unsigned int ColumnDimensions = M.Cols();

  D.SetSize(RowDimensions, RowDimensions);
  V.SetSize(RowDimensions, RowDimensions);

  D.Fill(0.0);
  vnl_matrix<double> v(RowDimensions, ColumnDimensions);
  for (unsigned int j = 0; j < ColumnDimensions; j++)
  {
    for (unsigned int i = 0; i < RowDimensions; i++)
    {
      v.put(i, j, M[i][j]);
    }
  }
  vnl_real_eigensystem eig(v);
  for (unsigned int j = 0; j < ColumnDimensions; j++)
  {
    for (unsigned int i = 0; i < RowDimensions; i++)
    {
      V[i][j] = static_cast<RealType>((eig.Vreal).get(i, j));
      if (i == j)
      {
        D[i][j] = static_cast<RealType>(std::real(eig.D(j)));
      }
    }
  }
}

template <typename TInput, typename TRealType, typename TOutput>
void
DecomposeTensorFunction<TInput, TRealType, TOutput>::EvaluateSymmetricEigenDecomposition(InputMatrixType &  M,
                                                                                         OutputMatrixType & D,
                                                                                         OutputMatrixType & V)
{
  //    The resulting eigenvectors and values are sorted in increasing order
  //    so <CODE> V.column(0) </CODE> is the eigenvector corresponding to the
  //    smallest eigenvalue.

  unsigned int RowDimensions = M.Rows();
  unsigned int ColumnDimensions = M.Cols();

  D.SetSize(RowDimensions, RowDimensions);
  V.SetSize(RowDimensions, RowDimensions);

  D.Fill(0.0);
  vnl_symmetric_eigensystem<RealType> eig(M.GetVnlMatrix());
  for (unsigned int j = 0; j < ColumnDimensions; j++)
  {
    for (unsigned int i = 0; i < RowDimensions; i++)
    {
      V[i][j] = (eig.V).get(i, j);
      if (i == j)
      {
        D[i][j] = eig.D(j);
      }
    }
  }
}

template <typename TInput, typename TRealType, typename TOutput>
void
DecomposeTensorFunction<TInput, TRealType, TOutput>::EvaluateRightPolarDecomposition(InputMatrixType &  M,
                                                                                     OutputMatrixType & R,
                                                                                     OutputMatrixType & S)
{
  OutputMatrixType U;
  OutputMatrixType W;
  OutputMatrixType V;

  this->EvaluateSVDDecomposition(M, U, W, V);

  R = U * V.GetTranspose();
  S = V * W * V.GetTranspose();
}

template <typename TInput, typename TRealType, typename TOutput>
void
DecomposeTensorFunction<TInput, TRealType, TOutput>::EvaluateLeftPolarDecomposition(InputMatrixType &  M,
                                                                                    OutputMatrixType & S,
                                                                                    OutputMatrixType & R)
{
  OutputMatrixType U;
  OutputMatrixType W;
  OutputMatrixType V;

  this->EvaluateSVDDecomposition(M, U, W, V);

  R = U * V.GetTranspose();
  S = U * W * U.GetTranspose();
}

template <typename TInput, typename TRealType, typename TOutput>
void
DecomposeTensorFunction<TInput, TRealType, TOutput>::EvaluateQRDecomposition(InputMatrixType &  M,
                                                                             OutputMatrixType & Q,
                                                                             OutputMatrixType & R)
{
  vnl_qr<RealType> qr(M.GetVnlMatrix());

  Q.SetSize(qr.Q().rows(), qr.Q().cols());
  R.SetSize(qr.R().rows(), qr.R().cols());
  for (unsigned int i = 0; i < Q.Rows(); i++)
  {
    for (unsigned int j = 0; j < Q.Cols(); j++)
    {
      Q[i][j] = qr.Q()(i, j);
    }
  }
  for (unsigned int i = 0; i < R.Rows(); i++)
  {
    for (unsigned int j = 0; j < R.Cols(); j++)
    {
      R[i][j] = qr.R()(i, j);
    }
  }
}

template <typename TInput, typename TRealType, typename TOutput>
void
DecomposeTensorFunction<TInput, TRealType, TOutput>::EvaluateSVDDecomposition(InputMatrixType &  M,
                                                                              OutputMatrixType & U,
                                                                              OutputMatrixType & W,
                                                                              OutputMatrixType & V)
{
  unsigned int RowDimensions = M.Rows();
  unsigned int ColumnDimensions = M.Cols();

  U.SetSize(RowDimensions, RowDimensions);
  V.SetSize(ColumnDimensions, ColumnDimensions);
  W.SetSize(RowDimensions, ColumnDimensions);

  vnl_svd<RealType> svd(M.GetVnlMatrix());
  for (unsigned int i = 0; i < RowDimensions; i++)
  {
    for (unsigned int j = 0; j < RowDimensions; j++)
    {
      U[i][j] = svd.U(i, j);
    }
  }
  for (unsigned int i = 0; i < ColumnDimensions; i++)
  {
    for (unsigned int j = 0; j < ColumnDimensions; j++)
    {
      V[i][j] = svd.V(i, j);
    }
  }

  W.Fill(0.0);
  unsigned int minDimensions = ColumnDimensions;
  if (static_cast<unsigned int>(RowDimensions) < static_cast<unsigned int>(ColumnDimensions))
  {
    minDimensions = RowDimensions;
  }
  for (unsigned int i = 0; i < minDimensions; i++)
  {
    W[i][i] = svd.W(i, i);
  }
}

template <typename TInput, typename TRealType, typename TOutput>
void
DecomposeTensorFunction<TInput, TRealType, TOutput>::EvaluateSVDEconomyDecomposition(InputMatrixType &  M,
                                                                                     OutputMatrixType & W,
                                                                                     OutputMatrixType & V)
{
  /**
   * Same as SVD except the routine does not return U ---
   * allows for faster computation.
   */

  unsigned int RowDimensions = M.Rows();
  unsigned int ColumnDimensions = M.Cols();

  V.SetSize(ColumnDimensions, ColumnDimensions);
  W.SetSize(RowDimensions, ColumnDimensions);

  vnl_svd_economy<RealType> svd(M.GetVnlMatrix());
  for (unsigned int i = 0; i < ColumnDimensions; i++)
  {
    for (unsigned int j = 0; j < ColumnDimensions; j++)
    {
      V[i][j] = svd.V()(i, j);
    }
  }

  W.Fill(0.0);
  unsigned int minDimensions = ColumnDimensions;
  if (static_cast<unsigned int>(RowDimensions) < static_cast<unsigned int>(ColumnDimensions))
  {
    minDimensions = RowDimensions;
  }
  for (unsigned int i = 0; i < minDimensions; i++)
  {
    W[i][i] = svd.lambdas()[i];
  }
}

template <typename TInput, typename TRealType, typename TOutput>
void
DecomposeTensorFunction<TInput, TRealType, TOutput>::EvaluateCholeskyDecomposition(InputMatrixType &  M,
                                                                                   OutputMatrixType & L)
{
  // Assumes symmetric tensor of type double
  vnl_matrix<double> m(M.Rows(), M.Cols());
  for (unsigned int j = 0; j < M.Cols(); j++)
  {
    for (unsigned int i = 0; i < M.Rows(); i++)
    {
      m.put(i, j, M[i][j]);
    }
  }

  vnl_cholesky cholesky(m, vnl_cholesky::quiet);

  L.SetSize(M.Rows(), M.Cols());
  for (unsigned int i = 0; i < L.Rows(); i++)
  {
    for (unsigned int j = 0; j < L.Cols(); j++)
    {
      L[i][j] = cholesky.L_badly_named_method().get(i, j);
    }
  }
}

template <typename TInput, typename TRealType, typename TOutput>
typename DecomposeTensorFunction<TInput, TRealType, TOutput>::RealType
DecomposeTensorFunction<TInput, TRealType, TOutput>::EvaluateDeterminant(InputMatrixType & M)
{
  RealType det = 0.0;

  if (M.Rows() == M.Cols() && M.Rows() >= 2 && M.Rows() <= 4)
  {
    vnl_matrix_fixed<RealType, 2, 2> m2;
    vnl_matrix_fixed<RealType, 3, 3> m3;
    vnl_matrix_fixed<RealType, 4, 4> m4;
    for (unsigned int i = 0; i < M.Rows(); i++)
    {
      for (unsigned int j = 0; j < M.Cols(); j++)
      {
        switch (M.Rows())
        {
          case 2:
          {
            m2.put(i, j, M[i][j]);
          }
          break;
          case 3:
          {
            m3.put(i, j, M[i][j]);
          }
          break;
          case 4:
          {
            m4.put(i, j, M[i][j]);
          }
          break;
        }
      }
    }

    switch (M.Rows())
    {
      case 2:
      {
        det = static_cast<RealType>(vnl_det(m2));
      }
      break;
      case 3:
      {
        det = static_cast<RealType>(vnl_det(m3));
      }
      break;
      case 4:
      {
        det = static_cast<RealType>(vnl_det(m4));
      }
      break;
    }
  }
  else
  {
    vnl_qr<RealType> qr(M.GetVnlMatrix());
    det = static_cast<RealType>(qr.determinant());
  }
  return det;
}

template <typename TInput, typename TRealType, typename TOutput>
void
DecomposeTensorFunction<TInput, TRealType, TOutput>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}
} // end namespace itk

#endif
