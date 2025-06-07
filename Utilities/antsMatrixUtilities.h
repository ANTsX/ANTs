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
#ifndef __antsMatrixUtilities_h
#define __antsMatrixUtilities_h
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_cholesky.h>
#include "itkImageToImageFilter.h"
namespace itk
{
namespace ants
{
template <typename TInputImage, typename TRealType = double>
class antsMatrixUtilities final : public ImageToImageFilter<TInputImage, TInputImage>
{
public:
  /** Standard class typdedefs. */
  typedef antsMatrixUtilities                          Self;
  typedef ImageToImageFilter<TInputImage, TInputImage> Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(antsMatrixUtilities);

  /** Dimension of the images. */
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;

  static constexpr unsigned int MatrixDimension = 2;

  /** Typedef support of input types. */
  typedef TInputImage                   ImageType;
  typedef typename ImageType::Pointer   ImagePointer;
  typedef typename ImageType::PixelType PixelType;
  typedef typename ImageType::IndexType IndexType;

  /** Some convenient typedefs. */
  typedef TRealType                                               RealType;
  typedef Image<RealType, Self::ImageDimension> RealImageType;

  /** note, eigen for pseudo-eigenvals  */
  typedef vnl_matrix<RealType>      MatrixType;
  typedef vnl_vector<RealType>      VectorType;
  typedef MatrixType                VariateType;
  typedef vnl_diag_matrix<RealType> DiagonalMatrixType;

  void
  NormalizeWeightsByCovariance();

  void
  SetPseudoInversePercentVariance(RealType p)
  {
    this->m_PercentVarianceForPseudoInverse = p;
  }

  MatrixType
  PseudoInverse(MatrixType p_in, bool take_sqrt = false)
  {
    return this->VNLPseudoInverse(p_in, take_sqrt);
  }

  MatrixType
  VNLPseudoInverse(MatrixType, bool take_sqrt = false);

  VectorType
  Orthogonalize(VectorType Mvec, VectorType V, MatrixType * projecterM = nullptr, MatrixType * projecterV = nullptr)
  {
    if (!projecterM && !projecterV)
    {
      double     ratio = inner_product(Mvec, V) / inner_product(V, V);
      VectorType ortho = Mvec - V * ratio;
      return ortho;
    }
    else if (!projecterM && projecterV)
    {
      double     ratio = inner_product(Mvec, *projecterV * V) / inner_product(*projecterV * V, *projecterV * V);
      VectorType ortho = Mvec - V * ratio;
      return ortho;
    }
    else if (!projecterV && projecterM)
    {
      double     ratio = inner_product(*projecterM * Mvec, V) / inner_product(V, V);
      VectorType ortho = (*projecterM * Mvec) - V * ratio;
      return ortho;
    }
    else
    {
      double ratio =
        inner_product(*projecterM * Mvec, *projecterV * V) / inner_product(*projecterV * V, *projecterV * V);
      VectorType ortho = Mvec - V * ratio;
      return ortho;
    }
  }

  MatrixType
  OrthogonalizeMatrix(MatrixType M, VectorType V)
  {
    for (unsigned int j = 0; j < M.cols(); j++)
    {
      VectorType Mvec = M.get_column(j);
      double     ratio = inner_product(Mvec, V) / inner_product(V, V);
      VectorType ortho = Mvec - V * ratio;
      M.set_column(j, ortho);
    }
    return M;
  }

  void
  SetMaskImageP(ImagePointer mask)
  {
    this->m_MaskImageP = mask;
  }

  void
  SetMatrixP(MatrixType matrix)
  {
    this->m_OriginalMatrixP.set_size(matrix.rows(), matrix.cols());
    this->m_MatrixP.set_size(matrix.rows(), matrix.cols());
    this->m_OriginalMatrixP.update(matrix);
    this->m_MatrixP.update(matrix);
  }

  itkSetMacro(FractionNonZeroQ, RealType);
  itkSetMacro(KeepPositiveQ, bool);
  void
  SetMaskImageQ(ImagePointer mask)
  {
    this->m_MaskImageQ = mask;
  }

  void
  SetMatrixQ(MatrixType matrix)
  {
    this->m_OriginalMatrixQ.set_size(matrix.rows(), matrix.cols());
    this->m_MatrixQ.set_size(matrix.rows(), matrix.cols());
    this->m_OriginalMatrixQ.update(matrix);
    this->m_MatrixQ.update(matrix);
  }

  itkSetMacro(FractionNonZeroR, RealType);
  itkSetMacro(KeepPositiveR, bool);
  void
  SetMaskImageR(ImagePointer mask)
  {
    this->m_MaskImageR = mask;
  }

  void
  SetMatrixR(MatrixType matrix)
  {
    this->m_OriginalMatrixR.set_size(matrix.rows(), matrix.cols());
    this->m_MatrixR.set_size(matrix.rows(), matrix.cols());
    this->m_OriginalMatrixR.update(matrix);
    this->m_MatrixR.update(matrix);
  }

  MatrixType
  GetMatrixP()
  {
    return this->m_MatrixP;
  }

  MatrixType
  GetMatrixQ()
  {
    return this->m_MatrixQ;
  }

  MatrixType
  GetMatrixR()
  {
    return this->m_MatrixR;
  }

  MatrixType
  GetOriginalMatrixP()
  {
    return this->m_OriginalMatrixP;
  }

  MatrixType
  GetOriginalMatrixQ()
  {
    return this->m_OriginalMatrixQ;
  }

  MatrixType
  GetOriginalMatrixR()
  {
    return this->m_OriginalMatrixR;
  }

  VectorType
  InitializeV(MatrixType p);

  MatrixType
  NormalizeMatrix(MatrixType p);

  MatrixType
  CovarianceMatrix(MatrixType p, RealType regularization = 1.e-2)
  {
    if (p.rows() < p.columns())
    {
      MatrixType invcov = p * p.transpose();
      invcov.set_identity();
      invcov = invcov * regularization + p * p.transpose();
      return invcov;
    }
    else
    {
      MatrixType invcov = p.transpose() * p;
      invcov.set_identity();
      invcov = invcov * regularization + p.transpose() * p;
      return invcov;
    }
  }

  VectorType
  GetCovMatEigenvector(MatrixType p, unsigned int evec);

  MatrixType
  GetCovMatEigenvectors(MatrixType p);

  VectorType
  AverageColumns(MatrixType p)
  {
    unsigned int ncol = p.columns();
    VectorType   v = p.get_column(0);

    v.fill(0);
    for (unsigned int i = 0; i < ncol; i++)
    {
      v = v + p.get_column(i);
    }
    return v / (RealType)ncol;
  }

  MatrixType
  WhitenMatrix(MatrixType p, RealType regularization = 1.e-2)
  {
    MatrixType invcov = this->CovarianceMatrix(p, regularization);

    invcov = this->PseudoInverse(invcov, true);
    if (p.rows() < p.columns())
    {
      return invcov * p;
    }
    else
    {
      return p * invcov;
    }
  }

  MatrixType
  WhitenMatrixByAnotherMatrix(MatrixType p, MatrixType op, RealType regularization = 1.e-2)
  {
    MatrixType invcov = this->CovarianceMatrix(op, regularization);

    invcov = this->PseudoInverse(invcov, true);
    if (p.rows() < p.columns())
    {
      return invcov * p;
    }
    else
    {
      return p * invcov;
    }
  }

  MatrixType
  ProjectionMatrix(MatrixType b)
  {
    b = this->NormalizeMatrix(b);
    b = this->WhitenMatrix(b);
    return b * b.transpose();
  }

  MatrixType
  DeleteCol(MatrixType p_in, unsigned int col)
  {
    unsigned int ncols = p_in.cols() - 1;

    if (col >= ncols)
    {
      ncols = p_in.cols();
    }
    MatrixType   p(p_in.rows(), ncols);
    unsigned int colct = 0;
    for (long i = 0; i < p.cols(); ++i) // loop over cols
    {
      if (i != col)
      {
        p.set_column(colct, p_in.get_column(i));
        colct++;
      }
    }
    return p;
  }

  RealType
  PearsonCorr(VectorType v1, VectorType v2)
  {
    double xysum = 0;

    for (unsigned int i = 0; i < v1.size(); i++)
    {
      xysum += v1(i) * v2(i);
    }
    double frac = 1.0 / (double)v1.size();
    double xsum = v1.sum(), ysum = v2.sum();
    double xsqr = v1.squared_magnitude();
    double ysqr = v2.squared_magnitude();
    double numer = xysum - frac * xsum * ysum;
    double denom = sqrt((xsqr - frac * xsum * xsum) * (ysqr - frac * ysum * ysum));
    if (denom <= 0)
    {
      return 0;
    }
    return numer / denom;
  }

  antsMatrixUtilities();
  ~antsMatrixUtilities() override = default;

  void
  PrintSelf(std::ostream & os, Indent indent) const override
  {
    os << indent;
  }

private:
  bool       m_Debug;
  MatrixType m_OriginalMatrixP;
  MatrixType m_OriginalMatrixQ;
  MatrixType m_OriginalMatrixR;

  antsMatrixUtilities(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented

  RealType m_PinvTolerance;
  RealType m_PercentVarianceForPseudoInverse;

  MatrixType   m_MatrixP;
  ImagePointer m_MaskImageP;
  RealType     m_FractionNonZeroP;
  bool         m_KeepPositiveP;

  MatrixType   m_MatrixQ;
  ImagePointer m_MaskImageQ;
  RealType     m_FractionNonZeroQ;
  bool         m_KeepPositiveQ;
  MatrixType   m_MatrixR;
  ImagePointer m_MaskImageR;
  RealType     m_FractionNonZeroR;
  bool         m_KeepPositiveR;
};
} // namespace ants
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "antsMatrixUtilities.hxx"
#endif

#endif
