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
#ifndef __antsSCCANObject_h
#  define __antsSCCANObject_h
// #define EIGEN_DEFAULT_TO_ROW_MAJOR
// #define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
// #include <Eigen/Dense>
// #include <Eigen/Sparse>
// #include <Eigen/SVD>
// #include "armadillo"
#  include <vnl/algo/vnl_matrix_inverse.h>
#  include <vnl/algo/vnl_cholesky.h>
#  include "itkImageToImageFilter.h"

#  include "itkMath.h"

/** Custom SCCA implemented with vnl and ITK: Flexible positivity constraints, image ops, permutation testing, etc. */
namespace itk
{
namespace ants
{
template <typename TInputImage, typename TRealType = double>
class antsSCCANObject final : public ImageToImageFilter<TInputImage, TInputImage>
{
public:
  /** Standard class typdedefs. */
  typedef antsSCCANObject                              Self;
  typedef ImageToImageFilter<TInputImage, TInputImage> Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(antsSCCANObject);

  /** Dimension of the images. */
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;

  static constexpr unsigned int MatrixDimension = 2;

  /** Typedef support of input types. */
  typedef TInputImage                   ImageType;
  typedef typename ImageType::Pointer   ImagePointer;
  typedef typename ImageType::PixelType PixelType;
  typedef typename ImageType::IndexType IndexType;

  /** Some convenient typedefs. */
  typedef TRealType                                                   RealType;
  typedef Image<RealType, Self::ImageDimension>     RealImageType;
  typedef Image<RealType, Self::ImageDimension - 1> RealImageTypeDminus1;

  /** Define eigen types */
  //  typedef Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic> eMatrix;
  //  typedef Eigen::Matrix<RealType, Eigen::Dynamic, 1> eVector;
  // typedef Eigen::DynamicSparseMatrix<RealType,Eigen::RowMajor> sMatrix;
  // typedef Eigen::FullPivHouseholderQR<eMatrix> svdobj2;
  // typedef Eigen::JacobiSVD<eMatrix> svdobj;

  /** note, eigen for pseudo-eigenvals  */
  typedef vnl_matrix<RealType>      MatrixType;
  typedef vnl_vector<RealType>      VectorType;
  typedef MatrixType                VariateType;
  typedef vnl_diag_matrix<RealType> DiagonalMatrixType;

  enum SCCANFormulationType
  {
    PQ,
    PminusRQ,
    PQminusR,
    PminusRQminusR,
    PQR
  };

  /** ivars Set/Get functionality */
  itkSetMacro(MaximumNumberOfIterations, unsigned int);
  itkGetConstMacro(MaximumNumberOfIterations, unsigned int);
  itkSetMacro(MinClusterSizeP, unsigned int);
  itkGetConstMacro(MinClusterSizeP, unsigned int);
  itkSetMacro(MinClusterSizeQ, unsigned int);
  itkGetConstMacro(MinClusterSizeQ, unsigned int);
  itkSetMacro(KeptClusterSize, unsigned int);
  itkGetConstMacro(KeptClusterSize, unsigned int);
  itkSetMacro(AlreadyWhitened, bool);
  itkGetConstMacro(AlreadyWhitened, bool);
  itkSetMacro(PriorWeight, RealType);
  itkGetConstMacro(PriorWeight, RealType);
  itkSetMacro(ConvergenceThreshold, RealType);
  itkGetConstMacro(ConvergenceThreshold, RealType);
  itkGetConstMacro(CurrentConvergenceMeasurement, RealType);
  itkGetConstMacro(ElapsedIterations, unsigned int);
  itkSetMacro(SCCANFormulation, SCCANFormulationType);
  itkGetConstMacro(SCCANFormulation, SCCANFormulationType);
  itkSetMacro(Silent, bool);
  itkGetMacro(Silent, bool);
  itkSetMacro(RowSparseness, RealType);
  itkGetMacro(RowSparseness, RealType);
  itkSetMacro(UseLongitudinalFormulation, RealType);
  itkGetMacro(UseLongitudinalFormulation, RealType);
  itkSetMacro(Smoother, RealType);
  itkGetMacro(Smoother, RealType);

  void
  NormalizeWeights(const unsigned int k);

  void
  NormalizeWeightsByCovariance(const unsigned int k, const TRealType taup = 0, const TRealType tauq = 0);

  void
  WhitenDataSetForRunSCCANMultiple(unsigned int nvecs = 0);

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

  void
  ZeroProduct(VectorType & v1, VectorType & v2)
  {
    for (unsigned int i = 0; i < v1.size(); i++)
    {
      if (fabs(v2(i)) > 1.e-9)
      {
        v1(i) = 0;
      }
    }
  }

  void
  DeleteRow(MatrixType &, unsigned int);

  void
  PosNegVector(VectorType & v1, bool pos)
  {
    for (unsigned int i = 0; i < v1.size(); i++)
    {
      if (v1(i) < 0 && pos)
      {
        v1(i) = 0;
      }
      else if (v1(i) > 0 && !pos)
      {
        v1(i) = 0;
      }
    }
  }

  RealType ReconstructionError(MatrixType, MatrixType);

  VectorType
  Orthogonalize(VectorType Mvec, VectorType V, MatrixType * projecterM = nullptr, MatrixType * projecterV = nullptr)
  {
    if ((!projecterM) && (!projecterV))
    {
      double ipv = inner_product(V, V);
      if (itk::Math::FloatAlmostEqual(ipv, 0.0))
      {
        return Mvec;
      }
      double     ratio = inner_product(Mvec, V) / ipv;
      VectorType ortho = Mvec - V * ratio;
      return ortho;
    }
    else if ((!projecterM) && (projecterV))
    {
      double     ratio = inner_product(Mvec, *projecterV * V) / inner_product(*projecterV * V, *projecterV * V);
      VectorType ortho = Mvec - V * ratio;
      return ortho;
    }
    else if ((!projecterV) && (projecterM))
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
      for (unsigned int i = 0; i < Mvec.size(); i++)
      {
        if (itk::Math::FloatAlmostEqual(Mvec(i), itk::NumericTraits<RealType>::ZeroValue()))
        {
          ortho(i) = 0;
        }
      }
      return ortho;
    }
  }

  MatrixType
  OrthogonalizeMatrix(MatrixType M, VectorType V)
  {
    for (unsigned int j = 0; j < M.cols(); j++)
    {
      VectorType Mvec = M.get_column(j);
      double     vnorm = inner_product(V, V);
      if (vnorm < this->m_Epsilon)
        vnorm = 1;
      double     ratio = inner_product(Mvec, V) / vnorm;
      VectorType ortho = Mvec - V * ratio;
      M.set_column(j, ortho);
    }
    return M;
  }

  MatrixType
  RankifyMatrixColumns(MatrixType M)
  {
    RealType rows = (RealType)M.rows();

    for (unsigned long j = 0; j < M.cols(); j++)
    {
      VectorType Mvec = M.get_column(j);
      VectorType rank = M.get_column(j);
      for (unsigned int i = 0; i < rows; i++)
      {
        double   rankval = 0;
        RealType xi = Mvec(i);
        for (unsigned int k = 0; k < rows; k++)
        {
          RealType yi = Mvec(k);
          RealType diff = fabs(xi - yi);
          if (diff > 0)
          {
            RealType val = (xi - yi) / diff;
            rankval += val;
          }
        }
        rank(i) = rankval / rows;
      }
      M.set_column(j, rank);
    }
    return M;
  }

  itkSetMacro(FractionNonZeroP, RealType);
  itkSetMacro(KeepPositiveP, bool);
  itkGetMacro(KeepPositiveP, bool);
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
  itkGetMacro(KeepPositiveQ, bool);

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

  itkSetMacro(Covering, unsigned int);
  itkSetMacro(GetSmall, bool);
  itkSetMacro(UseL1, bool);
  itkSetMacro(MaxBasedThresholding, bool);
  itkGetMacro(MaxBasedThresholding, bool);
  itkSetMacro(GradStep, RealType);
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
  GetMatrixU()
  {
    return this->m_MatrixU;
  }

  MatrixType
  GetOriginalMatrixP()
  {
    return this->m_OriginalMatrixP;
  }

  RealType
  GetLambda()
  {
    return this->m_lambda;
  }

  RealType
  SetLambda(RealType lambda)
  {
    return this->m_lambda = lambda;
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

  // Prior Constrained PCA
  MatrixType
  GetMatrixPriorROI()
  {
    return this->m_MatrixPriorROI;
  }

  // Prior Constrained PCA
  MatrixType
  GetMatrixPriorROI2()
  {
    return this->m_MatrixPriorROI2;
  }

  void
  SetFlagForSort()
  {
    this->flagForSort = true;
  }

  // Prior Constrained PCA
  MatrixType
  GetOriginalMatrixPriorROI()
  {
    return this->m_OriginalMatrixPriorROI;
  }

  RealType
  RunSCCAN2multiple(unsigned int n_vecs);

  RealType
  RunSCCAN2();

  RealType
  RunSCCAN3();

  RealType
  LineSearch(MatrixType & A,
             VectorType & x_k,
             VectorType & p_k,
             VectorType & b,
             RealType     minalph,
             RealType     maxalph,
             RealType);

  RealType
  EvaluateEnergy(MatrixType & A, VectorType & x_k, VectorType & p_k, VectorType & b, RealType minalph, RealType);

  RealType
  SparseConjGrad(VectorType &, VectorType, RealType, unsigned int);
  RealType
  ConjGrad(MatrixType & A, VectorType & x_k, VectorType b_in, RealType convcrit, unsigned int);

  RealType
  SparseConjGradRidgeRegression(MatrixType & A,
                                VectorType & x_k,
                                VectorType   b_in,
                                RealType     convcrit,
                                unsigned int,
                                bool);

  RealType
  IHTRegression(MatrixType & A,
                VectorType & x_k,
                VectorType   b_in,
                RealType     convcrit,
                unsigned int,
                RealType mu,
                bool     isp = true,
                bool     verbose = false);

  RealType
  MatchingPursuit(MatrixType & A, VectorType & x_k, RealType convcrit, unsigned int);

  RealType
  SparseNLConjGrad(MatrixType & A,
                   VectorType & x_k,
                   VectorType   b,
                   RealType,
                   unsigned int,
                   bool         makeprojsparse = false,
                   unsigned int loorth = 0,
                   unsigned int hiorth = 0);
  RealType
  SparseNLPreConjGrad(MatrixType & A, VectorType & x_k, VectorType b, RealType, unsigned int);
  RealType
  RidgeRegression(MatrixType & A,
                  VectorType & x_k,
                  VectorType   b,
                  RealType     lambda,
                  unsigned int,
                  bool makesparse = false);

  /** Return Rayleigh quotient */
  RealType
  PowerIteration(MatrixType & A, VectorType & x_k, unsigned int, bool);

  RealType
  IHTPowerIteration(MatrixType & A, VectorType & x_k, unsigned int, unsigned int);
  RealType
  IHTPowerIterationU(MatrixType & A, VectorType & x_k, unsigned int, unsigned int);

  RealType
  IHTPowerIterationPrior(MatrixType & A, VectorType & x_k, VectorType & x_k_1, unsigned int, unsigned int, double);

  void
  SoftClustThreshold(VectorType & v_in,
                     RealType     fractional_goal,
                     bool         allow_negative_weights,
                     unsigned int,
                     ImagePointer);
  void
  ReSoftThreshold(VectorType & v_in, RealType fractional_goal, bool allow_negative_weights);

  void
  ConstantProbabilityThreshold(VectorType & v_in, RealType probability_goal, bool allow_negative_weights);

  VectorType
  InitializeV(MatrixType p, unsigned long seed = 0);

  TRealType
  InitializeSCCA_simple(unsigned int n_vecs);

  TRealType
  InitializeSCCA(unsigned int n_vecs, unsigned int seeder);

  TRealType
  InitializeSPCA(unsigned int n_vecs, unsigned int seeder, MatrixType &, VectorType &);

  VectorType ComputeVectorLaplacian(VectorType, ImagePointer);
  VectorType ComputeVectorGradMag(VectorType, ImagePointer);
  VectorType
  SpatiallySmoothVector(VectorType, ImagePointer, bool surface = true);

  void
  SetSortFinalLocArray(VectorType locArray)
  {
    this->loc_Array = locArray;
  }

  VectorType
  GetSortFinalLocArray()
  {
    return this->loc_Array;
  }

  MatrixType
  NormalizeMatrix(MatrixType p, bool makepositive = true);

  /** needed for partial scca */
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

  MatrixType
  WhitenMatrix(MatrixType p, RealType regularization = 1.e-2)
  {
    double reg = 1.e-9;

    if (p.rows() < p.cols())
    {
      reg = regularization;
    }
    MatrixType cov = this->CovarianceMatrix(p, reg);
    MatrixType invcov = this->PseudoInverse(cov, true);
    bool       debug = false;
    if ((!this->m_Silent) && (debug))
    {
      std::cout << " cov " << std::endl;
      std::cout << cov << std::endl;
      std::cout << " invcov " << std::endl;
      std::cout << invcov << std::endl;
      std::cout << " id? " << std::endl;
      std::cout << cov * invcov << std::endl;
    }
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
  ProjectionMatrix(MatrixType b, double regularization = 0.001)
  {
    bool armadillo = false;

    b = this->NormalizeMatrix(b);
    MatrixType mat = b * b.transpose();
    if (!armadillo)
    {
      MatrixType cov(mat.rows(), mat.cols(), 0);
      cov.set_identity();
      mat = cov * regularization + mat;
      return vnl_svd<double>(mat).inverse();
      //      return vnl_svd<double>( mat ).pinverse( mc );
    }
    else
    { /*
  arma::mat amat( b.rows(), b.rows());
  for ( unsigned int i = 0 ; i < b.rows(); i++ )
for ( unsigned int j = 0 ; j < b.rows(); j++ )
{
amat( i , j ) = mat( i , j );
if ( i == j ) amat( i , j ) += regularization;
}
  arma::mat invamat = arma::inv( amat , 1.e-2 );

  for ( unsigned int i = 0 ; i < b.rows(); i++ )
for ( unsigned int j = 0 ; j < b.rows(); j++ ) mat( i , j ) = invamat( i , j );
   */
      return mat;
    }
  }

  VectorType
  TrueCCAPowerUpdate(RealType penaltyP, MatrixType p, VectorType w_q, MatrixType q, bool keep_pos, bool factorOutR);

  MatrixType PartialOutZ(MatrixType /*X*/, MatrixType /*Y*/, MatrixType /*Z*/)
  {
    if (!this->m_Silent)
      std::cout << "ERROR:  This function not yet implemented." << std::endl;
    /** compute the effect of Z and store it for later use */
  }

  // Prior Constrained PCA

  void
  SetMatrixPriorROI(MatrixType matrix)
  {
    this->m_OriginalMatrixPriorROI.set_size(matrix.rows(), matrix.cols());
    this->m_MatrixPriorROI.set_size(matrix.rows(), matrix.cols());
    this->m_OriginalMatrixPriorROI.update(matrix);
    this->m_MatrixPriorROI.update(matrix);
  }

  void
  SetMatrixPriorROI2(MatrixType matrix)
  {
    this->m_MatrixPriorROI2.set_size(matrix.rows(), matrix.cols());
    this->m_MatrixPriorROI2.update(matrix);
  }

  // itkSetMacro( priorScale, RealType );
  // itkGetMacro( priorScale, RealType );
  RealType
  GetPriorScaleMat()
  {
    return this->m_priorScaleMat;
  }

  void
  SetPriorScaleMat(MatrixType priorScaleMat)
  {
    this->m_priorScaleMat.set_size(priorScaleMat.rows(), priorScaleMat.cols());
    this->m_priorScaleMat.update(priorScaleMat);
  }

  VectorType
  GetPWeights()
  {
    return this->m_WeightsP;
  }

  VectorType
  GetQWeights()
  {
    return this->m_WeightsQ;
  }

  VectorType
  GetRWeights()
  {
    return this->m_WeightsR;
  }

  RealType
  GetCorrelationForSignificanceTest()
  {
    return this->CorrelationForSignificanceTest;
  }

  VectorType
  GetCanonicalCorrelations()
  {
    return this->m_CanonicalCorrelations;
  }

  VectorType
  GetVariateP(unsigned int i = 0)
  {
    return this->m_VariatesP.get_column(i);
  }

  VectorType
  GetVariateQ(unsigned int i = 0)
  {
    return this->m_VariatesQ.get_column(i);
  }

  MatrixType
  GetVariatesP()
  {
    return this->m_VariatesP;
  }

  MatrixType
  GetVariatesQ()
  {
    return this->m_VariatesQ;
  }

  VectorType
  FastOuterProductVectorMultiplication(VectorType & p, VectorType & v)
  { // computes  outer_product( p , p ) * v

    // std::cout << p.size() <<v.size() <<std::endl;
    if (p.size() != v.size())
    {
      if (!this->m_Silent)
        std::cout << "FastOuterProductVectorMultiplication Usage Error " << std::endl;
      if (!this->m_Silent)
        std::cout << "Size 1: " << p.size() << "Size 2: " << v.size() << std::endl;
      return v;
    }
    RealType   ip = inner_product(p, v);
    VectorType vout(p);
    return vout * ip;
  }

  RealType
  ComputeEnergySlope(std::vector<RealType> vexlist, unsigned int n)
  {
    unsigned int N = vexlist.size();
    unsigned int loline = N - n;

    if (N < n * 2)
    {
      return 1;
    }
    double s0 = (n + 1);
    double s1 = 0;
    double s2 = 0;
    double t0 = 0;
    double t1 = 0;
    for (unsigned int i = loline; i < N; ++i)
    {
      const double t = (i - loline);
      s1 += t;
      s2 += (t * t);
      const double e = vexlist[i] - vexlist[loline];
      t0 += e;
      t1 += (e * e);
    }
    double M = 1;
    double denom = (s0 * s2 - s1 * s1);
    if (denom > 0)
    {
      M = (s1 * t0 - s0 * t1) / denom;
    }
    return M;
    /*
    std::vector<RealType> sublist;
    for (int i=listsize-4; i<listsize; i++) sublist.push_back( vexlist[i] );
    bool allequal=true;
    for (int i = 4; i < listsize; ++i) {
      allequal=true;
      for (int j=0; j<4; j++) if ( sublist[j] != vexlist[i-j] ) allequal=false;
      if (allequal) return 1.e-7;
    }
    return sts/stt*(1);
*/
  }

  RealType
  SparseCCA(unsigned int nvecs);

  RealType
  RidgeCCA(unsigned int nvecs);

  RealType
  SparsePartialCCA(unsigned int nvecs);

  bool
  CCAUpdate(unsigned int nvecs, bool, bool, unsigned int);

  bool
  CCAUpdateLong(unsigned int nvecs, bool, bool);

  RealType
  SparsePartialArnoldiCCA(unsigned int nvecs);

  RealType
  IHTCCA(unsigned int nvecs);

  RealType
  SparseReconB(MatrixType &, VectorType &);

  RealType
  SparseRecon(unsigned int nvecs);

  RealType
  SparseReconPrior(unsigned int nvecs, bool prior);

  RealType
  IHTPowerIterationHome(MatrixType & A, VectorType & x_k, unsigned int, unsigned int);

  RealType
  SparseReconHome(unsigned int nvecs);

  RealType
  SparseArnoldiSVDGreedy(unsigned int nvecs);

  RealType
  SparseArnoldiSVD_Other(MatrixType & A);

  RealType
  SparseArnoldiSVD(unsigned int nvecs);

  RealType
  SparseArnoldiSVD_x(unsigned int nvecs);

  RealType
  SparseArnoldiSVD_z(unsigned int nvecs);

  RealType
  ComputeSPCAEigenvalues(unsigned int, RealType, bool);
  RealType
  BasicSVD();

  RealType
  CGSPCA(unsigned int);

  RealType
  NetworkDecomposition(unsigned int nvecs);

  RealType
  LASSO_Cross();

  RealType
  LASSO(unsigned int nvecs);

  void
  LASSO_alg(MatrixType & X, VectorType & y, VectorType & beta, RealType gamma, unsigned int its);

  inline RealType
  LASSOSoft(RealType beta, RealType gamma)
  {
    if (beta > 0 && gamma < beta)
    {
      return beta - gamma;
    }
    else if (beta > 0 && gamma > beta)
    {
      return 0;
    }
    else if (beta < 0 && gamma < (beta * (-1)))
    {
      return beta + gamma;
    }
    return 0;
  }

  inline RealType
  ComputeIntercept(MatrixType A, VectorType x, VectorType b)
  {
    RealType intercept = b.mean();

    for (unsigned int Acol = 0; Acol < A.cols(); Acol++)
    {
      intercept -= A.get_column(Acol).mean() * x(Acol);
    }
    return intercept;
  }

  inline RealType
  SimpleRegression(VectorType y, VectorType ypred)
  {
    RealType corr = this->PearsonCorr(y, ypred);
    double   sdy = sqrt((y - y.mean()).squared_magnitude() / (y.size() - 1));
    double   sdyp = sqrt((ypred - ypred.mean()).squared_magnitude() / (y.size() - 1));

    if (itk::Math::FloatAlmostEqual(sdyp, 0.0))
    {
      return 0;
    }
    return corr * sdy / sdyp;
  }

  MatrixType
  GetCovMatEigenvectors(MatrixType p);

  void
  MRFFilterVariateMatrix();

  ImagePointer
  ConvertVariateToSpatialImage(VectorType variate, ImagePointer mask, bool threshold_at_zero = false);

protected:
  void
  SortResults(unsigned int n_vecs);

  // for pscca
  void
  UpdatePandQbyR();

  void
  PositivePart(VectorType & x_k1, bool takemin = false)
  {
    if (takemin)
    {
      RealType minval = x_k1.min_value();
      x_k1 = x_k1 - minval;
      return;
    }
    for (unsigned int i = 0; i < x_k1.size(); i++)
    {
      if (x_k1[i] < 0)
      {
        x_k1[i] = itk::Math::abs(x_k1[i]);
        x_k1[i] = 0;
      }
    }
  }

  void
  SparsifyOther(VectorType & x_k1, bool doclassic = false)
  {
    RealType fnp = itk::Math::abs(this->m_RowSparseness);
    if (fnp < 1.e-11)
    {
      return;
    }
    if (doclassic)
    {
      if (this->m_RowSparseness < 0)
      {
        if (fnp > x_k1.max_value())
          fnp = x_k1.max_value() * 0.9;
        for (unsigned int i = 0; i < x_k1.size(); i++)
        {
          RealType delta = itk::Math::abs(x_k1(i)) - fnp;
          if (delta < 0)
            delta = 0;
          else if (x_k1(i) < 0)
            delta *= (-1);
          x_k1(i) = delta;
        }
        return;
      }
      if (this->m_RowSparseness > 0)
      {
        if (fnp > x_k1.max_value())
          fnp = x_k1.max_value() * 0.9;
        for (unsigned int i = 0; i < x_k1.size(); i++)
        {
          RealType delta = itk::Math::abs(x_k1(i)) - fnp;
          if (delta < 0)
            delta = 0;
          x_k1(i) = delta;
        }
        return;
      }
    } // doclassic

    bool usel1 = this->m_UseL1;
    this->m_UseL1 = true;
    bool keeppos = false;
    if (this->m_RowSparseness > 1.e-11)
    {
      keeppos = true;
    }
    this->Sparsify(x_k1, fnp, keeppos, 0, nullptr);
    this->m_UseL1 = usel1;
  }

  void
  SparsifyP(VectorType & x_k1)
  {
    RealType fnp = itk::Math::abs(this->m_FractionNonZeroP);
    this->Sparsify(x_k1, fnp, this->m_KeepPositiveP, this->m_MinClusterSizeP, this->m_MaskImageP);
  }

  void
  SparsifyQ(VectorType & x_k1)
  {
    RealType fnp = itk::Math::abs(this->m_FractionNonZeroQ);
    this->Sparsify(x_k1, fnp, this->m_KeepPositiveQ, this->m_MinClusterSizeQ, this->m_MaskImageQ);
  }


  void
  Sparsify(VectorType & x_k1, RealType fnp, bool keeppos, unsigned int clust, ImagePointer mask)
  {

    if (x_k1.size() <= 1)
    {
      return;
    }
    if (fnp >= 1 && keeppos)
    {
      this->PositivePart(x_k1);
      return;
    }
    if (fnp >= 1)
    {
      return;
    }
    bool negate = false;
    if (itk::Math::abs(x_k1.min_value()) > x_k1.max_value())
    {
      negate = true;
    }
    if (negate)
    {
      x_k1 = x_k1 * (-1);
    }
    RealType initmax = x_k1.max_value();
    x_k1 = x_k1 / initmax;
    if (!this->m_MaxBasedThresholding)
    {
      // remove this optimization and just use fnp * high value as threshold
      std::vector<RealType> x_k1sort(x_k1.size(), 0);
      for (unsigned long j = 0; j < x_k1.size(); ++j)
        x_k1sort[j] = (x_k1(j));
      sort(x_k1sort.begin(), x_k1sort.end(), std::less<RealType>());
      //    std::cout << "Sorted " << x_k1sort[1] << " " << x_k1sort[x_k1.size()-1] << std::endl;
      RealType      low = 0;
      RealType      maxj = static_cast<RealType>(x_k1.size());
      unsigned long maxjind = static_cast<unsigned long>((1 - fnp * 0.5) * maxj + 0.5);
      unsigned long minjind = static_cast<unsigned long>((fnp * 0.5) * maxj + 0.5);
      if (maxjind > (maxj - 1))
        maxjind = (maxj - 1);
      RealType maxval = x_k1sort[maxjind];
      RealType minval = x_k1sort[minjind];
      RealType high = maxval;
      if (itk::Math::abs(minval) > high)
        high = itk::Math::abs(minval);
      //    if ( high * 0.2  > 0 ) low = high * 0.2; // hack to speed convergence
      RealType     eng = fnp;
      RealType     mid = low + 0.5 * (high - low);
      unsigned int its = 0;
      RealType     fnm = 0;
      // NOT USED: RealType     lastfnm = 1;
      while (((eng > (fnp * 0.1)) && (itk::Math::abs(high - low) > this->m_Epsilon) && (its < 20)) || its < 3)
      {
        mid = low + 0.5 * (high - low);
        VectorType searcherm(x_k1);
        //      if ( its > 10 & fnm > 0.99 ) std::cout << " A " << searcherm << std::endl;
        this->SoftClustThreshold(searcherm, mid, keeppos, clust, mask);
        //      if ( its > 10 & fnm > 0.99 ) std::cout << " B " << searcherm << std::endl;
        searcherm = this->SpatiallySmoothVector(searcherm, mask);
        //      if ( its > 10 & fnm > 0.99 ) std::cout << " C " << searcherm << std::endl;
        //      if ( its > 10 & fnm > 0.99 ) exit(1);
        // NOT USED: lastfnm = fnm;
        fnm = this->CountNonZero(searcherm);
        if (fnm > fnp)
        {
          low = mid; // 0.5 * ( low + mid  ); // relax this b/c it may not be a strictly quadratic space
        }
        if (fnm < fnp)
        {
          high = mid; // 0.5 * ( high + mid  );
        }
        eng = itk::Math::abs(fnp - fnm);
        its++;
      }
      this->SoftClustThreshold(x_k1, mid, keeppos, clust, mask);
    }
    else
    {
      for (unsigned long j = 0; j < x_k1.size(); ++j)
        if (itk::Math::abs(x_k1[j]) < (1.0 - fnp))
          x_k1[j] = 0;
      this->SoftClustThreshold(x_k1, 0, keeppos, clust, mask);
    }
    x_k1 = this->SpatiallySmoothVector(x_k1, mask);
    if (negate)
    {
      x_k1 = x_k1 * (-1);
    }
    return;
  }


  void
  SparsifyOld(VectorType & x_k1, RealType fnp, bool keeppos, unsigned int clust, ImagePointer mask)
  {

    if (x_k1.size() <= 1)
    {
      return;
    }
    if (fnp >= 1 && keeppos)
    {
      this->PositivePart(x_k1);
      return;
    }
    if (fnp >= 1)
    {
      return;
    }
    bool negate = false;
    if (x_k1.mean() <= 0)
    {
      negate = true;
    }
    if (negate)
    {
      x_k1 = x_k1 * (-1);
    }
    RealType initmax = x_k1.max_value();
    x_k1 = x_k1 / initmax;
    RealType     low = 0;
    RealType     high = 1;
    RealType     eng = fnp;
    RealType     mid = low + 0.5 * (high - low);
    unsigned int its = 0;
    RealType     fnm = 0;
    RealType     lastfnm = 1;
    while (((eng > 5.e-3) && (itk::Math::abs(high - low) > this->m_Epsilon) && (its < 50)) || its < 5)
    {
      mid = low + 0.5 * (high - low);
      VectorType searcherm(x_k1);
      this->SoftClustThreshold(searcherm, mid, keeppos, clust, mask);
      searcherm = this->SpatiallySmoothVector(searcherm, mask);
      lastfnm = fnm;
      fnm = this->CountNonZero(searcherm);
      //      if ( mask ) std::cout <<" its " << its << " spar " << fnm << " low " << low << " mid " << mid << " high "
      //      << high << " eng " << eng << std::endl;
      if (fnm > fnp)
      {
        low = mid;
      }
      if (fnm < fnp)
      {
        high = mid;
      }
      eng = itk::Math::abs(fnp - fnm);
      its++;
    }

    this->SoftClustThreshold(x_k1, mid, keeppos, clust, mask);
    x_k1 = this->SpatiallySmoothVector(x_k1, mask);
    if (negate)
    {
      x_k1 = x_k1 * (-1);
    }
    return;
  }

  void
  SparsifyP(VectorType & x_k1, VectorType & refvec)
  {
    if (x_k1.size() != refvec.size())
    {
      if (!this->m_Silent)
      {
        std::cout << " sizes dont match " << std::endl;
      }
      std::exception();
    }
    for (unsigned int i = 0; i < x_k1.size(); i++)
    {
      if (refvec(i) == 0)
      {
        x_k1(i) = 0;
      }
    }
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
  CountNonZero(VectorType v)
  {
    unsigned long ct = 0;

    for (unsigned int i = 0; i < v.size(); i++)
    {
      if (itk::Math::abs(v[i]) > this->m_Epsilon)
      {
        ct++;
      }
    }
    return (RealType)ct / (RealType)v.size();
  }

  bool
  Close2Zero(RealType x)
  {
    RealType eps = this->m_Epsilon * 5.0;
    //    eps = 0.0001;
    if (itk::Math::abs(x - itk::NumericTraits<RealType>::ZeroValue()) < eps)
      return true;
    return false;
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
    return itk::Math::abs(numer / denom);
  }

  RealType
  RPearsonCorr(VectorType v1, VectorType v2)
  {
    std::vector<TRealType> zeromatch(v1.size(), 0);

    for (unsigned int zm = 0; zm < v1.size(); zm++)
    {
      if ((this->Close2Zero(v1(zm)) || this->Close2Zero(v2(zm))))
      {

        zeromatch[zm] = 1;
        v1(zm) = 0;
        v2(zm) = 0;
      }
    }

    double frac = 1.0 / (double)v1.size();
    double xysum = 0;
    double xsum = 0;
    double ysum = 0;
    double xsqr = 0;
    double ysqr = 0;
    for (unsigned int i = 0; i < v1.size(); i++)
    {
      if (zeromatch[i] == 0)
      {
        xysum += v1(i) * v2(i);
        xsum += v1(i);
        xsqr += v1(i) * v1(i);
        ysum += v2(i);
        ysqr += v2(i) * v2(i);
      }
    }
    double numer = xysum - frac * xsum * ysum;
    double denom = sqrt((xsqr - frac * xsum * xsum) * (ysqr - frac * ysum * ysum));
    if (denom <= 0)
    {
      return 0;
    }
    return itk::Math::abs(numer / denom);
  }


  RealType
  GoldenSection(MatrixType & A,
                VectorType & x_k,
                VectorType & p_k,
                VectorType & bsol,
                RealType     a,
                RealType     b,
                RealType     c,
                RealType     tau,
                RealType     lambda);

  //  VectorType vEtoV( eVector v ) {
  //   VectorType v_out( v.data() , v.size() );
  //  return v_out;
  // }

  // eVector vVtoE( VectorType v ) {
  //  eVector v_out( v.size() );
  //  for (unsigned int i=0; i < v.size() ; i++) v_out(i)=v(i);
  //  return v_out;
  // }
  /*
  MatrixType mEtoV( eMatrix m , unsigned int ncols = 0) {
    MatrixType m_out( m.data() , m.rows() , m.cols() );
    if (  m(0,1) != m_out(0,1) ) {
      std::cout << " WARNING!! in eigen to vnl coversion for matrices " << std::endl;
      std::cout <<" eigen " << m(0,1) << " vnl " << m_out(0,1) << std::endl;
    }
    //    std::cout <<" eigen at (0,1) " << m(0,1) << " vnl at (0,1) " << m_out(0,1) <<  " vnl at (1,0) " << m_out(1,0)
  << std::endl; if ( ncols == 0 ) return m_out; else return (m_out).get_n_columns(0,ncols);
    // use this if you dont set #define EIGEN_DEFAULT_TO_ROW_MAJOR (we do this)
    if ( ncols == 0 )
      return m_out.transpose();
    else return (m_out.transpose()).get_n_columns(0,ncols);
    }*/
  /*
  eMatrix mVtoE( MatrixType m ) {
// NOTE: Eigen matrices are the transpose of vnl matrices unless you set #define EIGEN_DEFAULT_TO_ROW_MAJOR which we do
    eMatrix m_out(m.rows(),m.cols());
    for ( long i=0; i<m.rows(); ++i)
      for ( long j=0; j<m.cols(); ++j)
    m_out(i,j)=m(i,j);
    return m_out;
   }
  */
  void
  GetSubMatrix(MatrixType & A, MatrixType & submatout)
  {
    VectorType   diag = this->m_Indicator.diagonal();
    unsigned int nzct = (unsigned int)diag.sum();
    MatrixType   submat(A.rows(), nzct, 0);

    nzct = 0;
    for (unsigned int i = 0; i < diag.size(); i++)
    {
      if (diag(i) > 0)
      {
        submat.set_column(nzct, A.get_column(i));
        nzct++;
      }
    }
    submatout = submat;
  }

  antsSCCANObject();
  ~antsSCCANObject() override = default;

  void
  PrintSelf(std::ostream &, /* os */ Indent /* indent */) const override
  {
    if (this->m_MaskImageP && this->m_MaskImageQ && this->m_MaskImageR)
    {
      if (!this->m_Silent)
        std::cout << " 3 matrices " << std::endl;
    }
    else if (this->m_MaskImageP && this->m_MaskImageQ)
    {
      if (!this->m_Silent)
        std::cout << " 2 matrices " << std::endl;
    }
    else
    {
      if (!this->m_Silent)
        std::cout << " fewer than 2 matrices " << std::endl;
    }
  }

  void
  RunDiagnostics(unsigned int);

  void
  AddColumnsToMatrix(MatrixType & mat_to_add_to, MatrixType & mat_to_take_from, unsigned int col0, unsigned int coln)
  {
    MatrixType outmat(mat_to_add_to.rows(), mat_to_add_to.cols() + (coln - col0) + 1, 0);

    for (unsigned int i = 0; i < mat_to_add_to.cols(); i++)
    {
      outmat.set_column(i, mat_to_add_to.get_column(i));
    }
    unsigned int ct = mat_to_add_to.cols();
    for (unsigned int i = col0; i <= coln; i++)
    {
      outmat.set_column(ct, mat_to_take_from.get_column(i));
      ct++;
    }
    mat_to_add_to = outmat;
  }

  RealType
  CurvatureSparseness(VectorType & x, RealType sparsenessgoal, unsigned int maxit, ImagePointer);

private:
  ImagePointer
  ConvertVariateToSpatialImage4D(VectorType variate, ImagePointer mask, bool threshold_at_zero = false);

  MatrixType m_OriginalMatrixPriorROI;
  VectorType
  ConvertImageToVariate(ImagePointer image, ImagePointer mask);
  VectorType
  ConvertImageToVariate4D(ImagePointer image, ImagePointer mask);
  VectorType
  ClusterThresholdVariate(VectorType &, ImagePointer mask, unsigned int);
  VectorType
  ClusterThresholdVariate4D(VectorType &, ImagePointer mask, unsigned int);

  bool       m_Debug;
  bool       m_Silent;
  bool       m_MaxBasedThresholding;
  MatrixType m_OriginalMatrixP;
  MatrixType m_OriginalMatrixQ;
  MatrixType m_OriginalMatrixR;
  RealType   m_RowSparseness;

  antsSCCANObject(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented

  unsigned int m_ElapsedIterations;
  unsigned int m_MaximumNumberOfIterations;
  RealType     m_CurrentConvergenceMeasurement;
  RealType     m_ConvergenceThreshold;

  SCCANFormulationType m_SCCANFormulation;
  RealType             m_PinvTolerance;
  RealType             m_PercentVarianceForPseudoInverse;
  RealType             m_Epsilon; /** used to prevent div by zero */

  // Prior constrained PCA --Refer to notation in the paper
  MatrixType m_MatrixPriorROI;
  MatrixType m_MatrixPriorROI2;
  MatrixType m_SortedIndicesAll;
  VectorType sortedIndicesLoop;
  MatrixType m_Ip;
  MatrixType m_Ik;
  MatrixType m_priorScaleMat;

  // VectorType m_WeightsP;
  VectorType loc_Array;
  bool       flagForSort;

  VectorType   m_WeightsP;
  MatrixType   m_MatrixP;
  ImagePointer m_MaskImageP;
  RealType     m_FractionNonZeroP;
  bool         m_KeepPositiveP;
  RealType     m_UseLongitudinalFormulation;
  RealType     m_Smoother;

  VectorType   m_WeightsQ;
  MatrixType   m_MatrixQ;
  ImagePointer m_MaskImageQ;
  RealType     m_FractionNonZeroQ;
  bool         m_KeepPositiveQ;

  MatrixType  m_Eigenvectors;
  VectorType  m_Eigenvalues;
  VectorType  m_CanonicalCorrelations;
  VariateType m_SparseVariatesP;
  VariateType m_VariatesP;
  VariateType m_VariatesQ;
  /** solution to   X - U V */
  MatrixType m_MatrixU;

  VectorType   m_WeightsR;
  MatrixType   m_MatrixR;
  ImagePointer m_MaskImageR;
  RealType     m_FractionNonZeroR;
  bool         m_KeepPositiveR;
  /** a special variable for pscca, holds R^T R */
  MatrixType m_MatrixRRt;
  MatrixType m_MatrixRp;
  MatrixType m_MatrixRq;

  /** softer = true will compute the update  : if ( beta > thresh )  beta <- beta - thresh
   *     rather than the default update      : if ( beta > thresh )  beta <- beta  */
  unsigned int m_Covering;
  unsigned int m_VecToMaskSize;
  bool         m_GetSmall;
  bool         m_UseL1;
  bool         m_AlreadyWhitened;
  bool         m_SpecializationForHBM2011;
  RealType     m_CorrelationForSignificanceTest;
  RealType     m_lambda;

  // this->ComputeIntercept( A, x, b );
  RealType                   m_Intercept;
  unsigned int               m_NTimeDimensions;
  unsigned int               m_MinClusterSizeP;
  unsigned int               m_MinClusterSizeQ;
  unsigned int               m_KeptClusterSize;
  unsigned int               m_GoldenSectionCounter;
  VectorType                 m_ClusterSizes;
  VectorType                 m_OriginalB;
  VectorType                 m_SparsenessP;
  VectorType                 m_SparsenessQ;
  vnl_diag_matrix<RealType>  m_Indicator;
  vnl_diag_matrix<TRealType> m_PreC; // preconditioning
  RealType                   m_GSBestSol;
  RealType                   m_GradStep;
  RealType                   m_GradStepP;
  RealType                   m_GradStepQ;
  RealType                   m_PriorWeight;
};
} // namespace ants
} // namespace itk

#  ifndef ITK_MANUAL_INSTANTIATION
#    include "antsSCCANObject.hxx"
#  endif

#endif

/*

  RealType SparseRayleighQuotientIteration( MatrixType& A, VectorType& x)
    {
    if ( x.two_norm() == 0 ) return;
    x = x / x.two_norm();
    for ( unsigned int i = 0 ; i < 5; i++)
      {
      VectorType Ax = A * x;
      RealType lambda = inner_product( Ax , Ax );
      vnl_diag_matrix<double> diag( A , lambda );
      this->RidgeRegression( A, x, 1.e2, 10, diag );
      }
    }

*/
