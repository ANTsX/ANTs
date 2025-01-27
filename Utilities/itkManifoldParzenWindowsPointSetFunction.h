/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkManifoldParzenWindowsPointSetFunction_h
#define __itkManifoldParzenWindowsPointSetFunction_h

#include "itkPointSetFunction.h"

#include "itkGaussianMembershipFunction.h"
#include "itkKdTreeGenerator.h"
#include "itkListSample.h"
#include "itkMatrix.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkMeshSource.h"
#include "itkPointSet.h"
#include "itkVector.h"
#include "itkWeightedCentroidKdTreeGenerator.h"

#include <vector>

namespace itk
{
/** \class ManifoldParzenWindowsPointSetFunction.h
 * \brief point set filter.
 */

template <typename TPointSet, typename TOutput = double, typename TCoordRep = double>
class ManifoldParzenWindowsPointSetFunction : public PointSetFunction<TPointSet, TOutput, TCoordRep>
{
public:
  typedef ManifoldParzenWindowsPointSetFunction           Self;
  typedef PointSetFunction<TPointSet, TOutput, TCoordRep> Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Extract dimension from output image. */
  static constexpr unsigned int Dimension = TPointSet::PointDimension;

  typedef typename Superclass::InputPointSetType InputPointSetType;
  typedef typename Superclass::InputPointType    InputPointType;

  /** Point set typedef support. */
  typedef TPointSet                                            PointSetType;
  typedef typename PointSetType::PointType                     PointType;
  typedef typename PointSetType ::PointsContainerConstIterator PointsContainerConstIterator;

  typedef Vector<typename PointSetType::CoordRepType, Self::Dimension> MeasurementVectorType;
  typedef typename Statistics::ListSample<MeasurementVectorType>                         SampleType;
  typedef typename Statistics ::WeightedCentroidKdTreeGenerator<SampleType>              TreeGeneratorType;
  typedef typename TreeGeneratorType::KdTreeType                                         KdTreeType;
  typedef typename KdTreeType ::InstanceIdentifierVectorType                             NeighborhoodIdentifierType;

  /** Other typedef */
  typedef TOutput RealType;
  typedef TOutput OutputType;

  typedef typename Statistics ::MersenneTwisterRandomVariateGenerator            RandomizerType;
  typedef typename Statistics::GaussianMembershipFunction<MeasurementVectorType> GaussianType;
  typedef std::vector<typename GaussianType::Pointer>                            GaussianContainerType;
  typedef typename GaussianType::MatrixType                                      CovarianceMatrixType;

  /** Helper functions */

  itkSetMacro(CovarianceKNeighborhood, unsigned int);
  itkGetConstMacro(CovarianceKNeighborhood, unsigned int);

  itkSetMacro(EvaluationKNeighborhood, unsigned int);
  itkGetConstMacro(EvaluationKNeighborhood, unsigned int);

  itkSetMacro(RegularizationSigma, RealType);
  itkGetConstMacro(RegularizationSigma, RealType);

  itkSetMacro(KernelSigma, RealType);
  itkGetConstMacro(KernelSigma, RealType);

  itkSetMacro(BucketSize, unsigned int);
  itkGetConstMacro(BucketSize, unsigned int);

  itkSetMacro(Normalize, bool);
  itkGetConstMacro(Normalize, bool);
  itkBooleanMacro(Normalize);

  itkSetMacro(UseAnisotropicCovariances, bool);
  itkGetConstMacro(UseAnisotropicCovariances, bool);
  itkBooleanMacro(UseAnisotropicCovariances);

  virtual void
  SetInputPointSet(const InputPointSetType * ptr);

  virtual TOutput
  Evaluate(const InputPointType & point) const;

  PointType
  GenerateRandomSample();

  typename GaussianType::Pointer
  GetGaussian(unsigned int i)
  {
    if (i < this->m_Gaussians.size())
    {
      return this->m_Gaussians[i].GetPointer();
    }
    else
    {
      return nullptr;
    }
    this->Modified();
  }

  void
  SetGaussian(unsigned int i, typename GaussianType::Pointer gaussian)
  {
    if (i >= this->m_Gaussians.size())
    {
      this->m_Gaussians.resize(i + 1);
    }
    this->m_Gaussians[i] = gaussian;
    this->Modified();
  }

  void
  GenerateKdTree();

  NeighborhoodIdentifierType
  GetNeighborhoodIdentifiers(MeasurementVectorType, unsigned int);
  NeighborhoodIdentifierType
  GetNeighborhoodIdentifiers(InputPointType, unsigned int);

protected:
  ManifoldParzenWindowsPointSetFunction();
  virtual ~ManifoldParzenWindowsPointSetFunction();
  void
  PrintSelf(std::ostream & os, Indent indent) const;

  void
  GenerateData();

private:
  // purposely not implemented
  ManifoldParzenWindowsPointSetFunction(const Self &);
  void
  operator=(const Self &);

  unsigned int m_CovarianceKNeighborhood;
  unsigned int m_EvaluationKNeighborhood;
  unsigned int m_BucketSize;
  RealType     m_RegularizationSigma;
  RealType     m_KernelSigma;

  typename TreeGeneratorType::Pointer m_KdTreeGenerator;
  typename SampleType::Pointer        m_SamplePoints;

  typename RandomizerType::Pointer m_Randomizer;
  GaussianContainerType            m_Gaussians;
  bool                             m_Normalize;
  bool                             m_UseAnisotropicCovariances;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkManifoldParzenWindowsPointSetFunction.hxx"
#endif

#endif
