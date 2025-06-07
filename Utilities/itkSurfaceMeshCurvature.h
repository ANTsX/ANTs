/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _SurfaceMeshCurvature_h
#define _SurfaceMeshCurvature_h

#include "itkSurfaceCurvatureBase.h"

namespace itk
{
/** \class SurfaceMeshCurvature
 *
 * This class takes a surface as input and creates a local
 * geometric frame for each surface point.
 *
 *
 */
template <typename TSurface, typename TSurfacePatch>
class SurfaceMeshCurvature : public SurfaceCurvatureBase<TSurface, 3>
{
public:
  /** Standard class typedefs. */
  typedef SurfaceMeshCurvature           Self;
  typedef SurfaceCurvatureBase<TSurface> Superclass;
  typedef SmartPointer<Self>             Pointer;
  typedef SmartPointer<const Self>       ConstPointer;

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(SurfaceMeshCurvature);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  typedef typename Superclass::RealType   RealType;
  typedef typename Superclass::PointType  VectorType;
  typedef typename Superclass::PointType  FixedVectorType;
  typedef typename Superclass::PointType  PointType;
  typedef typename Superclass::MatrixType MatrixType;

  typedef TSurfacePatch                            SurfacePatch;
  typedef typename TSurfacePatch::Pointer          SurfacePatchPointer;
  typedef typename TSurfacePatch::NodeLocationType VertexType;
  typedef typename TSurfacePatch::PixelType        FunctionValueType;
  typedef TSurface *                               SurfacePointer; // probably vtkmesh

  /** Find all points within some distance of the origin.
   * The argument gives the number of times to apply the
   * mean shift algorithm to find the best neighborhood.
   */
  virtual void
  FindNeighborhood(unsigned int numMeanShifts = 0)
  {
    this->m_PointList.clear();
    this->m_FunctionValueList.clear();
    VertexType origin = m_SurfacePatch->GetLocation();
    PointType  pt;
    for (int j = 0; j < 3; j++)
    {
      pt(j) = origin[j];
    }
    this->m_Origin = pt;
    this->m_PointList.insert(this->m_PointList.begin(), pt);
    this->m_FunctionValueList.insert(this->m_FunctionValueList.begin(), m_SurfacePatch->GetValue());
    for (unsigned int i = 0; i < m_SurfacePatch->m_Neighbors.size(); i++)
    {
      VertexType neigh = m_SurfacePatch->m_Neighbors[i]->GetLocation();
      PointType  pti;
      for (int j = 0; j < 3; j++)
      {
        pti(j) = neigh[j];
      }
      this->m_PointList.insert(this->m_PointList.begin(), pti);
      this->m_FunctionValueList.insert(this->m_FunctionValueList.begin(),
                                       this->m_SurfacePatch->m_Neighbors[i]->GetValue(0));
    }
  }

  void
  SetSurface(SurfacePointer s)
  {
    m_Surface = s;
  }

  void
  SetSurfacePatch(SurfacePatchPointer s)
  {
    m_SurfacePatch = s;
  }

  SurfaceMeshCurvature()
  {
    this->m_Origin.fill(0.0);
    this->m_ArbitraryTangent.fill(0.);
    this->m_Normal.fill(0.);
    this->m_Tangent1.fill(0.);
    this->m_Tangent2.fill(0.);
    this->m_DirectionalKappa = 0.0;
    this->m_Kappa1 = 0.0;
    this->m_Kappa2 = 0.0;
    this->m_GaussianKappa = 0.0;
    this->m_A = 0.0;
    this->m_B = 0.0;
    this->m_C = 0.0;
    this->m_W1 = 3. / 8.;
    this->m_W2 = 1. / 8.;
    this->m_Eval0 = 0.0;
    this->m_Eval1 = 0.0;
    this->m_Eval2 = 0.0;
    this->m_CurrentNeighborhoodPointIndex = 0;
    this->m_ParameterFileName = "";
    this->m_Pi = 3.14159265358979323846;
    this->m_Debug = true;
    this->m_Debug = false;
    this->m_UseGeodesicNeighborhood = false;
    this->m_TotalArea = 0.0;
  }

  ~SurfaceMeshCurvature(){};

protected:
private:
  SurfacePointer      m_Surface;
  SurfacePatchPointer m_SurfacePatch;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
// #include "itkSurfaceMeshCurvature.hxx"
#endif

#endif
