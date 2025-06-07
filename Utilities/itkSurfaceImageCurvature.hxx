/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _SurfaceImageCurvature_hxx
#define _SurfaceImageCurvature_hxx
#include "antsAllocImage.h"
#include <vnl/algo/vnl_real_eigensystem.h>

// #include "itkLevelSetCurvatureFunction.h"
#include "itkCastImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include <queue>
#include <map>
#include <algorithm>

namespace itk
{
/** Quick node class to help get the geodesic neighborhood
 */

template <typename TSurface>
class GeodesicNode
{
public:
  /** Image related types. */
  typedef TSurface                      ImageType;
  typedef typename ImageType::IndexType IndexType;

  unsigned long neighborhoodindex;
  float         distance;
  bool          connected;
  IndexType     imageindex;

  GeodesicNode()
  {
    distance = 0.0;
    connected = false;
    neighborhoodindex = 0;
  }

  GeodesicNode(unsigned long i, float d, bool t, IndexType ind)
  {
    distance = d;
    connected = t;
    neighborhoodindex = i;
    imageindex = ind;
  }

  ~GeodesicNode() = default;
};

template <typename pclass>
class GeodesicNodePriority /* defines the comparison operator for the prioritiy queue */
{
public:
  bool
  operator()(pclass N1, pclass N2)
  {
    return N1.distance > N2.distance;
  }
};

template <typename TSurface>
SurfaceImageCurvature<TSurface>::SurfaceImageCurvature()
{
  // outputs are curvature image and normal image
  this->ProcessObject::SetNumberOfRequiredOutputs(2);
  m_SurfaceLabel = 1;

  m_GradientImage = nullptr;

  m_UseLabel = false;
  m_kSign = -1.0;
  m_FunctionImage = nullptr;
  this->m_Vinterp = nullptr;
  this->m_MinSpacing = itk::NumericTraits<RealType>::max();
}

template <typename TSurface>
void
SurfaceImageCurvature<TSurface>::ProcessLabelImage()
{
  ImageType * image = GetInput();

  if (!image)
  {
    return;
  }

  IndexType index;

  typename ImageType::RegionType requestedRegion;
  m_ImageSize = image->GetLargestPossibleRegion().GetSize();
  ImageIteratorType ti(image, image->GetLargestPossibleRegion());

  ti.GoToBegin();
  while (!ti.IsAtEnd())
  {
    PixelType pix = ti.Get();
    index = ti.GetIndex();
    if (ti.Get() == 2)
    {
      ti.Set(m_SurfaceLabel);
    }
    else
    {
      ti.Set(0);
    }

    ti.Set(0);

    //    if (index[0] == 20 && index[1] < 100 && index[2] < 100) ti.Set(m_SurfaceLabel);

    float     rad = 23, d = 0;
    IndexType cind = { { 120, 120, 80 } };
    for (unsigned int j = 0; j < ImageDimension; j++)
    {
      d += (cind[j] - index[j]) * (cind[j] - index[j]);
    }

    d = sqrt(d);

    if (std::fabs(d - rad) <= 0.5f)
    {
      ti.Set(m_SurfaceLabel);
    }

    rad = 12, d = 0;
    IndexType cind2 = { { 20, 170, 45 } };
    for (unsigned int j = 0; j < ImageDimension; j++)
    {
      d += (cind2[j] - index[j]) * (cind2[j] - index[j]);
    }

    d = sqrt(d);

    //    if (fabs(d-rad) <= 0.5) ti.Set(m_SurfaceLabel);

    ++ti;
  }
}

template <typename TSurface>
void
SurfaceImageCurvature<TSurface>::FindEuclideanNeighborhood(
  typename SurfaceImageCurvature<TSurface>::PointType rootpoint)
{
  this->m_AveragePoint = this->m_Origin;
  this->m_PointList.insert(this->m_PointList.begin(), this->m_Origin);
  typename ImageType::SizeType  rad;
  IndexType                     oindex, index;
  typename ImageType::PointType tempp;
  tempp[0] = rootpoint[0];
  tempp[1] = rootpoint[1];
  tempp[2] = rootpoint[2];
  oindex = this->m_FunctionImage->TransformPhysicalPointToIndex(tempp);
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    rad[i] = (long)(m_NeighborhoodRadius/this->m_Spacing[i]+0.5);
  }
  m_ti.SetLocation(oindex);
  unsigned int temp = 0;
  for (temp = 0; temp < m_ti.Size(); temp++)
  {
    index = m_ti.GetIndex(temp);
    if (this->IsValidIndex(index))
    {
      float pix = m_ti.GetPixel(temp);
      if (this->IsValidSurface(pix, index))
      {
        PointType                     p;
        typename ImageType::PointType ipt;
        this->m_FunctionImage->TransformIndexToPhysicalPoint(index, ipt);
        float dist = 0.0;
        bool  isorigin = true;
        for (unsigned int k = 0; k < ImageDimension; k++)
        {
          if (index[k] != oindex[k])
          {
            isorigin = false;
          }
          p[k] = ipt[k];
          RealType delt = rootpoint[k] - ipt[k];
//          RealType delt = static_cast<RealType>(oindex[k]) - static_cast<RealType>(index[k]);
          dist += delt * delt;
        }
        dist = sqrt(dist);
        if (!isorigin && dist <= (m_NeighborhoodRadius))
        {
          this->m_AveragePoint = this->m_AveragePoint + p;
          this->m_PointList.insert(this->m_PointList.begin(), p);
        }
      }
    } // test valid index and valid surface
  }

  unsigned int npts = this->m_PointList.size();
  if (npts > 0)
  {
    this->m_AveragePoint /= (RealType)npts;
  }
  else
  {
    this->m_AveragePoint = this->m_Origin;
  }
  if (this->m_Debug)
  {
    std::cout << " point list size " << this->m_PointList.size() << std::endl;
    for (unsigned int i = 0; i < this->m_PointList.size(); i++)
    {
      std::cout << " point  " << this->m_PointList[i];
    }
    std::cout << std::endl;
  }
}

template <typename TSurface>
void
SurfaceImageCurvature<TSurface>::FindGeodesicNeighborhood()
{
  typedef std::priority_queue<GeodesicNode<ImageType>,
                              std::vector<GeodesicNode<ImageType>>,
                              GeodesicNodePriority<GeodesicNode<ImageType>>>
    QType;

  QType                                           nodeq;
  std::map<unsigned int, GeodesicNode<ImageType>> nodes;

  this->m_AveragePoint = this->m_Origin;

  this->m_PointList.insert(this->m_PointList.begin(), this->m_Origin);

  typename ImageType::SizeType rad;
  IndexType                    oindex, index;

  float         dist = 0.0;
  unsigned int  k = 0;
  unsigned long longindex = 0;
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    rad[i] = static_cast<long>(m_NeighborhoodRadius);
    oindex[i] = static_cast<long>(static_cast<double>(this->m_Origin[i]) + 0.5);
  }
  index = oindex;
  for (k = 0; k < ImageDimension; k++)
  {
    if (k == 0)
    {
      longindex = index[0];
    }
    if (k == 1)
    {
      longindex = index[1] + longindex + index[0] * m_ImageSize[0];
    }
    if (k == 2)
    {
      longindex = index[2] + longindex + index[2] * m_ImageSize[0] * m_ImageSize[1];
    }
  }
  GeodesicNode<ImageType> gnode(longindex, 0.0, true, oindex);
  nodes[longindex] = gnode;
  nodeq.push(gnode);

  float lastdist = 0.0;

  // if ( this->m_Origin[1]==146 && this->m_Origin[0] > 167 )
  // if ( this->m_Origin[1]==146 && this->m_Origin[0] == 168 && this->m_Origin[2]==215)
  // {
  // std::cout << " origin " << this->m_Origin << std::endl;
  // }

  while (!nodeq.empty() && lastdist <= m_NeighborhoodRadius)
  {
    GeodesicNode<ImageType> g = nodeq.top();

    lastdist = g.distance;

    PointType q;

    //    if ( g.distance < 2.0)
    //    {
    //      this->m_PointList.insert(this->m_PointList.begin(),q);
    //      nodes[g.neighborhoodindex].connected=true;
    //    }
    //    else
    if (lastdist <= m_NeighborhoodRadius)
    {
      m_ti2.SetLocation(g.imageindex);
      for (unsigned int jj = 0; jj < m_ti2.Size(); jj++)
      {
        index = m_ti2.GetIndex(jj);

        if ( // m_ti2.GetPixel(jj) == m_SurfaceLabel &&
          this->IsValidSurface(m_ti2.GetPixel(jj), index) && index[0] < m_ImageSize[0] - m_NeighborhoodRadius &&
          index[0] > m_NeighborhoodRadius && index[1] < m_ImageSize[1] - m_NeighborhoodRadius &&
          index[1] > m_NeighborhoodRadius && index[2] < m_ImageSize[2] - m_NeighborhoodRadius &&
          index[2] > m_NeighborhoodRadius)
        {
          longindex = 0;
          dist = 0;
          for (k = 0; k < ImageDimension; k++)
          {
            if (k == 0)
            {
              longindex = index[0];
            }
            if (k == 1)
            {
              longindex = index[1] + longindex + index[0] * m_ImageSize[0];
            }
            if (k == 2)
            {
              longindex = index[2] + longindex + index[2] * m_ImageSize[0] * m_ImageSize[1];
            }
            q[k] = (RealType)index[k];
            //            dist+=(g.imageindex[k]-oindex[k])*(g.imageindex[k]-oindex[k]);
            //            dist+=(index[k]-oindex[k])*(index[k]-oindex[k]);
            dist += (float)(g.imageindex[k] - index[k]) * (g.imageindex[k] - index[k]);
          }
          dist = sqrt(dist);

          // if ( this->m_Origin[1]==146 && this->m_Origin[0] == 168 && this->m_Origin[2]==215)
          // {
          // std::cout << " testing point " << index << " longind " << longindex << " dist " << dist <<
          // " bool " << nodes[longindex].connected << std::endl;
          // }
          //          if (!nodes[longindex].connected ) //&& !nodes[g.neighborhoodindex].connected)
          if (!nodes[longindex].connected && (dist + lastdist) <= m_NeighborhoodRadius)
          {
            GeodesicNode<ImageType> _gnode(longindex, dist + lastdist, true, index);
            //            GeodesicNode<ImageType> _gnode(longindex,dist,true,index);
            nodes[longindex] = _gnode;
            nodeq.push(_gnode);
            // if ( this->m_Origin[1]==146 && this->m_Origin[0] == 168 && this->m_Origin[2]==215)
            // /{
            // std::cout << " inserting point " << index << std::endl;
            // }
            this->m_PointList.insert(this->m_PointList.begin(), q);
            this->m_AveragePoint = this->m_AveragePoint + q;
          }
          //          else if ( dist > 0 && dist+lastdist < nodes[longindex].distance )
          //          {
          //            nodes[longindex].distance=dist+lastdist;
          //          }
        }
      }
    }
    nodeq.pop();
  }

  this->m_AveragePoint = this->m_AveragePoint / ((float)this->m_PointList.size());
}

template <typename TSurface>
void
SurfaceImageCurvature<TSurface>::FindNeighborhood(unsigned int numMeanShifts)
{
  if (this->m_UseGeodesicNeighborhood)
  {
    this->FindGeodesicNeighborhood();
  }
  else
  {
    this->FindEuclideanNeighborhood(this->GetOrigin());
    for (unsigned int dd = 0; dd < numMeanShifts; dd++)
    {
      this->m_PointList.clear();
      this->FindEuclideanNeighborhood(this->GetAveragePoint());
    }
  }
}

template <typename TSurface>
void
SurfaceImageCurvature<TSurface>::LevelSetMeanCurvature()
{
  /*
    ImageType* image=GetInput();

    if (!image) return;

    IndexType index;

    typename ImageType::RegionType requestedRegion;
    this->m_ImageSize=image->GetLargestPossibleRegion().GetSize();
    ImageIteratorType ti( image, image->GetLargestPossibleRegion() );

    // Define a level set curvature calculator
    typedef LevelSetCurvatureFunction<ImageType> CurvatureType;
    typename CurvatureType::Pointer inCurvature = CurvatureType::New();
    inCurvature->SetInputImage( image );

    ti.GoToBegin();
    while(!ti.IsAtEnd()  )
    {
      index=ti.GetIndex();
        //if (ti.Get() == this->m_SurfaceLabel)
      if(this->IsValidSurface(ti.Get(),index))
      {
        double curvature = inCurvature->EvaluateAtIndex( index );
        this->m_FunctionImage->SetPixel(index,fabs(curvature));
        }
      ++ti;
    }
  */
}

template <typename TSurface>
void
SurfaceImageCurvature<TSurface>::EstimateNormalsFromGradient()
{
  typename ImageType::Pointer image = GetInput();
  this->m_Spacing = image->GetSpacing();

  if (!image)
  {
    return;
  }

  typename ImageType::SizeType rad;
  typename ImageType::SizeType rad2;
  for (unsigned int t = 0; t < ImageDimension; t++)
  {
    rad[t] = (unsigned long)(this->m_NeighborhoodRadius);
    rad2[t] = 1;
  }
  this->m_ti.Initialize(rad, image, image->GetLargestPossibleRegion());
  this->m_ti2.Initialize(rad2, image, image->GetLargestPossibleRegion());

  typedef itk::ImageRegionIteratorWithIndex<TSurface> IteratorType;
  IteratorType                                        Iterator(image, image->GetLargestPossibleRegion().GetSize());
  bool                                                wmgmcurv = true;
  Iterator.GoToBegin();
  while (!Iterator.IsAtEnd())
  {
    float pix = Iterator.Get();
    if (!itk::Math::FloatAlmostEqual(pix, 0.0f) && !itk::Math::FloatAlmostEqual(pix, 1.0f) &&
        !itk::Math::FloatAlmostEqual(pix, 2.0f))
    {
      wmgmcurv = false;
    }
    ++Iterator;
  }

  if (wmgmcurv)
  {
    typename OutputImageType::Pointer laplacian = AllocImage<OutputImageType>(image->GetLargestPossibleRegion());
    laplacian->SetSpacing(image->GetSpacing());
    laplacian->SetDirection(image->GetDirection());
    laplacian->SetOrigin(image->GetOrigin());
    Iterator.GoToBegin();
    while (!Iterator.IsAtEnd())
    {
      IndexType ind = Iterator.GetIndex();
      if (itk::Math::FloatAlmostEqual(image->GetPixel(ind), static_cast<PixelType>(2)))
      {
        laplacian->SetPixel(ind, 1);
      }
      else if (itk::Math::FloatAlmostEqual(image->GetPixel(ind), itk::NumericTraits<PixelType>::OneValue()))
      {
        laplacian->SetPixel(ind, 0.);
      }
      else
      {
        laplacian->SetPixel(ind, 0.);
      }
      ++Iterator;
    }

    // smooth and then reset the values
    unsigned int totit = 50;
    for (unsigned int iterations = 0; iterations < totit; iterations++)
    {
      while (!Iterator.IsAtEnd())
      {
        IndexType ind = Iterator.GetIndex();
        if (itk::Math::FloatAlmostEqual(image->GetPixel(ind), static_cast<PixelType>(2)))
        {
          laplacian->SetPixel(ind, 1);
        }
        else if (itk::Math::FloatAlmostEqual(image->GetPixel(ind), itk::NumericTraits<PixelType>::ZeroValue()))
        {
          laplacian->SetPixel(ind, 0.);
        }
        ++Iterator;
      }

      typedef itk::DiscreteGaussianImageFilter<TSurface, TSurface> dgf;
      typename dgf::Pointer                                        filter = dgf::New();
      filter->SetVariance(0.5);
      filter->SetUseImageSpacing(true);
      filter->SetMaximumError(.01f);
      filter->SetInput(laplacian);
      filter->Update();
      laplacian = filter->GetOutput();
      Iterator.GoToBegin();
    }
    //    ANTs::WriteImage<TSurface>(laplacian,"lap.hdr");
    GradientImageFilterPointer filter = GradientImageFilterType::New();
    filter->SetInput(laplacian);
    RealType sigma = this->m_Sigma;
    filter->SetSigma(sigma);
    // Execute the filter
    filter->Update();
    this->m_GradientImage = filter->GetOutput();
    GradientPixelType zero;
    zero.Fill(0);
    Iterator.GoToBegin();
    while (!Iterator.IsAtEnd())
    {
      IndexType ind = Iterator.GetIndex();
      if (!itk::Math::FloatAlmostEqual(image->GetPixel(ind), itk::NumericTraits<PixelType>::OneValue()))
      {
        this->m_GradientImage->SetPixel(ind, zero);
      }
      ++Iterator;
    }
  }
  else
  {
    GradientImageFilterPointer filter = GradientImageFilterType::New();
    filter->SetInput(image);
    RealType sigma = this->m_Sigma;
    filter->SetSigma(sigma);
    // Execute the filter
    filter->Update();
    this->m_GradientImage = filter->GetOutput();
  }

  this->m_Vinterp = VectorInterpolatorType::New();
  this->m_Vinterp->SetInputImage(this->m_GradientImage);
}

template <typename TSurface>
void
SurfaceImageCurvature<TSurface>::WeingartenMap()
{
  typename ImageType::Pointer image = this->GetInput();
  if (!image)
  {
    return;
  }

  MatrixType   D;
  unsigned int j = 0;
  unsigned int i = 0;
  unsigned int npts = this->m_PointList.size();
  unsigned int vars = 3;
  if (npts < (vars + 1))
  {
    this->m_MeanKappa = 0;
    this->m_GaussianKappa = 0;
    this->m_Kappa1 = 0;
    this->m_Kappa2 = 0;
    this->m_Area = 0;
    return;
  }
  D.set_size(npts, vars); // each row contains [u^2 , uv, v^2, u, v, 1] for point p
  D.fill(0.0);

  MatrixType W(2, 2);
  W.fill(0.0);

  vnl_vector<double> xdists(npts);
  xdists.fill(0.0);
  vnl_vector<double> ydists(npts);
  ydists.fill(0.0);
  vnl_vector<double> zdists(npts);
  zdists.fill(0.0);
  vnl_vector<double> f_uvs(npts);
  f_uvs.fill(0.0);

  // go through all the points
  //  compute weight
  //  compute dist of unit dif and tangents
  //  compute dif of normal with grad at point

  this->m_Area = 0;
  PointType Q = this->m_Origin;
  RealType  areaelt = static_cast<RealType>(1000.0) / static_cast<RealType>(npts);
  for (j = 0; j < npts; j++)
  {
    typename ImageType::PointType pt;
    pt[0] = this->m_PointList[j][0];
    pt[1] = this->m_PointList[j][1];
    pt[2] = this->m_PointList[j][2];
    GradientPixelType norm = this->m_Vinterp->Evaluate(pt);
    PointType         Grad;
    for (i = 0; i < SurfaceDimension; i++)
    {
      Grad[i] = norm[i];
    }
    PointType PN = Grad / (Grad.magnitude());

    // get the surface parameterization ...
    float u1 = 0.0;
    float u2 = 0.0;
    float f_uv = 0.0;
    //    PointType Dif = Q + this->m_Normal - this->m_PointList[j] - PN;
    PointType Dif = Q - this->m_PointList[j];
    // this is the predefined u_i parameter
    u1 = this->innerProduct(Dif, this->m_Tangent1);
    // this is the predefined v_i parameter
    u2 = this->innerProduct(Dif, this->m_Tangent2);
    // now the inner product of PN and the normal is f_uv ...
    f_uv = this->innerProduct(PN, this->m_Normal);
    // the point is therefore defined as:
    //    PointType surfacePoint = this->m_Tangent1 * u1 + this->m_Tangent2 * u2 + PN;
    xdists[j] = (PN[0]);
    ydists[j] = (PN[1]);
    zdists[j] = (PN[2]);
    f_uvs[j] = f_uv;

    if (vars == 6)
    {
      // each row contains [u^2 , uv, v^2, u, v, 1] for point p
      D(j, 5) = u2 * u2; // (0   , 2*u2)
      D(j, 4) = u1 * u1; // (2*u1, 0)
      D(j, 3) = u1 * u2; // (u2  , u1)
    }
    // each row contains [ u, v, 1] for point p
    D(j, 2) = u2; // (1   , 0)
    D(j, 1) = u1; // (0   , 1)
    D(j, 0) = 1.0;
    // NOT USED: RealType dfuv_u = 0;
    // NOT USED: RealType dfuv_v = 0;
    // NOT USED: if ( ( itk::Math::abs (u1) > 0 ) && ( itk::Math::abs (u2) < 1.e-6 ) )
    // NOT USED:  dfuv_u = itk::Math::abs ( f_uv - 1.0 ) / itk::Math::abs (u1) * 100.0;
    // NOT USED:if ( ( itk::Math::abs (u2) > 0 ) && ( itk::Math::abs (u1) < 1.e-6 ) )
    // NOT USED:  dfuv_v = itk::Math::abs ( f_uv - 1.0 ) / itk::Math::abs (u2) * 100.0;
    // this->m_Area += sqrt( 1.0 + dfuv_u*dfuv_u + dfuv_v*dfuv_v );
    this->m_Area += itk::Math::abs(f_uv - 1.0f);
  }
  this->m_Area *= areaelt;
  vnl_svd<double>    svd(D);
  vnl_vector<double> ax = svd.solve(xdists);    // /totwt);
  vnl_vector<double> ay = svd.solve(ydists);    // /totwt);
  vnl_vector<double> az = svd.solve(zdists);    // /totwt);
  vnl_vector<double> df_uvs = svd.solve(f_uvs); // /totwt);

  // now get the first partials of each of these terms w.r.t. u and v

  // dN/du = (dN/du \dot T_1) T_1+ (dNdu dot T_2) T_2

  PointType dNdu;
  dNdu[0] = ax[1];
  dNdu[1] = ay[1];
  dNdu[2] = az[1];
  PointType dNdv;
  dNdv[0] = ax[2];
  dNdv[1] = ay[2];
  dNdv[2] = az[2];

  //  df_uvs = df_uvs * 1.e5; // scale up
  //  this->m_Area = sqrt( 1.0 + df_uvs[1] * df_uvs[1] + df_uvs[2] * df_uvs[2] );

  float a = 0;
  float b = 0;
  float c = 0;
  float d = 0;
  for (i = 0; i < SurfaceDimension; i++)
  {
    a += dNdu[i] * this->m_Tangent1[i];
    b += dNdv[i] * this->m_Tangent1[i];

    c += dNdu[i] * this->m_Tangent2[i];
    d += dNdv[i] * this->m_Tangent2[i];
  }
  W(0, 0) = a;
  W(0, 1) = b;
  W(1, 0) = c;
  W(1, 1) = d;

  // Compute estimated frame using eigensystem of D'*D
  {
    vnl_real_eigensystem                  eig(W);
    vnl_diag_matrix<std::complex<double>> DD(eig.D.rows()); //
    this->m_Kappa1 = std::real(eig.D(1, 1));
    this->m_Kappa2 = std::real(eig.D(0, 0));
    this->m_MeanKappa = (this->m_Kappa1 + this->m_Kappa2) * 0.5f;
    this->m_GaussianKappa = (this->m_Kappa1 * this->m_Kappa2);
  }
}


template <typename TSurface>
void
SurfaceImageCurvature<TSurface>::WeingartenMapGradients()
{
  typename ImageType::Pointer image = this->GetInput();
  if (!image)
  {
    return;
  }

  // rebuild the point list from gradients defined at u and v
  this->m_PointList.clear();
  PointType p;
  RealType  paramdelt = this->m_MinSpacing * static_cast<RealType>(0.5);
  RealType  eps = 1.e-6;
  for (RealType zi = -paramdelt; zi <= paramdelt + eps; zi = zi + paramdelt)
    for (RealType ui = -paramdelt; ui <= paramdelt + eps; ui = ui + paramdelt)
    {
      for (RealType vi = -paramdelt; vi <= paramdelt + eps; vi = vi + paramdelt)
      {
        p = this->m_Origin + this->m_Tangent1 * ui + this->m_Tangent2 * vi + this->m_Normal * zi;
        this->m_PointList.insert(this->m_PointList.begin(), p);
      }
    }
  this->WeingartenMap();
  return;
}


template <typename TSurface>
void
SurfaceImageCurvature<TSurface>::ComputeSurfaceArea()
{
  ImageType *                     image = GetInput();
  typename ImageType::SpacingType ispacing = image->GetSpacing();
  if (!image)
  {
    return;
  }
  FixedVectorType spacing;
  spacing[0] = ispacing[0];
  spacing[1] = ispacing[1];
  spacing[2] = ispacing[2];
  // BUG FIXME
  if (!this->m_GradientImage)
  {
    this->EstimateNormalsFromGradient();
  }

  IndexType index;

  typename ImageType::RegionType requestedRegion;
  this->m_ImageSize = image->GetLargestPossibleRegion().GetSize();
  ImageIteratorType ti(image, image->GetLargestPossibleRegion());

  RealType area = 0.0;
  this->m_TotalArea = 0.0;
  unsigned long ct = 0;

  ti.GoToBegin();
  while (!ti.IsAtEnd())
  {
    index = ti.GetIndex();

    if ( // ti.Get() == this->m_SurfaceLabel &&
      (this->IsValidSurface(ti.Get(), index)) && index[0] < this->m_ImageSize[0] - this->m_NeighborhoodRadius &&
      index[0] > this->m_NeighborhoodRadius && index[1] < this->m_ImageSize[1] - this->m_NeighborhoodRadius &&
      index[1] > this->m_NeighborhoodRadius && index[2] < this->m_ImageSize[2] - this->m_NeighborhoodRadius &&
      index[2] > this->m_NeighborhoodRadius)
    {
      ct++;
      // BUG FIXME
      area = this->ComputeLocalArea(spacing);
      //      area = 1.0;
      this->m_FunctionImage->SetPixel(index, area);
      this->m_TotalArea += area;
      this->m_PointList.clear();
      if (ct % 1000 == 0)
      {
        // std::cout << " ind " << index << " area " << area << std::endl;
      }
    }
    ++ti;
  }

  // std::cout << " surface area " << this->m_TotalArea << std::endl;
  return;
}


template <typename TSurface>
void
SurfaceImageCurvature<TSurface>::EstimateFrameFromGradient(IndexType index)
{
  GradientPixelType g = this->m_GradientImage->GetPixel(index);
  RealType          mag = 0.0;
  for (int i = 0; i < ImageDimension; i++)
  {
    this->m_Normal(i) = (RealType)g[i];
    mag += g[i] * g[i];
  }
  mag = sqrt(mag);
  if (mag <= 1.e-9)
  {
    this->m_Normal.fill(0.);
  }
  else
  {
    this->m_Normal /= sqrt(mag);
  }
  this->SetFrameFromNormal(this->m_Normal);
}

template <typename TSurface>
void
SurfaceImageCurvature<TSurface>::EstimateFrameFromGradient(ImagePointType ipt)
{
  GradientPixelType g = this->m_Vinterp->Evaluate(ipt);
  RealType          mag = itk::NumericTraits<RealType>::ZeroValue();
  for (int i = 0; i < ImageDimension; i++)
  {
    this->m_Normal(i) = (RealType)g[i];
    mag += g[i] * g[i];
  }
  mag = sqrt(mag);
  if (mag <= static_cast<RealType>(1.e-9))
  {
    this->m_Normal.fill(0.);
  }
  else
  {
    this->m_Normal /= sqrt(mag);
  }
  this->SetFrameFromNormal(this->m_Normal);
}

template <typename TSurface>
typename SurfaceImageCurvature<TSurface>::RealType
SurfaceImageCurvature<TSurface>::IntegrateFunctionOverSurface(bool norm)
{
  typename OutputImageType::Pointer image = this->m_FunctionImage;
  typename OutputImageType::Pointer tempimage = OutputImageType::New();
  tempimage->SetLargestPossibleRegion(image->GetLargestPossibleRegion());
  tempimage->SetBufferedRegion(image->GetLargestPossibleRegion());
  tempimage->AllocateInitialized();

  typename ImageType::SizeType rad;
  typename ImageType::SizeType rad2;
  for (unsigned int t = 0; t < ImageDimension; t++)
  {
    rad[t] = (unsigned long)(this->m_NeighborhoodRadius);
    rad2[t] = 1;
  }

  this->m_ti.Initialize(rad, this->GetInput(), image->GetLargestPossibleRegion());
  this->m_ti2.Initialize(rad2, this->GetInput(), image->GetLargestPossibleRegion());
  // unsigned long ct = 0;

  IndexType index;

  ImageIteratorType ti(this->GetInput(), this->GetInput()->GetLargestPossibleRegion());

  //  std::cout << " begin integrate ";

  ti.GoToBegin();
  while (!ti.IsAtEnd())
  {
    index = ti.GetIndex();
    tempimage->SetPixel(index, 0);
    if ( // ti.Get() == this->m_SurfaceLabel &&
      (this->IsValidSurface(ti.Get(), index)) && index[0] < this->m_ImageSize[0] - this->m_NeighborhoodRadius &&
      index[0] > this->m_NeighborhoodRadius && index[1] < this->m_ImageSize[1] - this->m_NeighborhoodRadius &&
      index[1] > this->m_NeighborhoodRadius && index[2] < this->m_ImageSize[2] - this->m_NeighborhoodRadius &&
      index[2] > this->m_NeighborhoodRadius)
    {
      PointType p;
      // ct++;
      for (unsigned int k = 0; k < ImageDimension; k++)
      {
        p[k] = (RealType)index[k];
      }
      this->SetOrigin(p);
      //          std::cout << " find nhood ";
      this->FindNeighborhood();
      RealType area = this->IntegrateFunctionOverNeighborhood(norm);
      tempimage->SetPixel(index, area);
    }
    ++ti;
  }

  this->CopyImageToFunctionImage(tempimage, this->m_FunctionImage);

  return 0;
}

template <typename TSurface>
void
SurfaceImageCurvature<TSurface>::CopyImageToFunctionImage(OutputImagePointer i1, OutputImagePointer i2)
{
  if (!i1 || !i2)
  {
    return;
  }

  OutputImageIteratorType ti1(i1, i1->GetLargestPossibleRegion());
  OutputImageIteratorType ti2(i2, i2->GetLargestPossibleRegion());

  ti1.GoToBegin();
  ti2.GoToBegin();
  while (!ti1.IsAtEnd())
  {
    ti2.Set(ti1.Get());
    ++ti1;
    ++ti2;
  }
}

template <typename TSurface>
typename SurfaceImageCurvature<TSurface>::RealType
SurfaceImageCurvature<TSurface>::IntegrateFunctionOverNeighborhood(bool norm)
{
  unsigned int npts = this->m_PointList.size();
  double       curvature = 0.0, tw = 0;
  for (unsigned int pp = 0; pp < npts; pp++)
  {
    IndexType localindex;
    for (unsigned int k = 0; k < ImageDimension; k++)
    {
      localindex[k] = (long)this->m_PointList[pp][k];
    }
    PointType dd = this->m_Origin - this->m_PointList[pp];
    double    wi = dd.magnitude();
    if (!itk::Math::FloatAlmostEqual(wi, 0.0))
    {
      wi = 1. / wi;
    }
    tw += wi;
    RealType func = this->m_FunctionImage->GetPixel(localindex);
    if (norm)
    {
      curvature += wi * static_cast<double>(func);
    }
    else
    {
      curvature += static_cast<double>(func);
    }
    //    curvature*=this->ComputeLocalArea(spacing);
  }
  // if (norm ) curvature/=tw;
  // SD sometimes tw is zero making curvature = NaN
  if (norm && !itk::Math::FloatAlmostEqual(tw, 0.0))
  {
    curvature /= tw;
  }
  this->m_PointList.clear();

  return curvature;
}

template <typename TSurface>
void
SurfaceImageCurvature<TSurface>::PostProcessGeometry()
{
  typename ImageType::Pointer image = GetInput();

  if (!image)
  {
    return;
  }

  IndexType index;

  typename ImageType::RegionType requestedRegion;
  this->m_ImageSize = image->GetLargestPossibleRegion().GetSize();
  ImageIteratorType ti(image, image->GetLargestPossibleRegion());

  std::vector<double> kvec;
  ti.GoToBegin();
  while (!ti.IsAtEnd())
  {
    PixelType pix = ti.Get();
    index = ti.GetIndex();
    if ( // ti.Get() == this->m_SurfaceLabel &&
      (this->IsValidSurface(ti.Get(), index)) && index[0] < this->m_ImageSize[0] - this->m_NeighborhoodRadius &&
      index[0] > this->m_NeighborhoodRadius && index[1] < this->m_ImageSize[1] - this->m_NeighborhoodRadius &&
      index[1] > this->m_NeighborhoodRadius && index[2] < this->m_ImageSize[2] - this->m_NeighborhoodRadius &&
      index[2] > this->m_NeighborhoodRadius) //
    {
      PointType p;
      for (unsigned int k = 0; k < ImageDimension; k++)
      {
        p[k] = (RealType)index[k];
      }
      this->SetOrigin(p);
      this->FindNeighborhood();
      int    npts = this->m_PointList.size() - 1;
      double curvature = 0.0;
      for (int pp = 0; pp < npts; pp++)
      {
        IndexType localindex;
        for (unsigned int k = 0; k < ImageDimension; k++)
        {
          localindex[k] = (long)this->m_PointList[pp][k];
        }
        //          PointType dd=this->m_Origin-this->m_PointList[pp];
        //          double wi=dd.magnitude();
        //          if (wi!=0.0) wi=1./wi;
        //          tw+=wi;vector<int> vec;

        curvature = this->m_FunctionImage->GetPixel(localindex);
        kvec.push_back(curvature);
      }

      std::sort(kvec.begin(), kvec.end()); // Sort the vector

      this->m_PointList.clear();
      //        curvature/=tw;
      this->m_FunctionImage->SetPixel(index, kvec[kvec.size() / 2]);
      kvec.clear();
    }
    ++ti;
  }
}

template <typename TSurface>
void
SurfaceImageCurvature<TSurface>::ComputeFrameOverDomain(unsigned int which)
{
  ImageType * image = this->GetInput();
  if (!image)
  {
    return;
  }
  for (unsigned int d = 0; d < ImageDimension; d++)
    if (static_cast<RealType>(image->GetSpacing()[d]) < this->m_MinSpacing)
      this->m_MinSpacing = image->GetSpacing()[d];
  IndexType index;
  this->m_ImageSize = image->GetLargestPossibleRegion().GetSize();
  ImageIteratorType ti(image, image->GetLargestPossibleRegion());

  // std::exception();
  // Get Normals First!
  this->EstimateNormalsFromGradient();

  // unsigned int  ct = 1;
  // unsigned long ct2 = 0;
  RealType      kpix = 0;
  // double thresh = 0.0;


  ti.GoToBegin();
  while (!ti.IsAtEnd())
  {
    index = ti.GetIndex();
    typename ImageType::PointType pt;
    image->TransformIndexToPhysicalPoint(index, pt);
    kpix = 0.0;
    if ( // ti.Get() == this->m_SurfaceLabel &&
      this->IsValidSurface(ti.Get(), index) && index[0] < this->m_ImageSize[0] - 2 * this->m_NeighborhoodRadius &&
      index[0] > 2 * this->m_NeighborhoodRadius && index[1] < this->m_ImageSize[1] - 2 * this->m_NeighborhoodRadius &&
      index[1] > 2 * this->m_NeighborhoodRadius && index[2] < this->m_ImageSize[2] - 2 * this->m_NeighborhoodRadius &&
      index[2] > 2 * this->m_NeighborhoodRadius) //
    {
      PointType p;
      for (unsigned int k = 0; k < ImageDimension; k++)
      {
        p[k] = pt[k];
      }
      this->SetOrigin(p);
      this->EstimateFrameFromGradient(pt);
      this->FindNeighborhood();
      switch (which)
      {
        case (0):
        {
          this->ComputeJoshiFrame(this->m_Origin);
        }
        break;
        case (1):
        {
          this->JainMeanAndGaussianCurvature(this->m_Origin);
        }
        break;
        case (2):
        {
          this->ShimshoniFrame(this->m_Origin);
        }
        break;
        case (3):
        {
          this->WeingartenMapGradients();
        }
        break;
        case (4):
        {
          kpix = this->ComputeMeanEuclideanDistance();
        }
        break;
        default:
        {
          this->WeingartenMapGradients();
        }
      }

      //      this->PrintFrame();

      //       bool geterror = false;
      //       if( geterror )
      //         {
      //         float error = 0.0;
      //         float temp1 = this->ErrorEstimate(this->GetOrigin() );
      //         float temp2 = this->ErrorEstimate(this->GetOrigin(), -1);
      //         if( temp1 < temp2 )
      //           {
      //           error = temp1;
      //           }
      //         else
      //           {
      //           error = temp2;
      //           this->SwitchNormalSign();
      // //         this->ComputeWeightsAndDirectionalKappaAndAngles(this->GetOrigin());
      // //         this->EstimateCurvature(this->m_A,this->m_B,this->m_B,this->m_C);
      // //         this->EstimateCurvature();
      //           }
      //         // std::cout << " best error " << error << std::endl;
      //         }

      kpix = 0;
      float fval = this->m_GaussianKappa;
      fval = this->m_MeanKappa;
      //      if( fabs(fval) > 1 )
      //        {
      //        fval = 0;
      //        }
      kpix = this->m_kSign * fval; // sulci
      if (std::isnan(kpix) || std::isinf(kpix))
      {
        this->m_Kappa1 = 0.0;
        this->m_Kappa2 = 0.0;
        this->m_MeanKappa = 0.0;
        this->m_GaussianKappa = 0.0;
        this->m_Area = 0.0;
        kpix = 0.0;
      }
      if (which == 5)
      {
        kpix = this->CharacterizeSurface();
      }
      if (which == 6)
      {
        kpix = this->m_GaussianKappa;
      }
      if (which == 7)
      {
        kpix = this->m_Area;
      }
      // ct++;
      this->m_PointList.clear();
    }
    // thresh += static_cast<double>(kpix);
    float offset = 0;
    if (which == 5)
    {
      offset = 0;
    }
    //    if (  ()  && ( )  ) kpix=0;
    this->m_FunctionImage->SetPixel(index, offset + kpix);
    // ct2++;
    ++ti;
  }
}

template <typename TSurface>
typename SurfaceImageCurvature<TSurface>::ImageType *
SurfaceImageCurvature<TSurface>::GetInput(void)
{
  if (this->GetNumberOfInputs() < 1)
  {
    return nullptr;
  }

  return static_cast<ImageType *>(this->ProcessObject::GetInput(0));
}

/**
 *
 */
template <typename TSurface>
typename SurfaceImageCurvature<TSurface>::OutputImageType *
SurfaceImageCurvature<TSurface>::GetOutput()
{
  return static_cast<OutputImageType *>(this->ProcessObject::GetOutput(0));
}

template <typename TSurface>
void
SurfaceImageCurvature<TSurface>::SetInputImage(typename ImageType::Pointer & input)
{
  this->ProcessObject::SetNthInput(0, input);
  this->m_ImageSize = input->GetLargestPossibleRegion().GetSize();
  typename OutputImageType::RegionType region = input->GetLargestPossibleRegion();
  if (!this->m_FunctionImage)
  {
    this->m_FunctionImage = OutputImageType::New();
    this->m_FunctionImage->CopyInformation(input);
    this->m_FunctionImage->SetRegions(region);
    this->m_FunctionImage->AllocateInitialized();
  }
}
} // namespace itk

#endif
