/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: itkSurfaceImageCurvature.hxx,v $
  Language:  C++
  Date:      $Date: 2008/11/15 23:46:06 $
  Version:   $Revision: 1.12 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _SurfaceImageCurvature_hxx
#define _SurfaceImageCurvature_hxx
#include "antsAllocImage.h"
#include <vnl/algo/vnl_real_eigensystem.h>

#include "itkSurfaceImageCurvature.h"
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

  ~GeodesicNode()
  {
  }
};

template <class pclass>
class GeodesicNodePriority /* defines the comparison operator for the prioritiy queue */
{
public:
  bool operator()( pclass N1, pclass N2)
  {
    return N1.distance > N2.distance;
  }
};

template <typename TSurface>
SurfaceImageCurvature<TSurface>
::SurfaceImageCurvature()
{
// outputs are curvature image and normal image
  this->ProcessObject::SetNumberOfRequiredOutputs( 2 );
  m_SurfaceLabel = 1;

  m_GradientImage = NULL;

  m_UseLabel = false;
  m_kSign = -1.0;
  m_FunctionImage = NULL;
  m_Sigma = 1.0;
}

template <typename TSurface>
void  SurfaceImageCurvature<TSurface>::ProcessLabelImage()
{
  ImageType* image = GetInput();

  if( !image )
    {
    return;
    }

  IndexType index;

  typename ImageType::RegionType requestedRegion;
  m_ImageSize = image->GetLargestPossibleRegion().GetSize();
  ImageIteratorType ti( image, image->GetLargestPossibleRegion() );

  ti.GoToBegin();
  while( !ti.IsAtEnd()  )
    {
    PixelType pix = ti.Get();
    index = ti.GetIndex();
    if( ti.Get() == 2 )
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
    IndexType cind = {{120, 120, 80}};
    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      d += (cind[j] - index[j]) * (cind[j] - index[j]);
      }

    d = sqrt(d);

    if( fabs(d - rad) <= 0.5 )
      {
      ti.Set(m_SurfaceLabel);
      }

    rad = 12, d = 0;
    IndexType cind2 = {{20, 170, 45}};
    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      d += (cind2[j] - index[j]) * (cind2[j] - index[j]);
      }

    d = sqrt(d);

//    if (fabs(d-rad) <= 0.5) ti.Set(m_SurfaceLabel);

    ++ti;
    }
}

template <typename TSurface>
void  SurfaceImageCurvature<TSurface>::FindEuclideanNeighborhood
  (typename SurfaceImageCurvature<TSurface>::PointType rootpoint)
{
  // Superclass::FindNeighborhood();
  this->m_AveragePoint = this->m_Origin;

  this->m_PointList.insert(this->m_PointList.begin(), this->m_Origin);

  typename ImageType::SizeType rad;
  IndexType oindex, index;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    rad[i] = (long)(m_NeighborhoodRadius);
    oindex[i] = (long) (rootpoint[i] + 0.5);
    }

  m_ti.SetLocation(oindex);

  unsigned int temp = 0;
  for( temp = 0; temp < m_ti.Size(); temp++ )
    {
    index = m_ti.GetIndex(temp);
    float pix = m_ti.GetPixel(temp);

//    if (pix == m_SurfaceLabel)
    if( this->IsValidSurface(pix, index) )
      {
      PointType p;
      float     dist = 0.0;
      bool      isorigin = true;
      for( unsigned int k = 0; k < ImageDimension; k++ )
        {
        if( index[k] != oindex[k] )
          {
          isorigin = false;
          }
        p[k] = (RealType) index[k];
        dist += (p(k) - this->m_Origin[k]) * (p(k) - this->m_Origin[k]);
        }
      dist = sqrt(dist);
      if( !isorigin && dist <= (m_NeighborhoodRadius) )
        {
        this->m_AveragePoint = this->m_AveragePoint + p;
        this->m_PointList.insert(this->m_PointList.begin(), p);
        }
      }
    }

  unsigned int npts = this->m_PointList.size();
  if( npts > 0 )
    {
    this->m_AveragePoint /= (RealType)npts;
    }
  else
    {
    this->m_AveragePoint = this->m_Origin;
    }

  if( this->m_Debug )
    {
    std::cout << " point list size " << this->m_PointList.size() << std::endl;
    //  for(int i = 0; i < this->m_PointList.size(); i++) {
    //    std:: cout << " point  " << this->m_PointList[i];
    //  }
    std::cout << std::endl;
    }
}

template <typename TSurface>
void  SurfaceImageCurvature<TSurface>::FindGeodesicNeighborhood()
{
  typedef std::priority_queue
    <GeodesicNode<ImageType>, std::vector<GeodesicNode<ImageType> >,
     GeodesicNodePriority<GeodesicNode<ImageType> > >  QType;

  QType                                            nodeq;
  std::map<unsigned int, GeodesicNode<ImageType> > nodes;

  this->m_AveragePoint = this->m_Origin;

  this->m_PointList.insert(this->m_PointList.begin(), this->m_Origin);

  typename ImageType::SizeType rad;
  IndexType oindex, index;

  float         dist = 0.0;
  unsigned int  k = 0;
  unsigned long longindex = 0;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    rad[i] = (long)(m_NeighborhoodRadius);
    oindex[i] = (long) (this->m_Origin[i] + 0.5);
    }
  index = oindex;
  for( k = 0; k < ImageDimension; k++ )
    {
    if( k == 0 )
      {
      longindex = index[0];
      }
    if( k == 1 )
      {
      longindex = index[1] + longindex + index[0] * m_ImageSize[0];
      }
    if( k == 2 )
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

  while( !nodeq.empty() && lastdist <= m_NeighborhoodRadius )
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
    if( lastdist <= m_NeighborhoodRadius )
      {
      m_ti2.SetLocation(g.imageindex);
      for( unsigned int jj = 0; jj < m_ti2.Size(); jj++ )
        {
        index = m_ti2.GetIndex(jj);

        if(  // m_ti2.GetPixel(jj) == m_SurfaceLabel &&
          this->IsValidSurface( m_ti2.GetPixel(jj), index) &&
          index[0] < m_ImageSize[0] - m_NeighborhoodRadius &&
          index[0] >  m_NeighborhoodRadius &&
          index[1] < m_ImageSize[1] - m_NeighborhoodRadius &&
          index[1] >  m_NeighborhoodRadius &&
          index[2] < m_ImageSize[2] - m_NeighborhoodRadius &&
          index[2] >  m_NeighborhoodRadius )
          {
          longindex = 0;
          dist = 0;
          for( k = 0; k < ImageDimension; k++ )
            {
            if( k == 0 )
              {
              longindex = index[0];
              }
            if( k == 1 )
              {
              longindex = index[1] + longindex + index[0] * m_ImageSize[0];
              }
            if( k == 2 )
              {
              longindex = index[2] + longindex + index[2] * m_ImageSize[0] * m_ImageSize[1];
              }
            q[k] = (RealType) index[k];
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
          if( !nodes[longindex].connected &&  (dist + lastdist) <= m_NeighborhoodRadius )
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

  this->m_AveragePoint = this->m_AveragePoint / ( (float)this->m_PointList.size() );
}

template <typename TSurface>
void  SurfaceImageCurvature<TSurface>::FindNeighborhood(unsigned int numMeanShifts)
{
  if( this->m_UseGeodesicNeighborhood )
    {
    this->FindGeodesicNeighborhood();
    }
  else
    {
    this->FindEuclideanNeighborhood(this->GetOrigin() );
    for( unsigned int dd = 0; dd < numMeanShifts; dd++ )
      {
      this->m_PointList.clear();
      this->FindEuclideanNeighborhood(this->GetAveragePoint() );
      }
    }

/*  if (this->m_Origin[0]==170 && this->m_Origin[1]==137 && this->m_Origin[2]==81)
  if ( this->m_Origin[1]==146 && this->m_Origin[0] > 167 )
  {
    std::cout << " origin " << this->m_Origin << std::endl;
    for (unsigned int tt=0; tt<this->m_PointList.size()-1; tt++)
    {
      PointType p=this->m_Origin-this->m_PointList[tt];
      float dist = p.magnitude();
      std::cout << " pt dist " << dist << " point " << this->m_PointList[tt] << std::endl;
    }

  }
*/
}

template <typename TSurface>
void  SurfaceImageCurvature<TSurface>
::LevelSetMeanCurvature()
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
void  SurfaceImageCurvature<TSurface>
::EstimateNormalsFromGradient()
{
  typename ImageType::Pointer image = GetInput();

  std::cout << " compute normals " << this->m_Sigma << " hood " << (this->m_NeighborhoodRadius)
                   << " spacing " << image->GetSpacing() <<  std::endl;

  if( !image )
    {
    return;
    }

  typename ImageType::SizeType rad;
  typename ImageType::SizeType rad2;
  for( unsigned int t = 0; t < ImageDimension; t++ )
    {
    rad[t] = (unsigned long) (this->m_NeighborhoodRadius);
    rad2[t] = 1;
    }
  this->m_ti.Initialize( rad, image, image->GetLargestPossibleRegion() );
  this->m_ti2.Initialize( rad2, image, image->GetLargestPossibleRegion() );

  typedef itk::ImageRegionIteratorWithIndex<TSurface> IteratorType;
  IteratorType Iterator( image, image->GetLargestPossibleRegion().GetSize() );
  bool         wmgmcurv = true;
  Iterator.GoToBegin();
  std::cout <<" ASS " << std::endl;
  while(  !Iterator.IsAtEnd()  )
    {
    float pix = Iterator.Get();
    if( pix != 0 && pix != 1 && pix != 2 )
      {
      wmgmcurv = false;
      }
    ++Iterator;
    }

  std::cout << " Using Binary Segmentation curv? " << wmgmcurv << std::endl;

  if( wmgmcurv )
    {
    typename OutputImageType::Pointer laplacian =
      AllocImage<OutputImageType>( image->GetLargestPossibleRegion() );
    laplacian->SetSpacing(image->GetSpacing() );
    Iterator.GoToBegin();
    while(  !Iterator.IsAtEnd()  )
      {
      IndexType ind = Iterator.GetIndex();
      if( image->GetPixel(ind) == 2 )
        {
        laplacian->SetPixel(ind, 1);
        }
      else if( image->GetPixel(ind) == 1 )
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
    for( unsigned int iterations = 0; iterations < totit; iterations++ )
      {
      std::cout << " % " << (float)iterations / (float)(totit) << std::endl;
      while(  !Iterator.IsAtEnd()  )
        {
        IndexType ind = Iterator.GetIndex();
        if( image->GetPixel(ind) == 2 )
          {
          laplacian->SetPixel(ind, 1);
          }
        else if( image->GetPixel(ind) == 0 )
          {
          laplacian->SetPixel(ind, 0.);
          }
        ++Iterator;
        }

      typedef itk::DiscreteGaussianImageFilter<TSurface, TSurface> dgf;
      typename dgf::Pointer filter = dgf::New();
      filter->SetVariance(0.5);
      filter->SetUseImageSpacingOn();
      filter->SetMaximumError(.01f);
      filter->SetInput(laplacian);
      filter->Update();
      laplacian = filter->GetOutput();
      Iterator.GoToBegin();
      }
//    WriteImage<TSurface>(laplacian,"lap.hdr");
//    std::cout << "Laplacian Solved " << std::endl;
    GradientImageFilterPointer filter = GradientImageFilterType::New();
    filter->SetInput( laplacian);
    RealType sigma = this->m_Sigma;
    filter->SetSigma( sigma );
    // Execute the filter
    filter->Update();
    this->m_GradientImage = filter->GetOutput();
    GradientPixelType zero;
    zero.Fill(0);
    Iterator.GoToBegin();
    while(  !Iterator.IsAtEnd()  )
      {
      IndexType ind = Iterator.GetIndex();
      if( image->GetPixel(ind) != 1 )
        {
        this->m_GradientImage->SetPixel(ind, zero);
        }
      ++Iterator;
      }
    }
  else
    {
    GradientImageFilterPointer filter = GradientImageFilterType::New();
    filter->SetInput(  image );
    RealType sigma = this->m_Sigma;
    filter->SetSigma( sigma );
    // Execute the filter
    filter->Update();
    this->m_GradientImage = filter->GetOutput();
    }
  std::cout << " compute normals done ";
}

template <typename TSurface>
void  SurfaceImageCurvature<TSurface>
::WeingartenMap()
{
  MatrixType D;

  float totwt = 0, wt = 0;

  unsigned int j = 0;
  unsigned int i = 0;
  unsigned int npts = this->m_PointList.size() - 1;

  if( npts < 4 )
    {
    this->m_MeanKappa = 0;
    this->m_GaussianKappa = 0;
    this->m_Kappa1 = 0;
    this->m_Kappa2 = 0;
    return;
    }

  unsigned int vars = 3;
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

// go through all the points
//  compute weight
//  compute dist of unit dif and tangents
//  compute dif of normal with grad at point

  PointType Q = this->m_Origin;
  for( j = 0; j < npts; j++ )
    {
    PointType Dif = Q - this->m_PointList[j];

    float u1 = 0.0;
    float u2 = 0.0;
    wt = 1.0; // /difmag;
    totwt += wt;

    IndexType index;
    for( i = 0; i < SurfaceDimension; i++ )
      {
      u1 += Dif[i] * this->m_Tangent1[i];
      u2 += Dif[i] * this->m_Tangent2[i];
      index[i] = (long) (this->m_PointList[j][i] + 0.5);
      }

//    if (fabs(u1) < 1.e-6) u1=1.e-6;
//    if (fabs(u2) < 1.e-6) u2=1.e-6;

    // get the normal at the point
    GradientPixelType norm = this->m_GradientImage->GetPixel(index);
    PointType         PN;

    float nmag = 0.0;
    for( i = 0; i < SurfaceDimension; i++ )
      {
      nmag += norm[i] * norm[i];
      PN[i] = norm[i];
      }

    nmag = sqrt(nmag);
    if( nmag >  1.e-9 )
      {
      PN /= (nmag);
      }
    else
      {
      PN *= 0.0;
      }

//    PointType dN=(PN-QN)*wt;

    xdists[j] = (PN[0]);
    ydists[j] = (PN[1]);
    zdists[j] = (PN[2]);
/*
    float a=0;
    float b=0;

    for (i=0; i<SurfaceDimension; i++)
    {
      a+=dN[i]*this->m_Tangent1[i];
      b+=dN[i]*this->m_Tangent2[i];
    }



    if ( u1*u1 > u2*u2 )
    {

      W(0,0)=W(0,0)+a;
      W(1,0)=W(1,0)+b;
    }
    else
    {
      W(0,1)=W(0,1)+a;
      W(1,1)=W(1,1)+b;
    }


   D(j,0)=u1*u1;
   D(j,1)=2.*u1*u2;
   D(j,2)=u2*u2;
   D(j,3)=u1;
   D(j,4)=u2;
   D(j,5)=1.0;
*/
// each row contains [u^2 , uv, v^2, u, v, 1] for point p

    if( vars == 6 )
      {
      D(j, 5) = u2 * u2; // (0   , 2*u2)
      D(j, 4) = u1 * u1; // (2*u1, 0)
      D(j, 3) = u1 * u2; // (u2  , u1)
      }
    D(j, 2) = u2; // (1   , 0)
    D(j, 1) = u1; // (0   , 1)
    D(j, 0) = 1.0;
    }

//  W=W/totwt;

  vnl_svd<double>    svd(D);
  vnl_vector<double> ax = svd.solve(xdists); // /totwt);
  vnl_vector<double> ay = svd.solve(ydists); // /totwt);
  vnl_vector<double> az = svd.solve(zdists); // /totwt);

// now get the first partials of each of these terms w.r.t. u and v

// dN/du = (dN/du \dot T_1) T_1+ (dNdu dot T_2) T_2

  PointType dNdu;
  dNdu[0] = ax[1];  dNdu[1] = ay[1];  dNdu[2] = az[1];
  PointType dNdv;
  dNdv[0] = ax[2];  dNdv[1] = ay[2];  dNdv[2] = az[2];

  float a = 0;
  float b = 0;

  float c = 0;
  float d = 0;
  for( i = 0; i < SurfaceDimension; i++ )
    {
    a += dNdu[i] * this->m_Tangent1[i];
    b += dNdv[i] * this->m_Tangent1[i];

    c += dNdu[i] * this->m_Tangent2[i];
    d += dNdv[i] * this->m_Tangent2[i];
    }
/*
    dNdu=a/(c+a)*this->m_Tangent1+c/(c+a)*this->m_Tangent2;
    dNdv=b/(b+d)*this->m_Tangent1+d/(b+d)*this->m_Tangent2;
    a=0; b=0; c=0; d=0;

    for (i=0; i<SurfaceDimension; i++)
    {
      a+=dNdu[i]*this->m_Tangent1[i];
      b+=dNdv[i]*this->m_Tangent1[i];

      c+=dNdu[i]*this->m_Tangent2[i];
      d+=dNdv[i]*this->m_Tangent2[i];
    }
*/
  W(0, 0) = a;
  W(0, 1) = b;
  W(1, 0) = c;
  W(1, 1) = d;

  // Compute estimated frame using eigensystem of D'*D
    {
    vnl_real_eigensystem eig(W);
//
    vnl_diag_matrix<vcl_complex<double> > DD(eig.D.rows() ); //
/*    for(i = 0; i < eig.D.n(); ++i) {
     vnl_test_assert("All real", vcl_imag(eig.D(i,i)) < 1e-15);
      DD(i,i) = vcl_real(eig.D(i,i));
    }

    vcl_cout << "D = " << eig.D << vcl_endl;
    vcl_cout << "V = " << eig.V << vcl_endl;
*/
    this->m_Kappa1 = vcl_real(eig.D(1, 1) );
    this->m_Kappa2 = vcl_real(eig.D(0, 0) );

// std::cout << " k1 " << this->m_Kappa1 << " k2 " << this->m_Kappa2 << " pt "<< this->m_Origin << std::endl;

    this->m_MeanKappa = (this->m_Kappa1 + this->m_Kappa2) * 0.5;
    this->m_GaussianKappa = (this->m_Kappa1 * this->m_Kappa2);
    }
}

template <typename TSurface>
void  SurfaceImageCurvature<TSurface>
::ComputeSurfaceArea()
{
  ImageType* image = GetInput();

  if( !image )
    {
    return;
    }

// BUG FIXME
  if( !this->m_GradientImage )
    {
    this->EstimateNormalsFromGradient();
    }

  IndexType index;

  typename ImageType::RegionType requestedRegion;
  this->m_ImageSize = image->GetLargestPossibleRegion().GetSize();
  ImageIteratorType ti( image, image->GetLargestPossibleRegion() );

  RealType area = 0.0;
  this->m_TotalArea = 0.0;

  ti.GoToBegin();
  unsigned int ct = 0;
  while( !ti.IsAtEnd()  )
    {
    index = ti.GetIndex();

    if(    // ti.Get() == this->m_SurfaceLabel &&
      (this->IsValidSurface(ti.Get(), index) ) &&
      index[0] < this->m_ImageSize[0] - this->m_NeighborhoodRadius &&
      index[0] >  this->m_NeighborhoodRadius &&
      index[1] < this->m_ImageSize[1] - this->m_NeighborhoodRadius &&
      index[1] >  this->m_NeighborhoodRadius &&
      index[2] < this->m_ImageSize[2] - this->m_NeighborhoodRadius &&
      index[2] >  this->m_NeighborhoodRadius )
      {
      ct++;
//        this->EstimateFrameFromGradient(index);
// BUG FIXME
//        area=this->ComputeLocalArea(spacing);
      area = 1.0;
      this->m_FunctionImage->SetPixel(index, area);
      this->m_TotalArea += area;
      this->m_PointList.clear();
      if( ct % 1000 == 0 )
        {
        std::cout << " ind " << index << " area " << area << std::endl;
        }
      }
    ++ti;
    }

  std::cout << " surface area " << this->m_TotalArea << std::endl;
  return;
}

template <typename TSurface>
void  SurfaceImageCurvature<TSurface>
::EstimateFrameFromGradient(IndexType index)
{
  GradientPixelType g = this->m_GradientImage->GetPixel(index);

  RealType mag = 0.0;

  for( int i = 0; i < ImageDimension; i++ )
    {
    this->m_Normal(i) = (RealType) g[i];
    mag += g[i] * g[i];
    }
  mag = sqrt(mag);
  if( mag <= 1.e-9 )
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
SurfaceImageCurvature<TSurface>
::IntegrateFunctionOverSurface(bool norm)
{
  typename OutputImageType::Pointer image = this->m_FunctionImage;

  if( !image )
    {
    std::cout << " no image " << std::endl; return 0;
    }

  //  std::cout << "  allocating temp image ";
  typename OutputImageType::Pointer tempimage = OutputImageType::New();
  tempimage->SetLargestPossibleRegion( image->GetLargestPossibleRegion() );
  tempimage->SetBufferedRegion( image->GetLargestPossibleRegion() );
  tempimage->Allocate();

  //  std::cout << "  done allocating  ";

  typename ImageType::SizeType rad;
  typename ImageType::SizeType rad2;
  for( unsigned int t = 0; t < ImageDimension; t++ )
    {
    rad[t] = (unsigned long) (this->m_NeighborhoodRadius);
    rad2[t] = 1;
    }

  this->m_ti.Initialize( rad, this->GetInput(), image->GetLargestPossibleRegion() );
  this->m_ti2.Initialize( rad2, this->GetInput(), image->GetLargestPossibleRegion() );

  IndexType index;

  typename ImageType::RegionType requestedRegion;
  ImageIteratorType ti( this->GetInput(), this->GetInput()->GetLargestPossibleRegion() );

  //  std::cout << " begin integrate ";

  ti.GoToBegin();
  unsigned int ct = 0;
  //  std::cout << " begin while " << std::endl;
  while( !ti.IsAtEnd()  )
    {
    index = ti.GetIndex();
    tempimage->SetPixel(index, 0);
    if(    // ti.Get() == this->m_SurfaceLabel &&
      (this->IsValidSurface(ti.Get(), index) ) &&
      index[0] < this->m_ImageSize[0] - this->m_NeighborhoodRadius &&
      index[0] >  this->m_NeighborhoodRadius &&
      index[1] < this->m_ImageSize[1] - this->m_NeighborhoodRadius &&
      index[1] >  this->m_NeighborhoodRadius &&
      index[2] < this->m_ImageSize[2] - this->m_NeighborhoodRadius &&
      index[2] >  this->m_NeighborhoodRadius )
      {
      PointType p;
      ct++;
      for( unsigned int k = 0; k < ImageDimension; k++ )
        {
        p[k] = (RealType) index[k];
        }
      this->SetOrigin(p);
//          std::cout << " find nhood ";
      this->FindNeighborhood();
//        std::cout << " get area ";
      RealType area = this->IntegrateFunctionOverNeighborhood(norm);
      tempimage->SetPixel(index, area);
      //      if( ct % 10000 == 0 )
      //        {
      //        std::cout << " area is : " << area << " ct " << ct << " pix " << ti.Get() << std::endl;
      //        }
//        if ( area > 1) std::cout << " ind " << index << " area " << area  << std::endl;
      // SD why sometimes a pixel is NaN ?
      //      if ( !(area > 0)) std::cout << " ind " << index << " area " << area << " pix " << ti.Get() <<
      // std::endl;
      }
    ++ti;
    }

  this->CopyImageToFunctionImage(tempimage, this->m_FunctionImage);

  return 0;
}

template <typename TSurface>
void
SurfaceImageCurvature<TSurface>
::CopyImageToFunctionImage(OutputImagePointer i1, OutputImagePointer i2)
{
  if( !i1 || !i2 )
    {
    return;
    }

  typename ImageType::RegionType requestedRegion;
  OutputImageIteratorType ti1( i1, i1->GetLargestPossibleRegion() );
  OutputImageIteratorType ti2( i2, i2->GetLargestPossibleRegion() );

  ti1.GoToBegin();
  ti2.GoToBegin();
  while( !ti1.IsAtEnd()  )
    {
    ti2.Set(ti1.Get() );
    ++ti1;
    ++ti2;
    }
}

template <typename TSurface>
typename SurfaceImageCurvature<TSurface>::RealType
SurfaceImageCurvature<TSurface>
::IntegrateFunctionOverNeighborhood(bool norm)
{
  unsigned int npts = this->m_PointList.size();
  double       curvature = 0.0, tw = 0;

  //  std::cout << " npts " << npts;

  for( unsigned int pp = 0; pp < npts; pp++ )
    {
    IndexType localindex;
    for( unsigned int k = 0; k < ImageDimension; k++ )
      {
      localindex[k] = (long) this->m_PointList[pp][k];
      }
    PointType dd = this->m_Origin - this->m_PointList[pp];
    double    wi = dd.magnitude();
    if( wi != 0.0 )
      {
      wi = 1. / wi;
      }
    tw += wi;
    RealType func = this->m_FunctionImage->GetPixel( localindex );
    //    std::cout << " pp " << pp << " func " << func << std::endl;
    if( norm )
      {
      curvature += wi * func;
      }
    else
      {
      curvature += func;
      }
//    curvature*=this->ComputeLocalArea(spacing);
    }
  // if (norm ) curvature/=tw;
  // SD sometimes tw is zero making curvature = NaN
  if( norm && tw != 0 )
    {
    curvature /= tw;
    }
  this->m_PointList.clear();

  return curvature;
}

template <typename TSurface>
void  SurfaceImageCurvature<TSurface>
::PostProcessGeometry()
{
  typename ImageType::Pointer  image = GetInput();

  if( !image )
    {
    return;
    }

  IndexType index;

  typename ImageType::RegionType requestedRegion;
  this->m_ImageSize = image->GetLargestPossibleRegion().GetSize();
  ImageIteratorType ti( image, image->GetLargestPossibleRegion() );

  std::vector<double> kvec;
  ti.GoToBegin();
  while( !ti.IsAtEnd()  )
    {
    PixelType pix = ti.Get();
    index = ti.GetIndex();
    if(  // ti.Get() == this->m_SurfaceLabel &&
      (this->IsValidSurface(ti.Get(), index) ) &&
      index[0] < this->m_ImageSize[0] - this->m_NeighborhoodRadius &&
      index[0] >  this->m_NeighborhoodRadius &&
      index[1] < this->m_ImageSize[1] - this->m_NeighborhoodRadius &&
      index[1] >  this->m_NeighborhoodRadius &&
      index[2] < this->m_ImageSize[2] - this->m_NeighborhoodRadius &&
      index[2] >  this->m_NeighborhoodRadius ) //
      {
      PointType p;
      for( unsigned int k = 0; k < ImageDimension; k++ )
        {
        p[k] = (RealType) index[k];
        }
      this->SetOrigin(p);
      this->FindNeighborhood();
      int    npts = this->m_PointList.size() - 1;
      double curvature = 0.0;
      for( int pp = 0; pp < npts; pp++ )
        {
        IndexType localindex;
        for( unsigned int k = 0; k < ImageDimension; k++ )
          {
          localindex[k] = (long) this->m_PointList[pp][k];
          }
//          PointType dd=this->m_Origin-this->m_PointList[pp];
//          double wi=dd.magnitude();
//          if (wi!=0.0) wi=1./wi;
//          tw+=wi;vector<int> vec;

        curvature = this->m_FunctionImage->GetPixel( localindex );
        kvec.push_back(curvature);
        }

      std::sort(kvec.begin(), kvec.end() ); // Sort the vector

      this->m_PointList.clear();
//        curvature/=tw;
      this->m_FunctionImage->SetPixel(index, kvec[kvec.size() / 2]);
      kvec.clear();
      }
    ++ti;
    }
}

template <typename TSurface>
void  SurfaceImageCurvature<TSurface>
::ComputeFrameOverDomain(unsigned int which)
{
  ImageType* image = GetInput();

  if( !image )
    {
    return;
    }

  IndexType index;

  typename ImageType::RegionType requestedRegion;
  this->m_ImageSize = image->GetLargestPossibleRegion().GetSize();
  ImageIteratorType ti( image, image->GetLargestPossibleRegion() );

  // std::exception();
// Get Normals First!
  this->EstimateNormalsFromGradient();

  unsigned int  ct = 1;
  unsigned long ct2 = 0;
  RealType      kpix = 0;

  double thresh = 0.0;

  ti.GoToBegin();
  while( !ti.IsAtEnd()  )
    {
    index = ti.GetIndex();
    if( ct2 % 200000 == 0 && ct2 > 0 )
      {
      std::cout << " ind " << index << " kp " << kpix << std::endl;
      }
    kpix = 0.0;
    if(  // ti.Get() == this->m_SurfaceLabel &&
      this->IsValidSurface(ti.Get(), index) &&
      index[0] < this->m_ImageSize[0] - 2 * this->m_NeighborhoodRadius &&
      index[0] >  2 * this->m_NeighborhoodRadius &&
      index[1] < this->m_ImageSize[1] - 2 * this->m_NeighborhoodRadius &&
      index[1] >  2 * this->m_NeighborhoodRadius &&
      index[2] < this->m_ImageSize[2] - 2 * this->m_NeighborhoodRadius &&
      index[2] >  2 * this->m_NeighborhoodRadius ) //
      {
      // std::cout << " val " << (RealType) ti.Get() << std::endl;
      PointType p;
      for( unsigned int k = 0; k < ImageDimension; k++ )
        {
        p[k] = (RealType) index[k];
        }
      this->SetOrigin(p);
      this->EstimateFrameFromGradient(index);
      this->FindNeighborhood();

      switch( which )
        {
        case ( 0 ):
          {
          this->ComputeJoshiFrame( this->m_Origin);
          }
          break;
        case ( 1 ):
          {
          this->JainMeanAndGaussianCurvature( this->m_Origin);
          }
          break;
        case ( 2 ):
          {
          this->ShimshoniFrame(this->m_Origin);
          }
          break;
        case ( 3 ):
          {
          this->WeingartenMap();
          }
          break;
        case ( 4 ):
          {
          kpix = this->ComputeMeanEuclideanDistance();
          }
          break;
        default:
          {
          this->WeingartenMap();
          }
        }

      //      this->PrintFrame();

      bool geterror = false;
      if( geterror )
        {
        float error = 0.0;
        float temp1 = this->ErrorEstimate(this->GetOrigin() );
        float temp2 = this->ErrorEstimate(this->GetOrigin(), -1);
        if( temp1 < temp2 )
          {
          error = temp1;
          }
        else
          {
          error = temp2;
          this->SwitchNormalSign();
//         this->ComputeWeightsAndDirectionalKappaAndAngles(this->GetOrigin());
//         this->EstimateCurvature(this->m_A,this->m_B,this->m_B,this->m_C);
//         this->EstimateCurvature();
          }
        std::cout << " best error " << error << std::endl;
        }

//    kpix=fabs(2.0/(3.1416)*atan((this->m_Kappa1+this->m_Kappa2)/(this->m_Kappa2-this->m_Kappa1)));
//    if (kpix >= 0)
//    kpix = (sqrt(this->m_Kappa1*this->m_Kappa1+this->m_Kappa2*this->m_Kappa2));
//      kpix=1.0;
//    if (this->m_MeanKappa > 0.)  kpix +=fabs(this->m_MeanKappa);
//      kpix =fabs(this->m_MeanKappa);
//    kpix = kpix+10.0;
//    if (this->m_Kappa1 > 0) kpix+=fabs(this->m_Kappa1);
//    if (this->m_Kappa2 > 0) kpix+=fabs(this->m_Kappa2);
//    else kpix = -1.0*(sqrt(this->m_Kappa1*this->m_Kappa1+this->m_Kappa2*this->m_Kappa2));

//    std::cout << " kpix " << kpix << " thresh " << thresh << std::endl;

//
//    if ( fabs(kpix) >  100 ) kpix=0.0;
//    else if (kpix < -1.0*thresh/ct ) kpix=2;
//    else kpix=0;

      kpix = 1.1;
      float fval = this->m_GaussianKappa;
      fval = this->m_MeanKappa;
//      fval=fval*(1.0+fabs(this->m_Kappa2-this->m_Kappa1));
      if( fabs(fval) > 1 )
        {
        fval /= fval;
        }
//      if (fval > 0.0) kpix+=(fval); // gyri
//      if (fval < 0.0) kpix-=(fval); // sulci bright
      kpix = this->m_kSign * fval; // sulci
      if( vnl_math_isnan(kpix)  || vnl_math_isinf(kpix) )
        {
        this->m_Kappa1 = 0.0;
        this->m_Kappa2 = 0.0;
        this->m_MeanKappa = 0.0;
        this->m_GaussianKappa = 0.0;
        kpix = 0.0;
        }
//      kpix=1.0+this->m_Kappa1*this->m_Kappa1+this->m_Kappa2*this->m_Kappa2;
//    kpix=fabs(this->m_Normal[0]);
      if( which == 5 )
        {
        kpix = this->CharacterizeSurface();
        }
      if( which == 6 )
        {
        kpix = this->m_GaussianKappa;
        }
      ct++;
      this->m_PointList.clear();
      }
    thresh += kpix;
    float offset = 0;
    if( fabs( image->GetPixel(index) - 0 ) > 1.e-6  )
      {
      offset = 128.0;
      }
    if( which == 5 )
      {
      offset = 0;
      }
    this->m_FunctionImage->SetPixel(index, offset + kpix);
    ct2++;
    ++ti;
    }

  std::cout << " average curvature " << thresh / (float)ct << " kSign " << this->m_kSign <<  std::endl;

/* now get s.d.
  float sd=0.0;
  float avgc=thresh/(float)ct;
  ti.GoToBegin();
  while(!ti.IsAtEnd()  )
  {
    PixelType pix=ti.Get();
    if ( //ti.Get() == this->m_SurfaceLabel
      this->IsValidSurface(ti.Get(),index)
    )
      sd+=(this->m_FunctionImage->GetPixel(ti.GetIndex())-avgc)*
          (this->m_FunctionImage->GetPixel(ti.GetIndex())-avgc);
    ++ti;
  }

  sd=sqrt(sd/(float)ct);
  float nsd=1.0;
  ti.GoToBegin();
  while(!ti.IsAtEnd()  )
  {
    PixelType pix=ti.Get();
    if ( //ti.Get() == this->m_SurfaceLabel
        this->IsValidSurface(ti.Get(),index) &&
        (this->m_FunctionImage->GetPixel(ti.GetIndex()) - avgc) > nsd*sd )
      this->m_FunctionImage->SetPixel(ti.GetIndex(),avgc+nsd*sd);

    ++ti;
  }
*/
}

template <typename TSurface>
typename SurfaceImageCurvature<TSurface>::ImageType
* SurfaceImageCurvature<TSurface>
::GetInput(void)
  {
  if( this->GetNumberOfInputs() < 1 )
    {
    return NULL;
    }

  return static_cast<ImageType *>
         (this->ProcessObject::GetInput(0) );
  }

/**
 *
 */
template <typename TSurface>
typename SurfaceImageCurvature<TSurface>::OutputImageType
* SurfaceImageCurvature<TSurface>
::GetOutput()
  {
  return static_cast<OutputImageType *>(this->ProcessObject::GetOutput(0) );
  }

template <typename TSurface>
void
SurfaceImageCurvature<TSurface>
::SetInputImage(typename  ImageType::Pointer & input)
{
  this->ProcessObject::SetNthInput(0,  input);

  this->m_ImageSize = input->GetLargestPossibleRegion().GetSize();

  typename OutputImageType::RegionType region;
  region.SetSize( this->m_ImageSize );

  if( !this->m_FunctionImage )
    {
    this->m_FunctionImage = OutputImageType::New();
    this->m_FunctionImage->SetLargestPossibleRegion( region );
    this->m_FunctionImage->SetBufferedRegion( region );
    this->m_FunctionImage->Allocate();
    this->m_FunctionImage->SetSpacing( input->GetSpacing() );
    this->m_FunctionImage->SetDirection( input->GetDirection() );
    this->m_FunctionImage->SetOrigin( input->GetOrigin() );
    }
  // this->ProcessLabelImage();

  // this->ProcessObject::SetNthOutput( 0, this->m_FunctionImage );
}
} // namespace itk

#endif
