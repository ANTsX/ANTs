/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for detailm_Solver.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _FEMDiscConformalMap_hxx
#define _FEMDiscConformalMap_hxx

#include <vcl_compiler.h>
#include <iostream>
#include <cmath>
#include <iostream>
#include <vnl/vnl_real_polynomial.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_vector_fixed.h>
#include <itkMath.h>
#include "vtkDelaunay2D.h"
#include "vtkSelectPolyData.h"
#include "vtkFloatArray.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkFEMDiscConformalMap.h"
#include "itkFEMLoadNode.h"
#include "itkSurfaceMeshCurvature.h"
#include "vtkClipPolyData.h"
#include "vtkContourFilter.h"
#include "vtkSmartPointer.h"

namespace itk
{
template <typename TSurface, typename TImage, unsigned int TDimension>
FEMDiscConformalMap<TSurface, TImage, TDimension>::FEMDiscConformalMap()
{
  m_Sigma = 2.0e-4;
  m_ParameterFileName = "";
  m_NorthPole = 0;
  m_SourceNodeNumber = 1;

  this->m_DistanceCostWeight = 1;
  this->m_LabelCostWeight = 0;

  m_Pi = 3.14159265358979323846;
  m_ReadFromFile = false;
  m_Debug = false;
  m_FindingRealSolution = true;
  this->m_SurfaceMesh = nullptr;
  for (int i = 0; i < 7; i++)
  {
    m_PoleElementsGN[i] = 0;
  }
  this->m_MapToCircle = true;
  this->m_MapToSquare = false;
  this->m_ParamWhileSearching = true;
  m_Smooth = 2.0;
  this->m_Label_to_Flatten = 0;
  this->m_FlatImage = nullptr;

  manifoldIntegrator = ManifoldIntegratorType::New();
}

template <typename TSurface, typename TImage, unsigned int TDimension>
bool
FEMDiscConformalMap<TSurface, TImage, TDimension>::InBorder(
  typename FEMDiscConformalMap<TSurface, TImage, TDimension>::GraphSearchNodePointer g)
{
  if (this->m_HelpFindLoop[g->GetIdentity()] > 0)
  {
    return true;
  }
  else
  {
    return false;
  }

  float dist = g->GetTotalCost();
  if (dist > m_MaxCost && dist < m_MaxCost + 1.0)
  {
    return true;
  }
  return false;
}

template <typename TSurface, typename TImage, unsigned int TDimension>
bool
FEMDiscConformalMap<TSurface, TImage, TDimension>::InDisc(
  typename FEMDiscConformalMap<TSurface, TImage, TDimension>::GraphSearchNodePointer g)
{
  //  if ( this->m_HelpFindLoop[ g->GetIdentity() ] != 0   ) return true;
  // else return false;
  //    if (  g->WasVisited() ) return true;

  float d = g->GetTotalCost();

  for (unsigned int i = 0; i < this->m_DiscBoundaryList.size(); i++)
  {
    if (d <= this->m_DiscBoundaryList[i]->GetTotalCost())
    {
      return true;
    }
  }
  return false;

  if (g->GetPredecessor())
  {
    return true;
  }
  else
  {
    return false;
  }
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMDiscConformalMap<TSurface, TImage, TDimension>::FindSource(IndexType index)
{
  vtkPoints * vtkpoints = this->m_SurfaceMesh->GetPoints();
  int         numPoints = vtkpoints->GetNumberOfPoints();
  float       mindist = 1.e9;

  for (int i = 0; i < numPoints; i++)
  {
    double * pt = vtkpoints->GetPoint(i);
    float    dist = 0.0;
    for (int j = 0; j < ImageDimension; j++)
    {
      dist += (pt[j] - (float)index[j]) * (pt[j] - (float)index[j]);
    }
    dist = sqrt(dist);
    if (dist < mindist)
    {
      mindist = dist;
      m_SourceNodeNumber = i;
    }
  }
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMDiscConformalMap<TSurface, TImage, TDimension>::FindMeanSourceInLabel(unsigned int label)
{
  vtkPoints *   vtkpoints = this->m_SurfaceMesh->GetPoints();
  int           numPoints = vtkpoints->GetNumberOfPoints();
  float         mindist = 1.e9;
  float         meanx = 0., meany = 0, meanz = 0;
  unsigned long ct = 0;

  vtkDataArray * labels = this->m_SurfaceMesh->GetPointData()->GetArray("Label");
  vtkDataArray * features = nullptr;

  if (this->m_SurfaceFeatureMesh)
  {
    if (this->m_SurfaceFeatureMesh->GetPointData()->GetArray("Feature"))
    {
      features = this->m_SurfaceFeatureMesh->GetPointData()->GetArray("Feature");
    }
    else
    {
      features = labels;
    }
  }

  if (!labels)
  {
    std::cout << " cant get Labels --- need an array named Label in the mesh ... " << std::endl;
    std::exception();
  }
  else
  {
    std::cout << " got Labels " << std::endl;
  }
  for (int i = 0; i < numPoints; i++)
  {
    double * pt = vtkpoints->GetPoint(i);
    manifoldIntegrator->GetGraphNode(i)->SetValue(labels->GetTuple1(i), 3);
    //    std::cout << " label " << labels->GetTuple1(i) <<  " & " <<  label <<std::endl;
    if (fabs(labels->GetTuple1(i) - label) < 0.5)
    {
      //      std::cout << " choose " <<  labels->GetTuple1(i) << " & " <<  label << std::endl;
      meanx += pt[0];
      meany += pt[1];
      meanz += pt[2];
      ct++;
      //      i=numPoints+1;
    }
    else
    {
      //      ct+=0;
      manifoldIntegrator->GetGraphNode(i)->SetUnVisitable();
    }
  }

  meanx /= (float)ct;
  meany /= (float)ct;
  meanz /= (float)ct;
  for (int i = 0; i < numPoints; i++)
  {
    double * pt = vtkpoints->GetPoint(i);
    float    dist = 0.0;
    dist += (pt[0] - meanx) * (pt[0] - meanx);
    dist += (pt[1] - meany) * (pt[1] - meany);
    dist += (pt[2] - meanz) * (pt[2] - meanz);
    dist = sqrt(dist);
    if (dist < mindist && fabs(label - labels->GetTuple1(i)) < 0.5)
    {
      mindist = dist;
      m_SourceNodeNumber = i;
      //    std::cout << "  label " << label  << " chose " << labels->GetTuple1(i)  << std::endl;
    }
  }
  if (this->FindLoopAroundNode(this->m_SourceNodeNumber) == 2)
  {
    std::cout << " found loop " << std::endl;
  }

  if (ct > 0)
  {
    std::cout << meanx << "  " << meany << "  " << meanz << std::endl;
  }
  else
  {
    std::cout << " no label " << label << " exiting " << std::endl;
    std::exception();
  }
}

template <typename TSurface, typename TImage, unsigned int TDimension>
float
FEMDiscConformalMap<TSurface, TImage, TDimension>::AssessNodeDistanceCost(unsigned int nodeid)
{
  return manifoldIntegrator->GetGraphNode(nodeid)->GetTotalCost();

  /*
  diff=G1->GetLocation()-G3->GetLocation();
  tangentlength=0;
  for (unsigned int d=0; d<diff.Size(); d++) tangentlength+=diff[d]*diff[d];
  elementlength+=sqrt(tangentlength);
  diff=G2->GetLocation()-G3->GetLocation();
  tangentlength=0;
  for (unsigned int d=0; d<diff.Size(); d++) tangentlength+=diff[d]*diff[d];
  elementlength+=sqrt(tangentlength);
  */
}

template <typename TSurface, typename TImage, unsigned int TDimension>
float
FEMDiscConformalMap<TSurface, TImage, TDimension>::GetBoundaryParameterForSquare(unsigned int nodeid,
                                                                                 unsigned int whichParam)
{
  // compute loop length and check validity
  float tangentlength = 0, totallength = 0, nodeparam = 0, x = 0, y = 0;

  typename GraphSearchNodeType::NodeLocationType diff;
  for (unsigned int i = 0; i < this->m_DiscBoundaryList.size(); i++)
  {
    unsigned int nxt = ((i + 1) % this->m_DiscBoundaryList.size());
    diff = this->m_DiscBoundaryList[i]->GetLocation() - this->m_DiscBoundaryList[nxt]->GetLocation();
    if (this->m_DiscBoundaryList[i]->GetIdentity() == nodeid)
    {
      nodeparam = totallength;
    }
    tangentlength = 0;
    for (unsigned int d = 0; d < diff.Size(); d++)
    {
      tangentlength += diff[d] * diff[d];
    }
    totallength += tangentlength;
  }

  float arclength = nodeparam / totallength;

  if (arclength <= 0.25)
  {
    x = 0;
    y = arclength / 0.25; /* 0 => 1 */
  }
  else if (arclength <= 0.5)
  {
    x = (0.5 - arclength) / 0.25; /* 0 => 1 */
    y = 1;
  }
  else if (arclength <= 0.75)
  {
    x = 1;
    y = (0.75 - arclength) / 0.25; /* 1 => 0 */
  }
  else if (arclength <= 1)
  {
    x = (1 - arclength) / 0.25; /* 1 => 0 */
    y = 0;
  }
  //  std::cout <<" x " << x << " y " << y << " al " << arclength << std::endl;

  if (whichParam == 0)
  {
    return x;
  }
  else if (whichParam == 1)
  {
    return y;
  }
  else
  {
    return arclength;
  }
}

template <typename TSurface, typename TImage, unsigned int TDimension>
float
FEMDiscConformalMap<TSurface, TImage, TDimension>::GetBoundaryParameterForCircle(unsigned int nodeid,
                                                                                 unsigned int whichParam)
{
  // compute loop length and check validity
  float tangentlength = 0, totallength = 0, nodeparam = 0;

  typename GraphSearchNodeType::NodeLocationType diff;
  for (unsigned int i = 0; i < this->m_DiscBoundaryList.size(); i++)
  {
    unsigned int nxt = ((i + 1) % this->m_DiscBoundaryList.size());
    diff = this->m_DiscBoundaryList[i]->GetLocation() - this->m_DiscBoundaryList[nxt]->GetLocation();
    if (this->m_DiscBoundaryList[i]->GetIdentity() == nodeid)
    {
      nodeparam = totallength;
    }
    tangentlength = 0;
    for (unsigned int d = 0; d < diff.Size(); d++)
    {
      tangentlength += diff[d] * diff[d];
    }
    totallength += tangentlength;
  }

  float arclength = nodeparam / totallength * M_PI * 2;

  if (whichParam == 0)
  {
    return cos(arclength);
  }
  else if (whichParam == 1)
  {
    return sin(arclength);
  }
  else
  {
    return arclength;
  }
}

template <typename TSurface, typename TImage, unsigned int TDimension>
unsigned int
FEMDiscConformalMap<TSurface, TImage, TDimension>::AddVertexToLoop()
{
  // compute loop length and check validity
  float tangentlength = 0, totallength = 0;
  bool  isvalid = true;

  typename GraphSearchNodeType::NodeLocationType diff;
  for (unsigned int i = 0; i < this->m_DiscBoundaryList.size(); i++)
  {
    unsigned int nxt = ((i + 1) % this->m_DiscBoundaryList.size());
    //   std::cout << " this->m_DiscBoundaryList[i]->GetLocation() " <<
    // this->m_DiscBoundaryList[i]->GetLocation() <<
    // std::endl;
    diff = this->m_DiscBoundaryList[i]->GetLocation() - this->m_DiscBoundaryList[nxt]->GetLocation();
    tangentlength = 0;
    for (unsigned int d = 0; d < diff.Size(); d++)
    {
      tangentlength += diff[d] * diff[d];
    }
    totallength += tangentlength;
    // check if the next node is in the current node's neighborlist
    isvalid = false;
    for (unsigned int n = 0; n < this->m_DiscBoundaryList[i]->m_NumberOfNeighbors; n++)
    {
      if (this->m_DiscBoundaryList[i]->m_Neighbors[n] == this->m_DiscBoundaryList[nxt])
      {
        isvalid = true;
      }
    }
  }
  std::cout << " length " << totallength << " valid? " << isvalid << " entries " << this->m_DiscBoundaryList.size()
            << " gsz " << this->m_HelpFindLoop.size() << std::endl;

  /** now find a node with  HelpFindLoop value == 0 that minimally changes the length ...
      and that has a root neighbor and next neighbor consistent with the original edge */
  unsigned int           newnodeid = 0;
  unsigned int           nodeparam = 0;        // the position in the curve parameter
  float                  newedgelength = 1.e9; // minimize this
  GraphSearchNodePointer G1, G2, G3, BestG = nullptr;
  for (unsigned int i = 0; i < this->m_DiscBoundaryList.size(); i++)
  {
    unsigned int nxt = ((i + 1) % this->m_DiscBoundaryList.size());
    G1 = this->m_DiscBoundaryList[i];
    G2 = this->m_DiscBoundaryList[nxt];
    G3 = nullptr;
    /** 3 conditions : (1) HelpFindLoop == 0  (2) am neighbor of curnode (3) am neighbor of next node (4) minlength */
    for (unsigned int n = 0; n < G1->m_NumberOfNeighbors; n++)
    {
      for (unsigned int m = 0; m < G2->m_NumberOfNeighbors; m++)
      {
        long hfl = abs(this->m_HelpFindLoop[G1->m_Neighbors[n]->GetIdentity()]);
        hfl += abs(this->m_HelpFindLoop[G2->m_Neighbors[m]->GetIdentity()]);
        if (G1->m_Neighbors[n] == G2->m_Neighbors[m] && hfl == 0)
        {
          G3 = G1->m_Neighbors[n];
          n = G1->m_NumberOfNeighbors;
          m = G2->m_NumberOfNeighbors;

          float elementlength = this->AssessNodeDistanceCost(G3->GetIdentity());

          if (elementlength < newedgelength)
          {
            newedgelength = elementlength;
            BestG = G3;
            newnodeid = G3->GetIdentity();
            nodeparam = i;
          }
        }
      }
    }
  }
  if (!BestG)
  {
    std::cout << " does not exist " << std::endl;
    return false;
  }
  if (this->m_HelpFindLoop[BestG->GetIdentity()] != 0)
  {
    std::cout << " already done " << std::endl;
    return false;
  }
  std::vector<GraphSearchNodePointer> neighborlist;
  long                                paramct = 0;
  totallength = 0;
  for (unsigned int i = 0; i < this->m_DiscBoundaryList.size(); i++)
  {
    paramct++;
    neighborlist.push_back(this->m_DiscBoundaryList[i]);
    this->m_HelpFindLoop[this->m_DiscBoundaryList[i]->GetIdentity()] = paramct;
    if (i == nodeparam)
    {
      paramct++;
      neighborlist.push_back(BestG);
      this->m_HelpFindLoop[BestG->GetIdentity()] = paramct;
    }
  }
  //  std::cout << " len1 " << this->m_DiscBoundaryList.size()  <<  " len2 " <<  neighborlist.size() <<
  // std::endl;
  this->SetDiscBoundaryList(neighborlist);
  //  this->m_DiscBoundaryList.assign(neighborlist.begin(),neighborlist.end());
  /*
  for (unsigned int i=0; i<this->m_DiscBoundaryList.size(); i++) {
    unsigned int nxt=( ( i+1 ) % this->m_DiscBoundaryList.size() );
    diff=this->m_DiscBoundaryList[i]->GetLocation()-this->m_DiscBoundaryList[nxt]->GetLocation();
    tangentlength=0;
    for (unsigned int d=0; d<diff.Size(); d++) tangentlength+=diff[d]*diff[d];
    totallength+=tangentlength;
  }
  */
  /** find out if there is a "short-cut" that lets you traverse the loop without hitting all the nodes ... */
  /** at the same time, you want to maximize the enclosed area */
  neighborlist.clear();
  paramct = 0;
  for (unsigned int i = 0; i < this->m_DiscBoundaryList.size(); i++)
  {
    paramct++;
    neighborlist.push_back(this->m_DiscBoundaryList[i]);
    this->m_HelpFindLoop[this->m_DiscBoundaryList[i]->GetIdentity()] = paramct;
    if (this->m_DiscBoundaryList[i]->GetIdentity() == newnodeid)
    {
      // find the parameter of any of this fellows neighors ...
      // if it's greater than i+1 then skip ahead
      long bestparam = this->m_HelpFindLoop[this->m_DiscBoundaryList[i]->GetIdentity()] + 1;
      for (unsigned int n = 0; n < this->m_DiscBoundaryList[i]->m_NumberOfNeighbors; n++)
      {
        if (this->m_HelpFindLoop[this->m_DiscBoundaryList[i]->m_Neighbors[n]->GetIdentity()] > bestparam)
        {
          bestparam = this->m_HelpFindLoop[this->m_DiscBoundaryList[i]->m_Neighbors[n]->GetIdentity()];
          BestG = this->m_DiscBoundaryList[i]->m_Neighbors[n];
        }
      }
      //      neighborlist.push_back( BestG );
      for (unsigned int j = i + 1; j < (unsigned int)bestparam - 1; j++)
      {
        this->m_HelpFindLoop[this->m_DiscBoundaryList[j]->GetIdentity()] = -1;
      }
      i = (unsigned int)(bestparam - 2);
    }
  }
  this->SetDiscBoundaryList(neighborlist);
  //  this->m_DiscBoundaryList.assign(neighborlist.begin(),neighborlist.end());
  float newtotallength = 0;
  for (unsigned int i = 0; i < this->m_DiscBoundaryList.size(); i++)
  {
    unsigned int nxt = ((i + 1) % this->m_DiscBoundaryList.size());
    diff = this->m_DiscBoundaryList[i]->GetLocation() - this->m_DiscBoundaryList[nxt]->GetLocation();
    tangentlength = 0;
    for (unsigned int d = 0; d < diff.Size(); d++)
    {
      tangentlength += diff[d] * diff[d];
    }
    newtotallength += tangentlength;
  }
  std::cout << " total1 " << totallength << " total2 " << newtotallength << std::endl;
  // check for self-intersection
  for (unsigned int i = 0; i < this->m_DiscBoundaryList.size(); i++)
  {
    for (unsigned int j = 0; j < this->m_DiscBoundaryList.size(); j++)
    {
      if (i != j && this->m_DiscBoundaryList[i]->GetIdentity() == this->m_DiscBoundaryList[j]->GetIdentity())
      {
        isvalid = false;
      }
    }
  }
  for (unsigned long g = 0; g < this->m_HelpFindLoop.size(); g++)
  {
    if (this->m_HelpFindLoop[g] != 0)
    {
      this->m_HelpFindLoop[g] = -1;
    }
  }
  for (unsigned int j = 0; j < this->m_DiscBoundaryList.size(); j++)
  {
    m_HelpFindLoop[this->m_DiscBoundaryList[j]->GetIdentity()] = j + 1;
  }

  unsigned long MAXLIST = 200;
  if (this->m_DiscBoundaryList.size() > MAXLIST)
  {
    std::cout << this->m_RootNode->GetLocation() << std::endl;
    std::cout << "long list " << MAXLIST << std::endl;
    return 2;
  }
  return isvalid;
}

template <typename TSurface, typename TImage, unsigned int TDimension>
unsigned int
FEMDiscConformalMap<TSurface, TImage, TDimension>::FindLoopAroundNode(unsigned int j_in)
{
  vtkDataArray * labels = this->m_SurfaceMesh->GetPointData()->GetArray("Label");
  vtkDataArray * features = nullptr;

  if (this->m_SurfaceFeatureMesh)
  {
    if (this->m_SurfaceFeatureMesh->GetPointData()->GetArray("Feature"))
    {
      features = this->m_SurfaceFeatureMesh->GetPointData()->GetArray("Feature");
    }
    else
    {
      features = labels;
    }
  }

  this->m_DiscBoundaryList.clear();
  unsigned int gsz = manifoldIntegrator->GetGraphSize();
  // get distance from this node to all others
  // now measure the length distortion of the given solution
  this->m_RootNode = manifoldIntegrator->GetGraphNode(j_in);
  for (int i = 0; i < manifoldIntegrator->GetGraphSize(); i++)
  {
    manifoldIntegrator->GetGraphNode(i)->SetTotalCost(vnl_huge_val(manifoldIntegrator->GetMaxCost()));
    manifoldIntegrator->GetGraphNode(i)->SetUnVisited();
    manifoldIntegrator->GetGraphNode(i)->SetValue(0.0, 1);
    manifoldIntegrator->GetGraphNode(i)->SetValue(0.0, 2);
    float labval = labels->GetTuple1(i);
    manifoldIntegrator->GetGraphNode(i)->SetValue(labval, 3);
    manifoldIntegrator->GetGraphNode(i)->SetPredecessor(nullptr);
  }
  manifoldIntegrator->EmptyQ();
  manifoldIntegrator->SetSearchFinished(false);
  manifoldIntegrator->SetSource(this->m_RootNode);
  manifoldIntegrator->InitializeQueue();
  manifoldIntegrator->SetParamWhileSearching(this->m_ParamWhileSearching);
  manifoldIntegrator->SetWeights(this->m_MaxCost, this->m_DistanceCostWeight, this->m_LabelCostWeight);
  manifoldIntegrator->PrintWeights();
  manifoldIntegrator->FindPath();

  /** at this point, we should have extracted the disc
      now we want to find a boundary point and an instant
      parameterization via FEM */
  this->SetDiscBoundaryList(manifoldIntegrator->m_BoundaryList);
  //  this->m_DiscBoundaryList.assign(manifoldIntegrator->m_BoundaryList.begin(),
  //                  manifoldIntegrator->m_BoundaryList.end());

  this->m_HelpFindLoop.clear();
  this->m_HelpFindLoop.resize(gsz, 0);
  this->m_HelpFindLoop[j_in] = -1;
  unsigned int           furthestnode = 0;
  float                  maxdist = 0;
  GraphSearchNodePointer farnode1 = nullptr;
  for (unsigned int j = 0; j < this->m_DiscBoundaryList.size(); j++)
  {
    if (this->m_DiscBoundaryList[j]->GetTotalCost() > maxdist)
    {
      maxdist = this->m_DiscBoundaryList[j]->GetTotalCost();
      furthestnode = this->m_DiscBoundaryList[j]->GetIdentity();
      farnode1 = this->m_DiscBoundaryList[j];
    }
    m_HelpFindLoop[this->m_DiscBoundaryList[j]->GetIdentity()] = j + 1;
  }
  if (this->m_ParamWhileSearching)
  {
    return 2;
  }
  // butt
  ManifoldIntegratorTypePointer discParameterizer = ManifoldIntegratorType::New();
  discParameterizer->SetSurfaceMesh(this->m_SurfaceMesh);
  discParameterizer->InitializeGraph3();
  discParameterizer->SetMaxCost(1.e9);
  std::cout << " dcz " << discParameterizer->GetGraphSize() << std::endl;
  for (int i = 0; i < discParameterizer->GetGraphSize(); i++)
  {
    discParameterizer->GetGraphNode(i)->SetTotalCost(vnl_huge_val(discParameterizer->GetMaxCost()));
    discParameterizer->GetGraphNode(i)->SetUnVisitable();
    if (this->m_HelpFindLoop[i] > 0)
    {
      discParameterizer->GetGraphNode(i)->SetUnVisited();
    }
    discParameterizer->GetGraphNode(i)->SetValue(0.0, 1);
    discParameterizer->GetGraphNode(i)->SetValue(0.0, 2);
    discParameterizer->GetGraphNode(i)->SetPredecessor(nullptr);
  }
  discParameterizer->EmptyQ();
  discParameterizer->SetSearchFinished(false);
  discParameterizer->SetSource(discParameterizer->GetGraphNode(furthestnode));
  discParameterizer->InitializeQueue();
  discParameterizer->SetWeights(1.e9, 1, 0);
  discParameterizer->FindPath();
  float                  temp = 0;
  GraphSearchNodePointer farnode2 = nullptr;
  unsigned int           vct = 0;
  for (int i = 0; i < discParameterizer->GetGraphSize(); i++)
  {
    if (discParameterizer->GetGraphNode(i)->WasVisited())
    {
      float t = discParameterizer->GetGraphNode(i)->GetTotalCost();
      if (t > temp)
      {
        temp = t;
        farnode2 = discParameterizer->GetGraphNode(i);
      }
      vct++;
      // std::cout << " dist " << t << " vct " << vct << std::endl;
    }
  }
  discParameterizer->BackTrack(farnode2);
  // butt2
  ManifoldIntegratorTypePointer lastlooper = ManifoldIntegratorType::New();
  lastlooper->SetSurfaceMesh(this->m_SurfaceMesh);
  lastlooper->InitializeGraph3();
  lastlooper->SetMaxCost(1.e9);
  //     std::cout << " dcz "<< lastlooper->GetGraphSize() << std::endl;
  for (int i = 0; i < lastlooper->GetGraphSize(); i++)
  {
    lastlooper->GetGraphNode(i)->SetTotalCost(vnl_huge_val(lastlooper->GetMaxCost()));
    lastlooper->GetGraphNode(i)->SetUnVisitable();
    if (this->m_HelpFindLoop[i] > 0)
    {
      lastlooper->GetGraphNode(i)->SetUnVisited();
    }
    lastlooper->GetGraphNode(i)->SetValue(0.0, 1);
    lastlooper->GetGraphNode(i)->SetValue(0.0, 2);
    lastlooper->GetGraphNode(i)->SetPredecessor(nullptr);
  }
  lastlooper->EmptyQ();
  lastlooper->SetSearchFinished(false);
  for (unsigned int pp = 0; pp < discParameterizer->GetPathSize(); pp++)
  {
    unsigned int id = discParameterizer->GetPathAtIndex(pp)->GetIdentity();
    lastlooper->SetSource(lastlooper->GetGraphNode(id));
  }
  lastlooper->InitializeQueue();
  lastlooper->SetWeights(1.e9, 1, 0);
  lastlooper->FindPath();
  GraphSearchNodePointer farnode3 = nullptr;
  vct = 0;
  temp = 0;
  for (int i = 0; i < lastlooper->GetGraphSize(); i++)
  {
    if (lastlooper->GetGraphNode(i)->WasVisited())
    {
      float t = lastlooper->GetGraphNode(i)->GetTotalCost();
      if (t > temp)
      {
        temp = t;
        farnode3 = lastlooper->GetGraphNode(i);
      }
      vct++;
      //           std::cout << " dist " << t << " vct " << vct << std::endl;
    }
  }
  //     std::exception();
  // finally reuse lastlooper with farnode3 as root ...
  for (int i = 0; i < lastlooper->GetGraphSize(); i++)
  {
    lastlooper->GetGraphNode(i)->SetTotalCost(vnl_huge_val(lastlooper->GetMaxCost()));
    lastlooper->GetGraphNode(i)->SetUnVisitable();
    if (this->m_HelpFindLoop[i] > 0)
    {
      lastlooper->GetGraphNode(i)->SetUnVisited();
    }
    lastlooper->GetGraphNode(i)->SetValue(0.0, 1);
    lastlooper->GetGraphNode(i)->SetValue(0.0, 2);
    lastlooper->GetGraphNode(i)->SetPredecessor(nullptr);
  }
  lastlooper->EmptyQ();
  lastlooper->SetSearchFinished(false);
  lastlooper->SetSource(lastlooper->GetGraphNode(farnode3->GetIdentity()));
  lastlooper->InitializeQueue();
  lastlooper->SetWeights(1.e9, 1, 0);
  lastlooper->FindPath();
  // now assemble the parameterization...
  //
  this->m_DiscBoundaryList.clear();
  // add both the above to the list
  //     std::cout << " part 1 " << std::endl;
  for (unsigned int i = 0; i < discParameterizer->GetPathSize(); i++)
  {
    unsigned int id = discParameterizer->GetPathAtIndex(i)->GetIdentity();
    //  std::cout << manifoldIntegrator->GetGraphNode(id)->GetLocation() << std::endl;
    this->m_DiscBoundaryList.push_back(manifoldIntegrator->GetGraphNode(id));
  }
  //     std::cout << " part 2 " << std::endl;
  lastlooper->BackTrack(lastlooper->GetGraphNode(farnode1->GetIdentity()));
  for (unsigned int i = 0; i < lastlooper->GetPathSize(); i++)
  {
    unsigned int id = lastlooper->GetPathAtIndex(i)->GetIdentity();
    //         std::cout << manifoldIntegrator->GetGraphNode(id)->GetLocation() << std::endl;
    this->m_DiscBoundaryList.push_back(manifoldIntegrator->GetGraphNode(id));
  }

  //     std::cout << farnode1->GetLocation() << std::endl;
  // std::cout << farnode2->GetLocation() << std::endl;
  //     std::cout << farnode3->GetLocation() << std::endl;
  //     std::cout << " another idea --- get two points far apart  then solve a minimization problem across the
  // graph
  // that gives the average value in 0 => 1 ... " << std::endl;
  // std::cout << " path 1 sz " << discParameterizer->GetPathSize() << std::endl;

  // finally do but add in reverse order
  // std::cout << " part 3 " <<farnode2->GetIdentity() << " and " << farnode1->GetIdentity() << std::endl;
  lastlooper->EmptyPath();
  lastlooper->BackTrack(lastlooper->GetGraphNode(farnode2->GetIdentity()));
  for (unsigned int i = lastlooper->GetPathSize() - 1; i > 0; i--)
  {
    unsigned int id = lastlooper->GetPathAtIndex(i)->GetIdentity();
    //         std::cout << manifoldIntegrator->GetGraphNode(id)->GetLocation() << std::endl;
    this->m_DiscBoundaryList.push_back(manifoldIntegrator->GetGraphNode(id));
  }
  std::cout << " Almost ... " << std::endl;

  this->m_HelpFindLoop.clear();
  this->m_HelpFindLoop.resize(gsz, 0);
  for (unsigned int j = 0; j < this->m_DiscBoundaryList.size(); j++)
  {
    m_HelpFindLoop[this->m_DiscBoundaryList[j]->GetIdentity()] = j + 1;
  }

  std::cout << " Achievement!! " << std::endl;

  return 2;
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMDiscConformalMap<TSurface, TImage, TDimension>::LocateAndParameterizeDiscBoundary(unsigned int label, bool CheckCost)
{
  float effectivemaxcost = 0;

  this->m_DiscBoundaryList.clear();
  unsigned int      gsz = manifoldIntegrator->GetGraphSize();
  std::vector<bool> alreadyfound(gsz, false);
  this->m_DiscBoundarySorter.clear();
  this->m_DiscBoundarySorter.resize(gsz, 0);
  for (unsigned int j = 0; j < gsz; j++)
  {
    if (manifoldIntegrator->GetGraphNode(j))
    {
      float cost = manifoldIntegrator->GetGraphNode(j)->GetTotalCost();
      if (!CheckCost)
      {
        cost = 0;
      }
      if (fabs(manifoldIntegrator->GetGraphNode(j)->GetValue(3) - label) < 0.5 && cost <= this->m_MaxCost)
      {
        float inb = 0;
        for (unsigned int i = 0; i < manifoldIntegrator->GetGraphNode(j)->m_NumberOfNeighbors; i++)
        {
          if (fabs(manifoldIntegrator->GetGraphNode(j)->m_Neighbors[i]->GetValue(3) - label) > 0.5)
          {
            // CurrentNode is in the boundary
            inb = 1;
          }
        }            // neighborhood
        if (inb > 0) //      std::cout <<  " Node is in boundary " << std::endl;
        {
          inb = 0;
          for (unsigned int i = 0; i < manifoldIntegrator->GetGraphNode(j)->m_NumberOfNeighbors; i++)
          {
            if (fabs(manifoldIntegrator->GetGraphNode(j)->m_Neighbors[i]->GetValue(3) - label) < 0.5 &&
                cost <= this->m_MaxCost)
            {
              // CurrentNode is in the boundary
              inb += 1;
            }
          }                                         // neighborhood
          if (inb >= 2 && alreadyfound[j] == false) // need at least two neighbors with same label
          {
            alreadyfound[j] = true;
            this->m_DiscBoundaryList.push_back(manifoldIntegrator->GetGraphNode(j));
            if (cost > effectivemaxcost)
            {
              effectivemaxcost = cost;
            }
          }
        }
      } // less than max cost
    }   // if node exists
  }     // gsz
  std::cout << " Boundary has " << this->m_DiscBoundaryList.size() << " elements with eff. max cost "
            << effectivemaxcost << std::endl;
  if (CheckCost)
  {
    // very inefficient way to parameterize boundary ....
    unsigned int      bsz = this->m_DiscBoundaryList.size();
    std::vector<bool> alreadyfound(bsz, false);
    this->m_DiscBoundarySorter.resize(bsz, 0);
    unsigned int boundcount = 0;
    unsigned int rootind = 0, lastroot = 0;
    this->m_DiscBoundarySorter[rootind] = boundcount;
    alreadyfound[rootind] = true;
    bool paramdone = false;
    while (!paramdone)
    {
      std::cout << " start param " << bsz << std::endl;
      if (bsz == 0)
      {
        std::exception();
      }
      for (unsigned int myi = 0; myi < bsz; myi++)
      {
        if (this->m_DiscBoundaryList[myi] != this->m_DiscBoundaryList[rootind])
        {
          std::cout << " myi " << myi << " root " << rootind << " bc " << boundcount << std::endl;
          for (unsigned int n = 0; n < this->m_DiscBoundaryList[rootind]->m_NumberOfNeighbors; n++)
          {
            if (this->m_DiscBoundaryList[myi] == this->m_DiscBoundaryList[rootind]->m_Neighbors[n] &&
                !alreadyfound[myi]) //
                                    // its
                                    // in
                                    // the
                                    // bndry
            {                       // check that it's not an isolated bastard
              bool oknode = true;
              if (oknode)
              {
                boundcount++;
                alreadyfound[myi] = true;
                this->m_DiscBoundarySorter[myi] = boundcount;
                n = this->m_DiscBoundaryList[rootind]->m_NumberOfNeighbors + 1;
                std::cout << " cur " << this->m_DiscBoundaryList[rootind]->GetLocation() << " next "
                          << this->m_DiscBoundaryList[myi]->GetLocation() << " boundcount " << boundcount << " of "
                          << bsz << "  curroot " << rootind << std::endl;
                lastroot = rootind;
                rootind = myi;
                myi = bsz;
              }
            } // is in boundary
          }   // neighborhood loop
        }     // not currootnode
        if (boundcount >= bsz - 1)
        {
          paramdone = true;
        }
        else if (myi == (bsz - 1))
        {
          boundcount--;
          rootind = lastroot;
          this->m_DiscBoundarySorter[rootind] = -1;
          std::cout << " failure " << std::endl;
          std::exception();
        }
      } // all boundary nodes
    }   // while

    //    std::exception();
  } //  param boundary if
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMDiscConformalMap<TSurface, TImage, TDimension>::ExtractSurfaceDisc(unsigned int label)
{
  std::cout << " set surface mesh " << std::endl;

  manifoldIntegrator->SetSurfaceMesh(this->m_SurfaceMesh);
  std::cout << " begin initializing graph " << std::endl;
  manifoldIntegrator->InitializeGraph3();
  this->m_SurfaceMesh = manifoldIntegrator->GetSurfaceMesh();
  //   float frac=0;
  //  IndexType index;
  if (this->m_Label_to_Flatten == 0)
  {
    std::cout << " enter LabelToExtract ";
    std::cin >> m_Label_to_Flatten;
    float mc = 0;
    std::cout << " Enter max cost ";
    std::cin >> mc;
    m_MaxCost = mc;
  }
  manifoldIntegrator->SetMaxCost(this->m_MaxCost);
  this->FindMeanSourceInLabel(this->m_Label_to_Flatten);
  //  this->LocateAndParameterizeDiscBoundary(  m_Label_to_Flatten , false );
  //   this->LocateAndParameterizeDiscBoundary(  m_Label_to_Flatten , true );
  // std::cout << " findpath in extractsurfacedisk done ";

  // assign scalars to the original surface mesh
  //  typedef itk::SurfaceMeshCurvature<GraphSearchNodeType,GraphSearchNodeType> surfktype;
  //  typename surfktype::Pointer surfk=surfktype::New();

  vtkPoints *     vtkpoints = this->m_SurfaceMesh->GetPoints();
  int             numPoints = vtkpoints->GetNumberOfPoints();
  vtkFloatArray * param = vtkFloatArray::New();
  param->SetName("angle");
  for (int i = 0; i < numPoints; i++)
  {
    float temp = fabs(manifoldIntegrator->GetGraphNode(i)->GetTotalCost());
    if (temp > m_MaxCost)
    {
      temp = m_MaxCost;
    }
    param->InsertNextValue(temp * 255. / m_MaxCost);
  }

  std::cout << " extractsurfacedisk done ";
  //  m_SurfaceMesh->GetPointData()->SetScalars(param);
}

template <typename TSurface, typename TImage, unsigned int TDimension>
bool
FEMDiscConformalMap<TSurface, TImage, TDimension>::GenerateSystemFromSurfaceMesh()
{
  if (!this->m_SurfaceMesh)
  {
    std::cout << " NO MESH";
    return false;
  }
  std::cout << " Generate system from surface mesh " << std::endl;
  m_Smooth = m_Sigma;
  // create a material
  // Choose the material properties
  typename MaterialType::Pointer m;
  m = MaterialType::New();
  m->GN = 0;  // Global number of the material ///
  m->E = 4.0; // Young modulus -- used in the membrane ///
  m->A = 0.0; // Crossection area ///
  m->h = 0.0; // Crossection area ///
  m->I = 0.0; // Moment of inertia ///
  m->nu = 0.; // .0;    // poissons -- DONT CHOOSE 1.0!!///
  m->RhoC = 0.0;

  // Create the element type
  ElementType::Pointer e1 = ElementType::New();
  e1->m_mat = dynamic_cast<MaterialType *>(m);

  vtkPoints * vtkpoints = this->m_SurfaceMesh->GetPoints();
  int         numPoints = vtkpoints->GetNumberOfPoints();
  int         foundnum = 0;
  int         boundsz = 0;
  for (int i = 0; i < numPoints; i++)
  {
    double *                   pt = vtkpoints->GetPoint(i);
    typename NodeType::Pointer n;
    n = new NodeType(pt[0], pt[1], pt[2]);
    if (this->InDisc(manifoldIntegrator->GetGraphNode(i)))
    {
      n->GN = i;
      m_Solver.node.push_back(itk::fem::FEMP<NodeType>(n));
      foundnum++;
    }
    if (this->InBorder(manifoldIntegrator->GetGraphNode(i)))
    {
      boundsz++;
    }
  }

  typename NodeType::Pointer * narr = new NodeType::Pointer[numPoints];
  for (int i = 0; i < numPoints; i++)
  {
    if (this->InDisc(manifoldIntegrator->GetGraphNode(i)) || this->InBorder(manifoldIntegrator->GetGraphNode(i)))
    {
      narr[i] = m_Solver.node.Find(i);
    }
  }

  std::cout << " Found " << foundnum << " nodes " << std::endl;
  std::cout << " bound " << boundsz << std::endl;
  vtkCellArray * vtkcells = this->m_SurfaceMesh->GetPolys();

  vtkIdType     npts;
  vtkIdType *   pts;
  unsigned long i = 0;
  //   unsigned long toti = vtkcells->GetNumberOfCells();
  //   unsigned long rate = toti/50;
  //  std::cout << " progress ";
  for (vtkcells->InitTraversal(); vtkcells->GetNextCell(npts, pts);)
  {
    //  if ( i % rate == 0 && i > rate ) std::cout << "  " <<  (float) i / (float) toti << " ";
    // turn the cell into an element
    //      std::cout << " points ids a " << pts[0] << " b " << pts[1] << " c " << pts[2] << std::endl;
    bool                 eltok = true;
    ElementType::Pointer e;
    e = dynamic_cast<ElementType *>(e1->Clone());

    if (this->InDisc(manifoldIntegrator->GetGraphNode(pts[0])))
    {
      e->SetNode(2, narr[pts[0]]);
    }
    //          e->SetNode(2,m_Solver.node.Find( pts[0] ));
    else
    {
      eltok = false;
    }

    if (this->InDisc(manifoldIntegrator->GetGraphNode(pts[1])))
    {
      e->SetNode(1, narr[pts[1]]);
    }
    //          e->SetNode(1,m_Solver.node.Find( pts[1] ));
    else
    {
      eltok = false;
    }

    if (this->InDisc(manifoldIntegrator->GetGraphNode(pts[2])))
    {
      e->SetNode(0, narr[pts[2]]);
    }
    //          e->SetNode(0,m_Solver.node.Find( pts[2] ));
    else
    {
      eltok = false;
    }

    if (eltok)
    {
      e->GN = i;
      m_Solver.el.push_back(itk::fem::FEMP<itk::fem::Element>(e));
      i++;
    } // else std::cout <<" cannot find elt " << std::endl;
  }

  std::cout << " DONE: NUMBER OF CELLS " << i << std::endl;

  return true;
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMDiscConformalMap<TSurface, TImage, TDimension>::FixBoundaryPoints(unsigned int option)
{
  itk::fem::Element::ArrayType::iterator elt = m_Solver.el.begin();
  unsigned int                           dofs = (*elt)->GetNumberOfDegreesOfFreedomPerNode();

  int fixct = 0;
  int eltct = 0;
  while (elt != m_Solver.el.end())
  {
    for (unsigned int i = 0; i < dofs; i++)
    {
      int nodeid = (*elt)->GetNode(i)->GN;
      //      if( manifoldIntegrator->GetGraphNode(nodeid)->GetTotalCost() > 222 )
      if (this->InBorder(manifoldIntegrator->GetGraphNode(nodeid)))
      {
        itk::fem::LoadBC::Pointer l1;
        l1 = itk::fem::LoadBC::New();
        l1->m_element = (*elt);
        l1->m_dof = i;
        l1->m_value = vnl_vector<double>(1, 0.0); // for exp
        float fixvalue = 0.0;
        if (this->m_MapToCircle)
        {
          fixvalue = this->GetBoundaryParameterForCircle(nodeid, option);
        }
        else
        {
          fixvalue = this->GetBoundaryParameterForSquare(nodeid, option);
        }
        l1->m_value = vnl_vector<double>(1, fixvalue); // for direct rad
        m_Solver.load.push_back(itk::fem::FEMP<itk::fem::Load>(&*l1));
        /*
    itk::fem::LoadNode::Pointer ln1=itk::fem::LoadNode::New();
    ln1->m_pt=0;
    ln1->F.set_size(1);
    ln1->F.fill(fixvalue);
    ln1->m_element=(*elt);
    m_Solver.load.push_back( itk::fem::FEMP<itk::fem::Load>(&*ln1) );

    itk::fem::LoadNode::Pointer ln2=itk::fem::LoadNode::New();
    ln2->m_pt=1;
    ln2->F.set_size(1);
    ln2->F.fill(fixvalue);
    ln2->m_element=(*elt);
    m_Solver.load.push_back( itk::fem::FEMP<itk::fem::Load>(&*ln2) );

    itk::fem::LoadNode::Pointer ln3=itk::fem::LoadNode::New();
    ln3->m_pt=2;
    ln3->F.set_size(1);
    ln3->F.fill(fixvalue);
    ln3->m_element=(*elt);
    m_Solver.load.push_back( itk::fem::FEMP<itk::fem::Load>(&*ln3) );
        */
        if (i == 0)
        {
          fixct++;
        }
      }
    }
    ++el;
    eltct++;
  }

  std::cout << " Fixed elt number " << fixct << " of " << eltct << std::endl;
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMDiscConformalMap<TSurface, TImage, TDimension>::ApplyRealForces()
{
  /*
  itk::fem::Element::Pointer e = m_Solver.el[m_PoleElementsGN[0]];
    // load node 0
    {
      itk::fem::LoadNode::Pointer ln1=itk::fem::LoadNode::New();
      ln1->m_pt=0;
      ln1->F.set_size(1);
      ln1->F.fill(-4.0*m_Pi);
      ln1->m_element=(e);
      m_Solver.load.push_back( itk::fem::FEMP<itk::fem::Load>(&*ln1) );

      itk::fem::LoadNode::Pointer ln2=itk::fem::LoadNode::New();
      ln2->m_pt=1;
      ln2->F.set_size(1);
      ln2->F.fill(-4.0*m_Pi);
      ln2->m_element=(e);
      m_Solver.load.push_back( itk::fem::FEMP<itk::fem::Load>(&*ln2) );

      itk::fem::LoadNode::Pointer ln3=itk::fem::LoadNode::New();
      ln3->m_pt=2;
      ln3->F.set_size(1);
      ln3->F.fill(-4.0*m_Pi);
      ln3->m_element=(e);
      m_Solver.load.push_back( itk::fem::FEMP<itk::fem::Load>(&*ln3) );


    }
  */
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMDiscConformalMap<TSurface, TImage, TDimension>::ApplyImaginaryForces()
{
  itk::fem::Element::Pointer e = m_Solver.el[m_PoleElementsGN[0]];

  itk::fem::LoadBC::Pointer l1;

  l1 = itk::fem::LoadBC::New();
  l1->m_element = (e);
  l1->m_dof = 0;
  l1->m_value = vnl_vector<double>(1, -1. / 3.); // for exp
  m_Solver.load.push_back(itk::fem::FEMP<itk::fem::Load>(&*l1));

  itk::fem::LoadBC::Pointer l2;
  l2 = itk::fem::LoadBC::New();
  l2->m_element = (e);
  l2->m_dof = 1;
  l2->m_value = vnl_vector<double>(1, -1. / 3.); // for exp
  m_Solver.load.push_back(itk::fem::FEMP<itk::fem::Load>(&*l2));

  itk::fem::LoadBC::Pointer l3;
  l3 = itk::fem::LoadBC::New();
  l3->m_element = (e);
  l3->m_dof = 2;
  l3->m_value = vnl_vector<double>(1, 2. / 3.); // for exp
  m_Solver.load.push_back(itk::fem::FEMP<itk::fem::Load>(&*l3));
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMDiscConformalMap<TSurface, TImage, TDimension>::MakeFlatImage()
{
  // first declare the flat image
  typename FlatImageType::RegionType region;
  typename FlatImageType::SizeType   size;
  int                                sz = 256;
  size[0] = sz;
  size[1] = sz;
  region.SetSize(size);
  m_FlatImage = FlatImageType::New();
  m_FlatImage->SetRegions(region);
  m_FlatImage->AllocateInitialized();
  typename FlatImageType::IndexType index;

  std::cout << " Making flat image " << std::endl;
  int maxits = 100;
  for (int its = 0; its <= maxits; its++)
  {
    for (::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n != m_Solver.node.end(); ++n)
    {
      float temp = 255.0 - manifoldIntegrator->GetGraphNode((*n)->GN)->GetValue(3); // curvature
      //      float temp=255.0*manifoldIntegrator->GetGraphNode((*n)->GN)->GetValue(2); // extrinsic dist
      index[0] =
        (long int)(0.5 + (1.0 + manifoldIntegrator->GetGraphNode((*n)->GN)->GetValue(0)) * (float)(sz - 1) / 2.);
      index[1] =
        (long int)(0.5 + (1.0 + manifoldIntegrator->GetGraphNode((*n)->GN)->GetValue(1)) * (float)(sz - 1) / 2.);
      // std::cout << " ind " << index << std::endl;
      m_FlatImage->SetPixel(index, temp);
    }
    typedef itk::DiscreteGaussianImageFilter<FlatImageType, FlatImageType> dgf;
    typename dgf::Pointer                                                  filter = dgf::New();
    filter->SetVariance(1.0);
    filter->SetUseImageSpacing(false);
    filter->SetMaximumError(.01f);
    filter->SetInput(m_FlatImage);
    filter->Update();
    m_FlatImage = filter->GetOutput();

    if (its < maxits)
    {
      int                                              center = (int)sz / 2;
      itk::ImageRegionIteratorWithIndex<FlatImageType> it(m_FlatImage, region);
      it.GoToBegin();
      typename FlatImageType::IndexType index;
      while (!it.IsAtEnd())
      {
        index = it.GetIndex();
        float x = (float)index[0] - (float)sz / 2.0;
        float y = (float)index[1] - (float)sz / 2.0;
        float dist = sqrt(x * x + y * y);
        //        std::cout << "center " << center <<  " index " << index << "dist " << dist ;
        if (dist > center)
        {
          it.Set(0.0);
        }
        ++it;
      }
    }
  }
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMDiscConformalMap<TSurface, TImage, TDimension>::BuildOutputMeshes(float tval)
{
  std::cout << " build output mesh " << std::endl;

  typedef GraphSearchNodeType::NodeLocationType loctype;
  // Get the number of points in the mesh
  int numPoints = m_Solver.node.size();

  vtkDataArray * labels = this->m_SurfaceMesh->GetPointData()->GetArray("Label");
  vtkDataArray * features = nullptr;
  if (this->m_SurfaceFeatureMesh)
  {
    if (this->m_SurfaceFeatureMesh->GetPointData()->GetArray("Feature"))
    {
      features = this->m_SurfaceFeatureMesh->GetPointData()->GetArray("Feature");
    }
    else
    {
      features = labels;
    }
  }
  else
  {
    features = labels;
  }

  // Create vtk polydata
  vtkPolyData * polydata1 = vtkPolyData::New();
  vtkPolyData * polydata2 = vtkPolyData::New();

  // Create the vtkPoints object and set the number of points
  vtkPoints * vpoints1 = vtkPoints::New();
  vtkPoints * vpoints2 = vtkPoints::New();
  vpoints1->SetNumberOfPoints(numPoints);
  vpoints2->SetNumberOfPoints(numPoints);

  std::cout << " start pts ";
  int idx = 0;

  vtkFloatArray * param = vtkFloatArray::New();
  param->SetName("feature");

  vtkFloatArray * paramAngle = vtkFloatArray::New();
  paramAngle->SetName("angle");

  vtkFloatArray * paramDistance = vtkFloatArray::New();
  paramDistance->SetName("distance");

  vtkIdTypeArray * paramPoints = vtkIdTypeArray::New();
  paramPoints->SetName("points");
  for (::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n != m_Solver.node.end(); ++n)
  {
    loctype loc = manifoldIntegrator->GetGraphNode((*n)->GN)->GetLocation();
    float   pt1[3];
    pt1[0] = loc[0];
    pt1[1] = loc[1];
    pt1[2] = loc[2];
    float pt2[3];
    pt2[0] = manifoldIntegrator->GetGraphNode((*n)->GN)->GetValue(0) * 100. * (1.0 - tval) + tval * loc[0];
    pt2[1] = manifoldIntegrator->GetGraphNode((*n)->GN)->GetValue(1) * 100. * (1.0 - tval) + tval * loc[1];
    pt2[2] = 1. * (1.0 - tval) + tval * loc[2];
    vpoints1->SetPoint(idx, pt1);
    vpoints2->SetPoint(idx, pt2);

    float temp = features->GetTuple1((*n)->GN);
    float temp2 = manifoldIntegrator->GetGraphNode((*n)->GN)->GetValue(0) * 255; // for length
    param->InsertNextValue(temp);
    paramDistance->InsertNextValue(temp2);
    paramPoints->InsertNextValue((*n)->GN);
    paramAngle->InsertNextValue(temp); // curvature
    //    param->InsertNextValue(temp*255./m_MaxCost);

    (*n)->GN = idx;
    idx++;
  }

  std::cout << " done with pts " << std::endl;
  vtkCellArray * tris1 = vtkCellArray::New();
  vtkCellArray * tris2 = vtkCellArray::New();

  std::cout << " start with tris " << std::endl;
  for (::itk::fem::Solver::ElementArray::iterator n = m_Solver.el.begin(); n != m_Solver.el.end(); ++n)
  {
    tris1->InsertNextCell(3);
    tris2->InsertNextCell(3);
    for (unsigned int i = 0; i < (*n)->GetNumberOfNodes(); i++)
    {
      tris1->InsertCellPoint((*n)->GetNode(i)->GN);
      tris2->InsertCellPoint((*n)->GetNode(i)->GN);
    }
  }
  std::cout << " done with tris " << std::endl;
  // Assign points and cells
  polydata1->SetPoints(vpoints1);
  polydata2->SetPoints(vpoints2);
  polydata1->SetPolys(tris1);
  polydata2->SetPolys(tris2);
  polydata1->GetPointData()->SetScalars(param);
  polydata2->GetPointData()->SetScalars(param);
  polydata1->GetPointData()->AddArray(paramAngle);
  polydata2->GetPointData()->AddArray(paramAngle);
  polydata1->GetPointData()->AddArray(paramDistance);
  polydata2->GetPointData()->AddArray(paramDistance);
  polydata1->GetPointData()->AddArray(paramPoints);
  polydata2->GetPointData()->AddArray(paramPoints);

  m_ExtractedSurfaceMesh = polydata1; //

  vtkDelaunay2D * delny2 = vtkDelaunay2D::New();
  delny2->SetInput(polydata2);
  m_DiskSurfaceMesh = delny2->GetOutput();
  // m_DiskSurfaceMesh=polydata2;

  return;
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMDiscConformalMap<TSurface, TImage, TDimension>::ConformalMap()
{
  m_Solver.load.clear();
  m_Solver.node.clear();
  m_Solver.el.clear();

  /**
   * Open the file and assign it to stream object f
   */
  if (m_ReadFromFile)
  {
    const char * filename = m_ParameterFileName.c_str();
    std::cout << "Reading FEM problem from file: " << std::string(filename) << "\n";
    std::ifstream f;
    f.open(filename);
    if (!f)
    {
      std::cout << "File " << filename << " not found!\n";
      return;
    }

    try
    {
      m_Solver.Read(f);
    }
    catch (const ::itk::fem::FEMException & e)
    {
      std::cout << "Error reading FEM problem: " << filename << "!\n";
      e.Print(std::cout);
      return;
    }

    f.close();
  }
  else if (this->GenerateSystemFromSurfaceMesh())
  {
    ;
  }
  else
  {
    return;
  }

  /**
   * Assign a unique id (global freedom number - GFN)
   * to every degree of freedom (DOF) in a system.
   */
  m_Solver.GenerateGFN();
  m_ImagSolution.set_size(m_Solver.GetNumberOfDegreesOfFreedom());
  m_RealSolution.set_size(m_Solver.GetNumberOfDegreesOfFreedom());
  m_Radius.set_size(m_Solver.GetNumberOfDegreesOfFreedom());
  m_RealSolution.fill(0);
  m_ImagSolution.fill(0);
  m_Radius.fill(0);

  m_Debug = false;
  if (m_Debug)
  {
    for (::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n != m_Solver.node.end(); ++n)
    {
      std::cout << "Node#: " << (*n)->GN << ": ";
      std::cout << " coord " << (*n)->GetCoordinates() << " coord2 "
                << manifoldIntegrator->GetGraphNode((*n)->GN)->GetLocation() << std::endl;
    }
    for (::itk::fem::Solver::ElementArray::iterator n = m_Solver.el.begin(); n != m_Solver.el.end(); ++n)
    {
      std::cout << "Elt#: " << (*n)->GN << ": has " << (*n)->GetNumberOfNodes() << " nodes ";
      for (unsigned int i = 0; i < (*n)->GetNumberOfNodes(); i++)
      {
        std::cout << " coord " << (*n)->GetNode(i)->GetCoordinates() << std::endl;
      }
    }
  }

  unsigned int maxits = m_Solver.GetNumberOfDegreesOfFreedom(); // should be > twice ndofs
  // if (m_Debug)
  std::cout << " ndof " << maxits << std::endl;
  itpackWrapper.SetMaximumNumberIterations(maxits * 5);
  itpackWrapper.SetTolerance(1.e-4);
  itpackWrapper.SuccessiveOverrelaxation();
  // itpackWrapper.JacobianConjugateGradient();
  itpackWrapper.SetMaximumNonZeroValuesInMatrix(maxits * 50);
  m_Solver.SetLinearSystemWrapper(&itpackWrapper);

  this->FixBoundaryPoints(0);
  m_Solver.AssembleK();
  m_Solver.DecomposeK();
  //  this->ApplyRealForces();
  std::cout << " appl force ";
  m_Solver.AssembleF();
  //  for (int i=0; i<maxits; i++) if (m_Solver.GetVectorValue(i) != 0) m_Solver.SetVectorValue(i,1.0);
  std::cout << " b solve ";
  m_Solver.Solve();
  std::cout << " e solve ";
  m_Solver.UpdateDisplacements(); // copies solution to nodes
  unsigned long ct = 0;
  for (::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n != m_Solver.node.end(); ++n)
  {
    for (unsigned int d = 0, dof; (dof = (*n)->GetDegreeOfFreedom(d)) != ::itk::fem::Element::InvalidDegreeOfFreedomID;
         d++)
    {
      m_RealSolution[dof] = m_Solver.GetSolution(dof);
      //      if ( ct % 10 == 0) std::cout << " mrdof " <<  m_RealSolution[dof]  << " dof " << dof <<
      // std::endl;
    }
    ct++;
  }

  this->ConformalMap2();
  // this->ConformalMap3();
  this->ConjugateHarmonic();
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMDiscConformalMap<TSurface, TImage, TDimension>::ConformalMap2()
{
  m_Solver.load.clear();
  this->FixBoundaryPoints(1);
  m_Solver.AssembleK(); // need to reassemble b/c LoadBC's affect K
  m_Solver.AssembleF();
  m_Solver.Solve();
  m_Solver.UpdateDisplacements(); // copies solution to nodes
  unsigned long ct = 0;
  for (::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n != m_Solver.node.end(); ++n)
  {
    for (unsigned int d = 0, dof; (dof = (*n)->GetDegreeOfFreedom(d)) != ::itk::fem::Element::InvalidDegreeOfFreedomID;
         d++)
    {
      m_ImagSolution[dof] = m_Solver.GetSolution(dof);
      //      if (ct % 10 == 0) std::cout << " midof " <<  m_ImagSolution[dof]  <<  " dof " << dof <<
      // std::endl;
    }
    ct++;
  }
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMDiscConformalMap<TSurface, TImage, TDimension>::MeasureLengthDistortion()
{
  // now measure the length distortion of the given solution
  manifoldIntegrator->EmptyQ();
  manifoldIntegrator->SetSearchFinished(false);
  for (int i = 0; i < manifoldIntegrator->GetGraphSize(); i++)
  {
    manifoldIntegrator->GetGraphNode(i)->SetTotalCost(vnl_huge_val(manifoldIntegrator->GetMaxCost()));
    manifoldIntegrator->GetGraphNode(i)->SetUnVisited();
    manifoldIntegrator->GetGraphNode(i)->SetValue(0.0, 1);
    manifoldIntegrator->GetGraphNode(i)->SetValue(0.0, 2);
    manifoldIntegrator->GetGraphNode(i)->SetPredecessor(nullptr);
  }
  manifoldIntegrator->SetSource(manifoldIntegrator->GetGraphNode(this->m_SourceNodeNumber));
  manifoldIntegrator->InitializeQueue();
  float mchere = 1.2;
  manifoldIntegrator->SetMaxCost(mchere);
  manifoldIntegrator->SetWeights(this->m_MaxCost, this->m_DistanceCostWeight, this->m_LabelCostWeight);
  manifoldIntegrator->FindPath();

  // backtrack everywhere to set up forweard tracking
  float maxmanifolddist = 0;
  float distDistortion = 0;

  unsigned int ct = 0;
  for (int i = 0; i < manifoldIntegrator->GetGraphSize(); i++)
  {
    if (manifoldIntegrator->GetGraphNode(i))
    {
      if (manifoldIntegrator->GetGraphNode(i)->GetTotalCost() <= manifoldIntegrator->GetMaxCost())
      {
        float ttt = manifoldIntegrator->GetGraphNode(i)->GetValue(2);
        if (ttt > maxmanifolddist)
        {
          maxmanifolddist = ttt;
        }
        ct++;
      }
    }
  }
  ct = 0;
  for (int i = 0; i < manifoldIntegrator->GetGraphSize(); i++)
  {
    if (manifoldIntegrator->GetGraphNode(i))
    {
      if (manifoldIntegrator->GetGraphNode(i)->GetTotalCost() < manifoldIntegrator->GetMaxCost())
      {
        const float rad = manifoldIntegrator->GetGraphNode(i)->GetValue(0);
        const float manifolddist = manifoldIntegrator->GetGraphNode(i)->GetValue(2) / maxmanifolddist;
        manifoldIntegrator->GetGraphNode(i)->SetValue(manifolddist, 2);
        distDistortion += fabs(rad - manifolddist);
        ct++;
      }
    }
  }
  std::cout << "  distDistortion/ct " << distDistortion / (float)ct << " maxmfd " << maxmanifolddist << std::endl;
  return;

  /*    m_ExtractedSurfaceMesh=polydata1;//
  //   m_DiskSurfaceMesh=delny2->GetOutput();
  float total=0.0;
  for (int i=0; i<this->m_SurfaceMesh->GetPoints()->GetNumberOfPoints(); i++)
  {
  float length1=0;
  float length2=0;
  float pt1[3];

  for (int ne=0; ne<manifoldIntegrator->GetGraphNode(i)->GetNumberOfNeighbors(); ne++)
  {
  }

  }
*/
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMDiscConformalMap<TSurface, TImage, TDimension>::ConjugateHarmonic()
{
  std::cout << " Conformal coordinates " << std::endl;

  unsigned long ct = 0;

  for (::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n != m_Solver.node.end(); ++n)
  {
    ct++;
    unsigned long dof = (*n)->GetDegreeOfFreedom(0);
    if (dof < m_RealSolution.size())
    {
      const float U = m_RealSolution[dof];
      const float V = m_ImagSolution[dof];
      manifoldIntegrator->GetGraphNode((*n)->GN)->SetValue(U, 0);
      manifoldIntegrator->GetGraphNode((*n)->GN)->SetValue(V, 1);
    }
  }

  //  this->MakeFlatImage();
  this->BuildOutputMeshes();

  return;
}
} // namespace itk

/*
vtkSelectPolyData *loop = vtkSelectPolyData::New();
loop->SetInput(this->m_SurfaceMesh);
// put points inside ...
vtkPoints* points = vtkPoints::New();
points->SetNumberOfPoints( this->m_DiscBoundaryList.size());
unsigned int idx=0;
for ( unsigned int j = 0 ; j < this->m_DiscBoundaryList.size(); j ++ )
  {
    float pt1[3];
    typename GraphSearchNodeType::NodeLocationType loc=this->m_DiscBoundaryList[j]->GetLocation();
    pt1[0]=loc[0];
    pt1[1]=loc[1];
    pt1[2]=loc[2];
    //        unsigned int idx=this->m_DiscBoundaryList[j]->GetIdentity();
    points->SetPoint(idx,pt);
    idx++;
  }
loop->GenerateSelectionScalarsOff();
loop->SetSelectionModeToClosestPointRegion(); //negative scalars inside
loop->SetSelectionModeToSmallestRegion(); //negative scalars inside
loop->SetLoop(points);
loop->Modified();
vtkClipPolyData *clip = vtkClipPolyData::New(); //clips out positive region
clip->SetInput(loop->GetOutput());
vtkPolyDataMapper *clipMapper = vtkPolyDataMapper::New();
clipMapper->SetInput(clip->GetOutput());
vtkActor *clipActor = vtkActor::New();
clipActor->SetMapper(clipMapper);
clipActor->AddPosition(1, 0, 0);
//    clipActor->GetProperty()->SetColor(0, 0, 1); //Set colour blue


vtkRenderer *ren1 = vtkRenderer::New();
vtkRenderWindow* renWin = vtkRenderWindow::New();
renWin->AddRenderer(ren1);
vtkRenderWindowInteractor* inter = vtkRenderWindowInteractor::New();
inter->SetRenderWindow(renWin);
ren1->SetViewport(0.0, 0.0, 1.0, 1.0);
ren1->AddActor(clipActor);

renWin->Render();
inter->Start();
ren1->Delete();
renWin->Delete();

points-> Delete();
loop->Delete();
clip->Delete();
clipMapper->Delete();
clipActor->Delete();
*/

#endif
