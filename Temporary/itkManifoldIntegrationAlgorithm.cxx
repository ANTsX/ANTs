/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit (ITK)

=========================================================================*/
#ifndef _itkManifoldIntegrationAlgorithm_cxx_
#define _itkManifoldIntegrationAlgorithm_cxx_

#include "itkSurfaceMeshCurvature.h"
#include "vtkFeatureEdges.h"
#include "vtkPointLocator.h"
#include "vtkCellLocator.h"
#include "vtkTriangleFilter.h"
#include "vtkCleanPolyData.h"
#include "vtkPolyDataConnectivityFilter.h"

namespace itk
{
template <typename TGraphSearchNode>
ManifoldIntegrationAlgorithm<TGraphSearchNode>::ManifoldIntegrationAlgorithm()
{
  m_SurfaceMesh = nullptr;
  m_QS = DijkstrasAlgorithmQueue<TGraphSearchNode>::New();
  m_MaxCost = vnl_huge_val(m_MaxCost);
  m_PureDist = false;
  //  m_LabelCost=0;
  m_ParamWhileSearching = false;
  this->m_DistanceCostWeight = 1;
  this->m_LabelCostWeight = 0;
}

template <typename TGraphSearchNode>
float
ManifoldIntegrationAlgorithm<TGraphSearchNode>::dstarUestimate(typename TGraphSearchNode::Pointer G)
{
  typedef itk::SurfaceMeshCurvature<TGraphSearchNode, TGraphSearchNode> surfktype;
  typename surfktype::Pointer                                           surfk = surfktype::New();
  surfk->SetSurfacePatch(G);
  surfk->FindNeighborhood();
  float dsu = (float)surfk->dstarUestimate();
  return dsu;
}

template <typename TGraphSearchNode>
void
ManifoldIntegrationAlgorithm<TGraphSearchNode>::InitializeGraph3()
{
  if (!m_SurfaceMesh)
  {
    return;
  }

  // Construct simple triangles
  vtkTriangleFilter * fltTriangle = vtkTriangleFilter::New();
  fltTriangle->SetInput(m_SurfaceMesh);

  cout << "   converting mesh to triangles " << endl;
  fltTriangle->Update();

  /* Clean the data
  vtkCleanPolyData* fltCleaner = vtkCleanPolyData::New();
  fltCleaner->SetInput(fltTriangle->GetOutput());
  fltCleaner->SetTolerance(0);
  fltCleaner->ConvertPolysToLinesOn();

  cout << "   cleaning up triangle mesh " << endl;
  fltCleaner->Update();

  // Go through and delete the cells that are of the wrong type
  //m_SurfaceMesh
  vtkPolyData* clean= fltCleaner->GetOutput();
  for(vtkIdType i = clean->GetNumberOfCells();i > 0;i--)
    {
    if(clean->GetCellType(i-1) != VTK_TRIANGLE)
      clean->DeleteCell(i-1);
    }
  clean->BuildCells();
  m_SurfaceMesh=clean;*/
  m_SurfaceMesh = fltTriangle->GetOutput();

  typedef float                  labelType;
  typedef std::vector<labelType> LabelSetType;
  LabelSetType                   myLabelSet;

  vtkPoints *    vtkpoints = m_SurfaceMesh->GetPoints();
  vtkPointData * pd = m_SurfaceMesh->GetPointData();
  int            numPoints = vtkpoints->GetNumberOfPoints();
  vtkDataArray * scs = pd->GetScalars();
  m_GraphX.resize(numPoints);
  for (int i = 0; i < numPoints; i++)
  {
    NodeLocationType                                                       loc;
    double *                                                               pt = vtkpoints->GetPoint(i);
    typename GraphSearchNode<PixelType, CoordRep, GraphDimension>::Pointer G =
      GraphSearchNode<PixelType, CoordRep, GraphDimension>::New();
    G->SetUnVisited();
    G->SetTotalCost(m_MaxCost);
    G->SetValue(scs->GetTuple1(i), 3); /** here we put the label value */

    labelType label = scs->GetTuple1(i);
    if (find(myLabelSet.begin(), myLabelSet.end(), label) == myLabelSet.end())
    {
      myLabelSet.push_back(label);
    }
    //    std::cout << " label " <<      scs->GetTuple1(i) << std::endl;
    // std::cout << " set3 " <<      scs->GetTuple3(i) << std::endl;
    // std::cout << " set4 " <<      scs->GetTuple4(i) << std::endl;
    for (int j = 0; j < GraphDimension; j++)
    {
      loc[j] = pt[j];
    }
    G->SetLocation(loc);
    G->SetPredecessor(nullptr);
    G->m_NumberOfNeighbors = 0;
    G->SetIdentity(i);
    m_GraphX[i] = G;
  }
  std::cout << " you have " << myLabelSet.size() << " labels " << std::endl;
  for (unsigned int i = 0; i < myLabelSet.size(); i++)
  {
    std::cout << " label " << myLabelSet[i] << std::endl;
  }
  std::cout << " allocation of graph done ";

  // now loop through the cells to get triangles and also edges

  vtkCellArray * vtkcells = m_SurfaceMesh->GetPolys();

  vtkIdType   npts;
  vtkIdType * pts;
  /* count possible neighbors ... */
  for (vtkcells->InitTraversal(); vtkcells->GetNextCell(npts, pts);)
  {
    m_GraphX[pts[0]]->m_NumberOfNeighbors += 2;
    m_GraphX[pts[1]]->m_NumberOfNeighbors += 2;
    m_GraphX[pts[2]]->m_NumberOfNeighbors += 2;
  }
  for (int i = 0; i < numPoints; i++)
  {
    m_GraphX[i]->m_Neighbors.resize(m_GraphX[i]->m_NumberOfNeighbors);
    //    std::cout <<" Num Neigh " << i << " is " << m_GraphX[i]->m_NumberOfNeighbors << std::endl;
    m_GraphX[i]->m_NumberOfNeighbors = 0;
  }
  for (vtkcells->InitTraversal(); vtkcells->GetNextCell(npts, pts);)
  {
    m_GraphX[pts[0]]->m_Neighbors[m_GraphX[pts[0]]->m_NumberOfNeighbors] = m_GraphX[pts[1]];
    m_GraphX[pts[0]]->m_NumberOfNeighbors++;
    m_GraphX[pts[0]]->m_Neighbors[m_GraphX[pts[0]]->m_NumberOfNeighbors] = m_GraphX[pts[2]];
    m_GraphX[pts[0]]->m_NumberOfNeighbors++;

    m_GraphX[pts[1]]->m_Neighbors[m_GraphX[pts[1]]->m_NumberOfNeighbors] = m_GraphX[pts[0]];
    m_GraphX[pts[1]]->m_NumberOfNeighbors++;
    m_GraphX[pts[1]]->m_Neighbors[m_GraphX[pts[1]]->m_NumberOfNeighbors] = m_GraphX[pts[2]];
    m_GraphX[pts[1]]->m_NumberOfNeighbors++;

    m_GraphX[pts[2]]->m_Neighbors[m_GraphX[pts[2]]->m_NumberOfNeighbors] = m_GraphX[pts[0]];
    m_GraphX[pts[2]]->m_NumberOfNeighbors++;
    m_GraphX[pts[2]]->m_Neighbors[m_GraphX[pts[2]]->m_NumberOfNeighbors] = m_GraphX[pts[1]];
    m_GraphX[pts[2]]->m_NumberOfNeighbors++;
  }
  // now go through each node and make its list of neighbors unique
  // had to do this b/c it's easier than fixing the junk above ... too tired!
  for (int i = 0; i < numPoints; i++)
  {
    std::vector<unsigned int> neighlist;
    for (unsigned int n = 0; n < m_GraphX[i]->m_NumberOfNeighbors; n++)
    {
      neighlist.push_back(m_GraphX[i]->m_Neighbors[n]->GetIdentity());
    }
    std::sort(neighlist.begin(), neighlist.end());
    std::vector<unsigned int>::iterator new_end_pos;
    new_end_pos = std::unique(neighlist.begin(), neighlist.end());
    neighlist.erase(new_end_pos, neighlist.end());
    //      std::cout << " new leng " << neighlist.size() << " old " << len1 << std::endl;
  }
}

template <typename TGraphSearchNode>
void
ManifoldIntegrationAlgorithm<TGraphSearchNode>::InitializeGraph2()
{
  if (!m_SurfaceMesh)
  {
    return;
  }
  /*
    // Construct simple triangles
    vtkTriangleFilter* fltTriangle = vtkTriangleFilter::New();
    fltTriangle->SetInput(m_SurfaceMesh);

    cout << "   converting mesh to triangles " << endl;
    fltTriangle->Update();

    cout << "  this mesh has " << fltTriangle->GetOutput()->GetNumberOfPoints() << " points" << endl;
    cout << "  this mesh has " << fltTriangle->GetOutput()->GetNumberOfCells() << " cells" << endl;

    // Clean the data
    vtkCleanPolyData* fltCleaner = vtkCleanPolyData::New();
    fltCleaner->SetInput(fltTriangle->GetOutput());
    fltCleaner->SetTolerance(0);
    fltCleaner->ConvertPolysToLinesOn();

    cout << "   cleaning up triangle mesh " << endl;
    fltCleaner->Update();

    // Go through and delete the cells that are of the wrong type
    //m_SurfaceMesh
    vtkPolyData* clean= fltCleaner->GetOutput();
    for(vtkIdType i = clean->GetNumberOfCells();i > 0;i--)
      {
      if(clean->GetCellType(i-1) != VTK_TRIANGLE)
        clean->DeleteCell(i-1);
      }
    clean->BuildCells();
  */

  vtkFeatureEdges * fltEdge = vtkFeatureEdges::New();
  fltEdge->BoundaryEdgesOff();
  fltEdge->FeatureEdgesOff();
  fltEdge->NonManifoldEdgesOff();
  fltEdge->ManifoldEdgesOn();
  fltEdge->ColoringOff();
  fltEdge->SetInput(m_SurfaceMesh);

  cout << "   extracting edges from the mesh" << endl;
  fltEdge->Update();

  // Got the new poly data
  vtkPolyData * m_EdgePolys = fltEdge->GetOutput();
  m_EdgePolys->BuildCells();
  m_EdgePolys->BuildLinks();

  unsigned int nEdges = m_EdgePolys->GetNumberOfLines();
  cout << "      number of edges (lines) : " << nEdges << endl;
  cout << "      number of cells : " << m_EdgePolys->GetNumberOfCells() << endl;
  cout << "      number if points : " << m_EdgePolys->GetNumberOfPoints() << endl;

  vtkPoints * vtkpoints = m_EdgePolys->GetPoints();
  int         numPoints = vtkpoints->GetNumberOfPoints();
  m_GraphX.resize(numPoints);
  for (int i = 0; i < numPoints; i++)
  {
    NodeLocationType                                                       loc;
    double *                                                               pt = vtkpoints->GetPoint(i);
    typename GraphSearchNode<PixelType, CoordRep, GraphDimension>::Pointer G =
      GraphSearchNode<PixelType, CoordRep, GraphDimension>::New();
    G->SetUnVisited();
    G->SetTotalCost(m_MaxCost);
    for (int j = 0; j < GraphDimension; j++)
    {
      loc[j] = pt[j];
    }
    G->SetLocation(loc);
    G->SetPredecessor(nullptr);
    G->m_NumberOfNeighbors = 0;
    m_GraphX[i] = G;
  }

  std::cout << " allocation of graph done ";

  vtkIdType   nPoints = 0;
  vtkIdType * xPoints = nullptr;
  for (unsigned int i = 0; i < nEdges; i++)
  {
    // Get the next edge
    m_EdgePolys->GetCellPoints(i, nPoints, xPoints);

    // Place the edge into the Edge structure
    assert(nPoints == 2);
    // Place the edge into the Edge structure
    //    std::cout << " nPoints " << nPoints << std::endl;
    //    std::cout << " pt " << xPoints[0] << " connects " << xPoints[1] << std::endl;
    assert(nPoints == 2);
    m_GraphX[xPoints[0]]->m_NumberOfNeighbors++;
  }

  std::cout << " counting nhood done ";
  // second, resize the vector for each G
  for (int i = 0; i < numPoints; i++)
  {
    m_GraphX[i]->m_Neighbors.resize(m_GraphX[i]->m_NumberOfNeighbors);
    m_GraphX[i]->m_NumberOfNeighbors = 0;
  }
  for (unsigned int i = 0; i < nEdges; i++)
  {
    // Get the next edge
    m_EdgePolys->GetCellPoints(i, nPoints, xPoints);
    // Place the edge into the Edge structure
    assert(nPoints == 2);
    m_GraphX[xPoints[0]]->m_Neighbors[m_GraphX[xPoints[0]]->m_NumberOfNeighbors] = m_GraphX[xPoints[1]];
    m_GraphX[xPoints[0]]->m_NumberOfNeighbors++;
  }

  //    vtkPolyDataConnectivityFilter* con = vtkPolyDataConnectivityFilter::New();
  //    con->SetExtractionModeToLargestRegion();
  //    con->SetInput(m_EdgePolys);
  //    m_SurfaceMesh=con->GetOutput();
  m_SurfaceMesh = m_EdgePolys;
}

template <typename TGraphSearchNode>
void
ManifoldIntegrationAlgorithm<TGraphSearchNode>::InitializeGraph()
{
  if (!m_SurfaceMesh)
  {
    return;
  }
  std::cout << " Generate graph from surface mesh " << std::endl;
  // get size of the surface mesh

  vtkExtractEdges * edgeex = vtkExtractEdges::New();
  edgeex->SetInput(m_SurfaceMesh);
  edgeex->Update();
  vtkPolyData * edg1 = edgeex->GetOutput();
  vtkIdType     nedg = edg1->GetNumberOfCells();
  vtkIdType     vers = m_SurfaceMesh->GetNumberOfPoints();
  int           nfac = m_SurfaceMesh->GetNumberOfPolys();
  float         g = 0.5 * (2.0 - vers + nedg - nfac);
  std::cout << " Genus " << g << std::endl;
  edg1->BuildCells();

  // now cruise through all edges and add to each node's neighbor list
  // first, count the num of edges for each node
  //    m_SurfaceMesh=edg1;

  vtkPoints * vtkpoints = edg1->GetPoints();
  int         numPoints = vtkpoints->GetNumberOfPoints();
  m_GraphX.resize(numPoints);
  for (int i = 0; i < numPoints; i++)
  {
    NodeLocationType                                                       loc;
    double *                                                               pt = vtkpoints->GetPoint(i);
    typename GraphSearchNode<PixelType, CoordRep, GraphDimension>::Pointer G =
      GraphSearchNode<PixelType, CoordRep, GraphDimension>::New();
    G->SetUnVisited();
    G->SetTotalCost(m_MaxCost);
    for (int j = 0; j < GraphDimension; j++)
    {
      loc[j] = pt[j];
    }
    G->SetLocation(loc);
    G->SetPredecessor(nullptr);
    G->m_NumberOfNeighbors = 0;
    m_GraphX[i] = G;
  }
  std::cout << " allocation of graph done ";

  std::cout << " begin edg iter ";
  vtkIdType   nPoints = 0;
  vtkIdType * xPoints = nullptr;
  for (unsigned int i = 0; i < nedg; i++)
  {
    // Get the next edge
    edg1->GetCellPoints(i, nPoints, xPoints);
    // Place the edge into the Edge structure
    //    std::cout << " nPoints " << nPoints << std::endl;
    //    std::cout << " pt " << xPoints[0] << " connects " << xPoints[1] << std::endl;
    assert(nPoints == 2);
    m_GraphX[xPoints[0]]->m_NumberOfNeighbors++;
  }

  std::cout << " counting nhood done ";
  // second, resize the vector for each G
  for (int i = 0; i < numPoints; i++)
  {
    m_GraphX[i]->m_Neighbors.resize(m_GraphX[i]->m_NumberOfNeighbors);
    m_GraphX[i]->m_NumberOfNeighbors = 0;
  }
  for (unsigned int i = 0; i < nedg; i++)
  {
    // Get the next edge
    edg1->GetCellPoints(i, nPoints, xPoints);
    // Place the edge into the Edge structure
    assert(nPoints == 2);
    m_GraphX[xPoints[0]]->m_Neighbors[m_GraphX[xPoints[0]]->m_NumberOfNeighbors] = m_GraphX[xPoints[1]];
    m_GraphX[xPoints[0]]->m_NumberOfNeighbors++;
  }

  m_SurfaceMesh = edg1;

  return;
}

template <typename TGraphSearchNode>
void
ManifoldIntegrationAlgorithm<TGraphSearchNode>::ConvertGraphBackToMesh()
{ // this is a sanity check
}

template <typename TGraphSearchNode>
void
ManifoldIntegrationAlgorithm<TGraphSearchNode>::InitializeQueue()
{
  int n = m_QS->m_SourceNodes.size();
  //  GraphIteratorType GraphIterator( m_Graph, m_GraphRegion );
  //  GraphIterator.GoToBegin();
  //  m_GraphIndex = GraphIterator.GetIndex();
  NodeLocationType loc;

  // make sure the graph contains the right pointers
  for (int i = 0; i < n; i++)
  {
    typename GraphSearchNode<PixelType, CoordRep, GraphDimension>::Pointer G = m_QS->m_SourceNodes[i];
    G->SetPredecessor(G);
    m_QS->m_Q.push(G);
    loc = G->GetLocation();
    //      for (int d=0;d<GraphDimension;d++) m_GraphIndex[d]=loc[d];
    //        m_Graph->SetPixel(m_GraphIndex,G);
  }
  for (unsigned int i = 0; i < m_QS->m_SinkNodes.size(); i++)
  {
    typename GraphSearchNode<PixelType, CoordRep, GraphDimension>::Pointer G = m_QS->m_SinkNodes[i];
    G->SetPredecessor(nullptr);
    loc = G->GetLocation();
    //      for (int d=0;d<GraphDimension;d++) m_GraphIndex[d]=(long)loc[d];
    //        m_Graph->SetPixel(m_GraphIndex,G);
  }
  m_SearchFinished = false;
}

/**
 *  parameterize the boundary --- an estimate
 */

template <typename TGraphSearchNode>
bool
ManifoldIntegrationAlgorithm<TGraphSearchNode>::ParameterizeBoundary(
  ManifoldIntegrationAlgorithm<TGraphSearchNode>::SearchNodePointer rootNode)
{
  std::vector<SearchNodePointer> neighborlist;
  bool                           I_Am_A_Neighbor = false;
  SearchNodePointer              neighbor = nullptr;
  SearchNodePointer              curNode = rootNode;
  //  unsigned int rootnn=rootNode>m_NumberOfNeighbors;
  unsigned int ct = 0;
  bool         canparam = false;
  unsigned int qsz = this->m_QS->m_Q.size();

  while (!I_Am_A_Neighbor && ct <= qsz * 3)
  {
    unsigned int limit = curNode->m_NumberOfNeighbors;
    for (unsigned int i = 0; i < limit; i++)
    {
      neighbor = curNode->m_Neighbors[i];
      bool inb = false;
      for (unsigned int q = 0; q < neighborlist.size(); q++)
      {
        if (neighbor == neighborlist[q])
        {
          inb = true;
        }
      }
      if (neighbor == rootNode && !inb && ct > 2)
      {
        I_Am_A_Neighbor = true;
        canparam = true;
      }
      if (neighbor->IsInQueue() && !inb) // add to border list
      {
        neighborlist.push_back(neighbor);
        curNode = neighbor;
        i = limit;
      } // add to border
    }   // neighborhood
    ct++;
  } // while

  if (neighborlist.size() >= this->m_BoundaryList.size() && canparam)
  {
    neighborlist.push_back(rootNode);
    this->m_BoundaryList.clear();
    this->m_BoundaryList.assign(neighborlist.begin(), neighborlist.end());
  }

  // if ( ct > 0 && canparam)  std::cout <<" qfrac " << this->m_BoundaryList.size()  << " canp  "<< canparam <<
  // " qsz "
  // << qsz << " cost " << m_CurrentCost << std::endl;
  return canparam;
}

/**
 *  Compute the local cost using Manhattan distance.
 */
template <typename TGraphSearchNode>
typename ManifoldIntegrationAlgorithm<TGraphSearchNode>::PixelType
ManifoldIntegrationAlgorithm<TGraphSearchNode>::MyLocalCost()
{
  NodeLocationType dif = m_CurrentNode->GetLocation() - m_NeighborNode->GetLocation();
  float            mag = 0.0;

  for (int jj = 0; jj < GraphDimension; jj++)
  {
    mag += dif[jj] * dif[jj];
  }
  mag = sqrt(mag);
  if (m_PureDist)
  {
    return mag;
  }
  else
  {
    float dL = fabs(m_CurrentNode->GetValue(3) - m_NeighborNode->GetValue(3));
    if (dL > 0.5)
    {
      dL = this->m_MaxCost * this->m_LabelCostWeight;
    }
    else
    {
      dL = 0;
    }
    return mag * this->m_DistanceCostWeight + dL;
  }
}

template <typename TGraphSearchNode>
bool
ManifoldIntegrationAlgorithm<TGraphSearchNode>::TerminationCondition()
{
  if (!m_QS->m_SinkNodes.empty())
  {
    if (m_NeighborNode == m_QS->m_SinkNodes[0] && !m_SearchFinished)
    {
      //      std::cout << " FOUND SINK ";
      m_SearchFinished = true;
      m_NeighborNode->SetTotalCost(m_CurrentCost + MyLocalCost());
      m_NeighborNode->SetPredecessor(m_CurrentNode);
    }
  }
  if (m_CurrentCost >= m_MaxCost)
  {
    m_SearchFinished = true;
  }
  return m_SearchFinished;
}

template <typename TGraphSearchNode>
void
ManifoldIntegrationAlgorithm<TGraphSearchNode>::SearchEdgeSet()
{
  int i = 0; // ,j=0;

  for (i = 0; i < m_CurrentNode->m_NumberOfNeighbors; i++)
  {
    m_NeighborNode = m_CurrentNode->m_Neighbors[i];
    //        std::cout << " i " << i << " position " << m_NeighborNode->GetLocation() << endl;
    //        std::cout << " i " << i << " position " << m_NeighborNode->GetLocation() << " label " <<
    // m_CurrentNode->GetValue() << endl;
    TerminationCondition();
    if (!m_SearchFinished && m_CurrentNode != m_NeighborNode && !m_NeighborNode->GetDelivered())
    {
      m_NewCost = m_CurrentCost + MyLocalCost();
      CheckNodeStatus();
    }
  }
}

template <typename TGraphSearchNode>
void
ManifoldIntegrationAlgorithm<TGraphSearchNode>::GetSearchBoundary()
{
  unsigned int gsz = this->GetGraphSize();

  for (unsigned int j = 0; j < gsz; j++)
  {
    this->m_CurrentNode = this->m_GraphX[j];
    if (this->m_CurrentNode)
    {
      const float cost = m_CurrentNode->GetTotalCost();
      if (cost <= this->m_MaxCost && (cost >= this->m_MaxCost - 4))
      {
        this->m_BoundaryList.push_back(this->m_CurrentNode);
      }
    }
  }
}

template <typename TGraphSearchNode>
void
ManifoldIntegrationAlgorithm<TGraphSearchNode>::CheckNodeStatus()
// checks a graph neighbor's status
{
  NodeLocationType dif = m_CurrentNode->GetLocation() - m_NeighborNode->GetLocation();

  //    std::cout << " visited? " << m_NeighborNode->GetVisited() <<
  //        " old cost " << m_NeighborNode->GetTotalCost() << " new cost " <<m_NewCost << std::endl;
  if (!m_NeighborNode->GetVisited() && !m_NeighborNode->GetUnVisitable())
  {
    // set the cost and put into the queue
    m_NeighborNode->SetTotalCost(m_NewCost);
    float delt = fabs(m_CurrentNode->GetValue() - m_NeighborNode->GetValue()); // *dif.magnitude();
    m_NeighborNode->SetValue(m_CurrentNode->GetValue() + delt);
    m_NeighborNode->SetPredecessor(m_CurrentNode);
    m_NeighborNode->SetVisited();

    float mag = 0.0;
    for (int jj = 0; jj < GraphDimension; jj++)
    {
      mag += dif[jj] * dif[jj];
    }
    mag = sqrt(mag);
    m_NeighborNode->SetValue(m_CurrentNode->GetValue(2) + mag, 2); // the actual manifold distance travelled
    //    if  (
    m_QS->m_Q.push(m_NeighborNode);
    // }
    // else {
    //  m_NeighborNode->SetUnVisitable();
    // }

    //    std::cout << " Pushing new node on " << m_NewCost << std::endl;
  }
  else if (m_NewCost < m_NeighborNode->GetTotalCost() && !m_NeighborNode->GetUnVisitable())
  {
    //      std::cout << " Updating " << std::endl;
    float delt = fabs(m_CurrentNode->GetValue() - m_NeighborNode->GetValue()); // *dif.magnitude();
    m_NeighborNode->SetValue(m_CurrentNode->GetValue() + delt);
    m_NeighborNode->SetTotalCost(m_NewCost);
    m_NeighborNode->SetPredecessor(m_CurrentNode);

    float mag = 0.0;
    for (int jj = 0; jj < GraphDimension; jj++)
    {
      mag += dif[jj] * dif[jj];
    }
    mag = sqrt(mag);
    m_NeighborNode->SetValue(m_CurrentNode->GetValue(2) + mag, 2); // the actual manifold distance travelled
    m_QS->m_Q.push(m_NeighborNode);
  }
}

template <typename TGraphSearchNode>
void
ManifoldIntegrationAlgorithm<TGraphSearchNode>::FindPath()
{
  if (m_QS->m_SourceNodes.empty())
  {
    std::cout << "ERROR !! DID NOT SET SOURCE!!\n";
    return;
  }

  //  std::cout << "MI start find path " << " Q size " << m_QS->m_Q.size() << " \n";

  while (!m_SearchFinished && !m_QS->m_Q.empty())
  {
    m_CurrentNode = m_QS->m_Q.top();
    m_CurrentCost = m_CurrentNode->GetTotalCost();
    if (this->m_ParamWhileSearching)
    {
      this->ParameterizeBoundary(this->m_CurrentNode);
    }
    m_QS->m_Q.pop();
    if (!m_CurrentNode->GetDelivered())
    {
      m_QS->IncrementTimer();
      // /std::cout << " searching " << m_CurrentNode->GetLocation()   << " \n";
      this->SearchEdgeSet();
      // if ( (m_CurrentNode->GetTimer() % 1.e5 ) == 0)
      // std::cout << " searched  " << m_CurrentNode->GetTimer()   << " \n";
    }
    m_CurrentNode->SetDelivered();
  } // end of while

  m_NumberSearched = (unsigned long)m_QS->GetTimer();
  //  std::cout << "Done with find path " << " Q size " << m_QS->m_Q.size() <<
  // " num searched " << m_NumberSearched << " \n";

  //  std::cout << " Max Distance " << m_CurrentCost << std::endl;
  if (!this->m_ParamWhileSearching)
  {
    this->GetSearchBoundary();
  }

  return;
}
} // namespace itk
#endif
