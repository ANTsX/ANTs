/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit (ITK)

Copyright (c) 2001 Insight Consortium
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

 * The name of the Insight Consortium, nor the names of any consortium members,
   nor of any contributors, may be used to endorse or promote products derived
   from this software without specific prior written permission.

  * Modified source versions must be plainly marked as such, and must not be
    misrepresented as being the original software.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
#ifndef _itkManifoldIntegrationAlgorithm_h_
#define _itkManifoldIntegrationAlgorithm_h_

namespace itk
{
/**
 * \class ManifoldIntegrationAlgorithm
 * \brief General shortest path / greedy dynamic programming solver.
 */
template <typename TGraphSearchNode>
class ManifoldIntegrationAlgorithm : public itk::LightObject
{
public:
  typedef ManifoldIntegrationAlgorithm Self;
  typedef LightObject                  Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;
  itkOverrideGetNameOfClassMacro(ManifoldIntegrationAlgorithm);
  itkNewMacro(Self);

  // Computation Data
  typedef TGraphSearchNode             SearchNode; /** dimension of the graph */
  typedef typename SearchNode::Pointer SearchNodePointer;
  enum
  {
    GraphDimension = SearchNode::GraphDimension
  };                                                                                  /** dimension of the graph */
  typedef typename SearchNode::PixelType                                   PixelType; /**  pixel type for the cost */
  typedef typename SearchNode::CoordRep                                    CoordRep;  /** coordinate type */
  typedef typename DijkstrasAlgorithmQueue<TGraphSearchNode>::Pointer      QType;
  typedef typename DijkstrasAlgorithmQueue<TGraphSearchNode>::NodeListType NodeListType;

  typedef typename TGraphSearchNode::NodeLocationType NodeLocationType;

  typedef vtkPolyData               TriangulationType;
  typedef vtkPolyData *             TriangulationTypePointer;
  typedef vector<SearchNodePointer> GraphType;

  // FUNCTIONS
  // estimate the metric tensor of the surface and also the (conjugate harmonic) function dstarU
  void
  GetSearchBoundary();

  float
  dstarUestimate(SearchNodePointer G);

  void
  InitializeGraph3();

  void
  InitializeGraph(); /** initializes all graph values appropriately */

  void
  InitializeGraph2(); /** initializes all graph values appropriately */

  void
  InitializeQueue(); /** initializes all queue values appropriately
                          call AFTER source and sink are set*/

  inline void
  EmptyQ()
  {
    m_QS->EmptyQ();
    m_QS->EmptyPath();
  }

  /* adds a source to the source set */
  void
  SetSource(typename TGraphSearchNode::Pointer G)
  {
    G->SetTotalCost(0.0);
    m_QS->m_SourceNodes.push_back(G);
  };

  // adds a sink to the sink set
  void
  SetSink(typename TGraphSearchNode::Pointer G)
  {
    m_QS->m_SinkNodes.push_back(G);
  }

  // Backtracks from the given node to its source node;
  float
  BackTrack(typename TGraphSearchNode::Pointer SinkNode)
  {
    m_QS->m_Path.clear();

    typename TGraphSearchNode::Pointer G = SinkNode;
    typename TGraphSearchNode::Pointer P = SinkNode->GetPredecessor();
    if (!P)
    {
      std::cout << " ERROR NO PRED TO SINK " << std::endl;
      return 1.;
    }
    m_QS->m_Path.push_back(G);
    //    if (P->GetAncestor() && G)
    //    {
    //      if (P->GetAncestor()->GetTotalCost() < G->GetTotalCost() )  P->SetAncestor(G);
    //    }
    //    else if (G)

    P->SetAncestor(G);
    //      if (P->GetValue(1) < G->GetValue(1) ) P->SetValue(G->GetValue(1),1);

    while (P && G != P)
    {
      //        std::cout << " Backtrack " << G->GetValue(0) << std::endl;
      G = P;
      P = G->GetPredecessor();
      //    if (P->GetValue(1) < G->GetValue(1) ) P->SetValue(G->GetValue(1),1);
      P->SetAncestor(G);
      if (P)
      {
        m_QS->m_Path.push_back(P);
      }
    }

    //    m_QS->m_Path.push_back(P);
    //    std::cout << " final cost " << P->GetTotalCost() << " high " << highcost << std::endl;
    if (!P)
    {
      cout << " null pred "; // else cout << " pred == self \n";
    }
    return P->GetValue(2);
  }

  // Backtracks from the given node to its source node performing integration
  void
  IntegrateBackward(typename TGraphSearchNode::Pointer SinkNode)
  {
    typename TGraphSearchNode::Pointer G = SinkNode;
    typename TGraphSearchNode::Pointer P = SinkNode->GetPredecessor();

    float intval = 0.0;
    G->SetTotalCost(intval);
    while (P && G != P)
    {
      //        NodeLocationType dif=P->GetLocation()-G->GetLocation();
      float dU = (P->GetValue() - G->GetValue());
      intval += dU;
      P->SetTotalCost(intval);
      G = P;
      P = G->GetPredecessor();
    }

    //    std::cout << " intval " << intval << " at " << G->GetLocation() << std::endl;
    if (!P)
    {
      cout << " null pred "; // else cout << " pred == self \n";
    }
    return;
  }

  // Backtracks from the given node to its source node performing integration
  void
  IntegrateForward(typename TGraphSearchNode::Pointer SinkNode)
  {
    typename TGraphSearchNode::Pointer G = SinkNode;
    typename TGraphSearchNode::Pointer P = SinkNode->GetAncestor();

    float intval = 0.0;
    G->SetTotalCost(intval);
    while (P && G != P)
    {
      //        NodeLocationType dif=P->GetLocation()-G->GetLocation();
      float dU = (P->GetValue() - G->GetValue());
      intval += dU;
      P->SetTotalCost(intval);
      G = P;
      P = G->GetAncestor();
    }

    //    std::cout << " intval " << intval << " at " << G->GetLocation() << std::endl;
    if (!P)
    {
      cout << " null pred "; // else cout << " pred == self \n";
    }
    return;
  }

  // Inverse of backtrack - from the given node to its sink node;
  float
  ForwardTrack(typename TGraphSearchNode::Pointer SinkNode)
  {
    typename TGraphSearchNode::Pointer G = SinkNode;
    typename TGraphSearchNode::Pointer P = SinkNode->GetAncestor();

    if (!P)
    {
      return G->GetValue(2);
    }
    //    float highcost=G->GetTotalCost();

    while (P && G != P && G)
    {
      G = P;
      P = G->GetAncestor();
      if (!P)
      {
        return G->GetValue(2);
      }
      //     if (P->GetTotalCost() > highcost) highcost=P->GetTotalCost();
    }

    //    SinkNode->SetValue(highcost);
    //    if (!P) cout << " null pred "; //else cout << " pred == self \n";
    return P->GetValue(2);
  }

  bool ParameterizeBoundary(SearchNodePointer);

  bool
  TerminationCondition(); /** decides when the algorithm stops */

  virtual void
  SearchEdgeSet(); /** loops over the neighbors in the graph */

  void
  CheckNodeStatus(); /** checks if the node has been explored already, its cost, etc. */

  virtual PixelType
  MyLocalCost(); /* computes the local cost */

  /* alternatively, we could pass the function as a template parameter
     or set a function pointer.  the latter method is used in dijkstrasegment. */

  virtual void
  FindPath(); /* runs the algorithm */

  inline unsigned int
  GetPathSize()
  {
    return m_QS->m_Path.size();
  }

  inline void
  EmptyPath()
  {
    m_QS->m_Path.clear();
  }

  inline typename TGraphSearchNode::Pointer
  GetPathAtIndex(unsigned int i)
  {
    return m_QS->m_Path[i];
  }

  inline typename TGraphSearchNode::Pointer
  GetNeighborNode()
  {
    return m_NeighborNode;
  }

  inline typename TGraphSearchNode::Pointer
  GetCurrentNode()
  {
    return m_CurrentNode;
  }

  PixelType
  GetMaxCost()
  {
    return this->m_MaxCost;
  }

  void
  SetMaxCost(PixelType m)
  {
    this->m_MaxCost = m;
  }

  void
  ResetMaxCost()
  {
    this->m_MaxCost = vnl_huge_val(this->m_MaxCost);
  }

  inline void
  SetDistanceCostWeight(float d)
  {
    this->m_DistanceCostWeight = d;
  }

  inline void
  SetLabelCostWeight(float d)
  {
    this->m_LabelCostWeight = d;
  }

  inline void
  SetWeights(float a, float b, float c)
  {
    this->m_MaxCost = a;
    this->m_DistanceCostWeight = b;
    this->m_LabelCostWeight = c;
  }

  inline void
  PrintWeights()
  {
    std::cout << this->m_MaxCost << " " << this->m_DistanceCostWeight << " " << this->m_LabelCostWeight << std::endl;
  }

  void
  SetSearchFinished(bool m)
  {
    m_SearchFinished = m;
  }

  /** sets the boolean that indicates if the algorithm is done */

  void
  SetSurfaceMesh(TriangulationTypePointer mesh)
  {
    m_SurfaceMesh = mesh;
  }

  TriangulationTypePointer
  GetSurfaceMesh()
  {
    return m_SurfaceMesh;
  }

  SearchNodePointer
  GetGraphNode(int i)
  {
    return m_GraphX[i];
  }

  int
  GetGraphSize()
  {
    return m_GraphX.size();
  }

  // sanity check to see if mesh to graph conversion is ok
  // see if genus is the same
  void
  ConvertGraphBackToMesh();

  void
  SetParamWhileSearching(bool b)
  {
    this->m_ParamWhileSearching = b;
  }

  std::vector<SearchNodePointer> m_BoundaryList;

protected:
  QType                m_QS;
  vector<unsigned int> m_EdgeTemplate; /** defines neighborhood connectivity */

  typename TGraphSearchNode::Pointer m_PredecessorNode; /** holds the predecessor node */
  typename TGraphSearchNode::Pointer m_CurrentNode;     /** holds the current node */
  typename TGraphSearchNode::Pointer m_NeighborNode;    /** holds the current neighbor node */

  bool      m_PureDist;
  bool      m_SearchFinished;
  PixelType m_NewCost;
  PixelType m_CurrentCost;
  float     m_MaxCost; // This creates an insurmountable barrier unless all costs are max
  float     m_DistanceCostWeight;
  float     m_LabelCostWeight;
  GraphType m_GraphX;

  unsigned long m_NumberSearched;

  TriangulationTypePointer m_SurfaceMesh;
  bool                     m_ParamWhileSearching;

  ManifoldIntegrationAlgorithm();
  ~ManifoldIntegrationAlgorithm(){};

private:
  ManifoldIntegrationAlgorithm(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkManifoldIntegrationAlgorithm.cxx"
#endif

#endif
