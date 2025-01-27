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
#ifndef _itkDijkstrasAlgorithm_h_
#define _itkDijkstrasAlgorithm_h_

#include <string>
#include <iostream>
#include <stack>
#include <vector>
#include <list>
#include <queue>
#include <map>

#include "itkMath.h"
// #include "vnl/vnl_matrix_fixed.h"
// #include "vnl/vnl_vector_fixed.h"
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNeighborhoodIterator.h"
#include "itkVector.h"
using namespace std;

namespace itk
{
/** The GraphSearchNode class defines a general shortest path graph node.
 *   The algorithm requires
 *   the node to have a pointer to itself and entry for the cumulative cost.
 *   We also define an index to its location and a couple of booleans
 *   for keeping track of the graph node's state.
 *   We assume the connectivity between nodes is defined externally.
 *   The class will also be useful for minimum spanning trees and
 *   other graph search algorithms.   Connectivity is defined externally.
 *   May be worthwhile to implement arbitrary connectivity e.g. for random graphs.
 *   One way to do this is to include a list of pointers which define
 *   the neighbors of this node, similar to how the predecessor is defined.
 */
template <typename TPixelType, typename TCoordRep = unsigned int, unsigned int NGraphDimension = 2>
typename GraphSearchNode : public itk::LightObject
{
public:
  /* Standard typedefs.*/
  typedef GraphSearchNode          Self;
  typedef LightObject              Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  itkOverrideGetNameOfClassMacro(GraphSearchNode);
  itkNewMacro(Self); /** Method for creation through the object factory.   */

  enum StateType
  {
    UnVisitedState,
    VisitedState,
    DeliveredState,
    UnVisitableState
  };
  enum
  {
    GraphDimension = NGraphDimension
  };
  typedef TPixelType                            PixelType; /** defines the cost data type */
  typedef TCoordRep                             CoordRep;  /** defines the location data type */
  typedef itk::Vector<CoordRep, GraphDimension> NodeLocationType;

  //  typedef typename itk::Image<float,GraphDimension>::IndexType  NodeLocationType;

  typedef vector<Pointer> NodeListType;

  //  typedef itk::Image<CoordRep,GraphDimension>::IndexType  NodeLocationType;

  inline void SetLocation(NodeLocationType loc) { m_Location = loc; }

  inline NodeLocationType GetLocation() { return m_Location; }

  inline void SetTotalCost(TPixelType cost) { m_TotalCost = cost; }

  inline void SetValue(TPixelType cost, int which = 0)
  {
    if (which <= 0)
    {
      m_Value1 = cost;
    }
    if (which == 1)
    {
      m_Value2 = cost;
    }
    if (which == 2)
    {
      m_Value3 = cost;
    }
    if (which >= 3)
    {
      m_Value4 = cost;
    }
  }

  inline TPixelType GetValue(int which = 0)
  {
    if (which <= 0)
    {
      return m_Value1;
    }
    if (which == 1)
    {
      return m_Value2;
    }
    if (which == 2)
    {
      return m_Value3;
    }
    if (which >= 3)
    {
      return m_Value4;
    }
    return m_Value1;
  }

  inline void SetUnVisited() { m_State = UnVisitedState; }

  inline void SetUnVisitable() { m_State = UnVisitableState; }

  inline void SetVisited() { m_State = VisitedState; }

  inline void SetDelivered() { m_State = DeliveredState; }

  inline bool IsInQueue()
  {
    if (m_State == VisitedState)
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  inline bool WasVisited()
  {
    if (m_State == VisitedState)
    {
      return true;
    }
    else if (m_State == DeliveredState)
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  inline TPixelType GetTotalCost() { return m_TotalCost; }

  inline void SetPredecessor(Pointer address) { m_PredecessorAddress = address; }

  inline Pointer GetPredecessor() { return m_PredecessorAddress; }

  inline void SetAncestor(Pointer address) { m_AncestorAddress = address; }

  inline Pointer GetAncestor() { return m_AncestorAddress; }

  inline bool GetUnVisited()
  {
    if (m_State == UnVisitedState)
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  inline bool GetUnVisitable()
  {
    if (m_State == UnVisitableState)
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  inline bool GetVisited()
  {
    if (m_State == VisitedState)
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  inline bool GetDelivered()
  {
    if (m_State == DeliveredState)
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  inline void SetState(StateType S) { m_State = S; }

  inline StateType GetState() { return m_State; }

  inline void SetIdentity(unsigned int i) { m_Identity = i; }

  inline unsigned int GetIdentity() { return m_Identity; }

  inline int GetNumberOfNeighbors() { return m_Neighbors.size(); }

  inline Pointer GetNeighbor(int i) { return m_Neighbors[i]; }

  void SetNeighborSize(int i) { m_Neighbors.resize(i); }

  NodeListType   m_Neighbors;
  unsigned short m_NumberOfNeighbors;
  unsigned int   m_Identity;

protected:
  GraphSearchNode()
  {
    m_TotalCost = 0.0;
    m_Value1 = 0.0;
    m_Value2 = 0.0;
    m_Value3 = 0.0;
    m_Value4 = 0.0;
    m_PredecessorAddress = nullptr;
    m_AncestorAddress = nullptr;
    m_State = UnVisitedState;
    m_NumberOfNeighbors = 0;
    m_Identity = 0;
  }

  ~GraphSearchNode() {}

private:
  TPixelType m_TotalCost; /** keeps track of the minimum accumulated cost. */
  TPixelType m_Value1;    /** keeps track of the minimum accumulated cost. */
  TPixelType m_Value2;    /** keeps track of the minimum accumulated cost. */
  TPixelType m_Value3;    /** keeps track of the minimum accumulated cost. */
  TPixelType m_Value4;    /** keeps track of the minimum accumulated cost. */

  StateType        m_State;
  NodeLocationType m_Location;           /** provides the location in the graph. */
  Pointer          m_PredecessorAddress; /** provides the best predecessor address */
  Pointer          m_AncestorAddress;    /** provides the best predecessor address */

  GraphSearchNode(const Self &); // purposely not implemented
  void operator=(const Self &);  // purposely not implemented
};

// Forward declaration of DijkstrasAlgorithm so it can be declared a friend
template <typename TGraphSearchNode>
class DijkstrasAlgorithm;

template <typename TGraphSearchNode>
class DijkstrasAlgorithmQueue : public itk::LightObject
/** \class DijkstrasAlgorithmQueue
the class containing the priority queue and associated data.
*/
{
private:
  template <typename G>
  typename GraphSearchNodePriority /* defines the comparison operator for the prioritiy queue */
    { public: bool operator()(typename G::Pointer N1,
                              typename G::Pointer N2){ return N1->GetTotalCost() > N2->GetTotalCost();
}
};      // namespace itk
public: /* Standard typedefs.*/
typedef DijkstrasAlgorithmQueue  Self;
typedef LightObject              Superclass;
typedef SmartPointer<Self>       Pointer;
typedef SmartPointer<const Self> ConstPointer;
itkOverrideGetNameOfClassMacro(DijkstrasAlgorithmQueue);
itkNewMacro(Self); /** Method for creation through the object factory.   */

typedef typename TGraphSearchNode::Pointer   TGraphSearchNodePointer;
typedef typename TGraphSearchNode::PixelType PixelType; /** pixel type for the cost */
typedef typename TGraphSearchNode::CoordRep  CoordRep;  /** type for coordinates */
typedef typename std::priority_queue<typename TGraphSearchNode::Pointer,
                                     std::vector<typename TGraphSearchNode::Pointer>,
                                     GraphSearchNodePriority<TGraphSearchNode>>
                                                   QType; /** the queue we are using */
typedef vector<typename TGraphSearchNode::Pointer> NodeListType;
inline QType
GetQ()
{
  return m_Q;
}

void
AddToPath(TGraphSearchNodePointer G)
{
  this->m_Path.push_back(G);
}

inline NodeListType
GetPath()
{
  return m_Path;
}

void
EmptyPath()
{
  m_Path.clear();
}

inline NodeListType
GetSourceNodes()
{
  return m_SourceNodes;
}

inline NodeListType
GetSinkNodes()
{
  return m_SinkNodes;
}

inline void
IncrementTimer()
{
  m_timer++;
}

inline long
GetTimer()
{
  return m_timer;
}

inline void
EmptyQ()
{
  while (!m_Q.empty())
  {
    m_Q.pop();
  }

  m_timer = 0;
  m_SourceNodes.clear();
  m_SinkNodes.clear();
}

NodeListType m_SinkNodes;
NodeListType m_SourceNodes;
QType        m_Q;
NodeListType m_Path;

protected:
friend class DijkstrasAlgorithm<TGraphSearchNode>; // so it can access this data easily

DijkstrasAlgorithmQueue()
{
  m_timer = 0;
}

~DijkstrasAlgorithmQueue() {}

private:
unsigned long m_timer;
DijkstrasAlgorithmQueue(const Self &); // purposely not implemented
void
operator=(const Self &); // purposely not implemented
}
;

/**
 * \class DijkstrasAlgorithm
 * \brief General shortest path / greedy dynamic programming solver.
 *
 *  This class uses Dijkstra's algorithm to solve the shortest path problem.
 *  We use the stl priority queue which is not optimal for this problem, but which
 *  works o.k.  It's implemented to be used for general regularly connected
 *  graphs with fixed costs, or for the variational solution of an integral
 *  curve matching energy.
 *  Note: we assume all edge weights are positive.
 *  The class is implemented as a an abstract base class, with virtual functions for
 *  LocalCost, SearchEdgeSet and FindPath.  LocalCost must be implemented by derived classes.
 *  The connectivity of the graph defined
 *  here is always regular and is controlled by a set of neighborhood indices.
 *  the default is a radius 1 neighborhood with all entries used.  However, the
 *  user may also supply her own regular connectivity by selecting the size of
 *  the neighborhood and a subset of the indices which define the edges.  If
 *  the GraphSearchNode contains its edges, they may be used with minor modification
 *  to the function SearchEdgeSet.
 *  Another improvement would be to make the LocalCost function a pointer
 *  to a function which could be set.
 */
template <typename TGraphSearchNode>
class DijkstrasAlgorithm : public itk::LightObject
{
public:
  typedef DijkstrasAlgorithm       Self;
  typedef LightObject              Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  itkOverrideGetNameOfClassMacro(DijkstrasAlgorithm);
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
  typedef Image<SearchNodePointer, GraphDimension>                         GraphType;
  typedef typename GraphType::SizeType                                     GraphSizeType;
  typedef ImageRegionIteratorWithIndex<GraphType>                          GraphIteratorType;
  typedef typename GraphType::RegionType                                   GraphRegionType;
  typedef typename DijkstrasAlgorithmQueue<TGraphSearchNode>::Pointer      QType;
  typedef typename DijkstrasAlgorithmQueue<TGraphSearchNode>::NodeListType NodeListType;
  typedef itk::NeighborhoodIterator<GraphType>                             GraphNeighborhoodIteratorType;
  typedef typename GraphNeighborhoodIteratorType::IndexType                GraphNeighborhoodIndexType;
  typedef typename GraphNeighborhoodIteratorType::RadiusType               RadiusType;
  typedef typename TGraphSearchNode::NodeLocationType                      NodeLocationType;
  typedef typename GraphType::IndexType                                    IndexType;
  // FUNCTIONS
  void
  InitializeGraph(); /** initializes all graph values appropriately */

  void
  InitializeQueue(); /** initializes all queue values appropriately
                          call AFTER source and sink are set*/

  void
  InitializeEdgeTemplate(); /** helper function initializes edge set appropriately */

  void
  InitializeEdgeTemplate(vector<unsigned int>, unsigned int); /** user supplied edge template */

  void
  SetGraphSize(typename GraphType::SizeType Sz); /** the rectangular size of the graph */

  inline void
  EmptyQ()
  {
    m_QS->EmptyQ();
    this->m_TotalCost = 0;
  }

  /* adds a source to the source set */
  void
  SetSource(typename TGraphSearchNode::Pointer G)
  {
    m_QS->m_SourceNodes.push_back(G);
    for (int i = 0; i < GraphDimension; i++)
    {
      m_GraphIndex[i] = (long int)(G->GetLocation()[i] + 0.5);
      //      ::std::cout << " mgi " << m_GraphIndex[i];
    }
    m_Graph->SetPixel(m_GraphIndex, G);
  };

  typename TGraphSearchNode::Pointer
  GetGraphNode(IndexType index)
  {
    //    ::std::cout << " get node "  << index << std::endl;
    return m_Graph->GetPixel(index);
  };

  // adds a sink to the sink set
  void
  SetSink(typename TGraphSearchNode::Pointer G)
  {
    m_QS->m_SinkNodes.push_back(G);
  }

  // Backtracks from the given node to its source node;
  void
  BackTrack(typename TGraphSearchNode::Pointer SinkNode)
  {
    m_QS->m_Path.clear();

    typename TGraphSearchNode::Pointer G = SinkNode;
    typename TGraphSearchNode::Pointer P = SinkNode->GetPredecessor();
    if (!P || !G)
    {
      return;
    }
    float highcost = G->GetValue();
    if (G->GetTotalCost() > P->GetValue())
    {
      P->SetAncestor(G);
      P->SetValue(G->GetTotalCost());
      highcost = G->GetTotalCost();
    }

    while (P && G != P)
    {
      m_QS->m_Path.push_back(G);
      G = P;
      P = G->GetPredecessor();
      if (P->GetValue() < highcost)
      {
        P->SetValue(highcost);
        P->SetAncestor(G);
      }
    }

    if (!P)
    {
      cout << " null pred "; // else cout << " pred == self \n";
    }
    return;
  }

  // Inverse of backtrack - from the given node to its sink node;
  void
  ForwardTrack(typename TGraphSearchNode::Pointer SinkNode)
  {
    typename TGraphSearchNode::Pointer G = SinkNode;
    typename TGraphSearchNode::Pointer P = SinkNode->GetAncestor();
    while (P && G != P && G)
    {
      if (P->GetValue() > G->GetValue())
      {
        G->SetValue(P->GetValue());
      }
      if (G->GetValue() > P->GetValue())
      {
        P->SetValue(G->GetValue());
      }
      G = P;
      P = G->GetAncestor();
    }

    return;
  }

  virtual bool
  TerminationCondition(); /** decides when the algorithm stops */

  virtual void
  SearchEdgeSet(); /** loops over the neighbors in the graph */

  void
  CheckNodeStatus(); /** checks if the node has been explored already, its cost, etc. */

  virtual PixelType
  LocalCost(); /* computes the local cost */

  /* alternatively, we could pass the function as a template parameter
     or set a function pointer.  the latter method is used in dijkstrasegment. */

  virtual void
  FindPath(); /* runs the algorithm */

  inline unsigned int
  GetPathSize()
  {
    return m_QS->m_Path.size();
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

  void
  SetMaxCost(PixelType m)
  {
    m_MaxCost = m;
  }

  double
  GetTotalCost()
  {
    return m_TotalCost;
  }

  void
  SetSearchFinished(bool m)
  {
    m_SearchFinished = m;
  }

  /** sets the boolean that indicates if the algorithm is done */
protected:
  QType                              m_QS;
  vector<unsigned int>               m_EdgeTemplate;    /** defines neighborhood connectivity */
  RadiusType                         m_Radius;          /** used by the neighborhood iterator */
  typename TGraphSearchNode::Pointer m_PredecessorNode; /** holds the predecessor node */
  typename TGraphSearchNode::Pointer m_CurrentNode;     /** holds the current node */
  typename TGraphSearchNode::Pointer m_NeighborNode;    /** holds the current neighbor node */
  typename GraphType::Pointer        m_Graph;           /** holds all the graph information */
  GraphRegionType                    m_GraphRegion;
  GraphSizeType                      m_GraphSize; /** rectangular size of graph */

  typename GraphType::IndexType m_GraphIndex;
  bool                          m_SearchFinished;
  PixelType                     m_NewCost;
  PixelType                     m_CurrentCost;
  PixelType                     m_MaxCost; // This creates an insurmountable barrier unless all costs are max

  double m_TotalCost;

  unsigned long m_NumberSearched;
  DijkstrasAlgorithm();
  ~DijkstrasAlgorithm(){};

private:
  DijkstrasAlgorithm(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkDijkstrasAlgorithm.cxx"
#endif

#endif
