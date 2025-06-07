/*=========================================================================

Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __antsFMarchingImageFilter_h
#define __antsFMarchingImageFilter_h

#include "itkArray.h"
#include "itkImageToImageFilter.h"
#include "itkIndex.h"
#include "itkLevelSet.h"
#include "itkNeighborhoodIterator.h"

#include "itkMath.h"

#include <functional>
#include <queue>

namespace itk
{
/** \class FMarchingImageFilter
 * \brief Solve an Eikonal equation using Fast Marching
 *
 * Fast marching solves an Eikonal equation where the speed is always
 * non-negative and depends on the position only. Starting from an
 * initial position on the front, fast marching systematically moves the
 * front forward one grid point at a time.
 *
 * Updates are preformed using an entropy satisfy scheme where only
 * "upwind" neighborhoods are used. This implementation of Fast Marching
 * uses a std::priority_queue to locate the next proper grid position to
 * update.
 *
 * Fast Marching sweeps through N grid points in (N log N) steps to obtain
 * the arrival time value as the front propagates through the grid.
 *
 * Implementation of this class is based on Chapter 8 of
 * "Level Set Methods and Fast Marching Methods", J.A. Sethian,
 * Cambridge Press, Second edition, 1999.
 *
 * This class is templated over the level set image type and the speed
 * image type. The initial front is specified by two containers: one
 * containing the known points and one containing the trial
 * points.  Alive points are those that are already part of the
 * object, and trial points are considered for inclusion.
 * In order for the filter to evolve, at least some trial
 * points must be specified.  These can for instance be specified as the layer of
 * pixels around the alive points.

 * The speed function can be specified as a speed image or a
 * speed constant. The speed image is set using the method
 * SetInput(). If the speed image is NULL, a constant speed function
 * is used and is specified using method the SetSpeedConstant().
 *
 * If the speed function is constant and of value one, fast marching results
 * in an approximate distance function from the initial alive points.
 * FMarchingImageFilter is used in the ReinitializeLevelSetImageFilter
 * object to create a signed distance function from the zero level set.
 *
 * The algorithm can be terminated early by setting an appropriate stopping
 * value. The algorithm terminates when the current arrival time being
 * processed is greater than the stopping value.
 *
 * There are two ways to specify the output image information
 * ( LargestPossibleRegion, Spacing, Origin): (a) it is copied directly from
 * the input speed image or (b) it is specified by the user. Default values
 * are used if the user does not specify all the information.
 *
 * The output information is computed as follows.
 * If the speed image is NULL or if the OverrideOutputInformation is set to
 * true, the output information is set from user specified parameters. These
 * parameters can be specified using methods SetOutputRegion(), SetOutputSpacing(), SetOutputDirection(),
 * and SetOutputOrigin(). Else if the speed image is not NULL, the output information
 * is copied from the input speed image.
 *
 * Possible Improvements:
 * In the current implemenation, std::priority_queue only allows
 * taking nodes out from the front and putting nodes in from the back.
 * To update a value already on the heap, a new node is added to the heap.
 * The defunct old node is left on the heap. When it is removed from the
 * top, it will be recognized as invalid and not used.
 * Future implementations can implement the heap in a different way
 * allowing the values to be updated. This will generally require
 * some sift-up and sift-down functions and
 * an image of back-pointers going from the image to heap in order
 * to locate the node which is to be updated.
 *
 * \sa LevelSetTypeDefault
 * \ingroup LevelSetSegmentation
 */
template <typename TLevelSet, typename TSpeedImage = Image<float, TLevelSet::ImageDimension>>
class FMarchingImageFilter final : public ImageToImageFilter<TSpeedImage, TLevelSet>
{
public:
  /** Standard class typdedefs. */
  typedef FMarchingImageFilter     Self;
  typedef ImageSource<TLevelSet>   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(FMarchingImageFilter);

  /** Typedef support of level set method types. */
  typedef LevelSetTypeDefault<TLevelSet>              LevelSetType;
  typedef typename LevelSetType::LevelSetImageType    LevelSetImageType;
  typedef typename LevelSetType::LevelSetPointer      LevelSetPointer;
  typedef typename LevelSetType::PixelType            PixelType;
  typedef typename LevelSetType::NodeType             NodeType;
  typedef typename LevelSetType::NodeContainer        NodeContainer;
  typedef typename LevelSetType::NodeContainerPointer NodeContainerPointer;
  typedef typename LevelSetImageType::SizeType        OutputSizeType;
  typedef typename LevelSetImageType::RegionType      OutputRegionType;
  typedef typename LevelSetImageType::SpacingType     OutputSpacingType;
  typedef typename LevelSetImageType::DirectionType   OutputDirectionType;
  typedef typename LevelSetImageType::PointType       OutputPointType;

  class AxisNodeType : public NodeType
  {
  public:
    AxisNodeType() { this->m_Axis = 0; }

    int
    GetAxis() const
    {
      return m_Axis;
    }

    void
    SetAxis(int axis)
    {
      m_Axis = axis;
    }

    const AxisNodeType &
    operator=(const NodeType & node)
    {
      this->NodeType::operator=(node);
      return *this;
    }

  private:
    int m_Axis;
  };

  /** SpeedImage typedef support. */
  typedef TSpeedImage SpeedImageType;

  /** SpeedImagePointer typedef support. */
  typedef typename SpeedImageType::Pointer      SpeedImagePointer;
  typedef typename SpeedImageType::ConstPointer SpeedImageConstPointer;

  /** Dimension of the level set and the speed image. */
  static constexpr unsigned int SetDimension = LevelSetType::SetDimension;
  static constexpr unsigned int SpeedImageDimension = SpeedImageType::ImageDimension;

  /** Index typedef support. */
  typedef Index<Self::SetDimension> IndexType;

  /** Enum of Fast Marching algorithm point types. FarPoints represent far
   * away points; TrialPoints represent points within a narrowband of the
   * propagating front; and AlivePoints represent points which have already
   * been processed.  Topology points were trial points but their inclusion
   * would have violated topology checks.
   */
  enum LabelType
  {
    FarPoint,
    AlivePoint,
    TrialPoint,
    InitialTrialPoint,
    TopologyPoint
  };

  /** LabelImage typedef support. */
  typedef Image<unsigned char, Self::SetDimension> LabelImageType;
  typedef NeighborhoodIterator<LabelImageType>                       NeighborhoodIteratorType;

  /** LabelImagePointer typedef support. */
  typedef typename LabelImageType::Pointer LabelImagePointer;

  /** ConnectedComponentImage typedef support. */
  typedef Image<unsigned int, Self::SetDimension> ConnectedComponentImageType;

  /** ConnectedComponentImagePointer typedef support. */
  typedef typename ConnectedComponentImageType::Pointer ConnectedComponentImagePointer;

  /** Set the container of Alive Points representing the initial front.
   * Alive points are represented as a VectorContainer of LevelSetNodes. */
  void
  SetAlivePoints(NodeContainer * points)
  {
    m_AlivePoints = points;
    this->Modified();
  }

  /** Get the container of Alive Points representing the initial front. */
  NodeContainerPointer
  GetAlivePoints()
  {
    return m_AlivePoints;
  }

  /** Set the container of Trial Points representing the initial front.
   * Trial points are represented as a VectorContainer of LevelSetNodes. */
  void
  SetTrialPoints(NodeContainer * points)
  {
    m_TrialPoints = points;
    this->Modified();
  }

  /** Get the container of Trial Points representing the initial front. */
  NodeContainerPointer
  GetTrialPoints()
  {
    return m_TrialPoints;
  }

  /** Get the point type label image. */
  LabelImagePointer
  GetLabelImage() const
  {
    return m_LabelImage;
  }

  /** Get the point type label image. */
  ConnectedComponentImagePointer
  GetConnectedComponentImage() const
  {
    return m_ConnectedComponentImage;
  }

  /** Set the Speed Constant. If the Speed Image is NULL,
   * the SpeedConstant value is used for the whole level set.
   * By default, the SpeedConstant is set to 1.0. */
  void
  SetSpeedConstant(double value)
  {
    m_SpeedConstant = value;
    m_InverseSpeed = -1.0 * itk::Math::sqr(1.0 / m_SpeedConstant);
    this->Modified();
  }

  /** Get the Speed Constant. */
  itkGetConstReferenceMacro(SpeedConstant, double);

  /** Set/Get the Normalization Factor for the Speed Image.
      The values in the Speed Image is divided by this
      factor. This allows the use of images with
      integer pixel types to represent the speed. */
  itkSetMacro(NormalizationFactor, double);
  itkGetConstMacro(NormalizationFactor, double);

  /** Set the Fast Marching algorithm Stopping Value. The Fast Marching
   * algorithm is terminated when the value of the smallest trial point
   * is greater than the stopping value. */
  itkSetMacro(StoppingValue, double);

  /** Get the Fast Marching algorithm Stopping Value. */
  itkGetConstReferenceMacro(StoppingValue, double);

  /** Set the Collect Points flag. Instrument the algorithm to collect
   * a container of all nodes which it has visited. Useful for
   * creating Narrowbands for level set algorithms that supports
   * narrow banding. */
  itkSetMacro(CollectPoints, bool);

  /** Get thConste Collect Points flag. */
  itkGetConstReferenceMacro(CollectPoints, bool);
  itkBooleanMacro(CollectPoints);

  using TopologyCheckType = itk::FastMarchingTraitsEnums::TopologyCheck;

  /** Set/Get boolean macro indicating whether the user wants to check topology. */
  itkSetMacro(TopologyCheck, TopologyCheckType);
  itkGetConstReferenceMacro(TopologyCheck, TopologyCheckType);

  /** Get the container of Processed Points. If the CollectPoints flag
   * is set, the algorithm collects a container of all processed nodes.
   * This is useful for defining creating Narrowbands for level
   * set algorithms that supports narrow banding. */
  NodeContainerPointer
  GetProcessedPoints() const
  {
    return m_ProcessedPoints;
  }

  /** The output largeset possible, spacing and origin is computed as follows.
   * If the speed image is NULL or if the OverrideOutputInformation is true,
   * the output information is set from user specified parameters. These
   * parameters can be specified using methods SetOutputRegion(), SetOutputSpacing(), SetOutputDirection(),
   * and SetOutputOrigin(). Else if the speed image is not NULL, the output information
   * is copied from the input speed image. */
  virtual void
  SetOutputSize(const OutputSizeType & size)
  {
    m_OutputRegion = size;
  }

  virtual OutputSizeType
  GetOutputSize() const
  {
    return m_OutputRegion.GetSize();
  }

  itkSetMacro(OutputRegion, OutputRegionType);
  itkGetConstReferenceMacro(OutputRegion, OutputRegionType);
  itkSetMacro(OutputSpacing, OutputSpacingType);
  itkGetConstReferenceMacro(OutputSpacing, OutputSpacingType);
  itkSetMacro(OutputDirection, OutputDirectionType);
  itkGetConstReferenceMacro(OutputDirection, OutputDirectionType);
  itkSetMacro(OutputOrigin, OutputPointType);
  itkGetConstReferenceMacro(OutputOrigin, OutputPointType);
  itkSetMacro(OverrideOutputInformation, bool);
  itkGetConstReferenceMacro(OverrideOutputInformation, bool);
  itkBooleanMacro(OverrideOutputInformation);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(SameDimensionCheck, (Concept::SameDimension<SetDimension, SpeedImageDimension>));
  itkConceptMacro(SpeedConvertibleToDoubleCheck, (Concept::Convertible<typename TSpeedImage::PixelType, double>));
  itkConceptMacro(DoubleConvertibleToLevelSetCheck, (Concept::Convertible<double, PixelType>));
  itkConceptMacro(LevelSetOStreamWritableCheck, (Concept::OStreamWritable<PixelType>));
  /** End concept checking */
#endif
protected:
  FMarchingImageFilter();
  ~FMarchingImageFilter() override{};
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  virtual void
  Initialize(LevelSetImageType *) final;

  virtual void
  UpdateNeighbors(const IndexType & index, const SpeedImageType *, LevelSetImageType *) final;

  virtual double
  UpdateValue(const IndexType & index, const SpeedImageType *, LevelSetImageType *) final;

  const AxisNodeType &
  GetNodeUsedInCalculation(unsigned int idx) const
  {
    return m_NodesUsed[idx];
  }

  void
  GenerateData() override;

  /** Generate the output image meta information. */
  void
  GenerateOutputInformation() override;

  void
  EnlargeOutputRequestedRegion(DataObject * output) override;

  /** Get Large Value. This value is used to
      represent the concept of infinity for the time assigned to pixels that
      have not been visited. This value is set by default to half the
      max() of the pixel type used to represent the time-crossing map. */
  itkGetConstReferenceMacro(LargeValue, PixelType);

  OutputRegionType                              m_BufferedRegion;
  typedef typename LevelSetImageType::IndexType LevelSetIndexType;
  LevelSetIndexType                             m_StartIndex;
  LevelSetIndexType                             m_LastIndex;

  itkGetConstReferenceMacro(StartIndex, LevelSetIndexType);
  itkGetConstReferenceMacro(LastIndex, LevelSetIndexType);

private:
  FMarchingImageFilter(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented

  NodeContainerPointer m_AlivePoints;
  NodeContainerPointer m_TrialPoints;

  LabelImagePointer              m_LabelImage;
  ConnectedComponentImagePointer m_ConnectedComponentImage;

  double m_SpeedConstant;
  double m_InverseSpeed;
  double m_StoppingValue;

  bool                 m_CollectPoints;
  NodeContainerPointer m_ProcessedPoints;

  OutputRegionType    m_OutputRegion;
  OutputPointType     m_OutputOrigin;
  OutputSpacingType   m_OutputSpacing;
  OutputDirectionType m_OutputDirection;
  bool                m_OverrideOutputInformation;

  typename LevelSetImageType::PixelType m_LargeValue;
  AxisNodeType                          m_NodesUsed[SetDimension];

  /** Trial points are stored in a min-heap. This allow efficient access
   * to the trial point with minimum value which is the next grid point
   * the algorithm processes. */
  typedef std::vector<AxisNodeType>                                      HeapContainer;
  typedef std::greater<AxisNodeType>                                     NodeComparer;
  typedef std::priority_queue<AxisNodeType, HeapContainer, NodeComparer> HeapType;

  HeapType m_TrialHeap;

  double m_NormalizationFactor;

  /**
   * Functions and variables to check for topology changes (2D/3D only).
   */
  TopologyCheckType m_TopologyCheck;

  // Functions/data for the 2-D case
  void
  InitializeIndices2D();

  bool IsChangeWellComposed2D(IndexType);
  bool
  IsCriticalC1Configuration2D(Array<short>);

  bool
  IsCriticalC2Configuration2D(Array<short>);

  bool
  IsCriticalC3Configuration2D(Array<short>);

  bool
  IsCriticalC4Configuration2D(Array<short>);

  bool IsSpecialCaseOfC4Configuration2D(PixelType, IndexType, IndexType, IndexType);

  Array<unsigned int> m_RotationIndices[4];
  Array<unsigned int> m_ReflectionIndices[2];

  // Functions/data for the 3-D case
  void
  InitializeIndices3D();

  bool
  IsCriticalC1Configuration3D(Array<short>);

  unsigned int
  IsCriticalC2Configuration3D(Array<short>);

  bool IsChangeWellComposed3D(IndexType);

  Array<unsigned int> m_C1Indices[12];
  Array<unsigned int> m_C2Indices[8];

  // Functions for both 2D/3D cases
  bool DoesVoxelChangeViolateWellComposedness(IndexType);
  bool DoesVoxelChangeViolateStrictTopology(IndexType);
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "antsFastMarchingImageFilter.hxx"
#endif

#endif
