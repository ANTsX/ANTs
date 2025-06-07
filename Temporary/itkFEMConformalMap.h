/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _FEMConformalMap_h
#define _FEMConformalMap_h

#include <vnl/vnl_matops.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

#include "vtkDataSetWriter.h"
#include "vtkDataSetMapper.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkActor.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkDataSetReader.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSet.h"
#include "vtkCellArray.h"
#include "vtkVolume16Reader.h"
#include "vtkImageReader2.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkOutlineFilter.h"
#include "vtkCamera.h"
#include "vtkProperty.h"
#include "vtkPolyData.h"
#include "vtkPolyVertex.h"
#include "vtkPointData.h"
#include "vtkExtractEdges.h"
#include "vtkPolyDataNormals.h"
#include "vtkMarchingCubes.h"
#include "vtkImageGaussianSmooth.h"
#include "vtkDecimatePro.h"
#include "vtkContourFilter.h"
#include "vtkPolyDataConnectivityFilter.h"
// #include "vtkKitwareContourFilter.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkSTLWriter.h"
#include "vtkUnstructuredGridToPolyDataFilter.h"
// #include "itkImageToVTKImageFilter.h"
#include "vtkDelaunay2D.h"
#include "vtkFloatArray.h"

#include "itkObject.h"
#include "itkProcessObject.h"

#include "itkVectorContainer.h"
#include "itkCastImageFilter.h"

#include "itkFEM.h"
#include "itkFEMLinearSystemWrapperItpack.h"
#include "itkFEMElement3DC0LinearTriangularLaplaceBeltrami.h"
#include "itkFEMElement3DC0LinearTriangularMembrane.h"

namespace itk
{
/** \class FEMConformalMap
 * Angenent, Haker conformal mapping algorithm using FEM.
 *
 * \note The origin of a neighborhood is always taken to be
 *       the first point entered into and the
 *       last point stored in the list.
 */
template <typename TSurface, typename TImage, unsigned int TDimension = 3>
class FEMConformalMap : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef FEMConformalMap          Self;
  typedef ProcessObject            Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(FEMConformalMap);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Surface (mesh) types. */
  typedef TImage                   ImageType;
  typedef typename TImage::Pointer ImageTypePointer;

  /** Surface (mesh) types. */
  typedef TSurface                                       SurfaceType;
  typedef typename SurfaceType::Pointer                  SurfaceTypePointer;
  typedef typename SurfaceType::PointType                PointType;
  typedef typename SurfaceType::CellsContainerPointer    InputCellsContainerPointer;
  typedef typename SurfaceType::CellsContainer::Iterator InputCellsContainerIterator;

  /** Image dimension. */
  //  static constexpr unsigned int ImageDimension = TImage::ImageDimension;
  static constexpr unsigned int ImageDimension = TDimension;
  static constexpr unsigned int SurfaceDimension = TDimension;

  typedef double                                                             RealType;
  typedef vnl_vector<RealType>                                               VectorType;
  typedef vnl_vector_fixed<RealType, Self::ImageDimension> FixedVectorType;
  typedef vnl_matrix<double>                                                 MatrixType;

  /** FEM types */
  typedef itk::fem::MaterialLinearElasticity                   MaterialType;
  typedef itk::fem::Node                                       NodeType;
  typedef itk::fem::LoadNode                                   LoadType;
  typedef itk::fem::Element3DC0LinearTriangularLaplaceBeltrami ElementType;
  typedef itk::fem::Element3DC0LinearTriangularMembrane        ElementType1;

  /** Set input parameter file */
  itkSetStringMacro(ParameterFileName);

  /** Set input parameter file */
  itkGetStringMacro(ParameterFileName);

  itkGetMacro(Sigma, RealType);
  itkSetMacro(Sigma, RealType);

  itkGetMacro(SurfaceMesh, SurfaceTypePointer);
  itkSetMacro(SurfaceMesh, SurfaceTypePointer);

  itkGetMacro(Image, ImageTypePointer);
  itkSetMacro(Image, ImageTypePointer);
  itkGetMacro(SphereImage, ImageTypePointer);
  itkSetMacro(SphereImage, ImageTypePointer);
  itkSetMacro(SouthPole, int);
  void
  SetNorthPole(int p)
  {
    if (p % 2 == 0)
    {
      m_NorthPole = p;
      m_SouthPole = p + 1;
    }
    else
    {
      m_NorthPole = p - 1;
      m_SouthPole = p - 1;
    }
  }

  void
  SetDebug(bool b)
  {
    m_Debug = b;
  }

  void
  SetReadFromFile(bool b)
  {
    m_ReadFromFile = b;
  }

  void
  FindPoles(int dim);

  void
  FixPoles(int dim);

  void
  ConformalParameterize();

  void
  ConformalMap();

  void
  ComputeStereographicCoordinates();

  void
  MapStereographicCoordinatesToImage(int dim);

  void
  MapImageToSphere(ImageType * img, float rad);

  void
  MapCheckerboardToImage(float increment);

  void
  BuildOutputMeshes(typename TImage::Pointer image);

  vtkPolyData * m_ExtractedSurfaceMesh;
  vtkPolyData * m_VtkSurfaceMesh;

protected:
  bool
  GenerateSystemFromSurfaceMesh();

  bool
  GenerateSystemFromVtkSurfaceMesh();

  void
  ApplyRealForces(int dim);

  void
  ApplyImaginaryForces(int dim);

  FEMConformalMap();
  virtual ~FEMConformalMap(){};

private:
  RealType    m_Sigma;
  RealType    m_Pi;
  std::string m_ParameterFileName;
  int         m_NorthPole;
  int         m_SouthPole;

  itk::fem::Solver m_Solver;

  bool m_ReadFromFile;
  bool m_Debug;
  bool m_FindingRealSolution;

  VectorType m_RealSolution;
  VectorType m_ImagSolution;

  ImageTypePointer                    m_Image;
  ImageTypePointer                    m_SphereImage;
  SurfaceTypePointer                  m_SurfaceMesh;
  itk::fem::LinearSystemWrapperItpack itpackWrapper;

  unsigned long m_PoleElementsGN[7];
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkFEMConformalMap.cxx"
#endif

#endif
