/********************************************************************************
 * Create a VTK mesh from the white matter image
 *******************************************************************************/
#ifndef __BinaryImageToMeshFilter_h_
#define __BinaryImageToMeshFilter_h_
#include <vtkSmartPointer.h>

#include "itkImage.h"
#include "itkResampleImageFilter.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkAntiAliasBinaryImageFilter.h"
#include "itkUnaryFunctorImageFilter.h"
// #include "itkBinaryWellComposed3DImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkCommand.h"

#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkImageMarchingCubes.h>
#include <vtkImageImport.h>
#include <itkVTKImageExport.h>
#include <vtkDecimatePro.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkStripper.h>
#include <vtkCleanPolyData.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkImageData.h>

#include <iostream>

using namespace std;

template <typename TImage>
void
ConnectITKToVTK(itk::VTKImageExport<TImage> * fltExport, vtkImageImport * fltImport)
{
  fltImport->SetUpdateInformationCallback(fltExport->GetUpdateInformationCallback());
  fltImport->SetPipelineModifiedCallback(fltExport->GetPipelineModifiedCallback());
  fltImport->SetWholeExtentCallback(fltExport->GetWholeExtentCallback());
  fltImport->SetSpacingCallback(fltExport->GetSpacingCallback());
  fltImport->SetOriginCallback(fltExport->GetOriginCallback());
  fltImport->SetScalarTypeCallback(fltExport->GetScalarTypeCallback());
  fltImport->SetNumberOfComponentsCallback(fltExport->GetNumberOfComponentsCallback());
  fltImport->SetPropagateUpdateExtentCallback(fltExport->GetPropagateUpdateExtentCallback());
  fltImport->SetUpdateDataCallback(fltExport->GetUpdateDataCallback());
  fltImport->SetDataExtentCallback(fltExport->GetDataExtentCallback());
  fltImport->SetBufferPointerCallback(fltExport->GetBufferPointerCallback());
  fltImport->SetCallbackUserData(fltExport->GetCallbackUserData());
}

class UnaryFunctorBinaryToFloat
{
public:
  UnaryFunctorBinaryToFloat() { m_InsideValue = 1.0f; }

  void
  SetInvertImage(bool invert)
  {
    m_InsideValue = invert ? -1.0f : 1.0f;
  }

  inline float
  operator()(short in)
  {
    return in == 0 ? -m_InsideValue : m_InsideValue;
  }

  inline bool
  operator!=(const UnaryFunctorBinaryToFloat & otro) const
  {
    return (!itk::Math::FloatAlmostEqual(this->m_InsideValue, otro.m_InsideValue));
  }

private:
  float m_InsideValue;
};

template <typename TImage>
class BinaryImageToMeshFilter final : public itk::ProcessObject
{
public:
  typedef BinaryImageToMeshFilter       Self;
  typedef itk::ProcessObject            Superclass;
  typedef itk::SmartPointer<Self>       Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;
  typedef itk::Image<float, 3>          FloatImageType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Image dimension. */
  static constexpr unsigned int ImageDimension = TImage::ImageDimension;

  /** Get the result mesh */
  vtkPolyData *
  GetMesh()
  {
    return fltTriangle->GetOutput();
  }

  /** Get the intermediate antialiased image */
  FloatImageType *
  GetAntiAliasImage()
  {
    return fltAlias->GetOutput();
  }

  /** Whether to invert the binary image */
  itkSetMacro(InvertInput, bool);
  itkGetMacro(InvertInput, bool);

  /** Whether to invert the binary image */
  itkSetMacro(ResampleScaleFactor, float);
  itkGetMacro(ResampleScaleFactor, float);

  /** Whether to decimate the mesh and how much, 0 for none*/
  itkSetMacro(DecimateFactor, double);
  itkGetMacro(DecimateFactor, double);

  /** Smoothing iterations, or 0 for none */
  itkSetMacro(SmoothingIterations, int);
  itkGetMacro(SmoothingIterations, int);

  /** Set the input */
  using Superclass::SetInput;
  virtual void
  SetInput(TImage * image)
  {
    this->SetNthInput(0, image);
  }

  /** Update method (why?) */
  void
  Update() override
  {
    this->GenerateData();
  }

  /** Set the anti-aliasing quality parameter */
  void
  SetAntiAliasMaxRMSError(double value)
  {
    fltAlias->SetMaximumRMSError(value);
  }

  /** Get the 'distance image' based on anti-aliasing the binary image */
  FloatImageType *
  GetDistanceImage()
  {
    return fltAlias->GetOutput();
  }

  /** Get the floating point 'resampled' image, if resampling was enabled */
  FloatImageType *
  GetResampledImage()
  {
    if (itk::Math::FloatAlmostEqual(m_ResampleScaleFactor, 1.0f))
    {
      return fltResample->GetOutput();
    }
    else
    {
      return fltToFloat->GetOutput();
    }
  }

  void
  PrintMeshStatistics(vtkPolyData * mesh)
  {
    vector<unsigned int> cellHist(100, 0);
    for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); i++)
    {
      cellHist[mesh->GetCellType(i)]++;
    }

    cout << "        mesh has " << mesh->GetNumberOfPoints() << " points." << endl;
    cout << "        mesh has " << mesh->GetNumberOfCells() << " cells. " << endl;
    cout << "        mesh has " << cellHist[VTK_VERTEX] << " vtk_vertex" << endl;
    cout << "        mesh has " << cellHist[VTK_LINE] << " vtk_line" << endl;
    cout << "        mesh has " << cellHist[VTK_POLY_LINE] << " vtk_poly_line" << endl;
    cout << "        mesh has " << cellHist[VTK_TRIANGLE] << " vtk_triangle" << endl;
    cout << "        mesh has " << cellHist[VTK_TRIANGLE_STRIP] << " vtk_triangle_strip" << endl;
  }

protected:
  BinaryImageToMeshFilter()
  {
    // Set the cardinality of the filter
    this->SetNumberOfIndexedInputs(1);
    this->SetNumberOfIndexedOutputs(1);

    // Begin with the well-connectedness filter
    // fltTopology = TopologyFilter::New();

    // Create the converter to float
    fltToFloat = ToFloatFilter::New();
    //    fltToFloat->SetInput(this->GetInput());//fltTopology->GetOutput());
    typename FloatImageType::Pointer imgPipeEnd = fltToFloat->GetOutput();

    // Initialize the interpolation function
    fnInterpolate = ResampleFunction::New();

    // Create a resampling image filter
    fltResample = ResampleFilter::New();
    fltResample->SetInput(imgPipeEnd);
    fltResample->SetTransform(itk::IdentityTransform<double, 3>::New());
    fltResample->SetInterpolator(fnInterpolate);
    imgPipeEnd = fltResample->GetOutput();

    // Create an anti-aliasing image filter
    fltAlias = AAFilter::New();
    fltAlias->SetMaximumRMSError(0.024);
    fltAlias->SetInput(imgPipeEnd);
    imgPipeEnd = fltAlias->GetOutput();

    // Cast the image to VTK
    fltExport = ExportFilter::New();
    fltImport = vtkImageImport::New();
    fltExport->SetInput(imgPipeEnd);
    ConnectITKToVTK(fltExport.GetPointer(), fltImport);

    // Compute marching cubes
    vtkImageData * importPipeEnd = fltImport->GetOutput();
    fltMarching = vtkImageMarchingCubes::New();
    fltMarching->SetInputData(importPipeEnd);
    fltMarching->ComputeScalarsOff();
    fltMarching->ComputeGradientsOff();
    fltMarching->SetNumberOfContours(1);
    fltMarching->SetValue(0, 0.0f);
    vtkPolyData * meshPipeEnd = fltMarching->GetOutput();

    // Keep the largest connected component
    fltConnect = vtkPolyDataConnectivityFilter::New();
    fltConnect->SetInputData(meshPipeEnd);
    fltConnect->SetExtractionModeToLargestRegion();
    meshPipeEnd = fltMarching->GetOutput();

    fltClean = vtkCleanPolyData::New();
    bool doclean = false;
    if (doclean)
    {
      // Clean up the data
      fltClean->SetInputData(meshPipeEnd);
      meshPipeEnd = fltClean->GetOutput();
      // Decimate the data
      fltDecimate = vtkDecimatePro::New();
      fltDecimate->SetInputData(meshPipeEnd);
      fltDecimate->PreserveTopologyOn();
      meshPipeEnd = fltDecimate->GetOutput();
    }

    // Smooth the data, keeping it on the contour
    vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smoothFilter->SetInputData(meshPipeEnd);
    smoothFilter->SetNumberOfIterations(5);
    meshPipeEnd = smoothFilter->GetOutput();

    // Compute triangle strips for faster display
    fltTriangle = vtkTriangleFilter::New();
    fltTriangle->SetInputData(meshPipeEnd);
    meshPipeEnd = fltTriangle->GetOutput();

    // Set up progress
    typedef itk::MemberCommand<Self> CommandType;
    typename CommandType::Pointer    cmd = CommandType::New();
    cmd->SetCallbackFunction(this, &Self::ProgressCommand);

    // Add progress to the two slow filters
    fltResample->AddObserver(itk::ProgressEvent(), cmd);
    fltAlias->AddObserver(itk::ProgressEvent(), cmd);

    // Invert - NO
    m_InvertInput = false;

    // Resample - NO
    m_ResampleScaleFactor = 1.0f;

    // Decimate - NO
    m_DecimateFactor = 0.0f;
  }

  ~BinaryImageToMeshFilter() override
  {
    // CLean up
    fltMarching->Delete();
    fltConnect->Delete();
    fltImport->Delete();
    fltTriangle->Delete();
    fltClean->Delete();
    fltDecimate->Delete();
    fltSmoothMesh->Delete();
  }

  /** Generate Data */
  void
  GenerateData(void) override
  {
    // Run the computation
    cout << "Computing mesh from binary image" << endl;
    // Get the input and output pointers
    typename TImage::ConstPointer inputImage = reinterpret_cast<TImage *>(this->GetInput(0));

    // Pass the input to the topology filter
    // fltTopology->SetInputData(inputImage);

    // Compute the max/min of the image to set fore/back
    typedef itk::MinimumMaximumImageCalculator<TImage> CalcType;
    typename CalcType::Pointer                         calc = CalcType::New();
    calc->SetImage(inputImage);
    calc->Compute();

    // Set the background and foreground in the topology filter
    // fltTopology->SetBackgroundValue(calc->GetMinimum());
    // fltTopology->SetForegroundValue(calc->GetMaximum());
    // fltTopology->Update();

    // Pass the input to the remapper
    fltToFloat->SetInput(inputImage);

    // Set the inversion if necessary
    m_ToFloatFunctor.SetInvertImage(m_InvertInput);
    fltToFloat->SetFunctor(m_ToFloatFunctor);

    // Convert the image to floats
    fltToFloat->Update();

    // Peform resampling only if necessary
    if (!itk::Math::FloatAlmostEqual(m_ResampleScaleFactor, 1.0f))
    {
      // Include the filter in the pipeline
      fltAlias->SetInput(fltResample->GetOutput());

      // Set the size parameter
      FloatImageType::SizeType szOutput = inputImage->GetBufferedRegion().GetSize();
      szOutput[0] = (unsigned long)(szOutput[0] * m_ResampleScaleFactor);
      szOutput[1] = (unsigned long)(szOutput[1] * m_ResampleScaleFactor);
      szOutput[2] = (unsigned long)(szOutput[2] * m_ResampleScaleFactor);
      fltResample->SetSize(szOutput);

      // Set the scale and origin
      FloatImageType::SpacingType xSpacing = inputImage->GetSpacing();
      xSpacing[0] /= static_cast<double>(m_ResampleScaleFactor);
      xSpacing[1] /= static_cast<double>(m_ResampleScaleFactor);
      xSpacing[2] /= static_cast<double>(m_ResampleScaleFactor);
      fltResample->SetOutputSpacing(xSpacing);
      fltResample->SetOutputOrigin(inputImage->GetOrigin());

      // Compute the resampling
      cout << "   resampling the image " << endl;
      fltResample->Update();
      cout << endl;
    }
    else
    {
      // Exclude the filter from the pipeline
      fltAlias->SetInput(fltToFloat->GetOutput());
    }

    // Run the filters
    if (fltAlias->GetMaximumRMSError() > 0.0)
    {
      cout << "   anti-aliasing the image " << endl;
      fltAlias->Update();
      fltExport->SetInput(fltAlias->GetOutput());
    }
    else
    {
      fltExport->SetInput(fltAlias->GetInput());
    }
    bool verbose = true;
    if (verbose)
    {
      cout << endl << "   converting image to VTK" << endl;
    }
    fltImport->Update();

    if (verbose)
    {
      cout << "   running marching cubes algorithm" << endl;
    }
    fltMarching->Update();

    if (verbose)
    {
      cout << "      mesh has " << fltMarching->GetOutput()->GetNumberOfCells() << " cells and "
           << fltMarching->GetOutput()->GetNumberOfPoints() << " points. " << endl;
    }

    if (verbose)
    {
      cout << "   extracting the largest component" << endl;
    }
    fltConnect->Update();

    if (verbose)
    {
      cout << "      mesh has " << fltConnect->GetOutput()->GetNumberOfCells() << " cells and "
           << fltConnect->GetOutput()->GetNumberOfPoints() << " points. " << endl;
    }

    bool doclean = false;
    if (doclean)
    {
      if (verbose)
        cout << "   cleaning the mesh " << endl;
      fltClean->Update();
      if (verbose)
        cout << "  after clean step mesh uses " << m_DecimateFactor << " decimation "
             << fltClean->GetOutput()->GetNumberOfCells() << " cells and " << fltClean->GetOutput()->GetNumberOfPoints()
             << " points. " << endl;
    }

    std::cout << "  to decimation with factor " << m_DecimateFactor << std::endl;
    // If decimation is on, run it
    if (m_DecimateFactor > 0.0)
    {
      if (verbose)
      {
        cout << "   decimating the mesh by factor of " << m_DecimateFactor << endl;
      }
      fltDecimate->SetTargetReduction(m_DecimateFactor);
      fltDecimate->SetInputData(fltConnect->GetOutput());
      fltDecimate->Update();
      fltTriangle->SetInputData(fltDecimate->GetOutput());
      //      if (verbose) cout << "      mesh has "
      // << fltClean->GetOutput()->GetNumberOfCells() << " cells and "
      //  << fltClean->GetOutput()->GetNumberOfPoints() << " points. " << endl;
    }
    else
    {
      fltTriangle->SetInputData(fltConnect->GetOutput());
    }

    if (verbose)
    {
      cout << "   converting mesh to triangles" << endl;
    }
    fltTriangle->Update();
    m_Result = fltTriangle->GetOutput();

    if (verbose)
    {
      cout << "      mesh has " << m_Result->GetNumberOfCells() << " cells and " << m_Result->GetNumberOfPoints()
           << " points. " << endl;
    }

    // If smoothing is on, run it
    if (m_SmoothingIterations > 0)
    {
      if (verbose)
      {
        cout << "   smoothing the mesh " << m_SmoothingIterations << endl;
      }
      fltSmoothMesh->SetNumberOfIterations(m_SmoothingIterations);
      fltSmoothMesh->SetInputData(m_Result);
      //      fltSmoothMesh->SetInputConnection(m_Result); // BA FIXME
      fltSmoothMesh->Update();
      m_Result = fltSmoothMesh->GetOutput();
      std::cout << " Done " << std::endl;
    }
  }

private:
  typedef itk::ImageRegionIterator<TImage>      IteratorType;
  typedef itk::ImageRegionConstIterator<TImage> ConstIteratorType;

  // Topology correction filter
  // typedef itk::BinaryWellComposed3DImageFilter<TImage> TopologyFilter;

  // Functor for remapping to float
  UnaryFunctorBinaryToFloat m_ToFloatFunctor;

  // Filter to remap image to floating point
  typedef itk::UnaryFunctorImageFilter<TImage, FloatImageType, UnaryFunctorBinaryToFloat> ToFloatFilter;

  // Windowed sinc for resampling
  typedef itk::Function::WelchWindowFunction<4>                                WindowFunction;
  typedef itk::NearestNeighborInterpolateImageFunction<FloatImageType, double> ResampleFunction;

  // Filter to resample image
  typedef itk::ResampleImageFilter<FloatImageType, FloatImageType> ResampleFilter;

  // Antialiasing filter
  typedef itk::AntiAliasBinaryImageFilter<FloatImageType, FloatImageType> AAFilter;

  // Export to VTK filter
  typedef itk::VTKImageExport<FloatImageType> ExportFilter;

  // typename TopologyFilter::Pointer fltTopology;
  typename AAFilter::Pointer         fltAlias;
  typename ExportFilter::Pointer     fltExport;
  typename ToFloatFilter::Pointer    fltToFloat;
  typename ResampleFilter::Pointer   fltResample;
  typename ResampleFunction::Pointer fnInterpolate;

  vtkImageImport *                fltImport;
  vtkImageMarchingCubes *         fltMarching;
  vtkPolyDataConnectivityFilter * fltConnect;
  vtkCleanPolyData *              fltClean;
  vtkTriangleFilter *             fltTriangle;
  vtkDecimatePro *                fltDecimate;
  vtkSmoothPolyDataFilter *       fltSmoothMesh;

  bool   m_InvertInput;
  float  m_ResampleScaleFactor;
  double m_DecimateFactor;
  int    m_SmoothingIterations;

  vtkPolyData * m_Result;

  void
  ProgressCommand(itk::Object * source, const itk::EventObject & /*evt*/)
  {
    // Get the elapsed progress
    itk::ProcessObject * po = reinterpret_cast<ProcessObject *>(source);
    float                progress = po->GetProgress();

    cout << setprecision(4) << (100 * progress) << "  " << flush;
  }
};

#endif
