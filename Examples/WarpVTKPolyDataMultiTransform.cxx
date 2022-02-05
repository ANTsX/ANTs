#include "antsUtilities.h"
#include "antsUtilities.h"
#include <vnl/vnl_inverse.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkTransformFactory.h"
#include "vtkPolyDataReader.h"
#include "itkDisplacementFieldFromMultiTransformFilter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"

#include "itkLabeledPointSetFileReader.h"
#include "itkLabeledPointSetFileWriter.h"
#include "itkMesh.h"
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPoints.h>

namespace ants
{
vnl_matrix_fixed<double, 4, 4>
ConstructNiftiSform(vnl_matrix<double> m_dir, vnl_vector<double> v_origin, vnl_vector<double> v_spacing)
{
  // Set the NIFTI/RAS transform
  vnl_matrix<double>      m_ras_matrix;
  vnl_diag_matrix<double> m_scale, m_lps_to_ras;
  vnl_vector<double>      v_ras_offset;

  // Compute the matrix
  m_scale.set(v_spacing);
  m_lps_to_ras.set(vnl_vector<double>(3, 1.0));
  m_lps_to_ras[0] = -1;
  m_lps_to_ras[1] = -1;
  m_ras_matrix = m_lps_to_ras * m_dir * m_scale;

  // Compute the vector
  v_ras_offset = m_lps_to_ras * v_origin;

  // Create the larger matrix
  vnl_vector<double> vcol(4, 1.0);
  vcol.update(v_ras_offset);

  vnl_matrix_fixed<double, 4, 4> m_sform;
  m_sform.set_identity();
  m_sform.update(m_ras_matrix);
  m_sform.set_column(3, vcol);
  return m_sform;
}

vnl_matrix_fixed<double, 4, 4>
ConstructVTKtoNiftiTransform(vnl_matrix<double> m_dir, vnl_vector<double> v_origin, vnl_vector<double> v_spacing)
{
  vnl_matrix_fixed<double, 4, 4> vox2nii = ConstructNiftiSform(m_dir, v_origin, v_spacing);
  vnl_matrix_fixed<double, 4, 4> vtk2vox;
  vtk2vox.set_identity();
  for (size_t i = 0; i < 3; i++)
  {
    vtk2vox(i, i) = 1.0 / v_spacing[i];
    vtk2vox(i, 3) = -v_origin[i] / v_spacing[i];
  }
  return vox2nii * vtk2vox;
}

static bool
WarpVTKPolyDataMultiTransform_ParseInput(int              argc,
                                         char **          argv,
                                         char *&          input_vtk_filename,
                                         char *&          output_vtk_filename,
                                         char *&          reference_image_filename,
                                         TRAN_OPT_QUEUE & opt_queue)
{
  opt_queue.clear();
  opt_queue.reserve(argc - 2);

  input_vtk_filename = argv[0];
  output_vtk_filename = argv[1];

  reference_image_filename = nullptr;

  int ind = 2;
  while (ind < argc)
  {
    if (strcmp(argv[ind], "-R") == 0)
    {
      ind++;
      if (ind >= argc)
      {
        return false;
      }
      reference_image_filename = argv[ind];
    }
    else if (strcmp(argv[ind], "-i") == 0)
    {
      ind++;
      if (ind >= argc)
      {
        return false;
      }
      TRAN_OPT opt;
      opt.filename = argv[ind];
      if (CheckFileType(opt.filename) != AFFINE_FILE)
      {
        std::cout << "file: " << opt.filename << " is not an affine .txt file. Invalid to use '-i' " << std::endl;
        return false;
      }
      opt.file_type = AFFINE_FILE;
      opt.do_affine_inv = true;
      opt_queue.push_back(opt);
    }
    else
    {
      TRAN_OPT opt;
      opt.filename = argv[ind];
      opt.file_type = CheckFileType(opt.filename);
      opt.do_affine_inv = false;
      opt_queue.push_back(opt);
    }
    ind++;
  }

  //    if (reference_image_filename == nullptr) {
  //        std::cout << "the reference image file (-R) must be given!!!"
  //        << std::endl;
  //        return false;
  //    }

  return true;
}

template <int ImageDimension>
void
WarpLabeledPointSetFileMultiTransform(char *           input_vtk_filename,
                                      char *           output_vtk_filename,
                                      char *           reference_image_filename,
                                      TRAN_OPT_QUEUE & opt_queue)
{
  using ImageType = itk::Image<float, ImageDimension>;
  using VectorType = itk::Vector<float, ImageDimension>;
  using DisplacementFieldType = itk::Image<VectorType, ImageDimension>;
  using AffineTransformType = itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension>;
  // typedef itk::WarpImageMultiTransformFilter<ImageType,ImageType, DisplacementFieldType, AffineTransformType>
  // WarperType;
  using WarperType =
    itk::DisplacementFieldFromMultiTransformFilter<DisplacementFieldType, DisplacementFieldType, AffineTransformType>;
  using FuncType = itk::LinearInterpolateImageFunction<ImageType>;

  itk::TransformFactory<AffineTransformType>::RegisterTransform();

  using ImageFileReaderType = itk::ImageFileReader<ImageType>;
  typename ImageFileReaderType::Pointer reader_img = ImageFileReaderType::New();
  typename ImageType::Pointer           img_ref;

  typename ImageFileReaderType::Pointer reader_img_ref = ImageFileReaderType::New();
  if (reference_image_filename)
  {
    reader_img_ref->SetFileName(reference_image_filename);
    reader_img_ref->Update();
    img_ref = reader_img_ref->GetOutput();
  }
  else
  {
    std::cout << "the reference image file (-R) must be given!!!" << std::endl;
    return;
  }

  typename WarperType::Pointer warper = WarperType::New();
  // warper->SetInput(img_mov);
  // warper->SetEdgePaddingValue( 0);
  VectorType pad;
  pad.Fill(0);
  // warper->SetEdgePaddingValue(pad);

  using TranReaderType = itk::TransformFileReader;

  using FieldReaderType = itk::ImageFileReader<DisplacementFieldType>;

  const int kOptQueueSize = opt_queue.size();
  for (int i = 0; i < kOptQueueSize; i++)
  {
    const TRAN_OPT & opt = opt_queue[i];

    switch (opt_queue[i].file_type)
    {
      case AFFINE_FILE:
      {
        typename TranReaderType::Pointer tran_reader = TranReaderType::New();
        tran_reader->SetFileName(opt.filename);
        tran_reader->Update();
        typename AffineTransformType::Pointer aff =
          dynamic_cast<AffineTransformType *>((tran_reader->GetTransformList())->front().GetPointer());
        if (opt_queue[i].do_affine_inv)
        {
          aff->GetInverse(aff);
        }
        // std::cout << aff << std::endl;
        warper->PushBackAffineTransform(aff);
        break;
      }
      case DEFORMATION_FILE:
      {
        typename FieldReaderType::Pointer field_reader = FieldReaderType::New();
        field_reader->SetFileName(opt.filename);
        field_reader->Update();
        typename DisplacementFieldType::Pointer field = field_reader->GetOutput();
        // std::cout << field << std::endl;
        warper->PushBackDisplacementFieldTransform(field);
        break;
      }
      default:
        std::cout << "Unknown file type!" << std::endl;
    }
  }

  warper->SetOutputParametersFromImage(img_ref);

  std::cout << "output size: " << warper->GetOutputSize() << std::endl;
  std::cout << "output spacing: " << warper->GetOutputSpacing() << std::endl;

  // warper->PrintTransformList();
  warper->DetermineFirstDeformNoInterp();
  warper->Update();

  typename DisplacementFieldType::Pointer field_output = DisplacementFieldType::New();
  field_output = warper->GetOutput();

  /**
   * Code to read the mesh
   */
  typename ImageType::PointType point;
  typename ImageType::PointType warpedPoint;

  using MeshType = itk::Mesh<float, ImageDimension>;
  //  typedef itk::LabeledPointSetFileReader<MeshType> VTKReaderType;
  vtkPolyDataReader * vtkreader = vtkPolyDataReader::New();
  vtkreader->SetFileName(input_vtk_filename);
  vtkreader->Update();
  vtkPolyData * mesh = vtkreader->GetOutput();

  /*
  typename VTKReaderType::Pointer vtkreader = VTKReaderType::New();
  vtkreader->SetFileName( input_vtk_filename );
  vtkreader->Update();
  typename MeshType::PointsContainerIterator It
    = vtkreader->GetOutput()->GetPoints()->Begin();

  while( It != vtkreader->GetOutput()->GetPoints()->End() )
    {
    point.CastFrom( It.Value() );
  */
  vnl_matrix_fixed<double, 4, 4> ijk2ras, vtk2ras, lps2ras;
  vtk2ras = ConstructVTKtoNiftiTransform(field_output->GetDirection().GetVnlMatrix(),
                                         field_output->GetOrigin().GetVnlVector(),
                                         field_output->GetSpacing().GetVnlVector());
  ijk2ras.set_identity();
  vtk2ras.set_identity();
  lps2ras.set_identity();
  lps2ras(0, 0) = -1;
  lps2ras(1, 1) = -1;

  // Set up the transforms
  ijk2ras = ConstructNiftiSform(field_output->GetDirection().GetVnlMatrix(),
                                field_output->GetOrigin().GetVnlVector(),
                                field_output->GetSpacing().GetVnlVector());

  vtk2ras = ConstructVTKtoNiftiTransform(field_output->GetDirection().GetVnlMatrix(),
                                         field_output->GetOrigin().GetVnlVector(),
                                         field_output->GetSpacing().GetVnlVector());

  vnl_matrix_fixed<double, 4, 4> ras2ijk = vnl_inverse(ijk2ras);
  // vnl_matrix_fixed<double, 4, 4> ras2vtk = vnl_inverse(vtk2ras);
  // vnl_matrix_fixed<double, 4, 4> ras2lps = vnl_inverse(lps2ras);
  // Update the coordinates
  for (int k = 0; k < mesh->GetNumberOfPoints(); k++)
  {
    // Get the point (in whatever format that it's stored)
    vnl_vector_fixed<double, 4> x_mesh, x_ras, x_ijk;
    x_mesh[0] = mesh->GetPoint(k)[0];
    x_mesh[1] = mesh->GetPoint(k)[1];
    x_mesh[2] = mesh->GetPoint(k)[2];
    x_mesh[3] = 1.0;

    // Map the point into RAS coordinates
    //    if(parm.mesh_coord == RAS)  x_ras = x_mesh;
    //    else if(parm.mesh_coord == LPS) x_ras = lps2ras * x_mesh;
    //    else if(parm.mesh_coord == IJKOS)
    x_ras = vtk2ras * x_mesh;
    //    else   x_ras = ijk2ras * x_mesh;

    // Map the point to IJK coordinates (continuous index)
    x_ijk = ras2ijk * x_ras;

    typename FuncType::ContinuousIndexType ind;
    ind[0] = x_ijk[0];
    ind[1] = x_ijk[1];
    ind[2] = x_ijk[2];
    field_output->TransformContinuousIndexToPhysicalPoint(ind, point);
    //      std::cout << " point " << point << std::endl;
    // std::cout << " point-t " << point << std::endl;
    bool isInside = warper->MultiTransformSinglePoint(point, warpedPoint);
    // if ( isInside ) std::cout << " point-w " << warpedPoint << std::endl;
    if (isInside)
    {
      typename MeshType::PointType newPoint;
      newPoint.CastFrom(warpedPoint);
      //      vtkreader->GetOutput()->SetPoint( It.Index(), newPoint );
      if (ImageDimension == 3)
      {
        mesh->GetPoints()->SetPoint(k, warpedPoint[0], warpedPoint[1], warpedPoint[2]);
      }
    }
    //    ++It;
  }

  std::string fn = std::string(output_vtk_filename);
  if (fn.rfind(".vtk") == fn.length() - 4)
  {
    vtkPolyDataWriter * writer = vtkPolyDataWriter::New();
    writer->SetFileName(output_vtk_filename);
    writer->SetInputData(mesh);
    writer->Update();
  }

  /*
  typedef itk::LabeledPointSetFileWriter<MeshType> VTKWriterType;
  typename VTKWriterType::Pointer vtkwriter = VTKWriterType::New();
  vtkwriter->SetFileName( output_vtk_filename );
  vtkwriter->SetInput( vtkreader->GetOutput() );
  vtkwriter->SetMultiComponentScalars( vtkreader->GetMultiComponentScalars() );
  vtkwriter->SetLines( vtkreader->GetLines() );
  vtkwriter->Update();
  */
  //    std::string filePrefix = output_vtk_filename;
  //    std::string::size_type pos = filePrefix.rfind(".");
  //    std::string extension = std::string(filePrefix, pos, filePrefix.length()
  //            - 1);
  //    filePrefix = std::string(filePrefix, 0, pos);
  //
  //    std::cout << "output extension is: " << extension << std::endl;
  //
  //    if (extension != std::string(".mha")) {
  //        typedef itk::VectorImageFileWriter<DisplacementFieldType, ImageType>
  //        WriterType;
  //        typename WriterType::Pointer writer = WriterType::New();
  //        writer->SetFileName(output_vtk_filename);
  //        writer->SetUseAvantsNamingConvention(true);
  //        writer->SetInput(field_output);
  //        writer->Update();
  //    } else {
  //        typedef itk::ImageFileWriter<DisplacementFieldType> WriterType;
  //        typename WriterType::Pointer writer = WriterType::New();
  //        writer->SetFileName(output_vtk_filename);
  //        writer->SetInput(field_output);
  //        writer->Update();
  //    }
}

template <int ImageDimension>
void
ComposeMultiAffine(char * /*input_affine_txt*/,
                   char *           output_affine_txt,
                   char *           reference_affine_txt,
                   TRAN_OPT_QUEUE & opt_queue)
{
  using ImageType = itk::Image<float, ImageDimension>;
  using VectorType = itk::Vector<float, ImageDimension>;
  using DisplacementFieldType = itk::Image<VectorType, ImageDimension>;
  using AffineTransformType = itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension>;
  using WarperType =
    itk::WarpImageMultiTransformFilter<ImageType, ImageType, DisplacementFieldType, AffineTransformType>;
  // typedef itk::DisplacementFieldFromMultiTransformFilter<DisplacementFieldType,
  // DisplacementFieldType, AffineTransformType> WarperType;

  itk::TransformFactory<AffineTransformType>::RegisterTransform();

  // typedef itk::ImageFileReader<ImageType> ImageFileReaderType;
  // typename ImageFileReaderType::Pointer reader_img = ImageFileReaderType::New();
  // typename ImageType::Pointer img_ref = ImageType::New();

  // typename ImageFileReaderType::Pointer reader_img_ref = ImageFileReaderType::New();

  typename WarperType::Pointer warper = WarperType::New();
  // warper->SetInput(img_mov);
  // warper->SetEdgePaddingValue( 0);
  VectorType pad;
  pad.Fill(0);
  // warper->SetEdgePaddingValue(pad);

  using TranReaderType = itk::TransformFileReader;

  int       cnt_affine = 0;
  const int kOptQueueSize = opt_queue.size();
  for (int i = 0; i < kOptQueueSize; i++)
  {
    const TRAN_OPT & opt = opt_queue[i];

    switch (opt_queue[i].file_type)
    {
      case AFFINE_FILE:
      {
        typename TranReaderType::Pointer tran_reader = TranReaderType::New();
        tran_reader->SetFileName(opt.filename);
        tran_reader->Update();
        typename AffineTransformType::Pointer aff =
          dynamic_cast<AffineTransformType *>((tran_reader->GetTransformList())->front().GetPointer());
        if (opt_queue[i].do_affine_inv)
        {
          aff->GetInverse(aff);
        }
        // std::cout << aff << std::endl;
        warper->PushBackAffineTransform(aff);
        cnt_affine++;
        break;
      }
      case DEFORMATION_FILE:
      {
        std::cout << "Compose affine only files: ignore " << opt.filename << std::endl;
        break;
      }
      default:
        std::cout << "Unknown file type!" << std::endl;
    }
  }

  using PointType = typename AffineTransformType::CenterType;
  PointType aff_center;

  typename AffineTransformType::Pointer aff_ref_tmp;
  if (reference_affine_txt)
  {
    typename TranReaderType::Pointer tran_reader = TranReaderType::New();
    tran_reader->SetFileName(reference_affine_txt);
    tran_reader->Update();
    aff_ref_tmp = dynamic_cast<AffineTransformType *>((tran_reader->GetTransformList())->front().GetPointer());
  }
  else
  {
    if (cnt_affine > 0)
    {
      std::cout << "the reference affine file for center is selected as the first affine!" << std::endl;
      aff_ref_tmp = ((warper->GetTransformList()).begin())->second.aex.aff;
    }
    else
    {
      std::cout << "No affine input is given. nothing to do ......" << std::endl;
      return;
    }
  }

  aff_center = aff_ref_tmp->GetCenter();
  std::cout << "new center is : " << aff_center << std::endl;

  // warper->PrintTransformList();

  // typename AffineTransformType::Pointer aff_output = warper->ComposeAffineOnlySequence(aff_center);
  typename AffineTransformType::Pointer aff_output = AffineTransformType::New();
  warper->ComposeAffineOnlySequence(aff_center, aff_output);
  using TranWriterType = itk::TransformFileWriter;
  typename TranWriterType::Pointer tran_writer = TranWriterType::New();
  tran_writer->SetFileName(output_affine_txt);
  tran_writer->SetInput(aff_output);
#if ITK_VERSION_MAJOR >= 5
  tran_writer->SetUseCompression(true);
#endif
  tran_writer->Update();

  std::cout << "wrote file to : " << output_affine_txt << std::endl;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
WarpVTKPolyDataMultiTransform(std::vector<std::string> args, std::ostream *)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "WarpVTKPolyDataMultiTransform");

  int     argc = args.size();
  char ** argv = new char *[args.size() + 1];
  for (unsigned int i = 0; i < args.size(); ++i)
  {
    // allocate space for the string plus a null character
    argv[i] = new char[args[i].length() + 1];
    std::strncpy(argv[i], args[i].c_str(), args[i].length());
    // place the null character in the end
    argv[i][args[i].length()] = '\0';
  }
  argv[argc] = nullptr;
  // class to automatically cleanup argv upon destruction
  class Cleanup_argv
  {
  public:
    Cleanup_argv(char ** argv_, int argc_plus_one_)
      : argv(argv_)
      , argc_plus_one(argc_plus_one_)
    {}

    ~Cleanup_argv()
    {
      for (unsigned int i = 0; i < argc_plus_one; ++i)
      {
        delete[] argv[i];
      }
      delete[] argv;
    }

  private:
    char **      argv;
    unsigned int argc_plus_one;
  };
  Cleanup_argv cleanup_argv(argv, argc + 1);

  // antscout->set_stream( out_stream );

  if (argc <= 4)
  {
    std::cout << "WarpLabeledPointSetFileMultiTransform ImageDimension inputVTKFile "
              << "outputVTKFile [-R reference_image] "
              << "{[deformation_field | [-i] affine_transform_txt ]}" << std::endl;
    return EXIT_FAILURE;
  }

  TRAN_OPT_QUEUE opt_queue;
  //    char *moving_image_filename = nullptr;
  char * input_vtk_filename = nullptr;
  char * output_vtk_filename = nullptr;
  char * reference_image_filename = nullptr;

  int kImageDim = std::stoi(argv[1]);

  const bool is_parsing_ok = WarpVTKPolyDataMultiTransform_ParseInput(
    argc - 2, argv + 2, input_vtk_filename, output_vtk_filename, reference_image_filename, opt_queue);

  if (is_parsing_ok)
  {
    switch (CheckFileType(output_vtk_filename))
    {
      case DEFORMATION_FILE:
      {
        if (reference_image_filename == nullptr)
        {
          std::cout << "the reference image file (-R) must be given!!!" << std::endl;
          return false;
        }

        std::cout << "output_vtk_filename: " << output_vtk_filename << std::endl;
        std::cout << "reference_image_filename: ";
        if (reference_image_filename)
        {
          std::cout << reference_image_filename << std::endl;
        }
        else
        {
          std::cout << "NULL" << std::endl;
        }
        DisplayOptQueue(opt_queue);

        switch (kImageDim)
        {
          case 2:
          {
            WarpLabeledPointSetFileMultiTransform<2>(
              input_vtk_filename, output_vtk_filename, reference_image_filename, opt_queue);
            break;
          }
          case 3:
          {
            WarpLabeledPointSetFileMultiTransform<3>(
              input_vtk_filename, output_vtk_filename, reference_image_filename, opt_queue);
            break;
          }
        }
        break;
      }

      case AFFINE_FILE:
      {
        std::cout << "output_affine_txt: " << output_vtk_filename << std::endl;
        std::cout << "reference_affine_txt: ";
        if (reference_image_filename)
        {
          std::cout << reference_image_filename << std::endl;
        }
        else
        {
          std::cout << "NULL" << std::endl;
        }
        DisplayOptQueue(opt_queue);

        switch (kImageDim)
        {
          case 2:
          {
            ComposeMultiAffine<2>(input_vtk_filename, output_vtk_filename, reference_image_filename, opt_queue);
            break;
          }
          case 3:
          {
            ComposeMultiAffine<3>(input_vtk_filename, output_vtk_filename, reference_image_filename, opt_queue);
            break;
          }
        }
        break;
      }

      default:
        std::cout << "Unknow output file format: " << output_vtk_filename << std::endl;
        break;
    }
  }
  else
  {
    std::cout << "Input error!" << std::endl;
  }

  return EXIT_FAILURE;
}
} // namespace ants
