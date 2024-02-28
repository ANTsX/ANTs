// compute the average of a list of affine transform

#include "antsUtilities.h"
#include "itkImageFileReader.h"

#include "itkImageFileWriter.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkTransformFactory.h"

#include "itkAverageAffineTransformFunction.h"

#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"

namespace ants
{
static bool
AverageAffineTransformNoRigid_ParseInput(int              argc,
                                         char **          argv,
                                         char *&          output_transform_filename,
                                         char *&          reference_transform_filename,
                                         TRAN_OPT_QUEUE & opt_queue)
{
  opt_queue.clear();
  opt_queue.reserve(argc);

  output_transform_filename = argv[0];

  reference_transform_filename = nullptr;

  int ind = 1;
  while (ind < argc)
  {
    if (strcmp(argv[ind], "-R") == 0)
    {
      ind++;
      if (ind >= argc)
      {
        return false;
      }
      reference_transform_filename = argv[ind];
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
        std::cerr << "file: " << opt.filename << " is not an affine .txt file. Invalid to use '-i' " << std::endl;
        return false;
      }
      opt.file_type = AFFINE_FILE;
      opt.do_affine_inv = true;

      opt.weight = 1.0;   // default value
      if (ind < argc - 1) // test if still has extra parameters
      {
        double weight;
        if (get_a_double_number(argv[ind + 1], weight))
        {
          ind++;
          opt.weight = weight;
        }
      }
      opt_queue.push_back(opt);
    }
    else
    {
      TRAN_OPT opt;
      opt.filename = argv[ind];
      if (CheckFileType(opt.filename) != AFFINE_FILE)
      {
        std::cerr << "file: " << opt.filename << " is not an affine .txt file." << std::endl;
        return false;
      }
      opt.file_type = CheckFileType(opt.filename);
      opt.do_affine_inv = false;

      opt.weight = 1.0;   // default value
      if (ind < argc - 1) // test if still has extra parameters
      {
        double weight;
        if (get_a_double_number(argv[ind + 1], weight))
        {
          ind++;
          opt.weight = weight;
        }
      }

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
AverageAffineTransformNoRigid(char * output_affine_txt, char * reference_affine_txt, TRAN_OPT_QUEUE & opt_queue)
{
  //    typedef itk::Image<float, ImageDimension> ImageType;
  //    typedef itk::Vector<float, ImageDimension> VectorType;
  //    typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;

  using AffineTransformType = itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension>;
  //    typedef itk::WarpImageMultiTransformFilter<ImageType, ImageType,
  //            DisplacementFieldType, AffineTransformType> WarperType;

  using WarperType = itk::AverageAffineTransformFunction<AffineTransformType>;

  itk::TransformFactory<AffineTransformType>::RegisterTransform();

  // typedef itk::ImageFileReader<ImageType> ImageFileReaderType;
  // typename ImageFileReaderType::Pointer reader_img = ImageFileReaderType::New();
  // typename ImageType::Pointer img_ref = ImageType::New();

  // typename ImageFileReaderType::Pointer reader_img_ref = ImageFileReaderType::New();

  WarperType average_func;
  average_func.verbose = true;
  average_func.useRigid = false;
  // warper->SetInput(img_mov);
  // warper->SetEdgePaddingValue( 0);
  //    VectorType pad;
  //    pad.Fill(0);
  // warper->SetEdgePaddingValue(pad);

  using TranReaderType = itk::TransformFileReader;

  //    typedef itk::ImageFileReader<DisplacementFieldType> FieldReaderType;

  int       cnt_affine = 0;
  const int kOptQueueSize = opt_queue.size();
  for (int i = 0; i < kOptQueueSize; i++)
  {
    const TRAN_OPT & opt = opt_queue[i];

    switch (opt.file_type)
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

        double weight = opt.weight;
        average_func.PushBackAffineTransform(aff, weight);
        cnt_affine++;
        break;
      }
      case DEFORMATION_FILE:
      {
        std::cout << "Average affine only files: ignore " << opt.filename << std::endl;
      }
      break;
      default:
        std::cout << "Unknown file type!" << std::endl;
    }
  }

  using PointType = typename WarperType::PointType;
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
      aff_ref_tmp = average_func.GetTransformList().begin()->aff;
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

  average_func.AverageMultipleAffineTransform(aff_center, aff_output);

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
AverageAffineTransformNoRigid(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "AverageAffineTransformNoRigid");
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

  if (argc <= 3)
  {
    std::cerr
      << "AverageAffineTransformNoRigid ImageDimension output_affine_transform [-R reference_affine_transform] "
      << "{[-i] affine_transform_txt [weight(=1)] ]}" << std::endl
      << std::endl
      << " Usage: Compute weighted average of input affine transforms. " << std::endl
      << "For 2D and 3D transform, the affine transform is first decomposed into "
         "scale x shearing x rotation. Then these parameters are averaged, using the weights if they provided. "
         "For 3D transform, the rotation component is the quaternion. After averaging, the quaternion will also "
         "be normalized to have unit norm. For 2D transform, the rotation component is the rotation angle. "
         "The weight for each transform is a non-negative number. The sum of all weights will be normalized to 1 "
         "before averaging. The default value for each weight is 1.0. "
      << std::endl
      << std::endl
      << "All affine transforms is a \"centerd\" transform, following ITK convention. A reference_affine_transform"
         " defines the center for the output transform. The first provided transform is the default reference "
         "transform"
      << std::endl
      << "Output affine transform is a MatrixOffsetBaseTransform." << std::endl
      << " -i option takes the inverse of the affine mapping." << std::endl
      << " For example: " << std::endl
      << " 2 output_affine.txt -R A.txt A1.txt 1.0 -i A2.txt 2.0 A3.txt A4.txt 6.0 A5.txt" << std::endl
      << "This computes: (1*A1 + 2*(A2)^-1 + A3 + A4*6 + A5 ) / (1+2+1+6+5)" << std::endl;
    return EXIT_SUCCESS;
  }

  TRAN_OPT_QUEUE opt_queue;

  char * output_transform_filename = nullptr;
  char * reference_transform_filename = nullptr;

  int kImageDim = std::stoi(argv[1]);

  const bool is_parsing_ok = AverageAffineTransformNoRigid_ParseInput(
    argc - 2, argv + 2, output_transform_filename, reference_transform_filename, opt_queue);

  if (is_parsing_ok)
  {
    std::cout << "output_transform_filename: " << output_transform_filename << std::endl;
    std::cout << "reference_transform_filename: ";

    if (reference_transform_filename)
    {
      std::cout << reference_transform_filename << std::endl;
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
        AverageAffineTransformNoRigid<2>(output_transform_filename, reference_transform_filename, opt_queue);
      }
      break;
      case 3:
      {
        AverageAffineTransformNoRigid<3>(output_transform_filename, reference_transform_filename, opt_queue);
      }
      break;
    }
  }

  else
  {
    std::cerr << "Input error!" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
} // namespace ants
