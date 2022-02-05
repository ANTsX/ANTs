
#include "antsUtilities.h"
#include "antsAllocImage.h"
#include "itkImageFileReader.h"
#include "itkVariableLengthVector.h"

#include "itkImageFileWriter.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkTransformFactory.h"
#include "itkWarpTensorImageMultiTransformFilter.h"
#include "itkTransformFileReader.h"
#include "itkVectorNearestNeighborInterpolateImageFunction.h"
#include "ReadWriteData.h"
#include "itkWarpImageMultiTransformFilter.h"
#include "itkExpTensorImageFilter.h"

namespace ants
{
static bool
WarpTensorImageMultiTransform_ParseInput(int              argc,
                                         char **          argv,
                                         char *&          moving_image_filename,
                                         char *&          output_image_filename,
                                         TRAN_OPT_QUEUE & opt_queue,
                                         MISC_OPT &       misc_opt)
{
  opt_queue.clear();
  opt_queue.reserve(argc - 2);

  misc_opt.reference_image_filename = nullptr;
  misc_opt.use_NN_interpolator = false;
  misc_opt.use_TightestBoundingBox = false;
  misc_opt.use_RotationHeader = false;

  moving_image_filename = argv[0];
  output_image_filename = argv[1];

  int  ind = 2;
  bool set_current_affine_inv = false;

  while (ind < argc)
  {
    if (strcmp(argv[ind], "--use-NN") == 0)
    {
      misc_opt.use_NN_interpolator = true;
    }
    else if (strcmp(argv[ind], "-R") == 0)
    {
      ind++;
      if (ind >= argc)
      {
        return false;
      }
      misc_opt.reference_image_filename = argv[ind];
    }
    else if ((strcmp(argv[ind], "--tightest-bounding-box") == 0) && (strcmp(argv[ind], "-R") != 0))
    {
      misc_opt.use_TightestBoundingBox = true;
    }
    else if (strcmp(argv[ind], "--reslice-by-header") == 0)
    {
      misc_opt.use_RotationHeader = true;
      TRAN_OPT opt;
      opt.file_type = IMAGE_AFFINE_HEADER;
      opt.do_affine_inv = false;
      opt_queue.push_back(opt);
    }
    else if (strcmp(argv[ind], "--Id") == 0)
    {
      TRAN_OPT opt;
      opt.filename = "--Id";
      opt.do_affine_inv = false;
      opt.file_type = IDENTITY_TRANSFORM;
      opt_queue.push_back(opt);
    }
    else if (strcmp(argv[ind], "--moving-image-header") == 0 || strcmp(argv[ind], "-mh") == 0)
    {
      TRAN_OPT opt;
      opt.file_type = IMAGE_AFFINE_HEADER;
      opt.filename = moving_image_filename;
      //            opt.do_affine_inv = false;
      SetAffineInvFlag(opt, set_current_affine_inv);
      opt_queue.push_back(opt);
    }
    else if (strcmp(argv[ind], "--reference-image-header") == 0 || strcmp(argv[ind], "-rh") == 0)
    {
      if (misc_opt.reference_image_filename == nullptr)
      {
        std::cout
          << "reference image filename is not given yet. Specify it with -R before --reference-image-header / -rh."
          << std::endl;
        return false;
      }

      TRAN_OPT opt;
      opt.file_type = IMAGE_AFFINE_HEADER;
      opt.filename = misc_opt.reference_image_filename;
      //            opt.do_affine_inv = false;
      SetAffineInvFlag(opt, set_current_affine_inv);
      opt_queue.push_back(opt);
    }
    else if (strcmp(argv[ind], "-i") == 0)
    {
      set_current_affine_inv = true;
    }

    else if (strcmp(argv[ind], "--ANTS-prefix") == 0)
    {
      ind++;
      std::string prefix = argv[ind];
      std::string path, name, ext;
      FilePartsWithgz(prefix, path, name, ext);
      if (ext.empty())
      {
        ext = ".nii.gz";
      }

      std::string deform_file_name, x_deform_name;
      deform_file_name = path + name + std::string("Warp") + ext;
      x_deform_name = path + name + std::string("Warpxvec") + ext;
      if (CheckFileExistence(x_deform_name.c_str()))
      {
        TRAN_OPT opt;
        opt.filename = deform_file_name.c_str();
        opt.file_type = CheckFileType(opt.filename.c_str());
        opt.do_affine_inv = false;
        opt_queue.push_back(opt);
        std::cout << "found deformation file: " << opt.filename << std::endl;
        DisplayOpt(opt);
      }

      std::string affine_file_name;
      affine_file_name = path + name + std::string("Affine.txt");
      if (CheckFileExistence(affine_file_name.c_str()))
      {
        TRAN_OPT opt;
        opt.filename = affine_file_name.c_str();
        opt.file_type = CheckFileType(opt.filename.c_str());
        opt.do_affine_inv = false;
        opt_queue.push_back(opt);
        std::cout << "found affine file: " << opt.filename << std::endl;
        DisplayOpt(opt);
      }
    }
    else if (strcmp(argv[ind], "--ANTS-prefix-invert") == 0)
    {
      ind++;
      std::string prefix = argv[ind];
      std::string path, name, ext;
      FilePartsWithgz(prefix, path, name, ext);
      if (ext.empty())
      {
        ext = ".nii.gz";
      }

      std::string affine_file_name;
      affine_file_name = path + name + std::string("Affine.txt");
      if (CheckFileExistence(affine_file_name.c_str()))
      {
        TRAN_OPT opt;
        opt.filename = affine_file_name.c_str();
        opt.file_type = CheckFileType(opt.filename.c_str());
        opt.do_affine_inv = true;
        opt_queue.push_back(opt);
        std::cout << "found affine file: " << opt.filename << std::endl;
        DisplayOpt(opt);
      }

      std::string deform_file_name, x_deform_name;
      deform_file_name = path + name + std::string("InverseWarp.nii.gz");
      x_deform_name = path + name + std::string("InverseWarpxvec.nii.gz");
      if (CheckFileExistence(x_deform_name.c_str()))
      {
        TRAN_OPT opt;
        opt.filename = deform_file_name.c_str();
        opt.file_type = CheckFileType(opt.filename.c_str());
        opt.do_affine_inv = false;
        opt_queue.push_back(opt);
        std::cout << "found deformation file: " << opt.filename << std::endl;
        DisplayOpt(opt);
      }
    }
    else
    {
      TRAN_OPT opt;
      opt.filename = argv[ind];
      opt.file_type = CheckFileType(opt.filename.c_str());
      opt.do_affine_inv = false;
      if (opt.file_type == AFFINE_FILE)
      {
        SetAffineInvFlag(opt, set_current_affine_inv);
      }
      else if (opt.file_type == DEFORMATION_FILE && set_current_affine_inv)
      {
        std::cout << "Ignore inversion of non-affine file type! " << std::endl;
        std::cout << "opt.do_affine_inv:" << opt.do_affine_inv << std::endl;
      }

      opt_queue.push_back(opt);
      DisplayOpt(opt);
    }
    ind++;
  }

  if (misc_opt.use_RotationHeader)
  {
    //                if (misc_opt.reference_image_filename) {
    //                    opt_queue[0].filename = misc_opt.reference_image_filename;
    //                } else {
    opt_queue[0].filename = "--Id";
    opt_queue[0].file_type = IDENTITY_TRANSFORM;
    opt_queue[0].do_affine_inv = false;
    //                }

    //               TRAN_OPT opt;
    //               opt.file_type = IMAGE_AFFINE_HEADER;
    //               opt.filename = moving_image_filename;
    //               opt.do_affine_inv = true;
    //               opt_queue.push_back(opt);
    //
    //               std::cout << "Use Rotation Header!" << std::endl;
  }

  return true;
}

template <typename TAffineTransform>
static void
GetIdentityTransform(typename TAffineTransform::Pointer & aff)
{
  aff = TAffineTransform::New();
  aff->SetIdentity();
}

#if 0
template <typename TensorImageType, typename ImageType>
static void
DirectionCorrect( typename TensorImageType::Pointer img_mov, typename ImageType::Pointer img_ref )
{
  itk::ImageRegionIteratorWithIndex<TensorImageType> it(img_mov, img_mov->GetLargestPossibleRegion() );

  typename TensorImageType::DirectionType::InternalMatrixType direction = img_mov->GetDirection().GetTranspose()
    * img_ref->GetDirection().GetVnlMatrix();

  if( !direction.is_identity( 0.00001 ) )
    {
    while( !it.IsAtEnd() )
      {
      typename TensorImageType::DirectionType::InternalMatrixType  dt;
      dt(0, 0) = it.Value()[0];
      dt(0, 1) = dt(1, 0) = it.Value()[1];
      dt(0, 2) = dt(2, 0) = it.Value()[2];
      dt(1, 1) = it.Value()[3];
      dt(1, 2) = dt(2, 1) = it.Value()[4];
      dt(2, 2) = it.Value()[5];

      dt = direction * dt * direction.transpose();

      typename TensorImageType::PixelType outDt;

      outDt[0] = dt(0, 0);
      outDt[1] = dt(0, 1);
      outDt[2] = dt(0, 2);
      outDt[3] = dt(1, 1);
      outDt[4] = dt(1, 2);
      outDt[5] = dt(2, 2);

      it.Set( outDt );

      ++it;
      }
    }
}

#endif

template <int ImageDimension>
static void
WarpImageMultiTransform(char *           moving_image_filename,
                        char *           output_image_filename,
                        TRAN_OPT_QUEUE & opt_queue,
                        MISC_OPT &       misc_opt)
{
  // typedef itk::Vector<float,6> PixelType;
  using PixelType = itk::SymmetricSecondRankTensor<float, 3>;
  using TensorImageType = itk::Image<PixelType, ImageDimension>;
  using ImageType = itk::Image<float, ImageDimension>;
  using VectorType = itk::Vector<float, ImageDimension>;
  using DisplacementFieldType = itk::Image<VectorType, ImageDimension>;
  using AffineTransformType = itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension>;
  using WarperType =
    itk::WarpImageMultiTransformFilter<ImageType, ImageType, DisplacementFieldType, AffineTransformType>;

  itk::TransformFactory<AffineTransformType>::RegisterTransform();

  using ImageFileReaderType = itk::ImageFileReader<ImageType>;
  // typename ImageFileReaderType::Pointer reader_img = ImageFileReaderType::New();
  // reader_img->SetFileName(moving_image_filename);
  // reader_img->Update();
  // typename ImageType::Pointer img_mov = ImageType::New();
  // img_mov = reader_img->GetOutput();
  typename TensorImageType::Pointer img_mov;
  ReadTensorImage<TensorImageType>(img_mov, moving_image_filename, true);

  typename ImageType::Pointer img_ref;

  typename ImageFileReaderType::Pointer reader_img_ref = ImageFileReaderType::New();
  if (misc_opt.reference_image_filename)
  {
    reader_img_ref->SetFileName(misc_opt.reference_image_filename);
    reader_img_ref->Update();
    img_ref = reader_img_ref->GetOutput();
  }

  // Convert to reference image tensor basis
  // DirectionCorrect<TensorImageType, ImageType>(img_mov, img_ref);

  typename TensorImageType::Pointer img_output = AllocImage<TensorImageType>(img_ref);
  for (unsigned int tensdim = 0; tensdim < 6; tensdim++)
  {
    using IndexSelectCasterType = itk::VectorIndexSelectionCastImageFilter<TensorImageType, ImageType>;
    typename IndexSelectCasterType::Pointer fieldCaster = IndexSelectCasterType::New();
    fieldCaster->SetInput(img_mov);
    fieldCaster->SetIndex(tensdim);
    fieldCaster->Update();
    typename ImageType::Pointer tenscomponent = fieldCaster->GetOutput();
    tenscomponent->SetSpacing(img_mov->GetSpacing());
    tenscomponent->SetOrigin(img_mov->GetOrigin());
    tenscomponent->SetDirection(img_mov->GetDirection());

    typename WarperType::Pointer warper = WarperType::New();
    warper->SetInput(tenscomponent);
    //      PixelType nullPix;
    // nullPix.Fill(0);
    warper->SetEdgePaddingValue(0);

    if (misc_opt.use_NN_interpolator)
    {
      using NNInterpolateType =
        typename itk::NearestNeighborInterpolateImageFunction<ImageType, typename WarperType::CoordRepType>;
      typename NNInterpolateType::Pointer interpolator_NN = NNInterpolateType::New();
      std::cout << "Haha" << std::endl;
      warper->SetInterpolator(interpolator_NN);
    }

    using TranReaderType = itk::TransformFileReader;
    using FieldReaderType = itk::ImageFileReader<DisplacementFieldType>;

    unsigned int transcount = 0;
    const int    kOptQueueSize = opt_queue.size();
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
          if (opt.do_affine_inv)
          {
            typename AffineTransformType::Pointer aff_inv = AffineTransformType::New();
            aff->GetInverse(aff_inv);
            aff = aff_inv;
          }
          // std::cout <<" aff " << transcount <<  std::endl;
          warper->PushBackAffineTransform(aff);
          if (transcount == 0)
          {
            warper->SetOutputParametersFromImage(img_mov);
          }
          transcount++;
        }
        break;
        case IDENTITY_TRANSFORM:
        {
          typename AffineTransformType::Pointer aff;
          GetIdentityTransform<AffineTransformType>(aff);
          // std::cout << " aff id" << transcount << std::endl;
          warper->PushBackAffineTransform(aff);
          transcount++;
        }
        break;
        case IMAGE_AFFINE_HEADER:
        {
          typename AffineTransformType::Pointer aff = AffineTransformType::New();
          typename ImageFileReaderType::Pointer reader_image_affine = ImageFileReaderType::New();
          reader_image_affine->SetFileName(opt.filename);
          reader_image_affine->Update();
          typename ImageType::Pointer img_affine = reader_image_affine->GetOutput();

          GetAffineTransformFromImage<ImageType, AffineTransformType>(img_affine, aff);

          if (opt.do_affine_inv)
          {
            typename AffineTransformType::Pointer aff_inv = AffineTransformType::New();
            aff->GetInverse(aff_inv);
            aff = aff_inv;
          }

          // std::cout <<" aff from image header " << transcount <<  std::endl;
          warper->PushBackAffineTransform(aff);

          //            if (transcount==0){
          //                warper->SetOutputParametersFromImage( img_mov );
          //            }

          transcount++;
        }
        break;
        case DEFORMATION_FILE:
        {
          typename FieldReaderType::Pointer field_reader = FieldReaderType::New();
          field_reader->SetFileName(opt.filename);
          field_reader->Update();
          typename DisplacementFieldType::Pointer field = field_reader->GetOutput();
          warper->PushBackDisplacementFieldTransform(field);
          warper->SetOutputParametersFromImage(field);
          transcount++;
        }
        break;
        default:
        {
          std::cout << "Unknown file type!" << std::endl;
        }
      }
    }

    // warper->PrintTransformList();

    if (img_ref.IsNotNull())
    {
      warper->SetOutputParametersFromImage(img_ref);
    }
    else
    {
      if (misc_opt.use_TightestBoundingBox == true)
      {
        // compute the desired spacking after inputting all the transform files using the

        typename ImageType::SizeType  largest_size;
        typename ImageType::PointType origin_warped;
        GetLargestSizeAfterWarp<WarperType, TensorImageType>(warper, img_mov, largest_size, origin_warped);
        warper->SetOutputParametersFromImage(img_mov);
        warper->SetOutputSize(largest_size);
        warper->SetOutputOrigin(origin_warped);
        {
          typename ImageType::DirectionType d;
          d.SetIdentity();
          warper->SetOutputDirection(d);
        }
      }
    }

    // std::cout << "output origin: " << warper->GetOutputOrigin() << std::endl;
    // std::cout << "output size: " << warper->GetOutputSize() << std::endl;
    // std::cout << "output spacing: " << warper->GetOutputSpacing() << std::endl;
    //    std::cout << "output direction: " << warper->GetOutputDirection() << std::endl;

    // warper->PrintTransformList();
    warper->DetermineFirstDeformNoInterp();
    warper->Update();

    using Iterator = itk::ImageRegionIteratorWithIndex<TensorImageType>;
    Iterator vfIter2(img_output, img_output->GetLargestPossibleRegion());
    for (vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2)
    {
      PixelType tens = vfIter2.Get();
      tens[tensdim] = warper->GetOutput()->GetPixel(vfIter2.GetIndex());
      vfIter2.Set(tens);
    }
  }

  // DirectionCorrect<TensorImageType>(img_output, img_mov);

  WriteTensorImage<TensorImageType>(img_output, output_image_filename, true);
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
WarpTensorImageMultiTransform(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "WarpTensorImageMultiTransform");

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
    std::cout
      << "WarpImageMultiTransform ImageDimension moving_image output_image [-R reference_image | "
         "--tightest-bounding-box] (--reslice-by-header) [--use-NN (use Nearest Neighbor Interpolator)]"
      << "[--ANTS-prefix prefix-name | --ANTS-prefix-invert prefix-name] {[deformation_field | [-i] "
         "affine_transform_txt | --Id | [-i] --moving-image-header / -mh  | [-i] --reference-image-header / -rh]}"
      << std::endl
      << "Example:" << std::endl
      << "Reslice the image: WarpImageMultiTransform 3 Imov.nii Iout.nii --tightest-bounding-box --reslice-by-header"
      << std::endl
      << "Reslice the image to a reference image: WarpImageMultiTransform 3 Imov.nii Iout.nii -R Iref.nii "
         "--tightest-bounding-box --reslice-by-header"
      << std::endl
      << "Note:" << std::endl
      << "-i will use the inversion of the following affine transform." << std::endl
      << "--tightest-bounding-box will be overrided by -R reference_image if given. It computes the tightest bounding "
         "box using all the affine transformations."
      << std::endl
      << "--Id uses the identity transform." << std::endl
      << "--moving-image-header or -mh in short will use the orientation header of the moving image file. This is "
         "typically not used with --reslice-by-header."
      << std::endl
      << "--reference-image-header or -rh in short will use the orientation header of the fixed image file. This is "
         "typically not used with --reslice-by-header."
      << std::endl
      << "--reslice-by-header uses the orientation matrix and origin encoded in the image file header. It can be used "
         "together with -R. "
      << "This is typically not used together with any other transforms. "
      << "--reslice-by-header is equvalient to -i -mh, or -fh -i -mh if used together with -R. " << std::endl;
    std::cout << std::endl
              << "For ANTS users:" << std::endl
              << "To use with the deformation field and the affine transform files generated from ANTS:" << std::endl
              << "--ANTS-prefix prefix-name" << std::endl
              << "--ANTS-prefix-invert prefix-name" << std::endl
              << "Example:" << std::endl
              << "3 moving_image output_image -R reference_image --ANTS-prefix abcd.nii.gz" << std::endl
              << "Applies abcdWarpxvec.nii.gz/abcdWarpyvec.nii.gz/abcdWarpzvec.nii.gz and then abcdAffine.txt. Use "
                 "this with ANTS to get the moving_image warped into the reference_image domain. "
              << std::endl
              << "3 reference_image output_image -R moving_image --ANTS-prefix-invert abcd.nii.gz --ANTS-invert"
              << std::endl
              << "Applies the inversion of abcdAffine.txt and then "
                 "abcdInverseWarpxvec.nii.gz/abcdInverseWarpyvec.nii.gz/abcdInverseWarpzvec.nii.gz. Use this with ANTS "
                 "to get the reference_image warped into the moving_image domain. "
              << std::endl
              << "Note: " << std::endl
              << R"(prefix name "abcd" without any extension will use ".nii.gz" by default)" << std::endl;
    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  TRAN_OPT_QUEUE opt_queue;
  char *         moving_image_filename = nullptr;
  char *         output_image_filename = nullptr;

  MISC_OPT misc_opt;

  const int  kImageDim = std::stoi(argv[1]);
  const bool is_parsing_ok = WarpTensorImageMultiTransform_ParseInput(
    argc - 2, argv + 2, moving_image_filename, output_image_filename, opt_queue, misc_opt);

  if (is_parsing_ok)
  {
    std::cout << "moving_image_filename: " << moving_image_filename << std::endl;
    std::cout << "output_image_filename: " << output_image_filename << std::endl;
    std::cout << "reference_image_filename: ";
    if (misc_opt.reference_image_filename)
    {
      std::cout << misc_opt.reference_image_filename << std::endl;
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
        WarpImageMultiTransform<2>(moving_image_filename, output_image_filename, opt_queue, misc_opt);
      }
      break;
      case 3:
      {
        WarpImageMultiTransform<3>(moving_image_filename, output_image_filename, opt_queue, misc_opt);
      }
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
