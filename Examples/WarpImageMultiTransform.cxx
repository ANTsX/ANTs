#include "antsUtilities.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkTransformFactory.h"
#include "itkWarpImageMultiTransformFilter.h"
#include "itkTransformFileReader.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "antsUtilities.h"
// Needed for the LabelImageGaussianInterpolateImageFunction to work on
// vector images
#include "itkLabelImageGaussianInterpolateImageFunction.h"

namespace ants
{
static bool IsInverseDeformation(const char *str)
{
  std::string            filename = str;
  std::string::size_type pos = filename.rfind( "Inverse" );

  if( pos == std::string::npos )
    {
    return false;
    }
  else
    {
    return true;
    }
}

static bool WarpImageMultiTransform_ParseInput(int argc, char * *argv, char *& moving_image_filename,
                                               char *& output_image_filename,
                                               TRAN_OPT_QUEUE & opt_queue, MISC_OPT & misc_opt,
                                               int NDimensions)
{
  opt_queue.clear();
  opt_queue.reserve(argc - 2);

  misc_opt.reference_image_filename = NULL;
  misc_opt.use_BSpline_interpolator = false;
  misc_opt.use_TightestBoundingBox = false;
  misc_opt.use_RotationHeader = false;

  misc_opt.use_NN_interpolator = false;
  misc_opt.use_MultiLabel_interpolator = false;
  misc_opt.use_BSpline_interpolator = false;

  moving_image_filename = argv[0];
  output_image_filename = argv[1];

  int  ind = 2;
  bool set_current_affine_inv = false;

  while( ind < argc )
    {
    if( strcmp(argv[ind], "--use-NN") == 0 )
      {
      misc_opt.use_NN_interpolator = true;
      }
    else if( strcmp(argv[ind], "--use-BSpline") == 0 )
      {
      misc_opt.use_BSpline_interpolator = true;
      }
    else if( strcmp(argv[ind], "--use-ML") == 0 )
      {
      misc_opt.use_MultiLabel_interpolator = true;
      ind++; if( ind >= argc )
        {
        return false;
        }

      char *s = argv[ind];
      if( strlen(s) > 3 && strcmp(s + strlen(s) - 3, "vox") == 0 )
        {
        misc_opt.opt_ML.physical_units = false;
        s[strlen(s) - 3] = 0;
        }
      else if( strlen(s) > 2 && strcmp(s + strlen(s) - 2, "mm") == 0 )
        {
        misc_opt.opt_ML.physical_units = true;
        s[strlen(s) - 2] = 0;
        }
      else
        {
        antscout << "Wrong specification of sigma in --use-ML. Must end with 'mm' or 'vox'" << std::endl;
        return false;
        }

      misc_opt.opt_ML.sigma.resize(NDimensions);
      if( strchr(s, 'x') )
        {
        char *tok = strtok(s, "x");
        int i = 0;
        while( tok!=NULL && i<NDimensions )
          {
          double x = atof(tok);
          if( x < 0 )
            {
            antscout << "Negative sigma specification:" << s << std::endl;
            }
          misc_opt.opt_ML.sigma[i] = x;
          tok = strtok(NULL, "x");
          ++i;
          }
        if( i!=NDimensions || tok!=NULL )
          {
          antscout << "Invalid sigma specification:" << s << std::endl;
          }
        }
      else
        {
        double x = atof(s);
        if( x < 0 )
          {
          antscout << "Negative sigma specification:" << s << std::endl;
          }
        misc_opt.opt_ML.sigma.resize(NDimensions);
        std::fill(misc_opt.opt_ML.sigma.begin(), misc_opt.opt_ML.sigma.end(), x);
        }
      }

    else if( strcmp(argv[ind], "-R") == 0 )
      {
      ind++; if( ind >= argc )
        {
        return false;
        }
      misc_opt.reference_image_filename = argv[ind];
      }
    else if( (strcmp(argv[ind], "--tightest-bounding-box") == 0) &&  (strcmp(argv[ind], "-R") != 0)  )
      {
      misc_opt.use_TightestBoundingBox = true;
      }
    else if( strcmp(argv[ind], "--reslice-by-header") == 0 )
      {
      misc_opt.use_RotationHeader = true;
      TRAN_OPT opt;
      opt.file_type = IMAGE_AFFINE_HEADER;
      opt.do_affine_inv = false;
      opt_queue.push_back(opt);
      }
    else if( strcmp(argv[ind], "--Id") == 0 )
      {
      TRAN_OPT opt;
      opt.filename = "--Id";
      opt.do_affine_inv = false;
      opt.file_type = IDENTITY_TRANSFORM;
      opt_queue.push_back(opt);
      }
    else if( strcmp(argv[ind], "--moving-image-header") == 0 || strcmp(argv[ind], "-mh") == 0 )
      {
      TRAN_OPT opt;
      opt.file_type = IMAGE_AFFINE_HEADER;
      opt.filename = moving_image_filename;
      //            opt.do_affine_inv = false;
      SetAffineInvFlag(opt, set_current_affine_inv);
      opt_queue.push_back(opt);
      }
    else if( strcmp(argv[ind], "--reference-image-header") == 0 || strcmp(argv[ind], "-rh") == 0 )
      {
      if( misc_opt.reference_image_filename == NULL )
        {
        antscout
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
    else if( strcmp(argv[ind], "-i") == 0 )
      {
      set_current_affine_inv = true;
      }

    else if( strcmp(argv[ind], "--ANTS-prefix") == 0 )
      {
      ind++;
      std::string prefix = argv[ind];
      std::string path, name, ext;
      FilePartsWithgz(prefix, path, name, ext);
      if( ext == "" )
        {
        ext = ".nii.gz";
        }

      std::string deform_file_name, x_deform_name;
      deform_file_name = path + name + std::string("Warp") + ext;
      x_deform_name = path + name + std::string("Warpxvec") + ext;
      if( CheckFileExistence(x_deform_name.c_str() ) )
        {
        TRAN_OPT opt;
        opt.filename = deform_file_name.c_str();
        opt.file_type = CheckFileType(opt.filename.c_str() );
        opt.do_affine_inv = false;
        opt_queue.push_back(opt);
        antscout << "found deformation file: " << opt.filename << std::endl;
        DisplayOpt(opt);
        }

      std::string affine_file_name;
      affine_file_name = path + name + std::string("Affine") + GetPreferredTransformFileType();
      if( CheckFileExistence(affine_file_name.c_str() ) )
        {
        TRAN_OPT opt;
        opt.filename = affine_file_name.c_str();
        opt.file_type = CheckFileType(opt.filename.c_str() );
        opt.do_affine_inv = false;
        opt_queue.push_back(opt);
        antscout << "found affine file: " << opt.filename << std::endl;
        DisplayOpt(opt);
        }
      }
    else if( strcmp(argv[ind], "--ANTS-prefix-invert") == 0 )
      {
      ind++;
      std::string prefix = argv[ind];
      std::string path, name, ext;
      FilePartsWithgz(prefix, path, name, ext);
      if( ext == "" )
        {
        ext = ".nii.gz";
        }

      std::string affine_file_name;
      affine_file_name = path + name + std::string("Affine") + GetPreferredTransformFileType();
      if( CheckFileExistence(affine_file_name.c_str() ) )
        {
        TRAN_OPT opt;
        opt.filename = affine_file_name.c_str();
        opt.file_type = CheckFileType(opt.filename.c_str() );
        opt.do_affine_inv = true;
        opt_queue.push_back(opt);
        antscout << "found affine file: " << opt.filename << std::endl;
        DisplayOpt(opt);
        }

      std::string deform_file_name, x_deform_name;
      deform_file_name = path + name + std::string("InverseWarp.nii.gz");
      x_deform_name = path + name + std::string("InverseWarpxvec.nii.gz");
      if( CheckFileExistence(x_deform_name.c_str() ) )
        {
        TRAN_OPT opt;
        opt.filename = deform_file_name.c_str();
        opt.file_type = CheckFileType(opt.filename.c_str() );
        opt.do_affine_inv = false;
        opt_queue.push_back(opt);
        antscout << "found deformation file: " << opt.filename << std::endl;
        DisplayOpt(opt);
        }
      }
    else
      {
      TRAN_OPT opt;
      opt.filename = argv[ind];
      opt.file_type = CheckFileType(opt.filename.c_str() );
      opt.do_affine_inv = false;
      if( opt.file_type == AFFINE_FILE )
        {
        SetAffineInvFlag(opt, set_current_affine_inv);
        }
      else if( opt.file_type == DEFORMATION_FILE && set_current_affine_inv )
        {
        antscout << "Ignore inversion of non-affine file type! " << std::endl;
        antscout << "opt.do_affine_inv:" << opt.do_affine_inv << std::endl;
        }

      opt_queue.push_back(opt);
      DisplayOpt(opt);
      }
    ind++;
    }

  if( misc_opt.use_RotationHeader )
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
    //               antscout << "Use Rotation Header!" << std::endl;
    }

  return true;
}

template <class AffineTransformPointer>
void GetIdentityTransform(AffineTransformPointer & aff)
{
  typedef typename AffineTransformPointer::ObjectType AffineTransform;
  aff = AffineTransform::New();
  aff->SetIdentity();
}

template <int ImageDimension, unsigned int NVectorComponents>
void WarpImageMultiTransform(char *moving_image_filename, char *output_image_filename,
                             TRAN_OPT_QUEUE & opt_queue, MISC_OPT & misc_opt)
{
  typedef float RealType;
  typedef itk::Vector<RealType,
                      NVectorComponents>                                                             PixelType;
  typedef itk::Image<PixelType,
                     ImageDimension>                                                                ImageType;
  typedef itk::VectorImage<RealType,
                           ImageDimension>                                                           RefImageType;
  typedef itk::Vector<RealType,
                      ImageDimension>                                                                VectorType;
  typedef itk::Image<VectorType,
                     ImageDimension>
    DisplacementFieldType;
  typedef itk::MatrixOffsetTransformBase<double, ImageDimension,
                                         ImageDimension>                               AffineTransformType;
  typedef itk::WarpImageMultiTransformFilter<ImageType, ImageType, DisplacementFieldType,
                                             AffineTransformType> WarperType;

  itk::TransformFactory<AffineTransformType>::RegisterTransform();

  typedef itk::ImageFileReader<ImageType>    ImageFileReaderType;
  typedef itk::ImageFileReader<RefImageType> VectorImageFileReaderType;
  typename ImageFileReaderType::Pointer reader_img = ImageFileReaderType::New();
  reader_img->SetFileName(moving_image_filename);
  reader_img->Update();
  typename ImageType::Pointer img_mov = reader_img->GetOutput();

  typename RefImageType::Pointer img_ref;

  typename VectorImageFileReaderType::Pointer reader_img_ref = VectorImageFileReaderType::New();
  if( misc_opt.reference_image_filename )
    {
    reader_img_ref->SetFileName(misc_opt.reference_image_filename);
    reader_img_ref->Update();
    img_ref = reader_img_ref->GetOutput();
    }
  // else
  //    img_ref = NULL;

  typename WarperType::Pointer  warper = WarperType::New();
  warper->SetInput(img_mov);
  PixelType zero; zero.Fill(0);
  warper->SetEdgePaddingValue( zero );

  if( misc_opt.use_NN_interpolator )
    {
    typedef typename itk::NearestNeighborInterpolateImageFunction<ImageType,
                                                                  typename WarperType::CoordRepType> NNInterpolateType;
    typename NNInterpolateType::Pointer interpolator_NN = NNInterpolateType::New();
    antscout << "User nearest neighbor interpolation (was Haha) " << std::endl;
    warper->SetInterpolator(interpolator_NN);
    }
  else if( misc_opt.use_MultiLabel_interpolator )
    {
    antscout << " Need to fix in main itk repository " << std::endl;
//      typedef VectorPixelCompare<RealType, NVectorComponents> CompareType;
//      typedef typename itk::LabelImageGaussianInterpolateImageFunction<ImageType,
//                                                                       typename WarperType::CoordRepType,
//                                                                       CompareType> MLInterpolateType;
//      typename MLInterpolateType::Pointer interpolator_ML = MLInterpolateType::New();
//
//
//      antscout << "Using multi-label anti-aliasing interpolation " << std::endl;
//      vnl_vector_fixed<double, ImageDimension> sigma;
//      for(size_t i = 0; i < ImageDimension; i++)
//        {
//        if(misc_opt.opt_ML.physical_units)
//          sigma[i] = misc_opt.opt_ML.sigma[i] / img_mov->GetSpacing()[i];
//        else
//          sigma[i] = misc_opt.opt_ML.sigma[i];
//        }
//
//      antscout << "  Sigma = " << sigma << " (voxel units)" << std::endl;
//
//      interpolator_ML->SetParameters(sigma.data_block(), 4.0);
//
//      warper->SetInterpolator(interpolator_ML);
    }

  else if( misc_opt.use_BSpline_interpolator )
    {
    antscout << " Not currently supported because of a lack of vector support " << std::endl;
    /*
      typedef typename itk::BSplineInterpolateImageFunction<ImageType, typename WarperType::CoordRepType> BSInterpolateType;
      typename BSInterpolateType::Pointer interpolator_BS = BSInterpolateType::New();
      interpolator_BS->SetSplineOrder(3);
      antscout << "User B-spline interpolation " << std::endl;
      warper->SetInterpolator(interpolator_BS);
    */
    }
  else
    {
    typedef typename itk::LinearInterpolateImageFunction<ImageType,
                                                         typename WarperType::CoordRepType> LinInterpolateType;
    typename LinInterpolateType::Pointer interpolator_LN = LinInterpolateType::New();
    antscout << "User Linear interpolation " << std::endl;
    warper->SetInterpolator(interpolator_LN);
    }

  typedef itk::TransformFileReader                    TranReaderType;
  typedef itk::ImageFileReader<DisplacementFieldType> FieldReaderType;
  bool         takeaffinv = false;
  unsigned int transcount = 0;
  const int    kOptQueueSize = opt_queue.size();
  for( int i = 0; i < kOptQueueSize; i++ )
    {
    const TRAN_OPT & opt = opt_queue[i];

    switch( opt.file_type )
      {
      case AFFINE_FILE:
        {
        typename TranReaderType::Pointer tran_reader = TranReaderType::New();
        tran_reader->SetFileName(opt.filename);
        tran_reader->Update();
        typename AffineTransformType::Pointer aff = dynamic_cast<AffineTransformType *>
          ( (tran_reader->GetTransformList() )->front().GetPointer() );
        if( opt.do_affine_inv )
          {
          typename AffineTransformType::Pointer aff_inv = AffineTransformType::New();
          aff->GetInverse(aff_inv);
          aff = aff_inv;
          takeaffinv = true;
          }
        // antscout <<" aff " << transcount <<  std::endl;
        warper->PushBackAffineTransform(aff);
        if( transcount == 0 )
          {
          warper->SetOutputParametersFromImage( img_mov );
          }
        transcount++;
        break;
        }

      case IDENTITY_TRANSFORM:
        {
        typename AffineTransformType::Pointer aff;
        GetIdentityTransform(aff);
        // antscout << " aff id" << transcount << std::endl;
        warper->PushBackAffineTransform(aff);
        transcount++;
        break;
        }

      case IMAGE_AFFINE_HEADER:
        {
        typename AffineTransformType::Pointer aff = AffineTransformType::New();
        typename ImageFileReaderType::Pointer reader_image_affine = ImageFileReaderType::New();
        reader_image_affine->SetFileName(opt.filename);
        reader_image_affine->Update();
        typename ImageType::Pointer img_affine = reader_image_affine->GetOutput();

        GetAffineTransformFromImage<ImageType, AffineTransformType>(img_affine, aff);

        if( opt.do_affine_inv )
          {
          typename AffineTransformType::Pointer aff_inv = AffineTransformType::New();
          aff->GetInverse(aff_inv);
          aff = aff_inv;
          takeaffinv = true;
          }

        // antscout <<" aff from image header " << transcount <<  std::endl;
        warper->PushBackAffineTransform(aff);

        //            if (transcount==0){
        //                warper->SetOutputParametersFromImage( img_mov);
        //            }

        transcount++;
        break;
        }

      case DEFORMATION_FILE:
        {
        typename FieldReaderType::Pointer field_reader = FieldReaderType::New();
        field_reader->SetFileName( opt.filename );
        field_reader->Update();
        typename DisplacementFieldType::Pointer field = field_reader->GetOutput();

        warper->PushBackDisplacementFieldTransform(field);
        warper->SetOutputParametersFromImage(field );

        transcount++;
        break;
        }
      default:
        antscout << "Unknown file type!" << std::endl;
      }
    }

  // antscout << " transcount " << transcount << std::endl; warper->PrintTransformList();
  if( transcount == 2 )
    {
    antscout << "  We check the syntax of your call .... " << std::endl;
    const TRAN_OPT & opt1 = opt_queue[0];
    const TRAN_OPT & opt2 = opt_queue[1];

    if( opt1.file_type == AFFINE_FILE  && opt2.file_type == DEFORMATION_FILE   )
      {
      bool defisinv = IsInverseDeformation(opt2.filename.c_str() );
      if( !takeaffinv )
        {
        antscout
          <<
          " Your 1st parameter should be an inverse affine map and the 2nd an InverseWarp  --- exiting without applying warp.  Check that , if using an inverse affine map, you pass the -i option before the Affine.txt."
          << std::endl;
        return;
        }
      if( !defisinv )
        {
        antscout
          <<
          " Your 2nd  parameter should be an InverseWarp when your 1st parameter is an inverse affine map  --- exiting without applying warp.  "
          << std::endl;
        return;
        }
      }
    if( opt2.file_type == AFFINE_FILE  && opt1.file_type == DEFORMATION_FILE   )
      {
      bool defisinv = IsInverseDeformation(opt1.filename.c_str() );
      if(  defisinv )
        {
        antscout
          <<
          " Your 1st parameter should be a Warp (not Inverse) when your 2nd parameter is an affine map --- exiting without applying warp.  "
          << std::endl;
        return;
        }
      if(  takeaffinv )
        {
        antscout
          <<
          " Your 2nd parameter should be a regular affine map (not inverted) if the 1st is a Warp --- exiting without applying warp. "
          << std::endl;
        return;
        }
      }
    antscout << " syntax probably ok. " << std::endl;
    }
  else
    {
    antscout << " You are doing something more complex -- we wont check syntax in this case " << std::endl;
    }

  if( img_ref.IsNotNull() )
    {
    warper->SetOutputParametersFromImage( img_ref );
    }
  else
    {
    if( misc_opt.use_TightestBoundingBox == true )
      {
      // compute the desired spacking after inputting all the transform files using the

      typename ImageType::SizeType largest_size;
      typename ImageType::PointType origin_warped;
      GetLargestSizeAfterWarp<WarperType, ImageType>(warper, img_mov, largest_size, origin_warped);
      // Use img_mov as initial template space, then overwrite individual components as desired
      warper->SetOutputParametersFromImage( img_mov );

      warper->SetOutputSize(largest_size);
      warper->SetOutputOrigin(origin_warped);
        {
        typename ImageType::DirectionType d;
        d.SetIdentity();
        warper->SetOutputDirection(d);
        }
      }
    }

  antscout << "output origin: " << warper->GetOutputOrigin() << std::endl;
  antscout << "output size: " << warper->GetOutputSize() << std::endl;
  antscout << "output spacing: " << warper->GetOutputSpacing() << std::endl;
  antscout << "output direction: " << warper->GetOutputDirection() << std::endl;

  // warper->PrintTransformList();
  warper->DetermineFirstDeformNoInterp();
  warper->Update();

  //    {
  //        typename ImageType::IndexType ind_orig, ind_warped;
  //        ind_orig[0] = 128;
  //        ind_orig[1] = 128;
  //        ind_orig[2] = 16;
  //        typename ImageType::PointType pt_orig, pt_warped;
  //        warper->GetOutput()->TransformIndexToPhysicalPoint(ind_orig, pt_orig);
  //        warper->MultiTransformSinglePoint(pt_orig, pt_warped);
  //        img_mov->TransformPhysicalPointToIndex(pt_warped, ind_warped);
  //        antscout << "Transform output index " << ind_orig << "("<<pt_orig<<")"
  //        << " from moving image index " << ind_warped << "("<<pt_warped<<")" << std::endl;
  //    }

  //    typename ImageType::PointType pt_in, pt_out;
  //    for(unsigned int i=0; i<ImageDimension; i++){
  //        pt_in[i] = warper->GetOutputSize()[i] * 0.5;
  //    }
  //    warper->MultiTransformSinglePoint(pt_in, pt_out);
  //    antscout << "pt_in=" << pt_in << " pt_out=" <<pt_out << std::endl;

  typename ImageType::Pointer img_output = warper->GetOutput();

  typedef itk::ImageFileWriter<ImageType> ImageFileWriterType;
  typename ImageFileWriterType::Pointer writer_img = ImageFileWriterType::New();
  if( img_ref )
    {
    img_output->SetDirection(img_ref->GetDirection() );
    }
  writer_img->SetFileName(output_image_filename);
  writer_img->SetInput(img_output);
  writer_img->Update();
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int WarpImageMultiTransform( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "WarpImageMultiTransform" );

  int     argc = args.size();
  char* * argv = new char *[args.size() + 1];
  for( unsigned int i = 0; i < args.size(); ++i )
    {
    // allocate space for the string plus a null character
    argv[i] = new char[args[i].length() + 1];
    std::strncpy( argv[i], args[i].c_str(), args[i].length() );
    // place the null character in the end
    argv[i][args[i].length()] = '\0';
    }
  argv[argc] = 0;
  // class to automatically cleanup argv upon destruction
  class Cleanup_argv
  {
public:
    Cleanup_argv( char* * argv_, int argc_plus_one_ ) : argv( argv_ ), argc_plus_one( argc_plus_one_ )
    {
    }

    ~Cleanup_argv()
    {
      for( unsigned int i = 0; i < argc_plus_one; ++i )
        {
        delete[] argv[i];
        }
      delete[] argv;
    }

private:
    char* *      argv;
    unsigned int argc_plus_one;
  };
  Cleanup_argv cleanup_argv( argv, argc + 1 );

  antscout->set_stream( out_stream );

  if( argc <= 3 )
    {
    antscout <<  " \n " << std::endl;
    antscout <<  "Usage: \n " << std::endl;

    //    antscout << argv[0] <<  " ImageDimension moving_image output_image [-R reference_image |
    // --tightest-bounding-box] (--reslice-by-header) [--use-NN]"
    // << "[--ANTS-prefix prefix-name | --ANTS-prefix-invert prefix-name] {[deformation_field | [-i]
    // InverseAffineTransform.txt | --Id | [-i] --moving-image-header / -mh  | [-i] --reference-image-header / -rh]} \n"
    // << std::endl;
    antscout << argv[0]
             <<
      " ImageDimension moving_image output_image  -R reference_image --use-NN   SeriesOfTransformations--(See Below) "
             << std::endl;
    antscout << " SeriesOfTransformations --- " << argv[0]
             <<  " can apply, via concatenation, an unlimited number of transformations to your data ." << std::endl;
    antscout
      <<
      " Thus, SeriesOfTransformations may be  an Affine transform followed by a warp  another affine and then another warp. "
      << std::endl;
    antscout << "  Inverse affine transformations are invoked by calling   -i MyAffine.txt " << std::endl;
    antscout
      << " InverseWarps are invoked by passing the InverseWarp.nii.gz  filename (see below for a note about this).  "
      << std::endl;
    antscout << std::endl;
    antscout
      <<
      " Example 1: Mapping a warped image into the reference_image domain by applying abcdWarp.nii.gz and then abcdAffine.txt\n"
      << std::endl;

    antscout << argv[0] <<  " 3 moving_image output_image -R reference_image abcdWarp.nii.gz abcdAffine.txt\n"
             << std::endl;

    antscout
      <<
      " Example 2: To map the fixed/reference_image warped into the moving_image domain by applying the inversion of abcdAffine.txt and then abcdInverseWarp.nii.gz .\n"
      << std::endl;

    antscout << argv[0]
             << " 3 reference_image output_image -R moving_image -i  abcdAffine.txt abcdInverseWarp.nii.gz \n \n"
             << std::endl;
    antscout
      <<
      "  Note that the inverse maps (Ex. 2) are passed to this program in the reverse order of the forward maps (Ex. 1). "
      << std::endl;
    antscout << " This makes sense, geometrically ... see ANTS.pdf for visualization of this syntax." << std::endl;
    antscout << std::endl;
    antscout << " Compulsory arguments:\n " << std::endl;

    antscout << " ImageDimension: 2 or 3 (for 2 or 3 Dimensional registration)\n " << std::endl;

    antscout << " moving_image: the image to apply the transformation to\n " << std::endl;

    antscout << " output_image: the resulting image\n \n " << std::endl;

    antscout << " Optional arguments:\n " << std::endl;

    antscout << " -R: reference_image space that you wish to warp INTO." << std::endl;
    antscout
      <<
      "       --tightest-bounding-box: Computes the tightest bounding box using all the affine transformations. It will be overrided by -R reference_image if given."
      << std::endl;
    antscout
      <<
      "       --reslice-by-header: equivalient to -i -mh, or -fh -i -mh if used together with -R. It uses the orientation matrix and origin encoded in the image file header. "
      << std::endl;
    antscout
      << "       It can be used together with -R. This is typically not used together with any other transforms.\n "
      << std::endl;

    antscout << " --use-NN: Use Nearest Neighbor Interpolation. \n " << std::endl;
    antscout << " --use-BSpline: Use 3rd order B-Spline Interpolation. \n " << std::endl;
    antscout
      <<
      " --use-ML sigma: Use anti-aliasing interpolation for multi-label images, with Gaussian smoothing with standard deviation sigma. \n "
      << std::endl;
    antscout
      <<
      "                 Sigma can be specified in physical or voxel units, as in Convert3D. It can be a scalar or a vector. \n "
      << std::endl;
    antscout << "                 Examples:  --use-ML 0.4mm    -use-ML 0.8x0.8x0.8vox    " << std::endl;

    //    antscout << " --ANTS-prefix prefix-name: followed by a deformation field filename. \n " << std::endl;

    //    antscout << " --ANTS-prefix-invert: . \n" << std::endl;

    antscout << " -i: will use the inversion of the following affine transform. \n " << std::endl;

    //    antscout << " --Id: use an identity transform. \n " << std::endl;

    // antscout << " --moving-image-header or -mh: will use the orientation header of the moving image file. This is
    // typically not used with --reslice-by-header.\n " << std::endl;

    //    antscout << " --reference-image-header or -rh: use the orientation matrix and origin encoded in the image
    // file header. It can be used together with -R.\n " << std::endl;
    antscout <<  " \n " << std::endl;

    //        antscout << " For ANTS users:" << std::endl;

    antscout << " Other Example Usages:" << std::endl;
    antscout
      <<
      " Reslice the image: WarpImageMultiTransform 3 Imov.nii.gz Iout.nii.gz --tightest-bounding-box --reslice-by-header"
      << std::endl;
    antscout
      <<
      " Reslice the image to a reference image: WarpImageMultiTransform 3 Imov.nii.gz Iout.nii.gz -R Iref.nii.gz --tightest-bounding-box --reslice-by-header\n"
      << std::endl;

    antscout << " Important Notes: " << std::endl;
    antscout << " Prefixname \"abcd\" without any extension will use \".nii.gz\" by default" << std::endl;
    antscout
      <<
      " The abcdWarp and abcdInverseWarp do not exist. They are formed on the basis of abcd(Inverse)Warp.nii.gz when calling "
      << argv[0] << ", yet you have to use them as if they exist." << std::endl;
    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }

  TRAN_OPT_QUEUE opt_queue;
  char *         moving_image_filename = NULL;
  char *         output_image_filename = NULL;

  MISC_OPT misc_opt;

  const int  kImageDim = atoi(argv[1]);
  const bool is_parsing_ok = WarpImageMultiTransform_ParseInput(argc - 2, argv + 2,
                                                     moving_image_filename, output_image_filename,
                                                     opt_queue, misc_opt, kImageDim);

  if( is_parsing_ok )
    {
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(moving_image_filename,
                                                                           itk::ImageIOFactory::ReadMode);
    imageIO->SetFileName(moving_image_filename);
    imageIO->ReadImageInformation();
    unsigned int ncomponents = imageIO->GetNumberOfComponents();

    antscout << "moving_image_filename: " << moving_image_filename << " components " << ncomponents << std::endl;
    antscout << "output_image_filename: " << output_image_filename << std::endl;
    antscout << "reference_image_filename: ";
    if( misc_opt.reference_image_filename )
      {
      antscout << misc_opt.reference_image_filename << std::endl;
      }
    else
      {
      antscout << "NULL" << std::endl;
      }
    DisplayOptQueue(opt_queue);

    try
      {
      switch( kImageDim )
        {
        case 2:

          switch( ncomponents )
            {
            case 2:
              WarpImageMultiTransform<2, 2>(moving_image_filename, output_image_filename, opt_queue, misc_opt);
              break;
            default:
              WarpImageMultiTransform<2, 1>(moving_image_filename, output_image_filename, opt_queue, misc_opt);
              break;
            }
          break;
        case 3:

          switch( ncomponents )
            {
            case 3:
              WarpImageMultiTransform<3, 3>(moving_image_filename, output_image_filename, opt_queue, misc_opt);
              break;
            case 6:
              WarpImageMultiTransform<3, 6>(moving_image_filename, output_image_filename, opt_queue, misc_opt);
              break;
            default:
              WarpImageMultiTransform<3, 1>(moving_image_filename, output_image_filename, opt_queue, misc_opt);
              break;
            }
          break;
        case 4:

          switch( ncomponents )
            {
            case 4:
              WarpImageMultiTransform<4, 4>(moving_image_filename, output_image_filename, opt_queue, misc_opt);
              break;
            default:
              WarpImageMultiTransform<4, 1>(moving_image_filename, output_image_filename, opt_queue, misc_opt);
              break;
            }
          break;
        default:
          antscout << " not supported " << kImageDim  << std::endl;
          return EXIT_FAILURE;
        }
      }
    catch( itk::ExceptionObject & e )
      {
      antscout << "Exception caught during WarpImageMultiTransform." << std::endl;
      antscout << e << std::endl;
      return EXIT_FAILURE;
      }
    //      WarpImageMultiTransform<2,2>(moving_image_filename, output_image_filename, opt_queue, misc_opt);
    }
  else
    {
    antscout << "Input error!" << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
