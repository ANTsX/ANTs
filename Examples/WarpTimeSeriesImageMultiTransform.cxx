
#include "antsUtilities.h"
#include "antsAllocImage.h"
#include "itkImageFileReader.h"
#include "itkVariableLengthVector.h"
#include "itkImageFileWriter.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkTransformFactory.h"
#include "itkTransformFileReader.h"
#include "itkVectorNearestNeighborInterpolateImageFunction.h"
#include "ReadWriteImage.h"
#include "itkWarpImageMultiTransformFilter.h"
#include "itkExtractImageFilter.h"

namespace ants
{
static bool WarpTimeSeriesImageMultiTransform_ParseInput(int argc, char * *argv, char *& moving_image_filename,
                                                         char *& output_image_filename,
                                                         TRAN_OPT_QUEUE & opt_queue, MISC_OPT & misc_opt)
{
  opt_queue.clear();
  opt_queue.reserve(argc - 2);

  misc_opt.reference_image_filename = NULL;
  misc_opt.use_NN_interpolator = false;
  misc_opt.use_TightestBoundingBox = false;
  misc_opt.use_RotationHeader = false;

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
      affine_file_name = path + name + std::string("Affine.txt");
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
      affine_file_name = path + name + std::string("Affine.txt");
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

template <int ImageDimension>
void WarpImageMultiTransformFourD(char *moving_image_filename, char *output_image_filename,
                                  TRAN_OPT_QUEUE & opt_queue, MISC_OPT & misc_opt)
{
  typedef itk::Image<float,
                     ImageDimension>                              VectorImageType; // 4D contains functional image
  typedef itk::Image<float, ImageDimension
                     - 1>                                         ImageType; // 3D image domain -R option
  typedef itk::Vector<float, ImageDimension
                      - 1>                                        VectorType; // 3D warp
  typedef itk::Image<VectorType, ImageDimension
                     - 1>                                         DisplacementFieldType; // 3D Field
  typedef itk::MatrixOffsetTransformBase<double, ImageDimension - 1, ImageDimension
                                         - 1>                     AffineTransformType;
  typedef itk::WarpImageMultiTransformFilter<ImageType, ImageType, DisplacementFieldType,
                                             AffineTransformType> WarperType;

  itk::TransformFactory<AffineTransformType>::RegisterTransform();

  typedef itk::ImageFileReader<ImageType> ImageFileReaderType;

  typename VectorImageType::Pointer img_mov;

  ReadImage<VectorImageType>(img_mov, moving_image_filename);
  antscout << " Four-D image size: " << img_mov->GetLargestPossibleRegion().GetSize() << std::endl;
  typename ImageType::Pointer img_ref;
  if( misc_opt.reference_image_filename )
    {
    ReadImage<ImageType>( img_ref, misc_opt.reference_image_filename );
    }

  typedef itk::ExtractImageFilter<VectorImageType, ImageType> ExtractFilterType;
  typedef itk::ImageRegionIteratorWithIndex<VectorImageType>  ImageIt;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>        SliceIt;

  // ORIENTATION ALERT -- the way this code sets up
  // transformedvecimage doesn't really make complete sense to me. In
  // particular, the 'grab upper dim-1 x dim-1 of directions' method
  // of setting lower-dimension dir cosines can lead to singular mattrices.
  // allocate output image
  typename VectorImageType::RegionType region = img_mov->GetLargestPossibleRegion();
  for( unsigned int i = 0; i < ImageDimension - 1; i++ )
    {
    region.SetSize( i, img_ref->GetLargestPossibleRegion().GetSize()[i] );
    }
  typename VectorImageType::Pointer transformedvecimage = AllocImage<VectorImageType>(region);

  typename VectorImageType::DirectionType direction = transformedvecimage->GetDirection();
  direction.Fill(0);
  typename VectorImageType::PointType origin;
  typename VectorImageType::SpacingType spc;
  for( unsigned int i = 0; i < ImageDimension - 1; i++ )
    {
    for( unsigned int j = 0; j < ImageDimension - 1; j++ )
      {
      direction[i][j] = img_ref->GetDirection()[i][j];
      }
    spc[i] = img_ref->GetSpacing()[i];
    origin[i] = img_ref->GetOrigin()[i];
    }
  direction[ImageDimension - 1][ImageDimension - 1] = 1;
  origin[ImageDimension - 1] = img_mov->GetOrigin()[ImageDimension - 1];
  spc[ImageDimension - 1] = img_mov->GetSpacing()[ImageDimension - 1];
  transformedvecimage->SetDirection(direction);
  transformedvecimage->SetSpacing(spc);
  transformedvecimage->SetOrigin(origin);

  antscout << " 4D-In-Spc " << img_mov->GetSpacing() << std::endl;
  antscout << " 4D-In-Org " << img_mov->GetOrigin() << std::endl;
  antscout << " 4D-In-Size " <<  img_mov->GetLargestPossibleRegion().GetSize() << std::endl;
  antscout << " 4D-In-Dir " << img_mov->GetDirection() << std::endl;
  antscout << " ...... " << std::endl;
  antscout << " 4D-Out-Spc " << transformedvecimage->GetSpacing() << std::endl;
  antscout << " 4D-Out-Org " << transformedvecimage->GetOrigin() << std::endl;
  antscout << " 4D-Out-Size " <<  transformedvecimage->GetLargestPossibleRegion().GetSize() << std::endl;
  antscout << " 4D-Out-Dir " << transformedvecimage->GetDirection() << std::endl;

  unsigned int timedims = img_mov->GetLargestPossibleRegion().GetSize()[ImageDimension - 1];
  for( unsigned int timedim = 0;  timedim < timedims;  timedim++ )
    {
    typename WarperType::Pointer  warper = WarperType::New();
    warper->SetEdgePaddingValue(0);

    if( misc_opt.use_NN_interpolator )
      {
      typedef typename itk::NearestNeighborInterpolateImageFunction<ImageType,
                                                                    typename WarperType::CoordRepType>
        NNInterpolateType;
      typename NNInterpolateType::Pointer interpolator_NN = NNInterpolateType::New();
      antscout <<  " Use Nearest Neighbor interpolation " << std::endl;
      warper->SetInterpolator(interpolator_NN);
      }

    typedef itk::TransformFileReader                    TranReaderType;
    typedef itk::ImageFileReader<DisplacementFieldType> FieldReaderType;

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
            }
          // antscout <<" aff " << transcount <<  std::endl;
          warper->PushBackAffineTransform(aff);
          if( transcount == 0 )
            {
            warper->SetOutputParametersFromImage( img_ref );
            }
          transcount++;
          }
          break;
        case IDENTITY_TRANSFORM:
          {
          typename AffineTransformType::Pointer aff;
          GetIdentityTransform(aff);
          // antscout << " aff id" << transcount << std::endl;
          warper->PushBackAffineTransform(aff);
          transcount++;
          }
          break;
        case DEFORMATION_FILE:
          {
          typename FieldReaderType::Pointer field_reader = FieldReaderType::New();
          field_reader->SetFileName( opt.filename );
          field_reader->Update();
          typename DisplacementFieldType::Pointer field = field_reader->GetOutput();

          warper->PushBackDisplacementFieldTransform(field);
          warper->SetOutputParametersFromImage( field );

          transcount++;
          }
          break;
        default:
          {
          antscout << "Unknown file type!" << std::endl;
          }
        }
      }

    // warper->PrintTransformList();
    if( img_ref.IsNotNull() )
      {
      warper->SetOutputParametersFromImage( img_ref );
      }
    else
      {
      if( misc_opt.use_TightestBoundingBox == true )
        {
        // compute the desired spacking after inputting all the transform files using the
        antscout << " not implemented " << std::endl;
        /*
          typename ImageType::SizeType largest_size;
          typename ImageType::PointType origin_warped;
          GetLaregstSizeAfterWarp(warper, warpthisimage , largest_size, origin_warped);
          warper->SetOutputParametersFromImage( warpthisimage );
          warper->SetOutputSize(largest_size);
          warper->SetOutputOrigin(origin_warped);
          {
          typename ImageType::DirectionType d;
          d.SetIdentity();
          warper->SetOutputDirection(d);
          }
          */
        }
      }

    if( timedim % vnl_math_max(timedims / 10, static_cast<unsigned int>(1) ) == 0 )
      {
      antscout << (float) timedim / (float)timedims * 100 << " % done ... " << std::flush;
      }
    typename VectorImageType::RegionType extractRegion = img_mov->GetLargestPossibleRegion();
    extractRegion.SetSize(ImageDimension - 1, 0);
    extractRegion.SetIndex(ImageDimension - 1, timedim );
    typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
    extractFilter->SetInput( img_mov );
    extractFilter->SetDirectionCollapseToSubmatrix();
    extractFilter->SetExtractionRegion( extractRegion );
    extractFilter->Update();
    typename ImageType::Pointer warpthisimage = extractFilter->GetOutput();
    typename ImageType::SpacingType qspc = warpthisimage->GetSpacing();
    typename ImageType::PointType qorg = warpthisimage->GetOrigin();
    typename ImageType::DirectionType qdir = warpthisimage->GetDirection();
    qdir.Fill(0);
    for( unsigned int qq = 0; qq < ImageDimension - 1; qq++ )
      {
      for( unsigned int pp = 0; pp < ImageDimension - 1; pp++ )
        {
        qdir[qq][pp] = img_mov->GetDirection()[qq][pp];
        }
      qspc[qq] = img_mov->GetSpacing()[qq];
      qorg[qq] = img_mov->GetOrigin()[qq];
      }
    warpthisimage->SetSpacing(qspc);
    warpthisimage->SetOrigin(qorg);
    warpthisimage->SetDirection(qdir);

    warper->SetInput( warpthisimage );
    warper->DetermineFirstDeformNoInterp();
    warper->Update();

    typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
    Iterator vfIter2(  warper->GetOutput(), warper->GetOutput()->GetLargestPossibleRegion() );
    for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
      {
      typename ImageType::PixelType  fval = vfIter2.Get();
      typename VectorImageType::IndexType ind;
      for( unsigned int xx = 0; xx < ImageDimension - 1; xx++ )
        {
        ind[xx] = vfIter2.GetIndex()[xx];
        }
      ind[ImageDimension - 1] = timedim;
      transformedvecimage->SetPixel(ind, fval);
      //    if ( ind[0] == 53 && ind[1] == 19 && ind[2] == 30 ) antscout << " fval " << fval << " td " << timedim <<
      // std::endl;
      }

    if( timedim == 0 )
      {
      antscout << warper->GetOutput()->GetDirection() << std::endl;
      }
    }
  antscout << " 100 % complete " << std::endl;
  WriteImage<VectorImageType>( transformedvecimage, output_image_filename);
}

template <int ImageDimension>
void WarpImageMultiTransform(char *moving_image_filename, char *output_image_filename,
                             TRAN_OPT_QUEUE & opt_queue, MISC_OPT & misc_opt)
{
  typedef itk::VectorImage<float,
                           ImageDimension>                        VectorImageType;
  typedef itk::Image<float,
                     ImageDimension>                              ImageType;
  typedef itk::Vector<float,
                      ImageDimension>                             VectorType;
  typedef itk::Image<VectorType,
                     ImageDimension>                              DisplacementFieldType;
  typedef itk::MatrixOffsetTransformBase<double, ImageDimension,
                                         ImageDimension>          AffineTransformType;
  typedef itk::WarpImageMultiTransformFilter<ImageType, ImageType, DisplacementFieldType,
                                             AffineTransformType> WarperType;

  itk::TransformFactory<AffineTransformType>::RegisterTransform();

  typename VectorImageType::Pointer img_mov;
  typename itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(moving_image_filename, itk::ImageIOFactory::ReadMode);
  imageIO->SetFileName(moving_image_filename);
  imageIO->ReadImageInformation();
  //    antscout << " Dimension " << imageIO->GetNumberOfDimensions()  << " Components "
  // <<imageIO->GetNumberOfComponents() << std::endl;
  unsigned int veclength = imageIO->GetNumberOfComponents();
  antscout << " read veclength as:: " << veclength << std::endl;
  ReadImage<VectorImageType>(img_mov, moving_image_filename);
  typename ImageType::Pointer img_ref;

  if( misc_opt.reference_image_filename )
    {
    ReadImage<ImageType>( img_ref, misc_opt.reference_image_filename );
    }

  typename VectorImageType::Pointer img_output =
    AllocImage<VectorImageType>(img_ref);
  img_output->SetNumberOfComponentsPerPixel(veclength);

  typename ImageType::IndexType index;
  index.Fill(0);
  typename VectorImageType::PixelType vec = img_mov->GetPixel(index);
  vec.Fill(0);
  img_output->FillBuffer( vec );
  for( unsigned int tensdim = 0;  tensdim < veclength;  tensdim++ )
    {
    typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, ImageType> IndexSelectCasterType;
    typename IndexSelectCasterType::Pointer fieldCaster = IndexSelectCasterType::New();
    fieldCaster->SetInput( img_mov );
    fieldCaster->SetIndex( tensdim );
    fieldCaster->Update();
    typename ImageType::Pointer tenscomponent = fieldCaster->GetOutput();
    tenscomponent->SetSpacing(img_mov->GetSpacing() );
    tenscomponent->SetOrigin(img_mov->GetOrigin() );
    tenscomponent->SetDirection(img_mov->GetDirection() );

    typename WarperType::Pointer  warper = WarperType::New();
    warper->SetInput(tenscomponent);
    //      PixelType nullPix;
    // nullPix.Fill(0);
    warper->SetEdgePaddingValue(0);

    if( misc_opt.use_NN_interpolator )
      {
      typedef typename itk::NearestNeighborInterpolateImageFunction<ImageType,
                                                                    typename WarperType::CoordRepType>
        NNInterpolateType;
      typename NNInterpolateType::Pointer interpolator_NN = NNInterpolateType::New();
      antscout << "Haha" << std::endl;
      warper->SetInterpolator(interpolator_NN);
      }

    typedef itk::TransformFileReader                    TranReaderType;
    typedef itk::ImageFileReader<DisplacementFieldType> FieldReaderType;

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
            }
          // antscout <<" aff " << transcount <<  std::endl;
          warper->PushBackAffineTransform(aff);
          if( transcount == 0 )
            {
            warper->SetOutputParametersFromImage( img_mov );
            }
          transcount++;
          }
          break;
        case IDENTITY_TRANSFORM:
          {
          typename AffineTransformType::Pointer aff;
          GetIdentityTransform(aff);
          // antscout << " aff id" << transcount << std::endl;
          warper->PushBackAffineTransform(aff);
          transcount++;
          }
          break;
        case DEFORMATION_FILE:
          {
          typename FieldReaderType::Pointer field_reader = FieldReaderType::New();
          field_reader->SetFileName( opt.filename );
          field_reader->Update();
          typename DisplacementFieldType::Pointer field = field_reader->GetOutput();

          warper->PushBackDisplacementFieldTransform(field);
          warper->SetOutputParametersFromImage( field );

          transcount++;
          }
          break;
        default:
          {
          antscout << "Unknown file type!" << std::endl;
          }
        }
      }

    // warper->PrintTransformList();

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
        GetLargestSizeAfterWarp<WarperType, VectorImageType>(warper, img_mov, largest_size, origin_warped);
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

    warper->DetermineFirstDeformNoInterp();
    warper->Update();

    typedef itk::ImageRegionIteratorWithIndex<VectorImageType> Iterator;
    Iterator vfIter2( img_output, img_output->GetLargestPossibleRegion() );
    for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
      {
      typename VectorImageType::PixelType  tens = vfIter2.Get();
      tens[tensdim] = warper->GetOutput()->GetPixel(vfIter2.GetIndex() );
      vfIter2.Set(tens);
      }
    }
  WriteImage<VectorImageType>(img_output, output_image_filename);
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int WarpTimeSeriesImageMultiTransform( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "WarpTimeSeriesImageMultiTransform" );

  std::remove( args.begin(), args.end(), std::string( "" ) );
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
    antscout << "\nUsage 1 (Forward warp): " << argv[0]
             <<
      " ImageDimension <moving_image.ext> <output_image.ext> -R <fixed_image.ext> <MyWarp.ext> <MyAffine.txt> [interpolation]"
             << std::endl;

    antscout << "\nUsage 2 (Inverse warp): " << argv[0]
             <<
      " ImageDimension <fixed_image.ext> <output_image.ext> -R <moving_image.ext> -i <MyAffine.txt> <MyInverseWarp.ext> [interpolation]"
             << std::endl;

    antscout << "\nUsage Information " << std::endl;
    antscout << " ImageDimension            : 3 or 4 (required argument)." << std::endl;
    antscout
      <<
      " <moving_image.ext>        : The image to apply the transformation to. The moving_image will be either a 3-D image with vector voxels or a 4D image with scalar voxels."
      << std::endl;
    antscout
      <<
      " <output_image.ext>        : The resulting image. Output will be of the same type as input, but will be resampled to the domain size defined by the -R image."
      << std::endl;
    antscout
      <<
      " <MyWarp.ext> <MyAffine.txt>    : Mappings can be stringed together, e.g.: MyAffine.txt MySecondAffine.txt MyWarp.nii.gz MySecondWarp.nii.gz -i MyInverseAffine.txt"
      << std::endl;

    antscout << "\nOptions:" << std::endl;
    antscout << " -i                : Will use the inversion of the following affine transform." << std::endl;
    antscout << " \n -R                : Reference image space that you wish to warp into." << std::endl;
    antscout
      <<
      " --reslice-by-header        : Equivalient to -i -mh, or -fh -i -mh if used together with -R. It uses the orientation matrix and origin encoded in the image file header. "
      << std::endl;
    antscout
      <<
      " --tightest-bounding-box    : Computes the tightest bounding box using all the affine transformations. It will be overrided by -R <reference_image.ext> if given."
      << std::endl;
    antscout
      << " These options can be used together with -R and are typically not used together with any other transforms."
      << std::endl;

    antscout << "\nInterpolation:" << std::endl;
    antscout << " --use-NN            : Use Nearest Neighbor Interpolator" << std::endl;
    antscout << " --use-BSpline            : Use 3rd order B-Spline Interpolation." << std::endl;

    antscout << "\n " << std::endl;
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

  int  kImageDim = atoi(argv[1]);

  const bool is_parsing_ok =
    WarpTimeSeriesImageMultiTransform_ParseInput(argc - 2, argv + 2, moving_image_filename, output_image_filename,
                                                 opt_queue,
                                                 misc_opt);

  if( is_parsing_ok )
    {
    antscout << "moving_image_filename: " << moving_image_filename << std::endl;
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

    switch( kImageDim )
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
      case 4:
        {
        WarpImageMultiTransformFourD<4>(moving_image_filename, output_image_filename, opt_queue, misc_opt);
        }
        break;
      }
    }
  else
    {
    antscout << "Input error!" << std::endl;
    }
  return EXIT_FAILURE;
}
} // namespace ants
