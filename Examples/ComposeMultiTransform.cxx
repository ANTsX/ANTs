
#include "antsUtilities.h"
#include "antsUtilities.h"
#include "itkImageFileReader.h"

#include "itkImageFileWriter.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkTransformFactory.h"

#include "itkDisplacementFieldFromMultiTransformFilter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"

namespace ants
{
static bool ComposeMultiTransform_ParseInput(int argc, char * *argv, char *& output_image_filename,
                                             char *& reference_image_filename, TRAN_OPT_QUEUE & opt_queue)
{
  opt_queue.clear();
  opt_queue.reserve(argc - 2);

  output_image_filename = argv[0];

  reference_image_filename = NULL;

  int ind = 1;
  while( ind < argc )
    {
    if( strcmp(argv[ind], "-R") == 0 )
      {
      ind++;
      if( ind >= argc )
        {
        return false;
        }
      reference_image_filename = argv[ind];
      }
    else if( strcmp(argv[ind], "-i") == 0 )
      {
      ind++;
      if( ind >= argc )
        {
        return false;
        }
      TRAN_OPT opt;
      opt.filename = argv[ind];
      if( CheckFileType(opt.filename) != AFFINE_FILE )
        {
        std::cout << "file: " << opt.filename
                 << " is not an affine .txt file. Invalid to use '-i' "
                 << std::endl;
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

//    if (reference_image_filename == NULL) {
//        std::cout << "the reference image file (-R) must be given!!!"
//        << std::endl;
//        return false;
//    }

  return true;
}

template <int ImageDimension>
void ComposeMultiTransform(char *output_image_filename,
                           char *reference_image_filename, TRAN_OPT_QUEUE & opt_queue)
{
  typedef itk::Image<float, ImageDimension>      ImageType;
  typedef itk::Vector<float, ImageDimension>     VectorType;
  typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;
  typedef itk::MatrixOffsetTransformBase<double, ImageDimension,
                                         ImageDimension> AffineTransformType;
  // typedef itk::WarpImageMultiTransformFilter<ImageType,ImageType, DisplacementFieldType, AffineTransformType>
  // WarperType;
  typedef itk::DisplacementFieldFromMultiTransformFilter<DisplacementFieldType,
                                                         DisplacementFieldType, AffineTransformType> WarperType;

  itk::TransformFactory<AffineTransformType>::RegisterTransform();

  typedef itk::ImageFileReader<ImageType> ImageFileReaderType;
  typename ImageFileReaderType::Pointer reader_img =
    ImageFileReaderType::New();
  typename ImageType::Pointer img_ref;

  typename ImageFileReaderType::Pointer reader_img_ref =
    ImageFileReaderType::New();
  if( reference_image_filename )
    {
    reader_img_ref->SetFileName(reference_image_filename);
    reader_img_ref->Update();
    img_ref = reader_img_ref->GetOutput();
    }
  else
    {
    std::cout << "the reference image file (-R) must be given!!!"
             << std::endl;
    return;
    }

  typename WarperType::Pointer warper = WarperType::New();
  // warper->SetInput(img_mov);
  // warper->SetEdgePaddingValue( 0);
  VectorType pad;
  pad.Fill(0);
  // warper->SetEdgePaddingValue(pad);

  typedef itk::TransformFileReader TranReaderType;

  typedef itk::ImageFileReader<DisplacementFieldType>
    FieldReaderType;

  const int kOptQueueSize = opt_queue.size();
  for( int i = 0; i < kOptQueueSize; i++ )
    {
    const TRAN_OPT & opt = opt_queue[i];

    switch( opt_queue[i].file_type )
      {
      case AFFINE_FILE:
        {
        typename TranReaderType::Pointer tran_reader = TranReaderType::New();
        tran_reader->SetFileName(opt.filename);
        tran_reader->Update();
        typename AffineTransformType::Pointer aff =
          dynamic_cast<AffineTransformType *>( (tran_reader->GetTransformList() )->front().GetPointer() );
        if( opt_queue[i].do_affine_inv )
          {
          aff->GetInverse(aff);
          }
        // std::cout << aff << std::endl;
        warper->PushBackAffineTransform(aff);
        }
        break;
      case DEFORMATION_FILE:
        {
        typename FieldReaderType::Pointer field_reader = FieldReaderType::New();
        field_reader->SetFileName(opt.filename);
        field_reader->Update();
        typename DisplacementFieldType::Pointer field = field_reader->GetOutput();
        // std::cout << field << std::endl;
        warper->PushBackDisplacementFieldTransform(field);
        }
        break;
      default:
        std::cout << "Unknown file type!" << std::endl;
      }
    }

  warper->SetOutputParametersFromImage( img_ref );
  std::cout << "output size: " << warper->GetOutputSize() << std::endl;
  std::cout << "output spacing: " << warper->GetOutputSpacing() << std::endl;

  // warper->PrintTransformList();
  warper->DetermineFirstDeformNoInterp();
  warper->Update();

  typename DisplacementFieldType::Pointer field_output =
    DisplacementFieldType::New();
  field_output = warper->GetOutput();

  std::string            filePrefix = output_image_filename;
  std::string::size_type pos = filePrefix.rfind(".");
  std::string            extension = std::string(filePrefix, pos, filePrefix.length()
                                                 - 1);
  filePrefix = std::string(filePrefix, 0, pos);

  std::cout << "output extension is: " << extension << std::endl;

    {
    typedef itk::ImageFileWriter<DisplacementFieldType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(output_image_filename);
    writer->SetInput(field_output);
    writer->Update();
    }
}

template <int ImageDimension>
void ComposeMultiAffine(char *output_affine_txt,
                        char *reference_affine_txt, TRAN_OPT_QUEUE & opt_queue)
{
  typedef itk::Image<float, ImageDimension>      ImageType;
  typedef itk::Vector<float, ImageDimension>     VectorType;
  typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;

  typedef itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension> AffineTransformType;
  // MatrixOffsetTransformBase is not usually registered, so register it here.
  // MatrixOffsetTransformBase should NOT be a valid transform type for writting,
  // but it is needed for historical reading purposes.
  itk::TransformFactory<AffineTransformType>::RegisterTransform();

  typedef itk::TransformFileReader TranReaderType;

  typedef itk::WarpImageMultiTransformFilter<ImageType, ImageType, DisplacementFieldType,
                                             AffineTransformType> WarperType;
  typename WarperType::Pointer warper = WarperType::New();
  bool      has_affine_tranform = false;
  const int kOptQueueSize = opt_queue.size();
  for( int i = 0; i < kOptQueueSize; i++ )
    {
    const TRAN_OPT & opt = opt_queue[i];

    switch( opt_queue[i].file_type )
      {
      case AFFINE_FILE:
        {
        typename TranReaderType::Pointer tran_reader = TranReaderType::New();
        tran_reader->SetFileName(opt.filename);
        tran_reader->Update();
        typename AffineTransformType::Pointer aff =
          dynamic_cast<AffineTransformType *>( (tran_reader->GetTransformList() )->front().GetPointer() );
        if( opt_queue[i].do_affine_inv )
          {
          aff->GetInverse(aff);
          }
        warper->PushBackAffineTransform(aff);
        has_affine_tranform = true;
        }
        break;
      case DEFORMATION_FILE:
        {
        std::cout << "Compose affine only files: ignore "
                 << opt.filename << std::endl;
        }
        break;
      default:
        {
        std::cout << "Unknown file type!" << std::endl;
        }
      }
    }

  typename AffineTransformType::Pointer aff_ref_tmp;
  if( reference_affine_txt )
    {
    typename TranReaderType::Pointer tran_reader = TranReaderType::New();
    tran_reader->SetFileName(reference_affine_txt);
    tran_reader->Update();
    aff_ref_tmp = dynamic_cast<AffineTransformType *>( (tran_reader->GetTransformList() )->front().GetPointer() );
    }
  else
    {
    if( has_affine_tranform == true )
      {
      std::cout << "the reference affine file for center is selected as the first affine!" << std::endl;
      aff_ref_tmp = ( (warper->GetTransformList() ).begin() )->second.aex.aff;
      }
    else
      {
      std::cout << "No affine input is given. nothing to do ......" << std::endl;
      return;
      }
    }
    {
    typedef typename AffineTransformType::CenterType PointType;
    const PointType aff_center = aff_ref_tmp->GetCenter();
    std::cout << "new center is : " << aff_center << std::endl;
      {
      typename AffineTransformType::Pointer aff_output = AffineTransformType::New();
      warper->ComposeAffineOnlySequence(aff_center, aff_output);
      typedef itk::TransformFileWriter TranWriterType;
      typename TranWriterType::Pointer tran_writer = TranWriterType::New();
      tran_writer->SetFileName(output_affine_txt);
      tran_writer->SetInput(aff_output);
      tran_writer->Update();
      }
    }
  std::cout << "wrote file to : " << output_affine_txt << std::endl;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int ComposeMultiTransform( std::vector<std::string> args, std::ostream* /*out_stream = NULL */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "ComposeMultiTransform" );
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

  // antscout->set_stream( out_stream );

  if( argc <= 3 )
    {
    std::cout
      << "ComposeMultiTransform ImageDimension output_field [-R reference_image] "
      << "{[deformation_field | [-i] affine_transform_txt ]}"
      << std::endl;
    std::cout << "  Usage has the same form as WarpImageMultiTransform " << std::endl;
    std::cout << " For Example: " << std::endl;
    std::cout << std::endl;
    std::cout <<   argv[0]  << " Dimension  outwarp.nii   -R template.nii   ExistingWarp.nii  ExistingAffine.nii "
             << std::endl;
    std::cout << " or for an inverse mapping : " << std::endl;
    std::cout << argv[0]
             << " Dimension  outwarp.nii   -R template.nii   -i ExistingAffine.nii ExistingInverseWarp.nii "
             << std::endl;
    std::cout << " recalling that the -i option takes the inverse of the affine mapping " << std::endl;
    std::cout << std::endl;
    std::cout << "Or: to compose multiple affine text file into one: "        << std::endl;
    std::cout << "ComposeMultiTransform ImageDimension output_affine_txt [-R reference_affine_txt] "
             << "{[-i] affine_transform_txt}" << std::endl
             << "This will be evoked if a text file is given as the second parameter. In this case "
             << "reference_affine_txt is used to define the center of the output affine.  "
             << "The default reference is the first given affine text file. "
             << "This ignores all non-txt files among the following parameters."
             << std::endl;
    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }

  TRAN_OPT_QUEUE opt_queue;
  //    char *moving_image_filename = NULL;
  char *output_image_filename = NULL;
  char *reference_image_filename = NULL;

  int  kImageDim = atoi(argv[1]);

  const bool is_parsing_ok = ComposeMultiTransform_ParseInput(argc - 2, argv + 2, output_image_filename,
                                                   reference_image_filename, opt_queue);

  if( is_parsing_ok )
    {
    switch( CheckFileType(output_image_filename) )
      {
      case DEFORMATION_FILE:
        {
        if( reference_image_filename == NULL )
          {
          std::cout << "the reference image file (-R) must be given!!!"
                   << std::endl;
          return false;
          }

        std::cout << "output_image_filename: " << output_image_filename
                 << std::endl;
        std::cout << "reference_image_filename: ";
        if( reference_image_filename )
          {
          std::cout << reference_image_filename << std::endl;
          }
        else
          {
          std::cout << "NULL" << std::endl;
          }
        DisplayOptQueue(opt_queue);

        switch( kImageDim )
          {
          case 2:
            {
            ComposeMultiTransform<2>(output_image_filename,
                                     reference_image_filename, opt_queue);
            }
            break;
          case 3:
            {
            ComposeMultiTransform<3>(output_image_filename,
                                     reference_image_filename, opt_queue);
            }
            break;
          }
        }
        break;
      case AFFINE_FILE:
        {
        std::cout << "output_affine_txt: " << output_image_filename
                 << std::endl;
        std::cout << "reference_affine_txt: ";
        if( reference_image_filename )
          {
          std::cout << reference_image_filename << std::endl;
          }
        else
          {
          std::cout << "NULL" << std::endl;
          }
        DisplayOptQueue(opt_queue);

        switch( kImageDim )
          {
          case 2:
            {
            ComposeMultiAffine<2>(output_image_filename,
                                  reference_image_filename, opt_queue);
            }
            break;
          case 3:
            {
            ComposeMultiAffine<3>(output_image_filename,
                                  reference_image_filename, opt_queue);
            }
            break;
          }
        }
        break;
      default:
        {
        std::cout << "Unknow output file format: " << output_image_filename << std::endl;
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
