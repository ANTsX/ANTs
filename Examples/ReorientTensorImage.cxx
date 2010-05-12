/*=========================================================================1

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: ReorientTensorImage.cxx,v $
  Language:  C++
  Date:      $Date: 2009/03/17 18:55:26 $
  Version:   $Revision: 1.2 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "ReadWriteImage.h"
#include "TensorFunctions.h"
#include "itkPreservationOfPrincipalDirectionTensorReorientationImageFilter.h"
#include "itkVectorImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVectorResampleImageFilter.h"
#include "itkWarpTensorImageMultiTransformFilter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"
#include "itkVectorNearestNeighborInterpolateImageFunction.h"

typedef enum { INVALID_FILE = 1, AFFINE_FILE, DEFORMATION_FILE, IMAGE_AFFINE_HEADER,
               IDENTITY_TRANSFORM } TRAN_FILE_TYPE;
typedef struct
  {
  std::string filename;
  TRAN_FILE_TYPE file_type;
  bool do_affine_inv;
  } TRAN_OPT;

typedef std::vector<TRAN_OPT> TRAN_OPT_QUEUE;

typedef struct
  {
  bool use_NN_interpolator;
  bool use_TightestBoundingBox;
  char * reference_image_filename;
  bool use_RotationHeader;
  } MISC_OPT;

void DisplayOptQueue(const TRAN_OPT_QUEUE & opt_queue);

void DisplayOpt(const TRAN_OPT & opt);

TRAN_FILE_TYPE CheckFileType(const char *str)
{
  std::string            filename = str;
  std::string::size_type pos = filename.rfind( "." );
  std::string            filepre = std::string( filename, 0, pos );

  if( pos != std::string::npos )
    {
    std::string extension = std::string( filename, pos, filename.length() - 1);
    if( extension == std::string(".gz") )
      {
      pos = filepre.rfind( "." );
      extension = std::string( filepre, pos, filepre.length() - 1 );
      }
    if( extension == ".txt" )
      {
      return AFFINE_FILE;
      }
    else
      {
      return DEFORMATION_FILE;
      }
    }
  else
    {
    return INVALID_FILE;
    }
  return AFFINE_FILE;
}

void FilePartsWithgz(const std::string & filename, std::string & path, std::string & name, std::string & ext)
{
  std::string            extension;
  std::string::size_type pos = filename.rfind( "." );
  std::string            filepre = std::string( filename, 0, pos );

  if( pos != std::string::npos )
    {
    extension = std::string( filename, pos, filename.length() - 1);
    if( extension == std::string(".gz") )
      {
      pos = filepre.rfind( "." );
      if( pos != std::string::npos )
        {
        extension = std::string( filepre, pos, filepre.length() - 1 ) + ".gz";
        filepre = std::string(filepre, 0, pos);
        }
      }
    }
  else
    {
    extension = std::string("");
    }

  ext = extension;

  pos = filepre.rfind('/');

  if( pos != std::string::npos )
    {
    path = std::string(filepre, 0, pos + 1);
    name = std::string(filepre, pos + 1, filepre.length() - 1);
    }
  else
    {
    path = std::string("");
    name = filepre;
    }

//    std::cout << "filename: " << filename << std::endl
//    << "path: " << path << std::endl
//    << "name: " << name << std::endl
//    << "ext: " << ext << std::endl;
}

bool CheckFileExistence(const char *str)
{
  std::ifstream myfile(str);
  bool          b = myfile.is_open();

  myfile.close();
  return b;
}

void SetAffineInvFlag(TRAN_OPT & opt, bool & set_current_affine_inv)
{
  opt.do_affine_inv = set_current_affine_inv;
  if( set_current_affine_inv )
    {
    set_current_affine_inv = false;
    }
}

bool ParseInput(int argc, char * *argv, char *& moving_image_filename,
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
        std::cout << "found deformation file: " << opt.filename << std::endl;
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
        std::cout << "found affine file: " << opt.filename << std::endl;
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
        std::cout << "found affine file: " << opt.filename << std::endl;
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
        std::cout << "found deformation file: " << opt.filename << std::endl;
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
        std::cout << "Ignore inversion of non-affine file type! " << std::endl;
        std::cout << "opt.do_affine_inv:" << opt.do_affine_inv << std::endl;
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
    //               std::cout << "Use Rotation Header!" << std::endl;
    }

  return true;
}

void DisplayOptQueue(const TRAN_OPT_QUEUE & opt_queue)
{
  const int kQueueSize = opt_queue.size();

  for( int i = 0; i < kQueueSize; i++ )
    {
    std::cout << "[" << i << "/" << kQueueSize << "]: ";

    switch( opt_queue[i].file_type )
      {
      case AFFINE_FILE:
        {
        std::cout << "AFFINE";
        }
        break;
      case DEFORMATION_FILE:
        {
        std::cout << "FIELD";
        }
        break;
      case IDENTITY_TRANSFORM:
        {
        std::cout << "IDENTITY";
        }
        break;
      case IMAGE_AFFINE_HEADER:
        {
        std::cout << "HEADER";
        }
        break;
      default:
        {
        std::cout << "Invalid Format!!!";
        }
        break;
      }
    if( opt_queue[i].do_affine_inv )
      {
      std::cout << "-INV";
      }
    std::cout << ": " << opt_queue[i].filename << std::endl;
    }
}

void DisplayOpt(const TRAN_OPT & opt)
{
  switch( opt.file_type )
    {
    case AFFINE_FILE:
      {
      std::cout << "AFFINE";
      }
      break;
    case DEFORMATION_FILE:
      {
      std::cout << "FIELD";
      }
      break;
    case IDENTITY_TRANSFORM:
      {
      std::cout << "IDENTITY";
      }
      break;
    case IMAGE_AFFINE_HEADER:
      {
      std::cout << "HEADER";
      }
      break;
    default:
      {
      std::cout << "Invalid Format!!!";
      }
      break;
    }
  if( opt.do_affine_inv )
    {
    std::cout << "-INV";
    }
  std::cout << ": " << opt.filename << std::endl;
}

template <class AffineTransformPointer>
void GetIdentityTransform(AffineTransformPointer & aff)
{
  typedef typename AffineTransformPointer::ObjectType AffineTransform;
  aff = AffineTransform::New();
  aff->SetIdentity();
}

template <class ImageTypePointer, class AffineTransformPointer>
void GetAffineTransformFromImage(const ImageTypePointer& img, AffineTransformPointer & aff)
{
  typedef typename ImageTypePointer::ObjectType                        ImageType;
  typedef typename ImageType::DirectionType                            DirectionType;
  typedef typename ImageType::PointType                                PointType;
  typedef typename ImageType::SpacingType                              SpacingType;
  typedef typename AffineTransformPointer::ObjectType::TranslationType VectorType;

  DirectionType direction = img->GetDirection();

  SpacingType spacing = img->GetSpacing();
  VectorType  translation;
  // translation.Fill(0);
  for( unsigned int i = 0; i < ImageType::GetImageDimension(); i++ )
    {
    translation[i] = img->GetOrigin()[i];
    }

  aff->SetMatrix(direction);
  // aff->SetCenter(pt);
  PointType pt; pt.Fill(0);
  aff->SetOffset(translation);
  aff->SetCenter(pt);

  std::cout << "aff from image:" << aff << std::endl;
}

template <class WarperPointerType, class ImagePointerType, class SizeType, class PointType>
void GetLaregstSizeAfterWarp(WarperPointerType & warper, ImagePointerType & img, SizeType & largest_size,
                             PointType & origin_warped)
{
  typedef typename ImagePointerType::ObjectType ImageType;
  const int ImageDimension = ImageType::GetImageDimension();

  // typedef typename ImageType::PointType PointType;
  typedef typename std::vector<PointType> PointList;

  typedef typename ImageType::IndexType IndexType;

  // PointList pts_orig;
  PointList pts_warped;

  typename ImageType::SizeType imgsz;
  imgsz = img->GetLargestPossibleRegion().GetSize();

  typename ImageType::SpacingType spacing;
  spacing = img->GetSpacing();

  pts_warped.clear();
  if( ImageDimension == 3 )
    {
    for( int i = 0; i < 8; i++ )
      {
      IndexType ind;

      switch( i )
        {
        case 0:
      { ind[0] = 0; ind[1] = 0; ind[2] = 0; }
                                            break;
        case 1:
      { ind[0] = imgsz[0] - 1; ind[1] = 0; ind[2] = 0; }
                                                       break;
        case 2:
      { ind[0] = 0; ind[1] = imgsz[1] - 1; ind[2] = 0; }
                                                       break;
        case 3:
      { ind[0] = imgsz[0] - 1; ind[1] = imgsz[1] - 1; ind[2] = 0; }
                                                                  break;
        case 4:
      { ind[0] = 0; ind[1] = 0; ind[2] = imgsz[2] - 1; }
                                                       break;
        case 5:
      { ind[0] = imgsz[0] - 1; ind[1] = 0; ind[2] = imgsz[2] - 1; }
                                                                  break;
        case 6:
      { ind[0] = 0; ind[1] = imgsz[1] - 1; ind[2] = imgsz[2] - 1; }
                                                                  break;
        case 7:
      { ind[0] = imgsz[0] - 1; ind[1] = imgsz[1] - 1; ind[2] = imgsz[2] - 1; }
                                                                             break;
        }
      PointType pt_orig, pt_warped;
      img->TransformIndexToPhysicalPoint(ind, pt_orig);
      if( warper->MultiInverseAffineOnlySinglePoint(pt_orig, pt_warped) == false )
        {
        std::cout << "ERROR: outside of numeric boundary with affine transform." << std::endl;
        exit(-1);
        }
      pts_warped.push_back(pt_warped);
      std::cout << '[' << i << ']' << ind << ',' << pt_orig << "->" << pt_warped << std::endl;
      }
    }
  else if( ImageDimension == 2 )
    {
    for( int i = 0; i < 4; i++ )
      {
      IndexType ind;

      switch( i )
        {
        case 0:
      { ind[0] = 0; ind[1] = 0; }
                                break;
        case 1:
      { ind[0] = imgsz[0] - 1; ind[1] = 0; }
                                           break;
        case 2:
      { ind[0] = 0; ind[1] = imgsz[1] - 1; }
                                           break;
        case 3:
      { ind[0] = imgsz[0] - 1; ind[1] = imgsz[1] - 1; }
                                                      break;
        }
      PointType pt_orig, pt_warped;
      img->TransformIndexToPhysicalPoint(ind, pt_orig);
      if( warper->MultiInverseAffineOnlySinglePoint(pt_orig, pt_warped) == false )
        {
        std::cout << "ERROR: outside of numeric boundary with affine transform." << std::endl;
        exit(-1);
        }
      pts_warped.push_back(pt_warped);
      std::cout << '[' << i << ']' << ind << ',' << pt_orig << "->" << pt_warped << std::endl;
      }
    }
  else
    {
    std::cout << "could not determine the dimension after warping for non 2D/3D volumes" << std::endl;
    exit(-1);
    }

  PointType pt_min, pt_max;
  pt_min = pts_warped[0];
  pt_max = pts_warped[0];
  for( unsigned int k = 0; k < pts_warped.size(); k++ )
    {
    for( int i = 0; i < ImageDimension; i++ )
      {
      pt_min[i] = (pt_min[i] < pts_warped[k][i]) ? (pt_min[i]) : (pts_warped[k][i]);
      pt_max[i] = (pt_max[i] > pts_warped[k][i]) ? (pt_max[i]) : (pts_warped[k][i]);
      }
    }
  for( int i = 0; i < ImageDimension; i++ )
    {
    largest_size[i] = (int) (ceil( (pt_max[i] - pt_min[i]) / spacing[i]) + 1);
    }

  origin_warped = pt_min;
  std::cout << "origin_warped: " << origin_warped << std::endl;
  std::cout << "pt_min: " << pt_min << " pt_max:" << pt_max << " largest_size:" << largest_size << std::endl;
}

template <int ImageDimension>
void ReorientTensorImage(char *moving_image_filename, char *output_image_filename,
                         TRAN_OPT_QUEUE & opt_queue, MISC_OPT & misc_opt)
{
  typedef itk::SymmetricSecondRankTensor<float, 3>                               TensorType;
  typedef itk::SymmetricSecondRankTensor<float, 3>                               PixelType;
  typedef itk::Image<PixelType, ImageDimension>                                  TensorImageType;
  typedef itk::Image<float, ImageDimension>                                      ImageType;
  typedef itk::Vector<float, ImageDimension>                                     VectorType;
  typedef itk::Image<VectorType, ImageDimension>                                 DeformationFieldType;
  typedef itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension> AffineTransformType;
  itk::TransformFactory<AffineTransformType>::RegisterTransform();

  typedef itk::ImageFileReader<ImageType> ImageFileReaderType;
  typename TensorImageType::Pointer img_mov;
  ReadTensorImage<TensorImageType>(img_mov, moving_image_filename, true);

  typename ImageType::Pointer img_ref;

  typename ImageFileReaderType::Pointer reader_img_ref = ImageFileReaderType::New();
  if( misc_opt.reference_image_filename )
    {
    reader_img_ref->SetFileName(misc_opt.reference_image_filename);
    reader_img_ref->Update();
    img_ref = reader_img_ref->GetOutput();
    }
  // else
  //    img_ref = NULL;

  typedef itk::TransformFileReader                                    TranReaderType;
  typedef itk::VectorImageFileReader<ImageType, DeformationFieldType> FieldReaderType;
  typename DeformationFieldType::Pointer field;

  unsigned int transcount = 0;
  const int    kOptQueueSize = opt_queue.size();
  for( int i = 0; i < kOptQueueSize; i++ )
    {
    const TRAN_OPT & opt = opt_queue[i];

    switch( opt.file_type )
      {
      case DEFORMATION_FILE:
        {
        typename FieldReaderType::Pointer field_reader = FieldReaderType::New();
        field_reader->SetFileName( opt.filename );
        field_reader->Update();
        field = field_reader->GetOutput();

        transcount++;
        break;
        }
      default:
        std::cout << "Unknown file type!" << std::endl;
      }
    }

  // warper->PrintTransformList();

  typedef itk::PreservationOfPrincipalDirectionTensorReorientationImageFilter<TensorImageType, DeformationFieldType>
    ReorientType;
  typename ReorientType::Pointer reo = ReorientType::New();
  reo->SetDeformationField( field );
  reo->SetInput( img_mov );
  reo->Update();

  typename TensorImageType::Pointer img_output = reo->GetOutput();
  WriteTensorImage<TensorImageType>(img_output, output_image_filename, true);
}

int main(int argc, char *argv[])
{
  std::cout << " Does not take into account reorientation needed when orientations change only in the header!! "
            << std::endl;
  std::cout << " consider the same DT image in 2 different orientations under an applied identity transform. "
            << std::endl;
  std::cout
    <<
  " the components will not be rotated correctly, but should be, b/c the header rotation is not accounted for in the reorientation filter."
    << std::endl;
  std::cout << " this is a bug we need to fix."  << std::endl;
  std::cout << " ... "   << std::endl;

  if( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " dimension infile.nii outfile.nii warp.nii " << std::endl;
    return 1;
    }

  TRAN_OPT_QUEUE opt_queue;
  char *         moving_image_filename = NULL;
  char *         output_image_filename = NULL;

  MISC_OPT misc_opt;

  bool is_parsing_ok = false;
  int  kImageDim = atoi(argv[1]);

  is_parsing_ok = ParseInput(argc - 2, argv + 2, moving_image_filename, output_image_filename, opt_queue, misc_opt);

  if( is_parsing_ok )
    {
    std::cout << "moving_image_filename: " << moving_image_filename << std::endl;
    std::cout << "output_image_filename: " << output_image_filename << std::endl;
    std::cout << "reference_image_filename: ";
    if( misc_opt.reference_image_filename )
      {
      std::cout << misc_opt.reference_image_filename << std::endl;
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
        // WarpImageMultiTransform<2>(moving_image_filename, output_image_filename, opt_queue, misc_opt);
        break;
        }
      case 3:
        {
        ReorientTensorImage<3>(moving_image_filename, output_image_filename, opt_queue, misc_opt);
        break;
        }
      }
    }
  else
    {
    std::cout << "Input error!" << std::endl;
    }

  exit(0);

  // ReorientTensorImage<3>(argc,argv);
//  WarpImageForward(argc,argv);
  return 0;
}
