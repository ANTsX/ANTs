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
#include "itkPreservationOfPrincipalDirectionTensorReorientationImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkWarpTensorImageMultiTransformFilter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"

typedef enum { INVALID_FILE = 1, AFFINE_FILE, DEFORMATION_FILE, IMAGE_AFFINE_HEADER,
               IDENTITY_TRANSFORM } TRAN_FILE_TYPE;
typedef struct
  {
  std::string filename;
  TRAN_FILE_TYPE file_type;
  bool do_affine_inv;
  } TRAN_OPT;

typedef std::vector<TRAN_OPT> TRAN_OPT_QUEUE;

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
                TRAN_OPT_QUEUE & opt_queue)
{
  opt_queue.clear();
  opt_queue.reserve(argc - 2);

  moving_image_filename = argv[0];
  output_image_filename = argv[1];

  int ind = 2;

  while( ind < argc )
    {
    TRAN_OPT opt;
    opt.filename = argv[ind];
    opt.file_type = CheckFileType(opt.filename.c_str() );
    opt.do_affine_inv = false;
    bool set_current_affine_inv = false;

    if( strcmp(argv[ind], "-i") == 0 )
      {
      std::cout << "ERROR - inverse transforms not yet supported\n" << std::endl;
      return false;
      }
    else
      {
      if( opt.file_type == AFFINE_FILE )
        {
        SetAffineInvFlag(opt, set_current_affine_inv);
        }
      else
        {
        if( opt.file_type == DEFORMATION_FILE && set_current_affine_inv )
          {
          std::cout << "Ignore inversion of non-affine file type! " << std::endl;
          std::cout << "opt.do_affine_inv:" << opt.do_affine_inv << std::endl;
          }
        }

      opt_queue.push_back(opt);
      DisplayOpt(opt);
      }
    ++ind;
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

template <int ImageDimension>
void ReorientTensorImage(char *moving_image_filename, char *output_image_filename, TRAN_OPT_QUEUE & opt_queue)
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

  // No reason to use log-euclidean space
  ReadTensorImage<TensorImageType>(img_mov, moving_image_filename, false);

  typename ImageType::Pointer img_ref = NULL;

  typename ImageFileReaderType::Pointer reader_img_ref = ImageFileReaderType::New();

  typedef itk::TransformFileReader                   TranReaderType;
  typedef itk::ImageFileReader<DeformationFieldType> FieldReaderType;
  typename DeformationFieldType::Pointer field = NULL;
  typename AffineTransformType::Pointer aff = NULL;

  const int kOptQueueSize = opt_queue.size();

  if( kOptQueueSize > 1 )
    {
    std::cout << "ERROR: Only 1 input transform is permitted" << std::endl;
    return;
    }

  typedef itk::PreservationOfPrincipalDirectionTensorReorientationImageFilter<TensorImageType,
                                                                              DeformationFieldType> PPDReorientType;
  typename PPDReorientType::Pointer reo = PPDReorientType::New();
  reo->SetInput( img_mov );

  const TRAN_OPT & opt = opt_queue[0];

  switch( opt.file_type )
    {
    case AFFINE_FILE:
      {
      typename TranReaderType::Pointer tran_reader = TranReaderType::New();
      tran_reader->SetFileName(opt.filename);
      tran_reader->Update();
      aff = dynamic_cast<AffineTransformType *>( (tran_reader->GetTransformList() )->front().GetPointer() );
      reo->SetAffineTransform( aff );

      std::cout << "Affine transform" << std::endl;
      break;
      }

    case DEFORMATION_FILE:
      {
      typename FieldReaderType::Pointer field_reader = FieldReaderType::New();
      field_reader->SetFileName( opt.filename );
      field_reader->Update();
      // field = field_reader->GetOutput();
      reo->SetDeformationField( field_reader->GetOutput() );
      std::cout << "Warp transform" << std::endl;
      break;
      }
    default:
      std::cout << "Unknown file type!" << std::endl;
    }

  reo->Update();

  typename TensorImageType::Pointer img_output = reo->GetOutput();

  // No reason to use log-euclidean space here
  WriteTensorImage<TensorImageType>(img_output, output_image_filename, false);
}

int main(int argc, char *argv[])
{
  if( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " Dimension infile.nii outfile.nii <warp.nii/affine.txt> " << std::endl;
    return 1;
    }

  TRAN_OPT_QUEUE opt_queue;
  char *         moving_image_filename = NULL;
  char *         output_image_filename = NULL;

  bool is_parsing_ok = false;
  int  dim = atoi(argv[1]);

  if( dim != 3 )
    {
    std::cerr << "ReorientTensorImage only supports 3D image volumes" << std::endl;
    exit(1);
    }

  is_parsing_ok = ParseInput(argc - 2, argv + 2, moving_image_filename, output_image_filename, opt_queue);

  if( is_parsing_ok )
    {
    std::cout << "moving_image_filename: " << moving_image_filename << std::endl;
    std::cout << "output_image_filename: " << output_image_filename << std::endl;
    DisplayOptQueue(opt_queue);

    ReorientTensorImage<3>(moving_image_filename, output_image_filename, opt_queue);
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
