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

#include "antsUtilities.h"
#include "antsUtilities.h"

#include "ReadWriteImage.h"
#include "itkPreservationOfPrincipalDirectionTensorReorientationImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkWarpTensorImageMultiTransformFilter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"

namespace ants
{
static bool ReorientTensorImage_ParseInput(int argc, char * *argv, char *& moving_image_filename,
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

    if( strcmp(argv[ind], "-i") == 0 )
      {
      antscout << "ERROR - inverse transforms not yet supported\n" << std::endl;
      return false;
      }
    else
      {
      bool set_current_affine_inv = false;
      if( opt.file_type == AFFINE_FILE )
        {
        SetAffineInvFlag(opt, set_current_affine_inv);
        }
      else
        {
        if( opt.file_type == DEFORMATION_FILE && set_current_affine_inv )
          {
          antscout << "Ignore inversion of non-affine file type! " << std::endl;
          antscout << "opt.do_affine_inv:" << opt.do_affine_inv << std::endl;
          }
        }

      opt_queue.push_back(opt);
      DisplayOpt(opt);
      }
    ++ind;
    }

  return true;
}

template <int ImageDimension>
void ReorientTensorImage(char *moving_image_filename, char *output_image_filename, TRAN_OPT_QUEUE & opt_queue)
{
  typedef itk::DiffusionTensor3D<double>                                         TensorType;
  typedef itk::DiffusionTensor3D<double>                                         PixelType;
  typedef itk::Image<PixelType, ImageDimension>                                  TensorImageType;
  typedef itk::Image<float, ImageDimension>                                      ImageType;
  typedef itk::Vector<double, ImageDimension>                                    VectorType;
  typedef itk::Image<VectorType, ImageDimension>                                 DisplacementFieldType;
  typedef itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension> AffineTransformType;
  itk::TransformFactory<AffineTransformType>::RegisterTransform();

  typedef itk::ImageFileReader<ImageType> ImageFileReaderType;
  typename TensorImageType::Pointer img_mov;

  // No reason to use log-euclidean space
  ReadTensorImage<TensorImageType>(img_mov, moving_image_filename, false);

  typename ImageType::Pointer img_ref = NULL;

  typename ImageFileReaderType::Pointer reader_img_ref = ImageFileReaderType::New();

  typedef itk::TransformFileReader                    TranReaderType;
  typedef itk::ImageFileReader<DisplacementFieldType> FieldReaderType;
  typename DisplacementFieldType::Pointer field = NULL;
  typename AffineTransformType::Pointer aff = NULL;

  const int kOptQueueSize = opt_queue.size();

  if( kOptQueueSize > 1 )
    {
    antscout << "ERROR: Only 1 input transform is permitted" << std::endl;
    return;
    }

  typedef itk::PreservationOfPrincipalDirectionTensorReorientationImageFilter<TensorImageType,
                                                                              DisplacementFieldType> PPDReorientType;
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
      antscout << "Affine transform" << std::endl;
      }
      break;
    case DEFORMATION_FILE:
      {
      typename FieldReaderType::Pointer field_reader = FieldReaderType::New();
      field_reader->SetFileName( opt.filename );
      field_reader->Update();
      // field = field_reader->GetOutput();
      reo->SetDisplacementField( field_reader->GetOutput() );
      antscout << "Warp transform" << std::endl;
      }
      break;
    default:
      {
      antscout << "Unknown file type!" << std::endl;
      }
    }
  reo->Update();

  typename TensorImageType::Pointer img_output = reo->GetOutput();
  // No reason to use log-euclidean space here
  WriteTensorImage<TensorImageType>(img_output, output_image_filename, false);
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int ReorientTensorImage( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "ReorientTensorImage" );

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

  if( argc < 4 )
    {
    antscout << "Usage: " << argv[0] << " Dimension infile.nii outfile.nii <warp.nii/affine.txt> " << std::endl;
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

  bool is_parsing_ok = false;
  int  dim = atoi(argv[1]);

  if( dim != 3 )
    {
    antscout << "ReorientTensorImage only supports 3D image volumes" << std::endl;
    return EXIT_FAILURE;
    }

  is_parsing_ok = ReorientTensorImage_ParseInput(argc - 2, argv + 2, moving_image_filename, output_image_filename,
                                                 opt_queue);

  if( is_parsing_ok )
    {
    antscout << "moving_image_filename: " << moving_image_filename << std::endl;
    antscout << "output_image_filename: " << output_image_filename << std::endl;
    DisplayOptQueue(opt_queue);

    ReorientTensorImage<3>(moving_image_filename, output_image_filename, opt_queue);
    }
  else
    {
    antscout << "Input error!" << std::endl;
    return EXIT_FAILURE;
    }
  // ReorientTensorImage<3>(argc,argv);
  // WarpImageForward(argc,argv);
  return EXIT_SUCCESS;
}
} // namespace ants
