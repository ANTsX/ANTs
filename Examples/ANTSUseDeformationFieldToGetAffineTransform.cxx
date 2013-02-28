/** ANTS Landmarks used to initialize an affine transform ... */

#include "antsUtilities.h"
#include <algorithm>

#include <stdio.h>
#include <time.h>

#include "itkLandmarkBasedTransformInitializer.h"
#include "itkImage.h"
#include "itkImageIOBase.h"
#include "itkImageIOFactory.h"
#include <math.h>
#include <iostream>
#include "ReadWriteImage.h"
#include "itkTransformFileWriter.h"

#include <vnl/vnl_matrix.h>

#include "vnl/algo/vnl_qr.h"
#include <algorithm>

namespace ants
{
template <class TransformType>
void WriteAffineTransformFile(typename TransformType::Pointer & transform,
                              const std::string & filename)
{
  itk::TransformFileWriter::Pointer transform_writer;

  transform_writer = itk::TransformFileWriter::New();
  transform_writer->SetFileName(filename);
  transform_writer->SetInput(transform);

  try
    {
    transform_writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    antscout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl
             << "Exception in writing tranform file: " << std::endl
             << filename << std::endl;
    return;
    }

  return;
}

// //////////////////////////////////////////////////////////////////////
// Stripped from ANTS_affine_registration2.h
template <class RunningAffineTransformType, class AffineTransformType>
inline void PostConversionInAffine(const typename RunningAffineTransformType::Pointer & transform_running,
                                   typename AffineTransformType::Pointer & transform)
{
  transform->SetCenter(*(reinterpret_cast<typename AffineTransformType::InputPointType *>
                         (const_cast<typename RunningAffineTransformType::InputPointType *>(&(transform_running->
                                                                                              GetCenter() ) ) ) ) );
  transform->SetTranslation(*(reinterpret_cast<typename AffineTransformType::OutputVectorType *>
                              (const_cast<typename RunningAffineTransformType::OutputVectorType *>(&(transform_running
                                                                                                     ->GetTranslation() ) ) ) ) );
  transform->SetMatrix(*(reinterpret_cast<typename AffineTransformType::MatrixType *>
                         (const_cast<typename RunningAffineTransformType::MatrixType *>(&(transform_running->GetMatrix() ) ) ) ) );

  // antscout << "transform_running" << transform_running << std::endl;
  // antscout << "transform" << transform << std::endl;
}

template <class TransformA>
void DumpTransformForANTS3D(const typename TransformA::Pointer & transform, const std::string & ANTS_prefix)
{
  const int ImageDimension = 3;

  // ANTS transform file type
  typedef itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension> AffineTransformType;
  AffineTransformType::Pointer transform_ANTS = AffineTransformType::New();

  //    typedef TransformAPointer::ObjectType TransformA;

  // antscout << " writing " << ANTS_prefix << " affine " << std::endl;
  // std::string ANTS_affine_filename = ANTS_prefix + std::string( "Affine.txt" );

  std::string ANTS_affine_filename = ANTS_prefix;

  antscout << " writing ANTS affine file:" << ANTS_affine_filename << std::endl;

  PostConversionInAffine<TransformA, AffineTransformType>(transform, transform_ANTS);

  WriteAffineTransformFile<AffineTransformType>(transform_ANTS,
                                                ANTS_affine_filename);
}

// ////////
// x: fixedLandmarks
// y: movingLandmarks
// (A,t,c) : affine transform, A:3*3, t: 3*1 c: 3*1 (c is the center of all points in x)
// y-c = A*(x-c) + t;
// steps:
// 1. c = average of points of x
// 2. let y1 = y-c; x1 = x - c; x11 = [x1; 1 ... 1] // extend x11
// 3. minimize (y1-A1*x11)^2, A1 is a 3*4 matrix
// 4. A = A1(1:3, 1:3), t = A1(1:3, 4);
// step 3:
//   A11 = (y1*x11')*(x11*x11')^(-1)
// type info:
//   assume PointContainerType is std::vector
//   assume TrnasformPointerType is MatrixOffsetTransformBase

template <class PointContainerType, class TransformType>
void GetAffineTransformFromTwoPointSets3D(PointContainerType & fixedLandmarks, PointContainerType & movingLandmarks,
                                          typename TransformType::Pointer & transform)
{
  const int Dim = 3;
  int       n = fixedLandmarks.size();

  vnl_matrix<double> y(Dim, n), x(Dim, n);
  for( int i = 0; i < n; i++ )
    {
    for( int j = 0; j < Dim; j++ )
      {
      x(j, i) = fixedLandmarks[i][j];
      y(j, i) = movingLandmarks[i][j];
      }
    }

  vnl_vector<double> c(Dim);
  for( int j = 0; j < Dim; j++ )
    {
    c[j] = x.get_row(j).mean();
    }

  vnl_matrix<double> y1(Dim, n), x11(Dim + 1, n);
  for( int i = 0; i < n; i++ )
    {
    y1.set_column(i, y.get_column(i) - c);

    vnl_vector<double> x_tmp(Dim), x1_tmp(Dim + 1);
    x_tmp = x.get_column(i) - c;
    for( int j = 0; j < Dim; j++ )
      {
      x1_tmp[j] = x_tmp[j];
      }
    x1_tmp[Dim] = 1;

    x11.set_column(i, x1_tmp);
    }

  vnl_matrix<double> A11(Dim, Dim + 1);
  vnl_matrix<double> x11t = x11.transpose();
  // vnl_matrix_inverse<double> tmp(x11 * x11t); // BA -- removed this -- not used?

  vnl_svd<double> qr( x11t );   // can use vnl_qr
  A11 = qr.inverse() * (y1.transpose() );
  A11 = A11.transpose();

  vnl_matrix<double> A(Dim, Dim);
  A = A11.extract(Dim, Dim, 0, 0);

  //    antscout << "y=" << y << std::endl;
  //    antscout << "x=" << x << std::endl;
  //
  //    antscout << "y1=" << y1 << std::endl;
  //    antscout << "x11=" << x11 << std::endl;
  antscout << "A11=" << A11 << std::endl;

  vnl_vector<double> t = A11.get_column(Dim);

  typedef typename TransformType::InputPointType   PointType;
  typedef typename TransformType::OutputVectorType VectorType;
  typedef typename TransformType::MatrixType       MatrixType;

  PointType center;
  for( int i = 0; i < Dim; i++ )
    {
    center[i] = c[i];
    }

  VectorType translation;
  for( int i = 0; i < Dim; i++ )
    {
    translation[i] = t[i];
    }

  MatrixType matrix(A);

  transform->SetCenter(center);
  transform->SetTranslation(translation);
  transform->SetMatrix(matrix);

  return;
}

template <class PointContainerType, class TTransform>
void GetRigidTransformFromTwoPointSets3D(PointContainerType & fixedLandmarks, PointContainerType & movingLandmarks,
                                         typename TTransform::Pointer & aff)
{
  // Set the transform type..
  typedef itk::VersorRigid3DTransform<double> TransformType;
  TransformType::Pointer transform = TransformType::New();

  typedef  float PixelType;
  const unsigned int Dimension = 3;
  typedef itk::Image<PixelType, Dimension> FixedImageType;
  typedef itk::Image<PixelType, Dimension> MovingImageType;
  typedef itk::Image<PixelType, Dimension> ImageType;

  typedef itk::LandmarkBasedTransformInitializer<TransformType,
                                                 FixedImageType, MovingImageType> TransformInitializerType;
  TransformInitializerType::Pointer initializer = TransformInitializerType::New();

  initializer->SetFixedLandmarks(fixedLandmarks);
  initializer->SetMovingLandmarks(movingLandmarks);
  initializer->SetTransform( transform );
  initializer->InitializeTransform();

  antscout << "rigid: " << transform << std::endl;

  // ANTS transform file type
  // typedef itk::MatrixOffsetTransformBase< double, Dimension, Dimension > AffineTransformType;
  // typename AffineTransformType::Pointer aff = AffineTransformType::New();

  PostConversionInAffine<TransformType, TTransform>(transform, aff);
}

template <class PointContainerType>
void FetchLandmarkMappingFromDisplacementField(const std::string& deformation_field_file_name,
                                               float load_ratio,
                                               PointContainerType & fixedLandmarks,
                                               PointContainerType & movingLandmarks,
                                               typename itk::Image<float, 3>::Pointer maskimg)
{
  const unsigned int ImageDimension = 3;

  typedef typename PointContainerType::value_type PointType;

  typedef itk::Image<float, ImageDimension>           ImageType;
  typedef itk::Vector<float, ImageDimension>          VectorType;
  typedef itk::Image<VectorType, ImageDimension>      DisplacementFieldType;
  typedef itk::ImageFileReader<DisplacementFieldType> FieldReaderType;

  typename FieldReaderType::Pointer field_reader = FieldReaderType::New();
  field_reader->SetFileName( deformation_field_file_name );
  field_reader->Update();
  typename DisplacementFieldType::Pointer field = field_reader->GetOutput();

  fixedLandmarks.clear();
  movingLandmarks.clear();

  unsigned int              nb_voxels = 1;
  itk::Size<ImageDimension> field_size = field->GetLargestPossibleRegion().GetSize();
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    nb_voxels *= field_size[i];
    }

  // float load_ratio = 0.01;
  unsigned int nb_try_to_load = (unsigned int) ( (float) nb_voxels * load_ratio);

  antscout << "trying to load " << nb_try_to_load << " from " <<  nb_voxels << " points." << std::endl;

  fixedLandmarks.reserve(nb_try_to_load );
  movingLandmarks.reserve(nb_try_to_load );

  typedef itk::ImageRegionIteratorWithIndex<DisplacementFieldType> FieldIteratorType;

  FieldIteratorType it(field, field->GetLargestPossibleRegion() );

  srand( time(NULL) );

  it.GoToBegin();
  unsigned int cnt = 0;
  for( ; (!it.IsAtEnd() ) & (cnt < nb_try_to_load); ++it )
    {
    bool getpoint = true;
    if( maskimg )
      {
      if( maskimg->GetPixel( it.GetIndex() ) < 0.5 )
        {
        getpoint = false;
        }
      }

    if( getpoint )
      {
      if( rand() % 32767 > load_ratio * 32767 )
        {
        continue;
        }

      PointType point1, point2;
      // get the output image index
      typename DisplacementFieldType::IndexType index = it.GetIndex();
      field->TransformIndexToPhysicalPoint(index, point1 );
      VectorType displacement = field->GetPixel(index);
      for( unsigned int j = 0; j < ImageDimension; j++ )
        {
        point2[j] = point1[j] + displacement[j];
        }

      fixedLandmarks.push_back(point1);
      movingLandmarks.push_back(point2);

      ++cnt;
      }
    }

  antscout << "total " << cnt << " points loaded from " << deformation_field_file_name << "." << std::endl;
  antscout << fixedLandmarks.size() << std::endl;
  antscout << movingLandmarks.size() << std::endl;

  return;
}

//
// The test specifies a bunch of fixed and moving landmarks and test if the
// fixed landmarks after transform by the computed transform coincides
// with the moving landmarks....

int DisplacementFieldBasedTransformInitializer3D(int argc, char * argv[])
{
  const unsigned int Dim = 3;

  typedef itk::Point<double, Dim> PointType;
  typedef itk::Image<float, Dim>  ImageType;
  typedef std::vector<PointType>  PointContainerType;
  const char *deformation_field_file_name = argv[1];
  float       load_ratio = atof(argv[2]);
  bool        bRigid = (strcmp(argv[3], "rigid") == 0);
  std::string ANTS_prefix(argv[4]);
  std::string maskfn = std::string("");
  if( argc > 5 )
    {
    maskfn = std::string(argv[5]);
    }
  antscout << " mask " << maskfn << std::endl;

  // input
  PointContainerType fixedLandmarks, movingLandmarks;
  // output
  typedef itk::MatrixOffsetTransformBase<double, 3, 3> AffineTransformType;
  AffineTransformType::Pointer aff = AffineTransformType::New();

  ImageType::Pointer maskimg = NULL;
  if( maskfn.length() > 4 )
    {
    ReadImage<ImageType>(maskimg, maskfn.c_str() );
    }

  FetchLandmarkMappingFromDisplacementField(deformation_field_file_name, load_ratio, fixedLandmarks, movingLandmarks,
                                            maskimg);

  if( bRigid )
    {
    GetRigidTransformFromTwoPointSets3D<PointContainerType, AffineTransformType>(fixedLandmarks, movingLandmarks, aff);
    }
  else
    {
    GetAffineTransformFromTwoPointSets3D<PointContainerType, AffineTransformType>(fixedLandmarks, movingLandmarks, aff);
    }

  antscout << "affine:" << aff;
  DumpTransformForANTS3D<AffineTransformType>(aff, ANTS_prefix);

  return EXIT_SUCCESS;

  //    initializer->SetFixedLandmarks(fixedLandmarks);
  //    initializer->SetMovingLandmarks(movingLandmarks);
  //    initializer->SetTransform( transform );
  //    initializer->InitializeTransform();
  //
  //    transform->Print(antscout);
  //
  //    // transform the transform to ANTS format
  //    std::string ANTS_prefix(argv[4]);
  //
  //
  //    typedef itk::MatrixOffsetTransformBase< double, 3, 3> AffineTransformType;
  //    AffineTransformType::Pointer aff = AffineTransformType::New();
  //    GetAffineTransformFromTwoPointSets3D(fixedLandmarks, movingLandmarks, aff);
  //    antscout << "affine:" << aff;
  //
  //
  //    if (bRigid)
  //        DumpTransformForANTS3D(transform, ANTS_prefix);
  //    else
  //        DumpTransformForANTS3D(aff, ANTS_prefix);
  //
  //
  //    return EXIT_SUCCESS;
}

// //////////////////////////////////////////////////////////////////////
// Stripped from ANTS_affine_registration2.h

int DisplacementFieldBasedTransformInitializer2D(int, char * [])
{
  antscout << " not implemented " << std::endl;
  return 1;

  /*
typedef  float PixelType;
const unsigned int Dimension = 2;
typedef itk::Image< PixelType, Dimension >  FixedImageType;
typedef itk::Image< PixelType, Dimension >  MovingImageType;
typedef itk::Image< PixelType, Dimension >  ImageType;
typename FixedImageType::Pointer fixedimage;
typename MovingImageType::Pointer movingimage;
ReadImage<ImageType>(fixedimage,argv[1]);
ReadImage<ImageType>(movingimage,argv[2]);

// Set the transform type..
typedef itk::Rigid2DTransform< double > TransformType;
   */

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int ANTSUseDeformationFieldToGetAffineTransform( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "ANTSUseDeformationFieldToGetAffineTransform" );
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

  if( argc < 3 )
    {
    antscout << "Usage:   " << argv[0]
             << " zzzWarp.nii.gz load_ratio(ex: 0.01) [rigid | affine] OutAffine.txt [mask.nii.gz]" << std::endl;
    antscout << " we expect the input deformation field in the same physical space as the images you want to "
             << std::endl;
    antscout << "load_ratio: ratio of points to be loaded from deformation field (to save memory) " << std::endl;
    antscout << " the mask gives the region from which points will be selected ... " << std::endl;
    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }

  // Get the image dimension
  // std::string fn = std::string(argv[1]) + "xvec.nii.gz";
  // itk::ImageIOBase::Pointer imageIO =
  //        itk::ImageIOFactory::CreateImageIO(fn.c_str(), itk::ImageIOFactory::ReadMode);
  // imageIO->SetFileName(fn.c_str());
  // imageIO->ReadImageInformation();

  int dim = 3;

  // switch ( imageIO->GetNumberOfDimensions() )
  switch( dim )
    {
    case 2:
      {
      DisplacementFieldBasedTransformInitializer2D(argc, argv);
      }
      break;
    case 3:
      {
      DisplacementFieldBasedTransformInitializer3D(argc, argv);
      }
      break;
    default:
      antscout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
