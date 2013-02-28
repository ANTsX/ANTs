/** ANTS Landmarks used to initialize an affine transform ... */

#include "antsUtilities.h"
#include <algorithm>

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
namespace ants
{
// //////////////////////////////////////////////////////////////////////
// Stripped from ANTS_affine_registration2.h
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
template <class RunningAffineTransformPointerType, class AffineTransformPointerType>
inline void PostConversionInAffine(RunningAffineTransformPointerType& transform_running,
                                   AffineTransformPointerType & transform)
{
  typedef typename RunningAffineTransformPointerType::ObjectType RunningAffineTransformType;
  typedef typename AffineTransformPointerType::ObjectType        AffineTransformType;

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
void DumpTransformForANTS3D(typename TransformA::Pointer & transform, const std::string & ANTS_prefix)
{
  const int ImageDimension = 3;

  // ANTS transform file type
  typedef itk::AffineTransform<double, ImageDimension> AffineTransformType;
  AffineTransformType::Pointer transform_ANTS = AffineTransformType::New();

  //    typedef TransformAPointer::ObjectType TransformA;

  // antscout << " writing " << ANTS_prefix << " affine " << std::endl;
  // std::string ANTS_affine_filename = ANTS_prefix + std::string( "Affine.txt" );

  std::string ANTS_affine_filename = ANTS_prefix;

  antscout << " writing ANTS affine file:" << ANTS_affine_filename << std::endl;
  PostConversionInAffine(transform, transform_ANTS);
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

  antscout << "y=" << y << std::endl;
  antscout << "x=" << x << std::endl;

  antscout << "y1=" << y1 << std::endl;
  antscout << "x11=" << x11 << std::endl;
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

//
// The test specifies a bunch of fixed and moving landmarks and test if the
// fixed landmarks after transform by the computed transform coincides
// with the moving landmarks....

int LandmarkBasedTransformInitializer3D(int, char * argv[])
{
  typedef  float PixelType;
  const unsigned int Dimension = 3;
  typedef itk::Image<PixelType, Dimension>             FixedImageType;
  typedef itk::Image<PixelType, Dimension>             MovingImageType;
  typedef itk::Image<PixelType, Dimension>             ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  ImageType::Pointer fixedimage;
  ImageType::Pointer movingimage;
  ReadImage<ImageType>(fixedimage, argv[1]);
  ReadImage<ImageType>(movingimage, argv[2]);

  bool bRigid = (strcmp(argv[3], "rigid") == 0);

  /** get all of the relevant labels in the fixed image and moving image */
  typedef std::vector<PixelType> LabelSetType;
  LabelSetType myFixLabelSet;
  LabelSetType myMovLabelSet;
  /** count the labels in the image */
  Iterator It( fixedimage, fixedimage->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    PixelType label = It.Get();
    if( fabs(label) > 0 )
      {
      if( find( myFixLabelSet.begin(), myFixLabelSet.end(), label )
          == myFixLabelSet.end() )
        {
        //          antscout <<" f-label " << label << std::endl;
        myFixLabelSet.push_back( label );
        }
      }
    }
  Iterator ItM( movingimage, movingimage->GetLargestPossibleRegion() );
  for( ItM.GoToBegin(); !ItM.IsAtEnd(); ++ItM )
    {
    PixelType label = ItM.Get();
    if( fabs(label) > 0 )
      {
      if( find( myMovLabelSet.begin(), myMovLabelSet.end(), label )
          == myMovLabelSet.end() )
        {
        //          antscout <<" m-label " << label << std::endl;
        myMovLabelSet.push_back( label );
        }
      }
    }

  std::sort(myFixLabelSet.begin(), myFixLabelSet.end() );
  std::sort(myMovLabelSet.begin(), myMovLabelSet.end() );

  LabelSetType::const_iterator fit;
  LabelSetType::const_iterator mit = myMovLabelSet.begin();
  for( fit = myFixLabelSet.begin(); fit != myFixLabelSet.end(); ++fit )
    {
    float fixlabel = *fit;
    float movlabel = *mit;
    antscout << " fix-label " << fixlabel << " movlabel " << movlabel << std::endl;
    if( movlabel != fixlabel )
      {
      antscout << " labels do not match -- exiting " << std::endl;
      exit(1);
      }
    ++mit;
    }

  // Set the transform type..
  typedef itk::VersorRigid3DTransform<double> TransformType;
  TransformType::Pointer transform = TransformType::New();
  typedef itk::LandmarkBasedTransformInitializer<TransformType,
                                                 FixedImageType, MovingImageType> TransformInitializerType;
  TransformInitializerType::Pointer initializer = TransformInitializerType::New();

  // Set fixed and moving landmarks
  typedef TransformInitializerType::LandmarkPointContainer PointsContainerType;
  PointsContainerType fixedLandmarks;
  PointsContainerType movingLandmarks;

  // compute the CoM's of all the landmarks
  ImageType::SpacingType spacing = fixedimage->GetSpacing();
  for( fit = myFixLabelSet.begin(); fit != myFixLabelSet.end(); ++fit )
    {
    float                                       currentlabel = *fit;
    float                                       totalct = 0;
    TransformInitializerType::LandmarkPointType myCenterOfMass;
    myCenterOfMass.Fill(0);
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      PixelType label = It.Get();
      if( fabs( label - currentlabel ) < 0.001  )
        {
        totalct++;
        // compute center of mass
        ImageType::PointType point;
        fixedimage->TransformIndexToPhysicalPoint(It.GetIndex(), point);
        for( unsigned int i = 0; i < spacing.Size(); i++ )
          {
          myCenterOfMass[i] += point[i];
          }
        antscout << " point " << point << std::endl;
        }
      }
    for( unsigned int i = 0; i < spacing.Size(); i++ )
      {
      myCenterOfMass[i] /= (float)totalct;
      }
    // antscout << " pushing-fix " <<  myCenterOfMass << std::endl;
    fixedLandmarks.push_back( myCenterOfMass );
    }

  // compute the CoM's of all the landmarks
  spacing = movingimage->GetSpacing();
  for( mit = myMovLabelSet.begin(); mit != myMovLabelSet.end(); ++mit )
    {
    float                                       currentlabel = *mit;
    float                                       totalct = 0;
    TransformInitializerType::LandmarkPointType myCenterOfMass;
    myCenterOfMass.Fill(0);
    for( ItM.GoToBegin(); !ItM.IsAtEnd(); ++ItM )
      {
      PixelType label = ItM.Get();
      if(  label == currentlabel  )
        {
        totalct++;
        // compute center of mass
        ImageType::PointType point;
        movingimage->TransformIndexToPhysicalPoint(ItM.GetIndex(), point);
        for( unsigned int i = 0; i < spacing.Size(); i++ )
          {
          myCenterOfMass[i] += point[i];
          }
        }
      }
    for( unsigned int i = 0; i < spacing.Size(); i++ )
      {
      myCenterOfMass[i] /= (float)totalct;
      }
    //    antscout << " pushing-mov " <<  myCenterOfMass << std::endl;
    movingLandmarks.push_back( myCenterOfMass );
    }

  TransformInitializerType::PointsContainerConstIterator
    fitr = fixedLandmarks.begin();
  TransformInitializerType::PointsContainerConstIterator
    mitr = movingLandmarks.begin();
  while( mitr != movingLandmarks.end() )
    {
    antscout << "  Fixed Landmark: " << *fitr << " Moving landmark " << *mitr << std::endl;
    ++fitr;
    ++mitr;
    }

  initializer->SetFixedLandmarks(fixedLandmarks);
  initializer->SetMovingLandmarks(movingLandmarks);

  //  initializer->SetFixedImage( fixedimage );
  //  initializer->SetMovingImage( movingimage );
  initializer->SetTransform( transform );
  initializer->InitializeTransform();

  transform->Print(antscout);
  // to transform a point
  //         transform->TransformPoint( *fitr ) << std::endl;

  // transform the transform to ANTS format
  std::string ANTS_prefix(argv[4]);

  typedef itk::AffineTransform<double, 3> AffineTransformType;
  AffineTransformType::Pointer aff = AffineTransformType::New();

  GetAffineTransformFromTwoPointSets3D<PointsContainerType, AffineTransformType>(fixedLandmarks, movingLandmarks, aff);

  antscout << "affine:" << aff;

  if( bRigid )
    {
    DumpTransformForANTS3D<TransformType>(transform, ANTS_prefix);
    }
  else
    {
    DumpTransformForANTS3D<AffineTransformType>(aff, ANTS_prefix);
    }

  return EXIT_SUCCESS;
}

int LandmarkBasedTransformInitializer2D(int, char * [])
{
  antscout << " not implemented " << std::endl;
  return EXIT_FAILURE;

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

  return EXIT_SUCCESS;
   */
}

int ANTSUseLandmarkImagesToGetAffineTransform( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "ANTSUseLandmarkImagesToGetAffineTransform" );
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
             <<
      " FixedImageWithLabeledLandmarks.nii.gz  MovingImageWithLabeledLandmarks.nii.gz [rigid | affine] OutAffine.txt "
             << std::endl;
    antscout
      << " we expect the input images to be (1) N-ary  (2) in the same physical space as the images you want to "
      << std::endl;
    antscout << " register and (3 ) to have the same landmark points defined within them ... " << std::endl;
    antscout << " landmarks will be defined from the center of mass of the labels in the input images . " << std::endl;
    antscout << " You can use ITK-snap to generate the label images. " << std::endl;
    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }

  // Get the image dimension
  std::string               fn = std::string(argv[1]);
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(fn.c_str(), itk::ImageIOFactory::ReadMode);
  imageIO->SetFileName(fn.c_str() );
  imageIO->ReadImageInformation();

  switch( imageIO->GetNumberOfDimensions() )
    {
    case 2:
      {
      LandmarkBasedTransformInitializer2D(argc, argv);
      }
      break;
    case 3:
      {
      LandmarkBasedTransformInitializer3D(argc, argv);
      }
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
  return EXIT_SUCCESS;
}
} // namespace ants
