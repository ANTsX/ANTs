/** ANTS Landmarks used to initialize an affine transform ... */

#include "itkLandmarkBasedTransformInitializer.h"
#include "itkImage.h"
#include "itkImageIOBase.h"
#include "itkImageIOFactory.h"
#include <math.h>
#include <iostream>
#include "ReadWriteImage.h"
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
        //	      std::cout <<" f-label " << label << std::endl;
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
        //	      std::cout <<" m-label " << label << std::endl;
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
    std::cout << " fix-label " << fixlabel << " movlabel " << movlabel << std::endl;
    if( movlabel != fixlabel )
      {
      std::cout << " labels do not match -- exiting " << std::endl;
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
  TransformInitializerType::LandmarkPointContainer fixedLandmarks;
  TransformInitializerType::LandmarkPointContainer movingLandmarks;
  TransformInitializerType::LandmarkPointType      point;
  TransformInitializerType::LandmarkPointType      tmp;

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
      if(  label == currentlabel  )
        {
        totalct++;
        // compute center of mass
        ImageType::PointType point;
        fixedimage->TransformIndexToPhysicalPoint(It.GetIndex(), point);
        for( unsigned int i = 0; i < spacing.Size(); i++ )
          {
          myCenterOfMass[i] += point[i];
          }
        // std::cout << " point " << point << std::endl;
        }
      }
    for( unsigned int i = 0; i < spacing.Size(); i++ )
      {
      myCenterOfMass[i] /= (float)totalct;
      }
    // std::cout << " pushing-fix " <<  myCenterOfMass << std::endl;
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
    //    std::cout << " pushing-mov " <<  myCenterOfMass << std::endl;
    movingLandmarks.push_back( myCenterOfMass );
    }

  TransformInitializerType::PointsContainerConstIterator
    fitr = fixedLandmarks.begin();
  TransformInitializerType::PointsContainerConstIterator
    mitr = movingLandmarks.begin();
  while( mitr != movingLandmarks.end() )
    {
    std::cout << "  Fixed Landmark: " << *fitr << " Moving landmark " << *mitr << std::endl;
    ++fitr;
    ++mitr;
    }

  initializer->SetFixedLandmarks(fixedLandmarks);
  initializer->SetMovingLandmarks(movingLandmarks);

  //  initializer->SetFixedImage( fixedimage );
  //  initializer->SetMovingImage( movingimage );
  initializer->SetTransform( transform );
  initializer->InitializeTransform();

  transform->Print(std::cout);
  // to transform a point
  //         transform->TransformPoint( *fitr ) << std::endl;

  return EXIT_SUCCESS;
}

int LandmarkBasedTransformInitializer2D(int, char * [])
{
  std::cout << " not implemented " << std::endl;

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

int main(int argc, char *argv[])
{
  if( argc < 3 )
    {
    std::cout << "Useage ex:   " << argv[0]
              << " FixedImageWithLabeledLandmarks.nii.gz  MovingImageWithLabeledLandmarks.nii.gz  OutAffine.txt "
              << std::endl;
    std::cout
      << " we expect the input images to be (1) N-ary  (2) in the same physical space as the images you want to "
      << std::endl;
    std::cout << " register and (3 ) to have the same landmark points defined within them ... " << std::endl;
    std::cout << " landmarks will be defined from the center of mass of the labels in the input images . " << std::endl;
    std::cout << " You can use ITK-snap to generate the label images. " << std::endl;
    return 1;
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

  return 0;
}
