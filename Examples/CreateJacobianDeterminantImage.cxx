#include "itkConstNeighborhoodIterator.h"
#include "itkDecomposeTensorFunction.h"
#include "itkImageRegionIteratorWithIndex.h"
// #include "itkVectorFieldGradientImageFunction.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkTimeProbe.h"
#include "itkVariableSizeMatrix.h"
#include "itkVectorImageFileReader.h"
#include "itkVector.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkANTSImageRegistrationOptimizer.h"

template <class TImage, class TDeformationField>
typename TDeformationField::PixelType
TransformVector(TDeformationField* field, typename TImage::IndexType index )
{
  enum { ImageDimension = TImage::ImageDimension };
  typename TDeformationField::PixelType vec = field->GetPixel(index);
  typename TDeformationField::PixelType newvec;
  newvec.Fill(0);
  for( unsigned int row = 0; row < ImageDimension; row++ )
    {
    for( unsigned int col = 0; col < ImageDimension; col++ )
      {
      newvec[row] += vec[col] * field->GetDirection()[row][col];
      }
    }

  return newvec;
}

template <unsigned int ImageDimension>
int CreateJacobianDeterminantImage( int argc, char *argv[] )
{
  typedef float                                  RealType;
  typedef itk::Image<RealType, ImageDimension>   ImageType;
  typedef itk::Vector<RealType, ImageDimension>  VectorType;
  typedef itk::Image<VectorType, ImageDimension> VectorImageType;

  /**
   * Read in vector field
   */
  typedef itk::VectorImageFileReader<ImageType, VectorImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->SetUseAvantsNamingConvention( true );
  reader->Update();
  typename VectorImageType::Pointer vecimg = reader->GetOutput();

  /** smooth before finite differencing */
  typedef itk::ANTSImageRegistrationOptimizer<ImageDimension, float> RegistrationOptimizerType;
  typedef typename RegistrationOptimizerType::Pointer                RegistrationOptimizerPointer;
  RegistrationOptimizerPointer reg = RegistrationOptimizerType::New();
  reg->SmoothDeformationFieldGauss(vecimg, 3);

  typename VectorImageType::SpacingType spacing
    = vecimg->GetSpacing();

  typename ImageType::Pointer jacobian = ImageType::New();
  jacobian->SetOrigin( vecimg->GetOrigin() );
  jacobian->SetSpacing( vecimg->GetSpacing() );
  jacobian->SetRegions( vecimg->GetLargestPossibleRegion() );
  jacobian->SetDirection(  vecimg->GetDirection() );
  jacobian->Allocate();

  itk::TimeProbe timer;
  timer.Start();

  bool calculateLogJacobian = false;
  if( argc > 4 )
    {
    calculateLogJacobian = static_cast<bool>( atoi( argv[4] ) );
    }

  typedef itk::ConstNeighborhoodIterator<VectorImageType>
    ConstNeighborhoodIteratorType;
  typename ConstNeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 2 );

  itk::ZeroFluxNeumannBoundaryCondition<VectorImageType> nbc;
  ConstNeighborhoodIteratorType                          bit;
  itk::ImageRegionIterator<ImageType>                    It;

  // Find the data-set boundary "faces"
  typename itk::NeighborhoodAlgorithm
  ::ImageBoundaryFacesCalculator<VectorImageType>::FaceListType faceList;
  typename itk::NeighborhoodAlgorithm
  ::ImageBoundaryFacesCalculator<VectorImageType> bC;
  faceList = bC( vecimg,
                 vecimg->GetLargestPossibleRegion(), radius );

  typedef itk::VariableSizeMatrix<RealType>                  MatrixType;
  typedef itk::DecomposeTensorFunction<MatrixType, RealType> DecomposerType;
  typename DecomposerType::Pointer decomposer = DecomposerType::New();

  typename itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator
  <VectorImageType>::FaceListType::iterator fit;
  for( fit = faceList.begin(); fit != faceList.end(); ++fit )
    {
    bit = ConstNeighborhoodIteratorType( radius, vecimg, *fit );
    bit.OverrideBoundaryCondition( &nbc );
    bit.GoToBegin();

    It = itk::ImageRegionIterator<ImageType>( jacobian, *fit );
    It.GoToBegin();

    while( !bit.IsAtEnd() )
      {
      MatrixType J;
      J.SetSize( ImageDimension, ImageDimension );
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        for( unsigned int j = 0; j < ImageDimension; j++ )
          {
          RealType x   = bit.GetCenterPixel()[j];
          RealType xp1 = bit.GetNext( i )[j];
          RealType xp2 = bit.GetNext( i, 2 )[j];
          RealType xm1 = bit.GetPrevious( i )[j];
          RealType xm2 = bit.GetPrevious( i, 2 )[j];

          RealType h = 0.5;
          xp1 = xp1 * h + x * (1.0 - h);
          xm1 = xm1 * h + x * (1.0 - h);
          xp2 = xp2 * h + xp1 * (1.0 - h);
          xp2 = xm2 * h + xm1 * (1.0 - h);

          J[i][j] = ( -xp2 + 8.0 * xp1 - 8.0 * xm1 + xm2 ) / ( 12.0 * spacing[i] );
          }
        J[i][i] += 1.0;
        }
      try
        {
        RealType jacDet = decomposer->EvaluateDeterminant( J );
        if( jacDet < 1.e-4 && calculateLogJacobian )
          {
          jacDet = 1.e-4;
          }
        if( vnl_math_isnan(jacDet) )
          {
          jacDet = 1;
          }
        It.Set( ( calculateLogJacobian ? vcl_log( jacDet ) : jacDet ) );
        }
      catch( ... )
        {
        It.Set( itk::NumericTraits<RealType>::max() );
        }
      ++bit;
      ++It;
      }
    }
  timer.Stop();
//  std::cout << "Elapsed time: " << timer.GetMeanTime() << std::endl;

  typedef itk::ImageFileWriter<ImageType> RealImageWriterType;
  typename RealImageWriterType::Pointer realwriter = RealImageWriterType::New();
  realwriter->SetFileName( argv[3] );
  realwriter->SetInput( jacobian );
  realwriter->Update();
  return 0;
}

int main( int argc, char *argv[] )
{
  if( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " ImageDimension deformationField outputImage log-jac?(default-false)"
              << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
    {
    case 2:
      CreateJacobianDeterminantImage<2>( argc, argv );
      break;
    case 3:
      CreateJacobianDeterminantImage<3>( argc, argv );
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
}
