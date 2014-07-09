/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RSfile: antsAffineInitializer.cxx,v $
  Language:  C++
  Date:      $Date: 2009/06/02 21:51:08 $
  Version:   $Revision: 1.103 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "antsUtilities.h"
#include <algorithm>
#include <vnl/vnl_inverse.h>
#include "antsAllocImage.h"
#include "itkImageMaskSpatialObject.h"
#include "itkANTSNeighborhoodCorrelationImageToImageMetricv4.h"
#include "itkArray.h"
#include "itkGradientImageFilter.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkBayesianClassifierImageFilter.h"
#include "itkBayesianClassifierInitializationImageFilter.h"
#include "itkBilateralImageFilter.h"
#include "itkCSVNumericObjectFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkCompositeValleyFunction.h"
#include "itkConjugateGradientLineSearchOptimizerv4.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkCorrelationImageToImageMetricv4.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkDistanceToCentroidMembershipFunction.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkDemonsImageToImageMetricv4.h"
#include "itkExpImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkGaussianImageSource.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkHistogram.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkImage.h"
#include "itkImageClassifierBase.h"
#include "itkImageDuplicator.h"
#include "itkImageFileWriter.h"
#include "itkImageGaussianModelEstimator.h"
#include "itkImageKmeansModelEstimator.h"
#include "itkImageMomentsCalculator.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkKdTree.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkLabelContourImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkLabeledPointSetFileReader.h"
#include "itkLabeledPointSetFileWriter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkListSample.h"
#include "itkMRFImageFilter.h"
#include "itkMRIBiasFieldCorrectionFilter.h"
#include "itkMaskImageFilter.h"
#include "itkMaximumImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkMultivariateLegendrePolynomial.h"
#include "itkMultiStartOptimizerv4.h"
#include "itkNeighborhood.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkNeighborhoodIterator.h"
#include "itkNormalVariateGenerator.h"
#include "itkOptimizerParameterScalesEstimator.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkRGBPixel.h"
#include "itkRegistrationParameterScalesFromPhysicalShift.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSampleToHistogramFilter.h"
#include "itkScalarImageKmeansImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkSize.h"
#include "itkSphereSpatialFunction.h"
#include "itkSTAPLEImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkTDistribution.h"
#include "itkTimeProbe.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkTranslationTransform.h"
#include "itkVariableSizeMatrix.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "vnl/vnl_matrix_fixed.h"
#include "itkTransformFactory.h"
#include "itkSurfaceImageCurvature.h"
#include "itkMultiScaleLaplacianBlobDetectorImageFilter.h"
#include "itkEuler2DTransform.h"
#include "itkEuler3DTransform.h"
#include "itkCenteredAffineTransform.h"
#include "itkCompositeTransform.h"

#include <fstream>
#include <iostream>
#include <map> // Here I'm using a map but you could choose even other containers
#include <sstream>
#include <string>

#include "ReadWriteData.h"
#include "TensorFunctions.h"
#include "antsMatrixUtilities.h"

namespace ants
{
template <unsigned int ImageDimension>
int antsAffineInitializerImp(int argc, char *argv[])
{
  typedef double RealType;
  typedef float  PixelType;

  /** Define All Parameters Here */
  double       pi = vnl_math::pi;                // probably a vnl alternative
  RealType     searchfactor = 10;                // in degrees, passed by user
  unsigned int mibins = 32;                      // for mattes MI metric
  RealType     degtorad = 0.0174532925;          // to convert degrees to radians
  RealType     localoptimizerlearningrate = 0.1; // for local search via conjgrad
  unsigned int localoptimizeriterations = 20;    // for local search via conjgrad
  // piover4 is (+/-) for cross-section of the sphere to multi-start search in increments
  // of searchfactor ( converted from degrees to radians ).
  // the search is centered +/- from the principal axis alignment of the images.
  RealType piover4 = pi / 4; // works in preliminary practical examples in 3D, in 2D use pi.
  bool     useprincaxis = false;

  typedef itk::TransformFileWriter                                        TransformWriterType;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;
  typedef typename itk::ImageMomentsCalculator<ImageType>                 ImageCalculatorType;
  typedef itk::AffineTransform<RealType, ImageDimension>                  AffineType;
  typedef itk::CompositeTransform<RealType, ImageDimension>               CompositeType;
  typedef itk::AffineTransform<RealType>                                  EulerTransformType;
  typedef typename ImageCalculatorType::MatrixType                        MatrixType;
  if( argc < 2 )
    {
    return 0;
    }
  int         argct = 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = std::string(argv[argct]);   argct++;
  std::string outname = std::string(argv[argct]); argct++;
  if(  argc > argct )
    {
    searchfactor = atof( argv[argct] );   argct++;
    }
  if(  argc > argct )
    {
    RealType temp = atof( argv[argct] );   argct++;
    if( temp > 1 )
      {
      temp = 1;
      }
    if( temp < 0.01 )
      {
      temp = 0.01;
      }
    piover4 = pi * temp;
    }
  if(  argc > argct )
    {
    useprincaxis = atoi( argv[argct] );   argct++;
    }
  if(  argc > argct )
    {
    localoptimizeriterations = atoi( argv[argct] );   argct++;
    }
  std::string maskfn = "";
  if(  argc > argct )
    {
    maskfn = std::string( argv[argct] );   argct++;
    }
  searchfactor *= degtorad; // convert degrees to radians
  typename ImageType::Pointer image1 = NULL;
  typename ImageType::Pointer image2 = NULL;
  ReadImage<ImageType>(image1, fn1.c_str() );
  ReadImage<ImageType>(image2, fn2.c_str() );
  VectorType ccg1;
  VectorType cpm1;
  MatrixType cpa1;
  VectorType ccg2;
  VectorType cpm2;
  MatrixType cpa2;

  typename ImageCalculatorType::Pointer calculator1 = ImageCalculatorType::New();
  typename ImageCalculatorType::Pointer calculator2 = ImageCalculatorType::New();
  calculator1->SetImage(  image1 );
  calculator2->SetImage(  image2 );
  typename ImageCalculatorType::VectorType fixed_center;
  fixed_center.Fill(0);
  typename ImageCalculatorType::VectorType moving_center;
  moving_center.Fill(0);
  try
    {
    calculator1->Compute();
    fixed_center = calculator1->GetCenterOfGravity();
    ccg1 = calculator1->GetCenterOfGravity();
    cpm1 = calculator1->GetPrincipalMoments();
    cpa1 = calculator1->GetPrincipalAxes();
    try
      {
      calculator2->Compute();
      moving_center = calculator2->GetCenterOfGravity();
      ccg2 = calculator2->GetCenterOfGravity();
      cpm2 = calculator2->GetPrincipalMoments();
      cpa2 = calculator2->GetPrincipalAxes();
      }
    catch( ... )
      {
      std::cerr << " zero image2 error ";
      fixed_center.Fill(0);
      }
    }
  catch( ... )
    {
    std::cerr << " zero image1 error ";
    }
  unsigned int eigind1 = 1;
  unsigned int eigind2 = 1;
  if( ImageDimension == 3 )
    {
    eigind1 = 2;
    }
  vnl_vector<RealType> evec1_2ndary =  cpa1.GetVnlMatrix().get_row( eigind2 );
  vnl_vector<RealType> evec1_primary = cpa1.GetVnlMatrix().get_row( eigind1 );
  vnl_vector<RealType> evec2_2ndary  = cpa2.GetVnlMatrix().get_row( eigind2 );
  vnl_vector<RealType> evec2_primary = cpa2.GetVnlMatrix().get_row( eigind1 );
  /** Solve Wahba's problem --- http://en.wikipedia.org/wiki/Wahba%27s_problem */
  vnl_matrix<RealType> B = outer_product( evec2_primary, evec1_primary );
  if( ImageDimension == 3 )
    {
    B = outer_product( evec2_2ndary, evec1_2ndary )
      + outer_product( evec2_primary, evec1_primary );
    }
  vnl_svd<RealType>    wahba( B );
  vnl_matrix<RealType> A_solution = wahba.V() * wahba.U().transpose();
  A_solution = vnl_inverse( A_solution );
  RealType det = vnl_determinant( A_solution  );
  if( det < 0 )
    {
    std::cerr << " bad det " << det << " v " <<  vnl_determinant( wahba.V() ) << " u "
             <<   vnl_determinant( wahba.U() )  << std::endl;
    vnl_matrix<RealType> id( A_solution );
    id.set_identity();
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      if( A_solution( i, i ) < 0 )
        {
        id( i, i ) = -1.0;
        }
      }
    A_solution =  A_solution * id.transpose();
    std::cerr << " bad det " << det << " v " <<  vnl_determinant( wahba.V() ) << " u "
             <<   vnl_determinant( wahba.U() )  << " new " << vnl_determinant( A_solution  ) << std::endl;
    }
  typename AffineType::Pointer affine1 = AffineType::New(); // translation to center
  typename AffineType::OffsetType trans = affine1->GetOffset();
  itk::Point<double, ImageDimension> trans2;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    trans[i] = moving_center[i] - fixed_center[i];
    trans2[i] =  fixed_center[i] * ( 1 );
    }
  affine1->SetIdentity();
  affine1->SetOffset( trans );
  if( useprincaxis )
    {
    affine1->SetMatrix( A_solution );
    }
  affine1->SetCenter( trans2 );
    {
    typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
    transformWriter->SetInput( affine1 );
    transformWriter->SetFileName( outname.c_str() );
    transformWriter->Update();
    }
  if( ImageDimension > 3  )
    {
    return EXIT_SUCCESS;
    }
  vnl_vector<RealType> evec_tert;
  if( ImageDimension == 3 )
    { // try to rotate around tertiary and secondary axis
    evec_tert = vnl_cross_3d( evec1_primary, evec1_2ndary );
    }
  if( ImageDimension == 2 )
    { // try to rotate around tertiary and secondary axis
    evec_tert = evec1_2ndary;
    evec1_2ndary = evec1_primary;
    }
  itk::Vector<RealType, ImageDimension> axis2;
  itk::Vector<RealType, ImageDimension> axis1;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    axis1[d] = evec_tert[d];
    axis2[d] = evec1_2ndary[d];
    }
  typename AffineType::Pointer affinesearch = AffineType::New();
  affinesearch->SetIdentity();
  affinesearch->SetCenter( trans2 );
  typedef  itk::MultiStartOptimizerv4         OptimizerType;
  typedef  typename OptimizerType::ScalesType ScalesType;
  typename OptimizerType::Pointer  mstartOptimizer = OptimizerType::New();
  typedef itk::MattesMutualInformationImageToImageMetricv4
    <ImageType, ImageType, ImageType> MetricType;
  typename MetricType::ParametersType newparams(  affine1->GetParameters() );
  typename MetricType::Pointer mimetric = MetricType::New();
  mimetric->SetNumberOfHistogramBins( mibins );
  mimetric->SetFixedImage( image1 );
  mimetric->SetMovingImage( image2 );
  mimetric->SetMovingTransform( affinesearch );
  mimetric->SetParameters( newparams );
  if( maskfn.length() > 3 )
    {
    typedef typename itk::ImageMaskSpatialObject<ImageDimension>::ImageType maskimagetype;
    typename maskimagetype::Pointer mask;
    ReadImage<maskimagetype>( mask, maskfn.c_str() );
    typename itk::ImageMaskSpatialObject<ImageDimension>::Pointer so =
      itk::ImageMaskSpatialObject<ImageDimension>::New();
    so->SetImage( const_cast<maskimagetype *>(mask.GetPointer() ) );
    mimetric->SetFixedImageMask( so );
    }
  mimetric->Initialize();
  typedef itk::RegistrationParameterScalesFromPhysicalShift<MetricType>
    RegistrationParameterScalesFromPhysicalShiftType;
  typename RegistrationParameterScalesFromPhysicalShiftType::Pointer shiftScaleEstimator =
    RegistrationParameterScalesFromPhysicalShiftType::New();
  shiftScaleEstimator->SetMetric( mimetric );
  shiftScaleEstimator->SetTransformForward( true ); // by default, scales for the moving transform
  typename RegistrationParameterScalesFromPhysicalShiftType::ScalesType
    movingScales( affinesearch->GetNumberOfParameters() );
  shiftScaleEstimator->EstimateScales( movingScales );
  mstartOptimizer->SetScales( movingScales );
  std::cout << " Scales: " << movingScales << std::endl;
  mstartOptimizer->SetMetric( mimetric );
  typename OptimizerType::ParametersListType parametersList = mstartOptimizer->GetParametersList();
  affinesearch->SetIdentity();
  affinesearch->SetCenter( trans2 );
  affinesearch->SetOffset( trans );
  parametersList.push_back( affinesearch->GetParameters() );
  for( double ang1 = ( piover4 * (-1) ); ang1 <= ( piover4 + searchfactor ); ang1 = ang1 + searchfactor )
    {
    if( ImageDimension == 3 )
      {
      for( double ang2 = ( piover4 * (-1) ); ang2 <= ( piover4 + searchfactor ); ang2 = ang2 + searchfactor )
        {
        affinesearch->SetIdentity();
        affinesearch->SetCenter( trans2 );
        affinesearch->SetOffset( trans );
        if( useprincaxis )
          {
          affinesearch->SetMatrix( A_solution );
          }
        affinesearch->Rotate3D(axis1, ang1, 1);
        affinesearch->Rotate3D(axis2, ang2, 1);
        parametersList.push_back( affinesearch->GetParameters() );
        }
      }
    if( ImageDimension == 2 )
      {
      affinesearch->SetIdentity();
      affinesearch->SetCenter( trans2 );
      affinesearch->SetOffset( trans );
      if( useprincaxis )
        {
        affinesearch->SetMatrix( A_solution );
        }
      affinesearch->Rotate2D( ang1, 1);
      parametersList.push_back( affinesearch->GetParameters() );
      }
    }
  mstartOptimizer->SetParametersList( parametersList );

  typedef  itk::ConjugateGradientLineSearchOptimizerv4 LocalOptimizerType;
  typename LocalOptimizerType::Pointer  localoptimizer = LocalOptimizerType::New();
  localoptimizer->SetMetric( mimetric );
  localoptimizer->SetScales( movingScales );
  localoptimizer->SetLearningRate( localoptimizerlearningrate );
  localoptimizer->SetMaximumStepSizeInPhysicalUnits( localoptimizerlearningrate ); // * sqrt( small_step )
  localoptimizer->SetNumberOfIterations( localoptimizeriterations );
  localoptimizer->SetLowerLimit( 0 );
  localoptimizer->SetUpperLimit( 2 );
  localoptimizer->SetEpsilon( 0.1 );
  localoptimizer->SetMaximumLineSearchIterations( 20 );
  localoptimizer->SetDoEstimateLearningRateOnce( true );
  localoptimizer->SetMinimumConvergenceValue( 1.e-6 );
  localoptimizer->SetConvergenceWindowSize( 3 );

  std::cout << "Begin MultiStart: " << parametersList.size() << " searches between -/+ " << piover4 / pi
           << " radians " << std::endl;
  if( localoptimizeriterations > 0 )
    {
    mstartOptimizer->SetLocalOptimizer( localoptimizer );
    }
  mstartOptimizer->StartOptimization();
  std::cout << "done" << std::endl;
  typename AffineType::Pointer bestaffine = AffineType::New();
  bestaffine->SetCenter( trans2 );
  bestaffine->SetParameters( mstartOptimizer->GetBestParameters() );

  typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
  transformWriter->SetInput( bestaffine );
  transformWriter->SetFileName( outname.c_str() );
  transformWriter->Update();
  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int antsAffineInitializer( std::vector<std::string> args, std::ostream* /*out_stream = NULL */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "antsAffineInitializer" );

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

  if( argc < 3 )
    {
    std::cerr << "\nUsage: " << argv[0]
             <<
      " ImageDimension <Image1.ext> <Image2.ext> TransformOutput.mat Optional-SearchFactor Optional-Radian-Fraction Optional-bool-UsePrincipalAxes Optional-uint-UseLocalSearch Optional-Image1Mask "
             << std::endl;
    std::cerr << " Optional-SearchFactor is in degrees --- e.g. 10 = search in 10 degree increments ." << std::endl;
    std::cerr << " Radian-Fraction should be between 0 and 1 --- will search this arc +/- around principal axis."
             << std::endl;
    std::cerr
      <<
      " Optional-bool-UsePrincipalAxes determines whether the rotation is searched around an initial principal axis alignment.  Default = false. "
      << std::endl;
    std::cerr
      <<
      " Optional-uint-UseLocalSearch determines if a local optimization is run at each search point for the set number of iterations. Default = 20."
      << std::endl;
    return 0;
    }

  switch( atoi(argv[1]) )
    {
    case 2:
      {
      return antsAffineInitializerImp<2>(argc, argv);
      }
    case 3:
      {
      return antsAffineInitializerImp<3>(argc, argv);
      }
    case 4:
      return antsAffineInitializerImp<4>(argc, argv);
    }
  return antsAffineInitializerImp<2>(argc, argv);
}
} // namespace ants
