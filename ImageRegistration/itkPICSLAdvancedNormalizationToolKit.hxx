/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkPICSLAdvancedNormalizationToolKit_hxx_
#define _itkPICSLAdvancedNormalizationToolKit_hxx_

// disable debug warnings in MS compiler
#ifdef _MSC_VER
#  pragma warning(disable : 4786)
#endif

#include "itkHistogramMatchingImageFilter.h"
#include "itkSpatialMutualInformationRegistrationFunction.h"
#include "itkIdentityTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLabeledPointSetFileReader.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkResampleImageFilter.h"

// Image similarity metrics
#include "itkAvantsMutualInformationRegistrationFunction.h"
#include "itkSyNDemonsRegistrationFunction.h"
#include "itkProbabilisticRegistrationFunction.h"
#include "itkCrossCorrelationRegistrationFunction.h"
#include "itkExpectationBasedPointSetRegistrationFunction.h"
// #include "itkRobustDemonsRegistrationFunction.h"
// #include "itkRobustOpticalFlow.h"
// #include "itkSectionMutualInformationRegistrationFunction.h"

// #include "itkJensenTsallisBSplineRegistrationFunction.h"

#include "itkMath.h"

#include "ANTS_affine_registration2.h"

namespace itk
{
template <unsigned int TDimension, typename TReal>
PICSLAdvancedNormalizationToolKit<TDimension, TReal>::PICSLAdvancedNormalizationToolKit()
{
  this->InitializeCommandLineOptions();
}

template <unsigned int TDimension, typename TReal>
void
PICSLAdvancedNormalizationToolKit<TDimension, TReal>::ParseCommandLine(int argc, char ** argv)
{
  this->m_Parser->Parse(argc, argv);

  typename ParserType::OptionListType unknownOptions = this->m_Parser->GetUnknownOptions();
  if (unknownOptions.size())
  {
    std::cout << std::endl << "WARNING:  Unknown options" << std::endl;
    typename ParserType::OptionListType::const_iterator its;
    for (its = unknownOptions.begin(); its != unknownOptions.end(); its++)
    {
      if ((*its)->GetShortName() != '\0')
      {
        std::cout << "   " << '-' << (*its)->GetShortName() << std::endl;
      }
      else
      {
        std::cout << "   "
                  << "--" << (*its)->GetLongName() << std::endl;
      }
    }
    std::cout << std::endl;
  }

  std::string  printhelp_long = this->m_Parser->GetOption("help")->GetFunction(0)->GetName();
  unsigned int help_long = this->m_Parser->template Convert<unsigned int>(printhelp_long);
  if (help_long)
  {
    this->m_Parser->PrintMenu(std::cout, 7, false);
    std::exception();
    exit(EXIT_SUCCESS);
  }

  std::string  printhelp_short = this->m_Parser->GetOption('h')->GetFunction(0)->GetName();
  unsigned int help_short = this->m_Parser->template Convert<unsigned int>(printhelp_short);
  if (help_short)
  {
    this->m_Parser->PrintMenu(std::cout, 7, true);
    std::exception();
    exit(EXIT_SUCCESS);
  }
}

template <unsigned int TDimension, typename TReal>
void
PICSLAdvancedNormalizationToolKit<TDimension, TReal>::InitializeTransformAndOptimizer()
{
  if (!this->m_TransformationModel)
  {
    this->m_TransformationModel = TransformationModelType::New();
  }
  if (!this->m_RegistrationOptimizer)
  {
    this->m_RegistrationOptimizer = RegistrationOptimizerType::New();
  }
  std::string temp = this->m_Parser->GetOption("output-naming")->GetFunction(0)->GetName();
  this->m_TransformationModel->SetNamingConvention(temp);
  this->m_TransformationModel->InitializeTransform();
  this->m_RegistrationOptimizer->SetParser(this->m_Parser);
  this->m_RegistrationOptimizer->SetSimilarityMetrics(this->m_SimilarityMetrics);
}

template <unsigned int TDimension, typename TReal>
void
PICSLAdvancedNormalizationToolKit<TDimension, TReal>::RunRegistration()
{
  /** parse the command line and get input objects */
  this->ReadImagesAndMetrics();
  //    std::exception();
  /** initializes the transformation model and the optimizer */
  this->InitializeTransformAndOptimizer();

  /** Get the mask if there is one */

  typename OptionType::Pointer option = this->m_Parser->GetOption("mask-image");
  if (option && option->GetNumberOfFunctions())
  {
    typedef ImageFileReader<ImageType> ReaderType;
    std::string                        maskfn = this->m_Parser->GetOption("mask-image")->GetFunction(0)->GetName();
    if (maskfn.length() > 4)
    {
      typename ReaderType::Pointer maskImageFileReader = ReaderType::New();
      maskImageFileReader->SetFileName(maskfn.c_str());
      maskImageFileReader->Update();
      ImagePointer maskImage = this->PreprocessImage(maskImageFileReader->GetOutput());
      this->m_RegistrationOptimizer->SetMaskImage(maskImage);
    }
  }

  typename TransformationModelType::AffineTransformPointer aff_init, aff, fixed_aff_init;

  // added by songgang
  // try initialize the affine transform

  typename OptionType::Pointer initialAffineOption = this->m_Parser->GetOption("initial-affine");
  if (initialAffineOption && initialAffineOption->GetNumberOfFunctions())
  {
    std::string initial_affine_filename = initialAffineOption->GetFunction(0)->GetName();

    std::cout << "Loading affine registration from: " << initial_affine_filename << std::endl;
    ReadAffineTransformFile<typename TransformationModelType::AffineTransformType>(initial_affine_filename, aff_init);
  }
  else
  {
    std::cout << "Use identity affine transform as initial affine para." << std::endl;
    std::cout << "aff_init.IsNull()==" << aff_init.IsNull() << std::endl;
  }

  std::string                  fixed_initial_affine_filename = std::string("");
  typename OptionType::Pointer fixedInitialAffineOption = this->m_Parser->GetOption("fixed-image-initial-affine");
  if (fixedInitialAffineOption && fixedInitialAffineOption->GetNumberOfFunctions())
  {
    fixed_initial_affine_filename = fixedInitialAffineOption->GetFunction(0)->GetName();
    std::cout << "Loading affine registration from: " << fixed_initial_affine_filename << std::endl;
    fixed_aff_init = TransformationModelType::AffineTransformType::New();
    ReadAffineTransformFile<typename TransformationModelType::AffineTransformType>(fixed_initial_affine_filename,
                                                                                   fixed_aff_init);
    std::cout << " FIXME!  currently, if one passes a fixed initial affine mapping, then NO affine mapping will be "
                 "performed subsequently! "
              << std::endl;

    typename OptionType::Pointer fixedImageInitialAffineRegImageOption =
      this->m_Parser->GetOption("fixed-image-initial-affine-ref-image");
    if (fixedImageInitialAffineRegImageOption && fixedImageInitialAffineRegImageOption->GetNumberOfFunctions())
    {
      std::string refheader = fixedImageInitialAffineRegImageOption->GetFunction(0)->GetName();
      std::cout << " Setting reference deformation space by " << refheader << std::endl;
      typedef ImageFileReader<ImageType> ReaderType;
      typename ReaderType::Pointer       fixedImageFileReader = ReaderType::New();
      fixedImageFileReader->SetFileName(refheader.c_str());
      fixedImageFileReader->Update();
      ImagePointer temp = fixedImageFileReader->GetOutput();
      this->m_RegistrationOptimizer->SetReferenceSpaceImage(temp);
    }
  }
  else
  {
    std::cout << "Use identity affine transform as initial fixed affine para." << std::endl;
    std::cout << "fixed_aff_init.IsNull()==" << fixed_aff_init.IsNull() << std::endl;
  }

  typename OptionType::Pointer useNNOption = this->m_Parser->GetOption("use-NN");
  if (useNNOption && useNNOption->GetNumberOfFunctions() &&
      this->m_Parser->template Convert<bool>(useNNOption->GetFunction(0)->GetName()))
  {
    this->m_RegistrationOptimizer->SetUseNearestNeighborInterpolation(true);
  }
  else
  {
    this->m_RegistrationOptimizer->SetUseNearestNeighborInterpolation(false);
  }

  std::string continue_affine = this->m_Parser->GetOption("continue-affine")->GetFunction(0)->GetName();
  if (fixed_initial_affine_filename != "")
  {
    continue_affine = std::string("false");
  }
  if (continue_affine == "true")
  {
    std::cout << "Continue affine registration from the input" << std::endl; // << aff_init << std::endl;

    OptAffineType affine_opt;
    // InitializeAffineOption()
    {
      affine_opt.transform_initial = aff_init;
      std::string temp = this->m_Parser->GetOption("number-of-affine-iterations")->GetFunction(0)->GetName();
      affine_opt.number_of_iteration_list = this->m_Parser->template ConvertVector<int>(temp);
      affine_opt.number_of_levels = affine_opt.number_of_iteration_list.size();
      temp = this->m_Parser->GetOption("affine-metric-type")->GetFunction(0)->GetName();
      if (temp == "MI")
      {
        affine_opt.metric_type = AffineWithMutualInformation;
      }
      if (temp == "MSQ")
      {
        affine_opt.metric_type = AffineWithMeanSquareDifference;
      }
      if (temp == "CCH")
      {
        affine_opt.metric_type = AffineWithHistogramCorrelation;
      }
      if (temp == "CC")
      {
        affine_opt.metric_type = AffineWithNormalizedCorrelation;
      }
      if (temp == "GD")
      {
        affine_opt.metric_type = AffineWithGradientDifference;
      }
      temp = this->m_Parser->GetOption("MI-option")->GetFunction(0)->GetName();
      std::vector<int> mi_option = this->m_Parser->template ConvertVector<int>(temp);
      affine_opt.MI_bins = mi_option[0];
      affine_opt.MI_samples = mi_option[1];
      temp = this->m_Parser->GetOption("rigid-affine")->GetFunction(0)->GetName();
      std::string temp2 = this->m_Parser->GetOption("do-rigid")->GetFunction(0)->GetName();
      affine_opt.is_rigid = ((temp == "true") || (temp2 == "true") || (temp == "1") || (temp2 == "1"));
      temp = this->m_Parser->GetOption("affine-gradient-descent-option")->GetFunction(0)->GetName();
      std::vector<double> gradient_option = this->m_Parser->template ConvertVector<double>(temp);
      affine_opt.maximum_step_length = gradient_option[0];
      affine_opt.relaxation_factor = gradient_option[1];
      affine_opt.minimum_step_length = gradient_option[2];
      affine_opt.translation_scales = gradient_option[3];
      // std::cout << affine_opt;
      temp = this->m_Parser->GetOption("use-rotation-header")->GetFunction(0)->GetName();
      affine_opt.use_rotation_header = (temp == "true");
      std::cout << "affine_opt.use_rotation_header = " << affine_opt.use_rotation_header << std::endl;

      temp = this->m_Parser->GetOption("ignore-void-origin")->GetFunction(0)->GetName();
      affine_opt.ignore_void_orgin = (temp == "true");
      std::cout << "affine_opt.ignore_void_orgin = " << affine_opt.ignore_void_orgin << std::endl;
    }

    aff = this->m_RegistrationOptimizer->AffineOptimization(affine_opt);
  }
  else
  {
    std::cout << "Use fixed initial affine para." << std::endl;
    if (aff_init.IsNull())
    {
      aff_init = TransformationModelType::AffineTransformType::New();
      aff_init->SetIdentity();
    }
    if (fixed_aff_init.IsNull())
    {
      fixed_aff_init = TransformationModelType::AffineTransformType::New();
      fixed_aff_init->SetIdentity();
    }
    aff = aff_init;
  }
  // std::cout << aff << std::endl;

  this->m_TransformationModel->SetAffineTransform(aff);
  this->m_TransformationModel->SetFixedImageAffineTransform(fixed_aff_init);
  this->m_RegistrationOptimizer->SetAffineTransform(aff);
  this->m_RegistrationOptimizer->SetFixedImageAffineTransform(
    this->m_TransformationModel->GetFixedImageAffineTransform());

  /** Second, optimize Diff */
  this->m_RegistrationOptimizer->DeformableOptimization();
  std::cout << " Registration Done " << std::endl;
  this->m_TransformationModel->SetDisplacementField(this->m_RegistrationOptimizer->GetDisplacementField());
  this->m_TransformationModel->SetInverseDisplacementField(
    this->m_RegistrationOptimizer->GetInverseDisplacementField());
}

template <unsigned int TDimension, typename TReal>
typename PICSLAdvancedNormalizationToolKit<TDimension, TReal>::ImagePointer
PICSLAdvancedNormalizationToolKit<TDimension, TReal>::PreprocessImage(ImagePointer image)
{
  image = this->ReplaceProblematicPixelValues(image, 0);
  return image;
}

template <unsigned int TDimension, typename TReal>
typename PICSLAdvancedNormalizationToolKit<TDimension, TReal>::ImagePointer
PICSLAdvancedNormalizationToolKit<TDimension, TReal>::ReplaceProblematicPixelValues(ImagePointer image, PixelType value)
{
  ImageRegionIterator<ImageType> It(image, image->GetRequestedRegion());
  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    PixelType pixel = It.Get();
    if (std::isinf(pixel) || std::isnan(pixel))
    {
      It.Set(value);
    }
  }

  typedef itk::MinimumMaximumImageFilter<ImageType> MinMaxFilterType;
  typename MinMaxFilterType::Pointer                minMaxFilter = MinMaxFilterType::New();

  minMaxFilter->SetInput(image);
  minMaxFilter->Update();

  double min = minMaxFilter->GetMinimum();
  double shift = -1.0 * static_cast<double>(min);
  double scale = static_cast<double>(minMaxFilter->GetMaximum());
  scale += shift;
  scale = 1.0 / scale;

  typedef itk::ShiftScaleImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer                             filter = FilterType::New();

  filter->SetInput(image);
  filter->SetShift(shift);
  filter->SetScale(scale);
  filter->Update();

  // when the image is a constant image, the scale filter will fail
  // replace again everything into $pixel
  ImagePointer                   image2 = filter->GetOutput();
  ImageRegionIterator<ImageType> It2(image2, image2->GetRequestedRegion());
  for (It2.GoToBegin(); !It2.IsAtEnd(); ++It2)
  {
    PixelType pixel = It2.Get();
    if (std::isinf(pixel) || std::isnan(pixel))
    {
      It2.Set(value);
    }
  }

  return image2;
}

/**
 * Standard "PrintSelf" method
 */
template <unsigned int TDimension, typename TReal>
void
PICSLAdvancedNormalizationToolKit<TDimension, TReal>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

template <unsigned int TDimension, typename TReal>
void
PICSLAdvancedNormalizationToolKit<TDimension, TReal>::ReadImagesAndMetrics()
{
  /**
   * Read in all images and preprocess them before
   * storing them in their corresponding image lists.
   */
  this->m_SimilarityMetrics.clear();

  typedef ImageFileReader<ImageType> ReaderType;
  bool                               useHistMatch = this->m_Parser->template Convert<bool>(
    this->m_Parser->GetOption("use-Histogram-Matching")->GetFunction()->GetName());

  /**
   * Read the metrics and image files
   */
  if (typename OptionType::Pointer option = this->m_Parser->GetOption("image-metric"))
  {
    std::cout << " values " << option->GetNumberOfFunctions() << std::endl;
    for (unsigned int i = 0; i < option->GetNumberOfFunctions(); i++)
    {
      SimilarityMetricPointer similarityMetric = SimilarityMetricType::New();
      RealType                similarityMetricScalarWeight = 1.0;
      similarityMetric->SetWeightScalar(similarityMetricScalarWeight);

      unsigned int parameterCount = 0;

      typename ReaderType::Pointer fixedImageFileReader = ReaderType::New();
      fixedImageFileReader->SetFileName(option->GetFunction(i)->GetParameter(parameterCount));
      fixedImageFileReader->Update();
      ImagePointer fixedImage = this->PreprocessImage(fixedImageFileReader->GetOutput());
      similarityMetric->SetFixedImage(fixedImage);
      parameterCount++;

      std::cout << "  Fixed image file: " << fixedImageFileReader->GetFileName() << std::endl;

      typename ReaderType::Pointer movingImageFileReader = ReaderType::New();
      movingImageFileReader->SetFileName(option->GetFunction(i)->GetParameter(parameterCount));
      movingImageFileReader->Update();
      ImagePointer movingImage = this->PreprocessImage(movingImageFileReader->GetOutput());
      similarityMetric->SetMovingImage(movingImage);
      typename SimilarityMetricType::RadiusType radius;
      radius.Fill(0);
      parameterCount++;

      std::cout << "  Moving image file: " << movingImageFileReader->GetFileName() << std::endl;

      /**
       * Check if similarity metric is image based or point-set based.
       *   Image metrics:
       *     > mean-squares/MeanSquares/MSQ
       *     > mutual-information/MutualInformation/MI
       *     > cross-correlation/CrossCorrelation/CC
       *     > probabilistic/Probabilistic/PR
       *
       *   Point-set metrics:
       *     > point-set-expectation/PointSetExpectation/PSE
       *     > jensen-tsallis-bspline/JensenTsallisBSpline/JTB
       *
       */

      bool isMetricPointSetBased = false;

      std::string whichMetric = option->GetFunction(i)->GetName();

      if (whichMetric == "point-set-expectation" || whichMetric == "PointSetExpectation" || whichMetric == "PSE"
          //                whichMetric == "jensen-tsallis-bspline" ||
          //                whichMetric == "JensenTsallisBSpline" ||
          //                whichMetric == "JTB"
      )
      {
        isMetricPointSetBased = true;
      }

      if (isMetricPointSetBased)
      {
        std::cout << "Metric " << i << ": "
                  << " Point-set " << whichMetric << " n-params " << option->GetFunction(i)->GetNumberOfParameters()
                  << std::endl;
        /**
         * Read in the point-set metric parameters
         */

        typedef LabeledPointSetFileReader<PointSetType> PointSetReaderType;

        typename PointSetReaderType::Pointer fixedPointSetReader = PointSetReaderType::New();
        fixedPointSetReader->SetFileName(option->GetFunction(i)->GetParameter(parameterCount));
        parameterCount++;

        typename PointSetReaderType::Pointer movingPointSetReader = PointSetReaderType::New();
        movingPointSetReader->SetFileName(option->GetFunction(i)->GetParameter(parameterCount));
        parameterCount++;

        if (option->GetFunction(i)->GetNumberOfParameters() > parameterCount)
        {
          similarityMetricScalarWeight =
            this->m_Parser->template Convert<RealType>(option->GetFunction(i)->GetParameter(parameterCount));
          parameterCount++;
        }
        similarityMetric->SetWeightScalar(similarityMetricScalarWeight);
        std::cout << "  similarity metric weight: " << similarityMetricScalarWeight << std::endl;

        TReal        pointSetPercent = 0.5;
        TReal        pointSetSigma = 5.0;
        bool         extractBoundaryPointsOnly = false;
        unsigned int kNeighborhood = 50;
        if (option->GetFunction(i)->GetNumberOfParameters() > parameterCount)
        {
          pointSetPercent =
            this->m_Parser->template Convert<TReal>(option->GetFunction(i)->GetParameter(parameterCount));
          parameterCount++;
        }
        if (option->GetFunction(i)->GetNumberOfParameters() > parameterCount)
        {
          pointSetSigma = this->m_Parser->template Convert<TReal>(option->GetFunction(i)->GetParameter(parameterCount));
          parameterCount++;
        }
        if (option->GetFunction(i)->GetNumberOfParameters() > parameterCount)
        {
          extractBoundaryPointsOnly =
            this->m_Parser->template Convert<bool>(option->GetFunction(i)->GetParameter(parameterCount));
          parameterCount++;
        }
        if (option->GetFunction(i)->GetNumberOfParameters() > parameterCount)
        {
          kNeighborhood =
            this->m_Parser->template Convert<unsigned int>(option->GetFunction(i)->GetParameter(parameterCount));
          parameterCount++;
        }
        std::cout << " point-set sigma = " << pointSetSigma << std::endl;
        std::cout << " percentage of points = " << pointSetPercent << std::endl;
        std::cout << " k-neighborhood = " << kNeighborhood << std::endl;
        if (extractBoundaryPointsOnly)
        {
          std::cout << " use only boundary points. " << pointSetPercent << std::endl;
        }

        fixedPointSetReader->SetRandomPercentage(pointSetPercent);
        fixedPointSetReader->SetExtractBoundaryPoints(extractBoundaryPointsOnly);
        fixedPointSetReader->Update();
        std::cout << "  Fixed point-set file: " << fixedPointSetReader->GetFileName() << std::endl;
        std::cout << "    Number of fixed labels: " << fixedPointSetReader->GetLabelSet()->size() << std::endl;
        std::cout << "    Distinct fixed labels: ";
        for (unsigned int n = 0; n < fixedPointSetReader->GetLabelSet()->size(); n++)
        {
          std::cout << fixedPointSetReader->GetLabelSet()->operator[](n) << " ";
        }
        std::cout << std::endl;

        movingPointSetReader->SetRandomPercentage(pointSetPercent);
        movingPointSetReader->SetExtractBoundaryPoints(extractBoundaryPointsOnly);
        movingPointSetReader->Update();

        movingPointSetReader->SetRandomPercentage(pointSetPercent);
        movingPointSetReader->SetExtractBoundaryPoints(extractBoundaryPointsOnly);
        movingPointSetReader->Update();
        std::cout << "  Moving point-set file: " << movingPointSetReader->GetFileName() << std::endl;
        std::cout << "    Number of moving labels: " << movingPointSetReader->GetLabelSet()->size() << std::endl;
        std::cout << "    Distinct moving labels: ";
        for (unsigned int n = 0; n < movingPointSetReader->GetLabelSet()->size(); n++)
        {
          std::cout << movingPointSetReader->GetLabelSet()->operator[](n) << " ";
        }
        std::cout << std::endl;

        if (whichMetric == "point-set-expectation" || whichMetric == "PointSetExpectation" || whichMetric == "PSE")
        {
          typedef itk::
            ExpectationBasedPointSetRegistrationFunction<ImageType, ImageType, DisplacementFieldType, PointSetType>
                                       MetricType;
          typename MetricType::Pointer metric = MetricType::New();
          metric->SetRadius(radius);
          metric->SetFixedPointSet(fixedPointSetReader->GetOutput());
          metric->SetMovingPointSet(movingPointSetReader->GetOutput());
          metric->SetFixedPointSetSigma(pointSetSigma);
          metric->SetMovingPointSetSigma(pointSetSigma);
          metric->SetKNeighborhood(kNeighborhood);

          if (option->GetFunction(i)->GetNumberOfParameters() > parameterCount)
          {
            unsigned int pm =
              this->m_Parser->template Convert<unsigned int>(option->GetFunction(i)->GetParameter(parameterCount));
            metric->SetUseSymmetricMatching(pm);
            std::cout << " Symmetric match iterations -- going Asymmeric for the rest " << pm << std::endl;
            parameterCount++;
          }

          similarityMetric->SetMetric(metric);
          similarityMetric->SetMaximizeMetric(true);
          similarityMetric->SetFixedPointSet(fixedPointSetReader->GetOutput());
          similarityMetric->SetMovingPointSet(movingPointSetReader->GetOutput());
          this->m_SimilarityMetrics.push_back(similarityMetric);
        }
        //              else if ( whichMetric == "jensen-tsallis-bspline" ||
        //                        whichMetric == "JensenTsallisBSpline" ||
        //                        whichMetric == "JTB" )
        //                {
        //                typedef itk::JensenTsallisBSplineRegistrationFunction
        //                        <ImageType, PointSetType, ImageType, PointSetType, DisplacementFieldType> MetricType;
        //                typename MetricType::Pointer metric = MetricType::New();
        //                metric->SetRadius( radius );
        //                metric->SetFixedPointSet( fixedPointSetReader->GetOutput() );
        //                metric->SetMovingPointSet( movingPointSetReader->GetOutput() );
        //                metric->SetFixedPointSetSigma( pointSetSigma );
        //                metric->SetMovingPointSetSigma( pointSetSigma );
        //                metric->SetFixedEvaluationKNeighborhood( kNeighborhood );
        //                metric->SetMovingEvaluationKNeighborhood( kNeighborhood );
        //
        //                if ( option->GetFunction( i )->GetNumberOfParameters() > parameterCount )
        //                  {
        //                  metric->SetAlpha( this->m_Parser->template
        //                  Convert<TReal>( option->GetFunction( i )->GetParameter(  parameterCount ) ) );
        //                  parameterCount++;
        //                  }
        //                if ( option->GetFunction( i )->GetNumberOfParameters() > parameterCount )
        //                  {
        //                  typename RegistrationOptimizerType::ArrayType meshResolution;
        //                  std::vector<TReal> resolution = this->m_Parser->template
        //                    ConvertVector<TReal>( option->GetFunction( i )->GetParameter(  parameterCount ) );
        //                  if ( resolution.size() != TDimension )
        //                    {
        //                    itkExceptionMacro( "Mesh resolution does not match image dimension." );
        //                    }
        //                  for ( unsigned int d = 0; d < TDimension; d++ )
        //                    {
        //                    meshResolution[d] = resolution[d];
        //                    }
        //                  metric->SetMeshResolution( meshResolution );
        //                  parameterCount++;
        //                  }
        //                if ( option->GetFunction( i )->GetNumberOfParameters() > parameterCount )
        //                  {
        //                  metric->SetSplineOrder( this->m_Parser->template
        //                  Convert<unsigned int>( option->GetFunction( i )->GetParameter(  parameterCount ) ) );
        //                  parameterCount++;
        //                  }
        //                if ( option->GetFunction( i )->GetNumberOfParameters() > parameterCount )
        //                  {
        //                  metric->SetNumberOfLevels( this->m_Parser->template
        //                  Convert<unsigned int>( option->GetFunction( i )->GetParameter(  parameterCount ) ) );
        //                  parameterCount++;
        //                  }
        //                if ( option->GetFunction( i )->GetNumberOfParameters() > parameterCount )
        //                  {
        //                  metric->SetUseAnisotropicCovariances(
        //                    this->m_Parser->template Convert<bool>(
        //                    option->GetFunction( i )->GetParameter(  parameterCount ) ) );
        //                  parameterCount++;
        //                  }
        //
        //
        //                std::cout << "  B-spline parameters " << std::endl;
        //                std::cout << "    mesh resolution: " << metric->GetMeshResolution() << std::endl;
        //                std::cout << "    spline order: " << metric->GetSplineOrder() << std::endl;
        //                std::cout << "    number of levels: " << metric->GetNumberOfLevels() << std::endl;
        //                std::cout << "  Alpha: " << metric->GetAlpha() << std::endl;
        //                if ( metric->GetUseAnisotropicCovariances() )
        //                  {
        //                  std::cout << "  using anisotropic covariances." << std::endl;
        //                  }
        //
        //                similarityMetric->SetMetric( metric );
        //                similarityMetric->SetMaximizeMetric( true );
        //                similarityMetric->SetFixedPointSet( fixedPointSetReader->GetOutput() );
        //                similarityMetric->SetMovingPointSet( movingPointSetReader->GetOutput() );
        //                this->m_SimilarityMetrics.push_back( similarityMetric );
        //                }
      }
      else // similarity metric is image-based
      {
        std::cout << "Metric " << i << ": "
                  << " Not a Point-set" << std::endl;
        std::cout << "  Fixed image file: " << fixedImageFileReader->GetFileName() << std::endl;
        std::cout << "  Moving image file: " << movingImageFileReader->GetFileName() << std::endl;

        similarityMetric->SetFixedPointSet(nullptr);
        similarityMetric->SetMovingPointSet(nullptr);

        if (option->GetFunction(i)->GetNumberOfParameters() > parameterCount)
        {
          similarityMetricScalarWeight =
            this->m_Parser->template Convert<RealType>(option->GetFunction(i)->GetParameter(parameterCount));
          parameterCount++;
        }
        similarityMetric->SetWeightScalar(similarityMetricScalarWeight);
        std::cout << "  similarity metric weight: " << similarityMetricScalarWeight << std::endl;

        radius.Fill(0);
        if (option->GetFunction(i)->GetNumberOfParameters() > parameterCount)
        {
          std::vector<unsigned int> rad =
            this->m_Parser->template ConvertVector<unsigned int>(option->GetFunction(i)->GetParameter(parameterCount));

          if (rad.size() == 1)
          {
            radius.Fill(rad[0]);
          }
          else if (rad.size() == TDimension)
          {
            for (unsigned int n = 0; n < TDimension; n++)
            {
              radius[n] = rad[n];
            }
          }
          else
          {
            std::cout << "Badly formed radius specification" << std::endl;
            std::exception();
          }
          parameterCount++;
        }
        std::cout << "  Radius: " << radius << std::endl;

        TReal extraparam = -1.e12;
        if (option->GetFunction(i)->GetNumberOfParameters() > parameterCount)
        {
          extraparam = this->m_Parser->template Convert<TReal>(option->GetFunction(i)->GetParameter(parameterCount));
          std::cout << " Setting Extra Param to :  " << extraparam
                    << " often used as a robustness parameter for longitudinal studies " << std::endl;
          parameterCount++;
        }
        std::cout << "  radius: " << radius << std::endl;

        unsigned int numberOfHistogramBins = 64;
        if (Dimension == 2)
        {
          numberOfHistogramBins = 32;
        }
        if (whichMetric == "mean-squares" || whichMetric == "MeanSquares" || whichMetric == "MSQ")
        {
          typedef itk::HistogramMatchingImageFilter<ImageType, ImageType> FilterType;
          typename FilterType::Pointer                                    filter = FilterType::New();
          filter->SetSourceImage(movingImage);
          filter->SetReferenceImage(fixedImage);
          filter->SetNumberOfHistogramLevels(256);
          filter->SetNumberOfMatchPoints(12);
          filter->ThresholdAtMeanIntensityOn();
          //  filter->ThresholdAtMeanIntensityOff();
          if (useHistMatch)
          {
            filter->Update();
            std::cout << " use Histogram Matching " << std::endl;
            movingImage = filter->GetOutput();
            movingImage = this->PreprocessImage(movingImage);
            similarityMetric->SetMovingImage(movingImage);
          }
          typedef SyNDemonsRegistrationFunction<ImageType, ImageType, DisplacementFieldType> MetricType;
          typename MetricType::Pointer                                                       metric = MetricType::New();
          if (radius[0] > 0)
          {
            metric->SetUseMovingImageGradient(true);
          }
          metric->SetIntensityDifferenceThreshold(extraparam);
          metric->SetRobustnessParameter(extraparam);
          //        metric->SetRobust( true );
          //        metric->SetSymmetric( false );
          // /        metric->SetNormalizeGradient( false );
          metric->SetRadius(radius);

          similarityMetric->SetMetric(metric);
          similarityMetric->SetMaximizeMetric(true);
          this->m_SimilarityMetrics.push_back(similarityMetric);
        }
        else if (whichMetric == "SSD")
        {
          if (useHistMatch)
          {
            movingImage = this->PreprocessImage(movingImage);
            typedef itk::HistogramMatchingImageFilter<ImageType, ImageType> FilterType;
            typename FilterType::Pointer                                    filter = FilterType::New();
            filter->SetSourceImage(movingImage);
            filter->SetReferenceImage(fixedImage);
            filter->SetNumberOfHistogramLevels(256);
            filter->SetNumberOfMatchPoints(12);
            filter->ThresholdAtMeanIntensityOn();
            filter->Update();
            std::cout << " use Histogram Matching " << std::endl;
            movingImage = filter->GetOutput();
          }
          similarityMetric->SetMovingImage(movingImage);
          typedef SyNDemonsRegistrationFunction<ImageType, ImageType, DisplacementFieldType> MetricType;
          typename MetricType::Pointer                                                       metric = MetricType::New();
          if (radius[0] > 0)
          {
            metric->SetUseMovingImageGradient(true);
          }
          metric->SetIntensityDifferenceThreshold(extraparam);
          metric->SetRobustnessParameter(extraparam);
          metric->SetUseSSD(true);
          metric->SetRadius(radius);
          similarityMetric->SetMetric(metric);
          similarityMetric->SetMaximizeMetric(true);
          this->m_SimilarityMetrics.push_back(similarityMetric);
        }
        else if (whichMetric == "mutual-information" || whichMetric == "MutualInformation" || whichMetric == "MI")
        {
          typedef itk::AvantsMutualInformationRegistrationFunction<ImageType, ImageType, DisplacementFieldType>
                                       MetricType;
          typename MetricType::Pointer metric = MetricType::New();
          metric->SetNumberOfHistogramBins(numberOfHistogramBins);
          metric->SetNormalizeGradient(false);
          metric->SetRobustnessParameter(extraparam);
          unsigned int histbins = radius[0];
          if (histbins < 8)
          {
            histbins = 8;
          }
          metric->SetNumberOfHistogramBins(histbins);
          radius.Fill(0);
          metric->SetRadius(radius);
          similarityMetric->SetMetric(metric);
          similarityMetric->SetMaximizeMetric(true);
          this->m_SimilarityMetrics.push_back(similarityMetric);
        }
        else if (whichMetric == "spatial-mutual-information" || whichMetric == "SpatialMutualInformation" ||
                 whichMetric == "SMI")
        {
          typedef itk::SpatialMutualInformationRegistrationFunction<ImageType, ImageType, DisplacementFieldType>
                                       MetricType;
          typename MetricType::Pointer metric = MetricType::New();
          metric->SetNumberOfHistogramBins(numberOfHistogramBins);
          metric->SetNormalizeGradient(false);
          metric->SetRobustnessParameter(extraparam);
          unsigned int histbins = radius[0];
          if (histbins < 8)
          {
            histbins = 8;
          }
          metric->SetNumberOfHistogramBins(histbins);
          radius.Fill(0);
          metric->SetRadius(radius);
          similarityMetric->SetMetric(metric);
          similarityMetric->SetMaximizeMetric(true);
          this->m_SimilarityMetrics.push_back(similarityMetric);
        }
        else if (whichMetric == "cross-correlation" || whichMetric == "CrossCorrelation" || whichMetric == "CC")
        {
          typedef itk::HistogramMatchingImageFilter<ImageType, ImageType> FilterType;
          typename FilterType::Pointer                                    filter = FilterType::New();
          filter->SetSourceImage(movingImage);
          filter->SetReferenceImage(fixedImage);
          filter->SetNumberOfHistogramLevels(256);
          filter->SetNumberOfMatchPoints(12);
          filter->ThresholdAtMeanIntensityOn();
          //  filter->ThresholdAtMeanIntensityOff();
          if (useHistMatch)
          {
            std::cout << " use Histogram Matching " << std::endl;
            filter->Update();
            movingImage = filter->GetOutput();
            std::cout << " prepro " << std::endl;
            movingImage = this->PreprocessImage(movingImage);
            std::cout << " set " << std::endl;
            similarityMetric->SetMovingImage(movingImage);
          }

          typedef itk::CrossCorrelationRegistrationFunction<ImageType, ImageType, DisplacementFieldType> MetricType;
          typename MetricType::Pointer metric = MetricType::New();
          metric->SetNormalizeGradient(false);
          metric->SetRadius(radius);
          metric->SetRobustnessParameter(extraparam);
          similarityMetric->SetMetric(metric);
          similarityMetric->SetMaximizeMetric(true);
          this->m_SimilarityMetrics.push_back(similarityMetric);
        }
        else if (whichMetric == "probabilistic" || whichMetric == "Probabilistic" || whichMetric == "PR")
        {
          typedef itk::HistogramMatchingImageFilter<ImageType, ImageType> FilterType;
          typename FilterType::Pointer                                    filter = FilterType::New();
          filter->SetSourceImage(movingImage);
          filter->SetReferenceImage(fixedImage);
          filter->SetNumberOfHistogramLevels(256);
          filter->SetNumberOfMatchPoints(12);
          filter->ThresholdAtMeanIntensityOn();
          //  filter->ThresholdAtMeanIntensityOff();
          if (useHistMatch)
          {
            filter->Update();
            std::cout << " use Histogram Matching " << std::endl;
            movingImage = filter->GetOutput();
            movingImage = this->PreprocessImage(movingImage);
            similarityMetric->SetMovingImage(movingImage);
          }
          typedef itk::ProbabilisticRegistrationFunction<ImageType, ImageType, DisplacementFieldType> MetricType;
          typename MetricType::Pointer metric = MetricType::New();
          metric->SetNormalizeGradient(false);
          metric->SetRadius(radius);
          metric->SetRobustnessParameter(extraparam);

          similarityMetric->SetMetric(metric);
          similarityMetric->SetMaximizeMetric(true);
          this->m_SimilarityMetrics.push_back(similarityMetric);
        }
        else
        {
          itkWarningMacro("Could not decipher image metric choice: " << whichMetric);
        }
      }
    }
  }
  else
  {
    itkExceptionMacro("No image metric specified on the command line.");
  }
}

template <unsigned int TDimension, typename TReal>
void
PICSLAdvancedNormalizationToolKit<TDimension, TReal>::InitializeCommandLineOptions()
{
  this->m_Parser = ParserType::New();
  //    this->m_Parser->SetCommandDescription( " PICSL Advanced Image Normalization Toolkit" );

  if (false)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("weight-image");
    option->SetShortName('w');
    option->SetDescription("blah-blah");

    this->m_Parser->AddOption(option);
  }

  if (true)
  {
    std::string description = std::string("this mask -- defined in the 'fixed' image space defines ") +
                              std::string("the region of interest over which the registration is ") +
                              std::string("computed ==> above 0.1 means inside mask ==> continuous ") +
                              std::string("values in range [0.1,1.0] effect optimization like a ") +
                              std::string("probability.  ==> values > 1 are treated as = 1.0 ");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("mask-image");
    option->SetShortName('x');
    //         option->SetDescription( "this mask -- defined in the 'fixed' image space defines the region of interest
    //         over
    // which the registration is computed ==> above 0.1 means inside mask \n\t\t==> continuous values in range [0.1
    // , 1.0] effect optimization like a probability.  \n\t\t==> values > 1 are treated as = 1.0 " );
    option->SetDescription(description);
    option->SetUsageOption(0, "maskFileName");
    this->m_Parser->AddOption(option);
  }

  if (false)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("mask-threshold");
    option->SetShortName('M');
    option->SetDescription("blah-blah");

    this->m_Parser->AddOption(option);
  }

  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("image-metric");
    option->SetShortName('m');

    std::string newLineTabs("\n\t\t");

    std::string intensityBasedDescription("Intensity-Based Metrics: ");
    std::string weightDescription(
      "The metric weights are relative to the weights on the N other metrics passed to ANTS --- N is unlimited.  So, "
      "the weight, w_i on the i^{th} metric will be  w_i=w_i/ ( sum_i  w_i ).");
    std::string intensityBasedOptions("[fixedImage,movingImage,weight,radius/OrForMI-#histogramBins]");
    std::string ccDescription("CC/cross-correlation/CrossCorrelation");
    std::string miDescription("MI/mutual-information/MutualInformation");
    std::string smiDescription("SMI/spatial-mutual-information/SpatialMutualInformation");
    std::string prDescription("PR/probabilistic/Probabilistic");
    std::string msqDescription(
      "MSQ/mean-squares/MeanSquares -- demons-like, radius > 0 uses moving image gradient in metric deriv.");
    std::string ssdDescription("SSD --- standard intensity difference.");
    intensityBasedDescription +=
      (newLineTabs + ccDescription + intensityBasedOptions + newLineTabs + miDescription + intensityBasedOptions +
       newLineTabs + smiDescription + intensityBasedOptions + newLineTabs + prDescription + intensityBasedOptions +
       newLineTabs + ssdDescription + intensityBasedOptions + newLineTabs + msqDescription + intensityBasedOptions);

    std::string pointBasedDescription("\n\t      Point-Set-Based Metrics: ");
    std::string pointBasedOptions = std::string("[fixedImage,movingImage,fixedPoints,movingPoints") +
                                    std::string(",weight,pointSetPercentage,pointSetSigma,boundaryPointsOnly") +
                                    std::string(",kNeighborhood");
    std::string pseDescription("PSE/point-set-expectation/PointSetExpectation");
    std::string pseOptions(", PartialMatchingIterations=100000]   \n the partial matching option assumes the complete "
                           "labeling is in the first set of label parameters ... more iterations leads to more "
                           "symmetry in the matching  - 0 iterations means full asymmetry ");
    //      std::string jtbDescription( "JTB/jensen-tsallis-bspline/JensenTsallisBSpline" );
    std::string jtbOptions =
      std::string(",alpha,meshResolution,splineOrder,numberOfLevels") + std::string(",useAnisotropicCovariances]");
    pointBasedDescription += (newLineTabs + pseDescription + pointBasedOptions + pseOptions + newLineTabs +
                              /*jtbDescription + */ pointBasedOptions + jtbOptions);
    std::string description = weightDescription + intensityBasedDescription + pointBasedDescription;

    option->SetDescription(description);
    this->m_Parser->AddOption(option);
  }

  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("output-naming");
    option->SetShortName('o');
    option->SetDescription(
      "The name for the output - a prefix or a name+type : e.g.  -o OUT or  -o OUT.nii  or -o OUT.mha ");
    this->m_Parser->AddOption(option);
  }

  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("roi");
    option->SetLongName("R");
    option->SetDescription("TODO/FIXME: the --R sets an ROI option -- it passes a vector of parameters that sets the "
                           "center and bounding box \n of the region of interest for a sub-field registration.   e.g. "
                           "in 3D  the option setting \n  -r   10x12x15x50x50x25 \n sets up a bounding box of size "
                           "50,50,25 with origin at   10,12,15 in voxel  (should this be physical?) coordinates. ");
    std::string roidefault = std::string("0");
    /** set up a default parameter */
    option->AddFunction(roidefault);
    this->m_Parser->AddOption(option);
  }

  if (false)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("number-of-levels");
    option->SetShortName('n');
    option->SetDescription("number of levels in multi-resolution optimization -- an integer :: 3 is a common choice ");
    std::string nlevdefault = std::string("3");
    /** set up a default parameter */
    option->AddFunction(nlevdefault);
    this->m_Parser->AddOption(option);
  }

  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("number-of-iterations");
    option->SetShortName('i');
    option->SetDescription("number of iterations per level -- a 'vector' e.g.  :  100x100x20 ");
    std::string nitdefault = std::string("10x10x5");
    /** set up a default parameter */
    option->AddFunction(nitdefault);
    this->m_Parser->AddOption(option);
  }

  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("Restrict-Deformation");
    option->SetDescription(
      "restrict the gradient that drives the deformation by scalar factors along specified dimensions -- a TReal "
      "'vector' of length ImageDimension to multiply against the similarity metric's gradient values ---  e.g. in 3D : "
      "0.1x1x0 --- will set the z gradient to zero and scale the x gradient by 0.1 and y by 1 (no change). Thus, you "
      "get a 2.5-Dimensional registration as there is still 3D continuity in the mapping. ");
    std::string nitdefault;
    if (TDimension == 2)
    {
      nitdefault = std::string("1x1");
    }
    if (TDimension == 3)
    {
      nitdefault = std::string("1x1x1");
    }
    /** set up a default parameter */
    option->AddFunction(nitdefault);
    this->m_Parser->AddOption(option);
  }

  if (false)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("number-to-interpolate");
    option->SetShortName('b');
    option->SetDescription("blah-blah");

    this->m_Parser->AddOption(option);
  }

  if (false)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("use-converge-criteria");
    option->SetShortName('a');
    option->SetDescription("blah-blah");

    this->m_Parser->AddOption(option);
  }

  if (false)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("thickness-limit");
    option->SetShortName('T');
    option->SetDescription("blah-blah");

    this->m_Parser->AddOption(option);
  }

  if (false)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("prior-image-list");
    option->SetShortName('P');
    option->SetDescription("blah-blah");

    this->m_Parser->AddOption(option);
  }

  if (false)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("use-laplacian");
    option->SetShortName('i');
    option->SetDescription("blah-blah");

    this->m_Parser->AddOption(option);
  }

  if (false)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("use-histogram-matching");
    option->SetShortName('e');
    option->SetDescription("blah-blah");

    this->m_Parser->AddOption(option);
  }

  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("verbose");
    option->SetShortName('v');
    option->SetDescription(" verbose output ");

    this->m_Parser->AddOption(option);
  }

  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("use-all-metrics-for-convergence");
    option->SetDescription(" enable to use weighted sum of all metric terms for convergence computation. By default, "
                           "only the first metric is used");
    std::string zero = std::string("0");
    option->AddFunction(zero);
    this->m_Parser->AddOption(option);
  }

  if (true)
  {
    std::string description = std::string("Print the help menu (short version).");

    OptionType::Pointer option = OptionType::New();
    option->SetShortName('h');
    option->SetDescription(description);
    std::string zero = std::string("0");
    option->AddFunction(zero);
    this->m_Parser->AddOption(option);
  }

  if (true)
  {
    std::string description = std::string("Print the help menu.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("help");
    option->SetDescription(description);
    std::string zero = std::string("0");
    option->AddFunction(zero);
    this->m_Parser->AddOption(option);
  }

  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("transformation-model");
    option->SetShortName('t');
    option->SetDescription(
      "TRANSFORMATION[gradient-step-length,number-of-time-steps,DeltaTime,symmetry-type].\n\t      Choose one of the "
      "following TRANSFORMATIONS:\n\t\tDiff = diffeomorphic\n\t\tElast = Elastic\n\t\tExp = exponential diff\n\t\t "
      "Greedy Exp = greedy exponential diff, like diffeomorphic demons. same parameters. \n\t\tSyN -- symmetric "
      "normalization \n \n DeltaTime is the integration time-discretization step - sub-voxel - n-time steps currently "
      "fixed at 2 ");
    std::string nitdefault = std::string("SyN[0.5]");
    /** set up a default parameter */
    option->AddFunction(nitdefault);
    this->m_Parser->AddOption(option);
  }

  //    if (true)
  //    {
  //        OptionType::Pointer option = OptionType::New();
  //        option->SetLongName( "def-field-sigma" );
  //        option->SetShortName( 's' );
  //        option->SetDescription( "smoothing of deformation field " );
  //
  //        this->m_Parser->AddOption( option );
  //    }
  //
  //    if (true)
  //    {
  //        OptionType::Pointer option = OptionType::New();
  //        option->SetLongName( "gradient-field-sigma" );
  //        option->SetShortName( 'g' );
  //        option->SetDescription( "this smooths the gradient update field" );
  //        option->AddFunction("0.0");
  //        this->m_Parser->AddOption( option );
  //    }
  //
  //    if (true)
  //    {
  //        OptionType::Pointer option = OptionType::New();
  //        option->SetLongName( "gradient-step-length" );
  //        option->SetShortName( 'l' );
  //        option->SetDescription( "gradient descent parameter - a TReal :: e.g. 1.0 " );
  //        option->AddFunction("1.0");
  //        this->m_Parser->AddOption( option );
  //    }

  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("regularization");
    option->SetShortName('r');
    option->SetDescription(
      "REGULARIZATION[gradient-field-sigma,def-field-sigma,truncation].\n\t      Choose one of the following "
      "REGULARIZATIONS:\n\t\tGauss = gaussian\n\t\tDMFFD = directly manipulated free form deformation");
    std::string nitdefault = std::string("Gauss[3,0.5]");
    /** set up a default parameter */
    option->AddFunction(nitdefault);
    this->m_Parser->AddOption(option);
  }

  // added by songgang
  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("initial-affine");
    option->SetShortName('a');
    option->SetDescription("use the input file as the initial affine parameter");
    this->m_Parser->AddOption(option);
  }
  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("fixed-image-initial-affine");
    option->SetShortName('F');
    option->SetDescription("use the input file as the initial affine parameter for the fixed image ");
    this->m_Parser->AddOption(option);
  }
  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("fixed-image-initial-affine-ref-image");
    option->SetDescription(
      "reference space for using the input file as the initial affine parameter for the fixed image ");
    this->m_Parser->AddOption(option);
  }
  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("geodesic");
    option->SetShortName('T');
    option->SetDescription(" = 0 / 1 / 2, 0 = not time-dependent, 1 = asymmetric , 2 = symmetric  ");
    // compute the Euclidean length of the diffeomorphism and write to an image -- OutputNamethick.nii.gz -- This is a
    // Beta version of thickness computation -- Not Full-Blown DiReCT , Das, 2009, Neuroimage ---  syn with time is the
    // only model that can be used with this option " );
    this->m_Parser->AddOption(option);
  }

  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("go-faster");
    option->SetShortName('G');
    option->SetDescription(" true / false -- if true, SyN is faster but loses some accuracy wrt inverse-identity "
                           "constraint, see Avants MIA 2008.");
    std::string nitdefault = std::string("false");
    option->AddFunction(nitdefault);
    this->m_Parser->AddOption(option);
  }

  // added by songgang
  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("continue-affine");
    option->SetDescription("true (default) | false, do (not) perform affine given the initial affine parameters");
    std::string nitdefault = std::string("true");
    option->AddFunction(nitdefault);
    this->m_Parser->AddOption(option);
  }

  // added by songgang

  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("number-of-affine-iterations");
    option->SetDescription("number of iterations per level -- a 'vector' e.g.  :  100x100x20 ");
    std::string nitdefault = std::string("10000x10000x10000");
    option->AddFunction(nitdefault);
    this->m_Parser->AddOption(option);
  }
  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("use-NN");
    option->SetDescription("use nearest neighbor interpolation ");
    std::string nitdefault = std::string("0");
    option->AddFunction(nitdefault);
    this->m_Parser->AddOption(option);
  }
  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("use-Histogram-Matching");
    option->SetDescription("use histogram matching of moving to fixed image ");
    std::string nitdefault = std::string("0");
    option->AddFunction(nitdefault);
    this->m_Parser->AddOption(option);
  }

  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("affine-metric-type");
    option->SetDescription(
      "MI: mutual information (default), MSQ: mean square error, SSD, CC: Normalized correlation, CCH: Histogram-based "
      "correlation coefficient (not recommended), GD: gradient difference (not recommended) ");
    std::string nitdefault = std::string("MI");
    option->AddFunction(nitdefault);
    this->m_Parser->AddOption(option);
  }

  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("MI-option");
    option->SetDescription("option of mutual information: MI_bins x MI_samples (default: 32x32000)");

    switch (TDimension)
    {
      case 2:
      {
        option->AddFunction(std::string("32x5000"));
      }
      break;
      case 3:
      {
        option->AddFunction(std::string("32x32000"));
      }
      break;
    }
    // std::string nitdefault=std::string("32x8000");
    // option->AddFunction(nitdefault);
    this->m_Parser->AddOption(option);
  }

  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("rigid-affine");
    option->SetDescription("use rigid transformation : true / false(default)");
    std::string nitdefault = std::string("false");
    option->AddFunction(nitdefault);
    this->m_Parser->AddOption(option);
  }
  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("do-rigid");
    option->SetDescription("use rigid transformation : true / false(default)");
    std::string nitdefault = std::string("false");
    option->AddFunction(nitdefault);
    this->m_Parser->AddOption(option);
  }

  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("affine-gradient-descent-option");
    option->SetDescription("option of gradient descent in affine transformation:  maximum_step_length x "
                           "relaxation_factor x minimum_step_length x translation_scales ");
    std::string nitdefault = std::string("0.1x0.5x1.e-4x1.e-4");
    option->AddFunction(nitdefault);
    this->m_Parser->AddOption(option);
  }

  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("use-rotation-header");
    option->SetDescription("use rotation matrix in image headers: true (default) / false");
    std::string nitdefault = std::string("false");
    option->AddFunction(nitdefault);
    this->m_Parser->AddOption(option);
  }

  if (true)
  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("ignore-void-origin");
    option->SetDescription("ignore the apparently unmatched origins (when use-rotation-header is false and the "
                           "rotation matrix is identity: true (default) / false");
    std::string nitdefault = std::string("false");
    option->AddFunction(nitdefault);
    this->m_Parser->AddOption(option);
  }

  if (true)
  {
    std::string description = std::string("At each resolution level the image is subsampled ") +
                              std::string("and smoothed by Gaussian convolution. This option ") +
                              std::string("allows the user to override the default smoothing ") +
                              std::string("by specifying sigma values (in mm) for smoothing ") +
                              std::string("both fixed and moving images for each resolution ") + std::string("level.");

    std::string defaultFunction = std::string("");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("gaussian-smoothing-sigmas");
    option->SetDescription(description);
    option->AddFunction(defaultFunction);
    this->m_Parser->AddOption(option);
  }

  if (true)
  {
    std::string description = std::string("At each resolution level the image is subsampled ") +
                              std::string("and smoothed by Gaussian convolution. This option ") +
                              std::string("allows the user to override the default subsampling ") +
                              std::string("by specifying the subsampling factor for ") +
                              std::string("both fixed and moving images for each resolution ") + std::string("level.");

    std::string defaultFunction = std::string("");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("subsampling-factors");
    option->SetDescription(description);
    option->AddFunction(defaultFunction);
    this->m_Parser->AddOption(option);
  }
}
} // end namespace itk
#endif
