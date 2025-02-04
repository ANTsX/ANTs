#ifndef antsDisplacementAndVelocityFieldRegistrationCommandIterationUpdate__h_
#define antsDisplacementAndVelocityFieldRegistrationCommandIterationUpdate__h_

#include "itkImageDuplicator.h"

namespace ants
{
/*
There are two types of registration that do not use generic "itkImageRegistrationMethodv4" filter and generic
optimization structures:
- DisplacementFieldRegistrationType
  including:
    * SyN registration
    * BSplineSyN registration

- VelocityFieldRegistrationType
  including:
    * TimeVaryingVelocityFeild
    * TimeVaryingBSplineVelocityField
As these registration types have their own specific optimization processes, a different observer is needed to watch
their internal optimization procedure.
*/
template <typename TFilter>
class antsDisplacementAndVelocityFieldRegistrationCommandIterationUpdate final : public itk::Command
{
public:
  typedef antsDisplacementAndVelocityFieldRegistrationCommandIterationUpdate Self;
  typedef itk::Command                                                       Superclass;
  typedef itk::SmartPointer<Self>                                            Pointer;
  itkNewMacro(Self);

  typedef typename TFilter::FixedImageType  FixedImageType;
  typedef typename TFilter::MovingImageType MovingImageType;

  /** ImageDimension constants */
  static constexpr unsigned int VImageDimension = FixedImageType::ImageDimension;

  typedef typename TFilter::OutputTransformType                                                OutputTransformType;
  typedef typename TFilter::OutputTransformType::ScalarType                                    RealType;
  typedef itk::ImageToImageMetricv4<FixedImageType, MovingImageType, FixedImageType, RealType> MetricType;
  typedef typename MetricType::MeasureType                                                     MeasureType;
  typedef typename MetricType::VirtualImageType                                                VirtualImageType;
  typedef itk::CompositeTransform<RealType, VImageDimension>                                   CompositeTransformType;
  typedef typename CompositeTransformType::TransformType                                       TransformBaseType;
  typedef itk::DisplacementFieldTransform<RealType, VImageDimension>     DisplacementFieldTransformType;
  typedef typename DisplacementFieldTransformType::DisplacementFieldType DisplacementFieldType;
  typedef typename DisplacementFieldType::PixelType                      DisplacementVectorType;
  typedef itk::ImageDuplicator<DisplacementFieldType>                    DisplacementFieldDuplicatorType;

protected:
  antsDisplacementAndVelocityFieldRegistrationCommandIterationUpdate()
  {
    m_clock.Start();
    m_clock.Stop();
    const itk::RealTimeClock::TimeStampType now = m_clock.GetTotal();
    this->m_lastTotalTime = now;
    m_clock.Start();
    this->m_LogStream = &std::cout;
    this->m_ComputeFullScaleCCInterval = 0;
    this->m_WriteIterationsOutputsInIntervals = 0;
    this->m_CurrentStageNumber = 0;
  }

public:
  void
  Execute(itk::Object * caller, const itk::EventObject & event) final
  {
    Execute((const itk::Object *)caller, event);
  }

  void
  Execute(const itk::Object * object, const itk::EventObject & event) final
  {
    TFilter const * const filter = dynamic_cast<const TFilter *>(object);

    if (typeid(event) == typeid(itk::InitializeEvent))
    {
      const unsigned int currentLevel = filter->GetCurrentLevel();

      typename TFilter::ShrinkFactorsPerDimensionContainerType shrinkFactors =
        filter->GetShrinkFactorsPerDimension(currentLevel);
      typename TFilter::SmoothingSigmasArrayType                 smoothingSigmas = filter->GetSmoothingSigmasPerLevel();
      typename TFilter::TransformParametersAdaptorsContainerType adaptors =
        filter->GetTransformParametersAdaptorsPerLevel();
      bool smoothingSigmasAreInPhysicalUnits = filter->GetSmoothingSigmasAreSpecifiedInPhysicalUnits();

      m_clock.Stop();
      const itk::RealTimeClock::TimeStampType now = m_clock.GetTotal();
      this->Logger() << "  Current level = " << currentLevel + 1 << " of " << this->m_NumberOfIterations.size()
                     << std::endl;
      this->Logger() << "    number of iterations = " << this->m_NumberOfIterations[currentLevel] << std::endl;
      this->Logger() << "    shrink factors = " << shrinkFactors << std::endl;
      this->Logger() << "    smoothing sigmas = " << smoothingSigmas[currentLevel];
      if (smoothingSigmasAreInPhysicalUnits)
      {
        this->Logger() << " mm" << std::endl;
      }
      else
      {
        this->Logger() << " vox" << std::endl;
      }
      this->Logger() << "    required fixed parameters = " << adaptors[currentLevel]->GetRequiredFixedParameters()
                     << std::flush << std::endl;
      // this->Logger() << "\n  LEVEL_TIME_INDEX: " << now << " SINCE_LAST: " << (now-this->m_lastTotalTime) <<
      // std::endl;
      this->m_lastTotalTime = now;
      m_clock.Start();

      typedef itk::GradientDescentOptimizerv4Template<RealType> GradientDescentOptimizerType;
      GradientDescentOptimizerType *                            optimizer =
        reinterpret_cast<GradientDescentOptimizerType *>(const_cast<TFilter *>(filter)->GetModifiableOptimizer());

      // TODO:  This looks very wrong.  There is a const_cast above, and then the change
      //       of the number of iterations here on what should be a const object.
      optimizer->SetNumberOfIterations(this->m_NumberOfIterations[currentLevel]);
    }
    else if (typeid(event) == typeid(itk::IterationEvent))
    {
      const unsigned int currentLevel = filter->GetCurrentLevel();
      const unsigned int lCurrentIteration = filter->GetCurrentIteration();
      if (lCurrentIteration == 1)
      {
        if (this->m_ComputeFullScaleCCInterval != 0)
        {
          // Print header line one time
          this->Logger() << "XXDIAGNOSTIC,Iteration,metricValue,convergenceValue,ITERATION_TIME_INDEX,SINCE_LAST,"
                            "FullScaleCCInterval="
                         << this->m_ComputeFullScaleCCInterval << std::flush << std::endl;
        }
        else
        {
          this->Logger() << "XXDIAGNOSTIC,Iteration,metricValue,convergenceValue,ITERATION_TIME_INDEX,SINCE_LAST"
                         << std::endl;
        }
      }
      m_clock.Stop();
      const itk::RealTimeClock::TimeStampType now = m_clock.GetTotal();

      MeasureType        metricValue = 0.0;
      const unsigned int lastIteration = this->m_NumberOfIterations[currentLevel];
      if ((this->m_ComputeFullScaleCCInterval != 0) &&
          (lCurrentIteration == 1 || (lCurrentIteration % this->m_ComputeFullScaleCCInterval == 0) ||
           lCurrentIteration == lastIteration))
      {
        // This function finds the similarity value between the original fixed image and the original moving images
        // using a CC metric type with radius 4.
        // The feature can be used to observe the progress of the registration process at each iteration.
        this->UpdateFullScaleMetricValue(filter, metricValue);
      }

      if ((this->m_WriteIterationsOutputsInIntervals != 0) &&
          (lCurrentIteration == 1 || (lCurrentIteration % this->m_WriteIterationsOutputsInIntervals == 0) ||
           lCurrentIteration == lastIteration))
      {
        // This function writes the output volume of each iteration to the disk.
        // The feature can be used to observe the progress of the registration process at each iteration,
        // and make a short movie from the the registration process.
        this->WriteIntervalVolumes(filter);
      }
      else
      {
        this->Logger() << " "; // if the output of current iteration is written to disk, and star
      }                        // will appear before line, else a free space will be printed to keep visual alignment.

      std::streamsize ss = std::cout.precision();

      this->Logger() << "1DIAGNOSTIC, " << std::setw(5) << lCurrentIteration << ", " << std::scientific
                     << std::setprecision(12) << filter->GetCurrentMetricValue() << ", " << std::scientific
                     << std::setprecision(12) << filter->GetCurrentConvergenceValue() << ", " << std::setprecision(4)
                     << now << ", " << std::setprecision(4) << (now - this->m_lastTotalTime) << ", ";
      if ((this->m_ComputeFullScaleCCInterval != 0) && fabs(metricValue) > 1e-7)
      {
        this->Logger() << std::scientific << std::setprecision(12) << metricValue << std::flush << std::endl;
      }
      else
      {
        this->Logger() << std::endl;
      }

      this->Logger() << std::setprecision(ss);
      this->Logger().unsetf(std::ios::fixed | std::ios::scientific);

      this->m_lastTotalTime = now;
      m_clock.Start();
    }
    else
    {
      // Invalid event type
      return;
    }
  }

  itkSetMacro(ComputeFullScaleCCInterval, unsigned int);

  itkSetMacro(WriteIterationsOutputsInIntervals, unsigned int);

  itkSetMacro(CurrentStageNumber, unsigned int);

  void
  SetOutputPrefix(const std::string & outputPrefix)
  {
    this->m_OutputPrefix = outputPrefix;
  }

  void
  SetNumberOfIterations(const std::vector<unsigned int> & iterations)
  {
    this->m_NumberOfIterations = iterations;
  }

  void
  SetLogStream(std::ostream & logStream)
  {
    this->m_LogStream = &logStream;
  }

  void
  SetOrigFixedImage(typename FixedImageType::Pointer origFixedImage)
  {
    this->m_origFixedImage = origFixedImage;
  }

  void
  SetOrigMovingImage(typename MovingImageType::Pointer origMovingImage)
  {
    this->m_origMovingImage = origMovingImage;
  }

  void
  UpdateFullScaleMetricValue(const TFilter * const filter, MeasureType & metricValue) const
  {
    // Get the registration metric from the filter, input metric is needed to find the type of input transform.
    typename MetricType::ConstPointer inputMetric(dynamic_cast<MetricType const *>(filter->GetMetric()));

    // //////////////////////////////////Define the CC Metric Type to Compute Similarity
    // Measure////////////////////////////
    // This metric type is used to measure the general similarity metric between the original input fixed and moving
    // images.
    typename MetricType::Pointer metric;
    typedef itk::
      ANTSNeighborhoodCorrelationImageToImageMetricv4<FixedImageType, MovingImageType, FixedImageType, MeasureType>
                                            CorrelationMetricType;
    typename CorrelationMetricType::Pointer correlationMetric = CorrelationMetricType::New();
    {
      typename CorrelationMetricType::RadiusType radius;
      radius.Fill(4); // NOTE: This is just a common reference for fine-tuning parameters, so perhaps a smaller window
                      // would be sufficient.
      correlationMetric->SetRadius(radius);
    }
    correlationMetric->SetUseMovingImageGradientFilter(false);
    correlationMetric->SetUseFixedImageGradientFilter(false);
    metric = correlationMetric;

    /*
     Below, the implementation is just provided for SyN registration filter.
     TODO: expand the similarity metric implementation for other registration types mentioned above.
     */
    if (strcmp(inputMetric->GetMovingTransform()->GetNameOfClass(), "DisplacementFieldTransform") == 0)
    {
      /*
      Filter returns the SyN internal trnasforms (MovingToMiddleTransform & FixedToMiddleTransform) at each iteration.
      These transforms are used to generate input transforms of full scale CC metric. NOTICE: Using const_cast for
      filter does not make any issue because the requested outputs are copied to another objects, and there will be no
      change to them at future.
      */
      // Copy the SyN internal transforms at each iteration
      typename DisplacementFieldTransformType::Pointer myFixedToMiddleTransform = DisplacementFieldTransformType::New();
      typename DisplacementFieldTransformType::Pointer myMovingToMiddleTransform =
        DisplacementFieldTransformType::New();

      // copy FixedToMiddleTransform
      typename DisplacementFieldDuplicatorType::Pointer FixedDisplacementDuplicator =
        DisplacementFieldDuplicatorType::New();
      FixedDisplacementDuplicator->SetInputImage(
        const_cast<DisplacementFieldTransformType *>(filter->GetFixedToMiddleTransform())->GetDisplacementField());
      FixedDisplacementDuplicator->Update();
      typename DisplacementFieldDuplicatorType::Pointer FixedInverseDisplacementDuplicator =
        DisplacementFieldDuplicatorType::New();
      FixedInverseDisplacementDuplicator->SetInputImage(
        const_cast<DisplacementFieldTransformType *>(filter->GetFixedToMiddleTransform())
          ->GetInverseDisplacementField());
      FixedInverseDisplacementDuplicator->Update();

      myFixedToMiddleTransform->SetDisplacementField(FixedDisplacementDuplicator->GetOutput());
      myFixedToMiddleTransform->SetInverseDisplacementField(FixedInverseDisplacementDuplicator->GetOutput());

      // copy MovingToMiddleTransform
      typename DisplacementFieldDuplicatorType::Pointer MovingDisplacementDuplicator =
        DisplacementFieldDuplicatorType::New();
      MovingDisplacementDuplicator->SetInputImage(
        const_cast<DisplacementFieldTransformType *>(filter->GetMovingToMiddleTransform())->GetDisplacementField());
      MovingDisplacementDuplicator->Update();
      typename DisplacementFieldDuplicatorType::Pointer MovingInverseDisplacementDuplicator =
        DisplacementFieldDuplicatorType::New();
      MovingInverseDisplacementDuplicator->SetInputImage(
        const_cast<DisplacementFieldTransformType *>(filter->GetMovingToMiddleTransform())
          ->GetInverseDisplacementField());
      MovingInverseDisplacementDuplicator->Update();

      myMovingToMiddleTransform->SetDisplacementField(MovingDisplacementDuplicator->GetOutput());
      myMovingToMiddleTransform->SetInverseDisplacementField(MovingInverseDisplacementDuplicator->GetOutput());

      // Based on SyN Registration implementation, fixed composite and moving composite transforms are generated to
      // compute the metric value at each iteration.
      typedef typename TFilter::InitialTransformType InitialTransformType;

      typename CompositeTransformType::Pointer fixedComposite = CompositeTransformType::New();
      typename CompositeTransformType::Pointer movingComposite = CompositeTransformType::New();

      fixedComposite->AddTransform(const_cast<InitialTransformType *>(filter->GetFixedInitialTransform()));
      fixedComposite->AddTransform(myFixedToMiddleTransform->GetInverseTransform());
      fixedComposite->FlattenTransformQueue();
      fixedComposite->SetOnlyMostRecentTransformToOptimizeOn();

      movingComposite->AddTransform(const_cast<InitialTransformType *>(filter->GetMovingInitialTransform()));
      movingComposite->AddTransform(myMovingToMiddleTransform->GetInverseTransform());
      movingComposite->FlattenTransformQueue();
      movingComposite->SetOnlyMostRecentTransformToOptimizeOn();

      // SyN uses the above composite transforms to compute the current metric value in two ways as follows:
      /*
       At the first method, the input images are downsampled by the fixed and moving transforms,
       and then, the output of resamplers are passed to the CC similarity metric with identity transforms.
      */
      if (filter->GetDownsampleImagesForMetricDerivatives())
      {
        typedef itk::ResampleImageFilter<FixedImageType, FixedImageType, RealType> FixedResamplerType;
        typename FixedResamplerType::Pointer fixedResampler = FixedResamplerType::New();
        fixedResampler->SetTransform(fixedComposite);
        fixedResampler->SetInput(this->m_origFixedImage);
        fixedResampler->SetOutputParametersFromImage(this->m_origFixedImage);
        fixedResampler->SetDefaultPixelValue(0);
        fixedResampler->Update();

        typedef itk::ResampleImageFilter<MovingImageType, MovingImageType, RealType> MovingResamplerType;
        typename MovingResamplerType::Pointer movingResampler = MovingResamplerType::New();
        movingResampler->SetTransform(movingComposite);
        movingResampler->SetInput(this->m_origMovingImage);
        movingResampler->SetOutputParametersFromImage(this->m_origFixedImage);
        movingResampler->SetDefaultPixelValue(0);
        movingResampler->Update();

        typedef typename itk::IdentityTransform<RealType, VImageDimension> IdentityTransformType;
        typename IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();

        typename DisplacementFieldType::Pointer identityField = DisplacementFieldType::New();
        identityField->CopyInformation(this->m_origFixedImage);
        identityField->SetRegions(this->m_origFixedImage->GetRequestedRegion());
        identityField->AllocateInitialized();

        typename DisplacementFieldTransformType::Pointer identityDisplacementFieldTransform =
          DisplacementFieldTransformType::New();
        identityDisplacementFieldTransform->SetDisplacementField(identityField);

        metric->SetFixedImage(fixedResampler->GetOutput());
        metric->SetFixedTransform(identityTransform);
        metric->SetMovingImage(movingResampler->GetOutput());
        metric->SetMovingTransform(identityDisplacementFieldTransform);
      }
      /*
       At the second method, the computed fixed and moving composite transforms are passed to the CC similarity metric
       directly with the full scale fixed and moving images.
      */
      else if (!(filter->GetDownsampleImagesForMetricDerivatives()))
      {
        metric->SetFixedImage(this->m_origFixedImage);
        metric->SetFixedTransform(fixedComposite);
        metric->SetMovingImage(this->m_origMovingImage);
        metric->SetMovingTransform(movingComposite);
      }
    }
    metric->SetVirtualDomainFromImage(this->m_origFixedImage);
    metric->Initialize();
    metricValue = metric->GetValue();
  }

  void
  WriteIntervalVolumes(TFilter const * const filter) const
  {
    // //////////////////////////
    // Get output transform from the registration filter at each iteration
    // It can be useful in some cases e.g. when we want to present the registration progress as a series of consequent
    // pictures.
    typename DisplacementFieldTransformType::Pointer OutputTransformAtCurrentIteration =
      DisplacementFieldTransformType::New();

    // Filter return the MovingToMiddleTransform and FixedToMiddleTransform of each iteration, so they should be
    // composed to generate the final output transform of current iteration

    // Notice that using const_cast for filter does not make any issue, because the inputMovingTransform will never be
    // used in any processing. It is only copied to another transform.

    typedef itk::ComposeDisplacementFieldsImageFilter<DisplacementFieldType, DisplacementFieldType> ComposerType;
    typename ComposerType::Pointer composer = ComposerType::New();
    composer->SetDisplacementField(const_cast<DisplacementFieldTransformType *>(filter->GetMovingToMiddleTransform())
                                     ->GetInverseDisplacementField());
    composer->SetWarpingField(
      const_cast<DisplacementFieldTransformType *>(filter->GetFixedToMiddleTransform())->GetDisplacementField());
    composer->Update();

    typename ComposerType::Pointer inverseComposer = ComposerType::New();
    inverseComposer->SetDisplacementField(
      const_cast<DisplacementFieldTransformType *>(filter->GetFixedToMiddleTransform())->GetInverseDisplacementField());

    inverseComposer->SetWarpingField(
      const_cast<DisplacementFieldTransformType *>(filter->GetMovingToMiddleTransform())->GetDisplacementField());
    inverseComposer->Update();

    OutputTransformAtCurrentIteration->SetDisplacementField(composer->GetOutput());
    OutputTransformAtCurrentIteration->SetInverseDisplacementField(inverseComposer->GetOutput());

    // Now this output transform is copied to another instance to prevent undesired changes.

    typename DisplacementFieldDuplicatorType::Pointer disDuplicator = DisplacementFieldDuplicatorType::New();
    disDuplicator->SetInputImage(OutputTransformAtCurrentIteration->GetDisplacementField());
    disDuplicator->Update();

    typename DisplacementFieldDuplicatorType::Pointer disInverseDuplicator = DisplacementFieldDuplicatorType::New();
    disInverseDuplicator->SetInputImage(OutputTransformAtCurrentIteration->GetInverseDisplacementField());
    disInverseDuplicator->Update();

    typename DisplacementFieldTransformType::Pointer outputTransformReadyToUse = DisplacementFieldTransformType::New();
    outputTransformReadyToUse->SetDisplacementField(disDuplicator->GetOutput());
    outputTransformReadyToUse->SetInverseDisplacementField(disInverseDuplicator->GetOutput());

    // Now add this updated transform to the composite transform including the initial trnasform
    typedef typename TFilter::InitialTransformType InitialTransformType;

    typename CompositeTransformType::Pointer outputCompositTransform = CompositeTransformType::New();
    if (filter->GetMovingInitialTransform())
    {
      outputCompositTransform->AddTransform(const_cast<InitialTransformType *>(filter->GetMovingInitialTransform()));
    }
    outputCompositTransform->AddTransform(outputTransformReadyToUse);
    outputCompositTransform->FlattenTransformQueue();
    outputCompositTransform->SetOnlyMostRecentTransformToOptimizeOn();

    // Now we use the output transform to get warped image using linear interpolation
    typedef itk::LinearInterpolateImageFunction<MovingImageType, RealType> LinearInterpolatorType;
    typename LinearInterpolatorType::Pointer linearInterpolator = LinearInterpolatorType::New();

    typedef itk::ResampleImageFilter<FixedImageType, MovingImageType, RealType> ResampleFilterType;
    typename ResampleFilterType::Pointer                                        resampler = ResampleFilterType::New();
    resampler->SetTransform(outputCompositTransform);
    resampler->SetInput(this->m_origMovingImage);
    resampler->SetOutputParametersFromImage(this->m_origFixedImage);
    resampler->SetInterpolator(linearInterpolator);
    resampler->SetDefaultPixelValue(0);
    resampler->Update();

    // write the results to the disk
    const unsigned int curLevel = filter->GetCurrentLevel();
    const unsigned int curIter = filter->GetCurrentIteration();
    std::stringstream  currentFileName;
    currentFileName << this->m_OutputPrefix << "Stage" << this->m_CurrentStageNumber + 1 << "_level" << curLevel + 1;
    /*
     The name arrangement of written files are important to us.
     To prevent: "Iter1 Iter10 Iter2 Iter20" we use the following style.
     Then the order is: "Iter1 Iter2 ... Iters10 ... Itert20"
    */
    if (curIter < 10)
    {
      currentFileName << "_Iter000" << curIter << ".nii.gz";
    }
    else if (curIter < 100)
    {
      currentFileName << "_Iter00" << curIter << ".nii.gz";
    }
    else if (curIter < 1000)
    {
      currentFileName << "_Iter0" << curIter << ".nii.gz";
    }
    else
    {
      currentFileName << "_Iter" << curIter << ".nii.gz";
    }
    std::cout << "*"; // The star befor each DIAGNOSTIC shows that its output is writtent out.
    std::cout << currentFileName.str()
              << std::endl; // The star befor each DIAGNOSTIC shows that its output is writtent out.

    typedef itk::ImageFileWriter<MovingImageType> WarpedImageWriterType;
    typename WarpedImageWriterType::Pointer       writer = WarpedImageWriterType::New();
    writer->SetFileName(currentFileName.str().c_str());
    writer->SetInput(resampler->GetOutput());
    try
    {
      writer->Update();
    }
    catch (const itk::ExceptionObject & err)
    {
      std::cout << "Can't write warped image " << currentFileName.str().c_str() << std::endl;
      std::cout << "Exception Object caught: " << std::endl;
      std::cout << err << std::endl;
    }
  }

private:
  std::ostream &
  Logger() const
  {
    return *m_LogStream;
  }

  /**
   *  WeakPointer to the Optimizer
   */
  // itk::WeakPointer<OptimizerType>   m_Optimizer;

  std::string                       m_OutputPrefix;
  std::vector<unsigned int>         m_NumberOfIterations;
  std::ostream *                    m_LogStream;
  itk::TimeProbe                    m_clock;
  itk::RealTimeClock::TimeStampType m_lastTotalTime;

  unsigned int m_ComputeFullScaleCCInterval;
  unsigned int m_WriteIterationsOutputsInIntervals;
  unsigned int m_CurrentStageNumber;

  typename FixedImageType::Pointer  m_origFixedImage;
  typename MovingImageType::Pointer m_origMovingImage;
};
} // end namespace ants
#endif // antsDisplacementAndVelocityFieldRegistrationCommandIterationUpdate__h_
