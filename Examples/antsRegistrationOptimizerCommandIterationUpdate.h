#ifndef antsRegistrationOptimizerCommandIterationUpdate__h_
#define antsRegistrationOptimizerCommandIterationUpdate__h_

namespace ants
{
/** \class antsRegistrationOptimizerCommandIterationUpdate
 *  \brief observe the optimizer for traditional registration methods
 */
template <class ParametersValueType, unsigned VImageDimension, class TOptimizer>
class antsRegistrationOptimizerCommandIterationUpdate : public itk::Command
{
public:
  typedef antsRegistrationOptimizerCommandIterationUpdate Self;
  typedef itk::Command                                    Superclass;
  typedef itk::SmartPointer<Self>                         Pointer;
  itkNewMacro( Self );

  typedef ParametersValueType     RealType;
  typedef ParametersValueType     PixelType;
  typedef typename itk::Image<PixelType, VImageDimension>              ImageType;
  typedef itk::ImageToImageMetricv4
                <ImageType, ImageType, ImageType, RealType>            MetricType;
  typedef typename MetricType::MeasureType                             MeasureType;
  typedef itk::CompositeTransform<RealType, VImageDimension>           CompositeTransformType;
  typedef typename CompositeTransformType::TransformType               TransformBaseType;
protected:
  antsRegistrationOptimizerCommandIterationUpdate()
  {
    m_clock.Start();
    m_clock.Stop();
    const itk::RealTimeClock::TimeStampType now = m_clock.GetTotal();
    this->m_lastTotalTime = now;
    m_clock.Start();
    this->m_LogStream = &::ants::antscout;
    this->m_origFixedImage = ImageType::New();
    this->m_origMovingImage = ImageType::New();
    this->m_ComputeFullScaleCCInterval = 0;
    this->m_WriteInterationsOutputsInIntervals = 0;
    this->m_CurrentStageNumber = 0;
    this->m_CurLevel = 0;
  }

public:

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *) caller, event);
  }

  void Execute(const itk::Object *, const itk::EventObject & event)
  {
#if 0
    if( typeid( event ) == typeid( itk::InitializeEvent ) )
      {
      const unsigned int currentLevel = this->m_Optimizer->GetCurrentLevel();

      typename TOptimizer::ShrinkFactorsPerDimensionContainerType shrinkFactors = this->m_Optimizer->GetShrinkFactorsPerDimension( currentLevel );
      typename TOptimizer::SmoothingSigmasArrayType smoothingSigmas = this->m_Optimizer->GetSmoothingSigmasPerLevel();
      typename TOptimizer::TransformParametersAdaptorsContainerType adaptors =
        this->m_Optimizer->GetTransformParametersAdaptorsPerLevel();
      bool smoothingSigmasAreInPhysicalUnits = this->m_Optimizer->GetSmoothingSigmasAreSpecifiedInPhysicalUnits();

      m_clock.Stop();
      const itk::RealTimeClock::TimeStampType now = m_clock.GetTotal();
      this->Logger() << "  Current level = " << currentLevel + 1 << " of " << this->m_NumberOfIterations.size()
                     << std::endl;
      this->Logger() << "    number of iterations = " << this->m_NumberOfIterations[currentLevel] << std::endl;
      this->Logger() << "    shrink factors = " << shrinkFactors << std::endl;
      this->Logger() << "    smoothing sigmas = " << smoothingSigmas[currentLevel];
      if( smoothingSigmasAreInPhysicalUnits )
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

      typedef itk::GradientDescentOptimizerv4<ParametersValueType> GradientDescentOptimizerType;
      GradientDescentOptimizerType * optimizer = reinterpret_cast<GradientDescentOptimizerType *>(
          const_cast<typename TOptimizer::OptimizerType *>( const_cast<TOptimizer *>( this->m_Optimizer )->GetOptimizer() ) );

      // TODO:  This looks very wrong.  There is a const_cast above, and then the change
      //       of the number of iterations here on what should be a const object.
      optimizer->SetNumberOfIterations( this->m_NumberOfIterations[currentLevel] );
      }
    else
#endif
    if( typeid( event ) == typeid( itk::IterationEvent ) )
      {
      const unsigned int lCurrentIteration = this->m_Optimizer->GetCurrentIteration() + 1;
      if( lCurrentIteration  == 1 )
        {
        if( this->m_ComputeFullScaleCCInterval != 0 )
          {
          // Print header line one time
          this->Logger()
            << "DIAGNOSTIC,Iteration,metricValue,convergenceValue,ITERATION_TIME_INDEX,SINCE_LAST,FullScaleCCInterval="
            << this->m_ComputeFullScaleCCInterval << std::flush << std::endl;
          }
        else
          {
          this->Logger() << "DIAGNOSTIC,Iteration,metricValue,convergenceValue,ITERATION_TIME_INDEX,SINCE_LAST"
                         << std::flush << std::endl;
          }
        }
      m_clock.Stop();
      const itk::RealTimeClock::TimeStampType now = m_clock.GetTotal();

      MeasureType        metricValue = 0.0;
      const unsigned int lastIteration = this->m_Optimizer->GetNumberOfIterations();
      if( ( this->m_ComputeFullScaleCCInterval != 0 ) &&
          ( lCurrentIteration == 1 || ( lCurrentIteration % this->m_ComputeFullScaleCCInterval == 0 ) ||
            lCurrentIteration == lastIteration) )
        {
        // This function finds the similarity value between the original fixed image and the original moving images
        // using a CC metric type with radius 4.
        // The feature can be used to observe the progress of the registration process at each iteration.
        this->UpdateFullScaleMetricValue(this->m_Optimizer, metricValue);
        }

      if( ( this->m_WriteInterationsOutputsInIntervals != 0 ) &&
          ( lCurrentIteration == 1 || (lCurrentIteration % this->m_WriteInterationsOutputsInIntervals == 0 ) ||
         lCurrentIteration == lastIteration) )
        {
        // This function writes the output volume of each iteration to the disk.
        // The feature can be used to observe the progress of the registration process at each iteration,
        // and make a short movie from the the registration process.
        this->WriteIntervalVolumes(this->m_Optimizer);
        }
      else
        {
        std::cout << " "; // if the output of current iteration is written to disk, and star
        }                 // will appear before line, else a free space will be printed to keep visual alignment.

      this->Logger() << "DIAGNOSTIC, "
                     << std::setw(5) << lCurrentIteration << ", "
                     << std::scientific << std::setprecision(12) << this->m_Optimizer->GetValue() << ", "
                     << std::scientific << std::setprecision(12) << this->m_Optimizer->GetConvergenceValue() << ", "
                     << std::setprecision(4) << now << ", "
                     << std::setprecision(4) << (now - this->m_lastTotalTime)  << ", ";
      if( ( this->m_ComputeFullScaleCCInterval != 0 ) &&  fabs(metricValue) > 1e-7 )
        {
        this->Logger() << std::scientific << std::setprecision(12) << metricValue
                       << std::flush << std::endl;
        }
      else
        {
        this->Logger() << std::flush << std::endl;
        }

      this->m_lastTotalTime = now;
      m_clock.Start();
      }
    else
      {
      // Unknown event type
      return;
      }
  }

  itkSetMacro( ComputeFullScaleCCInterval, unsigned int );

  itkSetMacro( WriteInterationsOutputsInIntervals, unsigned int );

  itkSetMacro( CurrentStageNumber, unsigned int );

  void SetNumberOfIterations( const std::vector<unsigned int> & iterations )
  {
    this->m_NumberOfIterations = iterations;
  }

  void SetLogStream(std::ostream & logStream)
  {
    this->m_LogStream = &logStream;
  }

  /**
   * Type defining the optimizer
   */
  typedef    TOptimizer OptimizerType;

  /**
   * Set Optimizer
   */
  void SetOptimizer( OptimizerType * optimizer )
  {
    this->m_Optimizer = optimizer;
    this->m_Optimizer->AddObserver( itk::IterationEvent(), this );
  }

  void SetOrigFixedImage(typename ImageType::Pointer origFixedImage)
  {
    this->m_origFixedImage = origFixedImage;
  }

  void SetOrigMovingImage(typename ImageType::Pointer origMovingImage)
  {
    this->m_origMovingImage = origMovingImage;
  }

  void UpdateFullScaleMetricValue(itk::WeakPointer<OptimizerType> myOptimizer,
                                  MeasureType & metricValue ) const
  {
    // Get the registration metric from the optimizer
    typename MetricType::ConstPointer inputMetric( dynamic_cast<MetricType const *>( myOptimizer->GetMetric() ) );

    // Define the CC metric type
    // This metric type is used to measure the general similarity metric between the original input fixed and moving
    // images.
    typedef itk::ANTSNeighborhoodCorrelationImageToImageMetricv4<ImageType, ImageType, ImageType, MeasureType> CorrelationMetricType;
    typename CorrelationMetricType::Pointer correlationMetric = CorrelationMetricType::New();
      {
      typename CorrelationMetricType::RadiusType radius;
      radius.Fill( 4 );  // NOTE: This is just a common reference for fine-tuning parameters, so perhaps a smaller
                         // window would be sufficient.
      correlationMetric->SetRadius( radius );
      }
    correlationMetric->SetUseMovingImageGradientFilter( false );
    correlationMetric->SetUseFixedImageGradientFilter( false );
    typename MetricType::Pointer metric = correlationMetric.GetPointer();

    // We need to create an exact copy from the composite fixed and moving transforms returned from the metric
    // We should roll off the composite transform and create a new instance from each of its sub transforms

    // For the fixed transform, first we should check that wether it is an identity transform or composite transform.
    typename TransformBaseType::Pointer fixedTransform;
    if( strcmp( inputMetric->GetFixedTransform()->GetNameOfClass(), "CompositeTransform" ) == 0 )
      {
      typename CompositeTransformType::Pointer myFixedTransform = CompositeTransformType::New();

      // Const_cast just makes it possible to cast the metric's transform to a composite transform, so we can copy each
      // of its sub transforms to a new instance.
      // Notice that the metric transform will not be changed inside this fuction.
      // NOTE:  This will not be needed when ITKv4 is updated to include const versions of all get functions.
      // TODO: Remove const_cast once ITKv4 is fixed to allow const Get Macro functions.
      typedef typename MetricType::FixedTransformType FixedTransformType;
      typename CompositeTransformType::ConstPointer inputFixedTransform =
        dynamic_cast<CompositeTransformType *>( const_cast<FixedTransformType *>( inputMetric->GetFixedTransform() ) );
      const unsigned int N = inputFixedTransform->GetNumberOfTransforms();
      for( unsigned int i = 0; i < N; i++ )
        {
        // Create a new instance from each sub transform.
        typename TransformBaseType::Pointer subTransform(
          dynamic_cast<TransformBaseType *>( inputFixedTransform->GetNthTransform(i)->CreateAnother().GetPointer() ) );
        // Copy the information to each sub transform and add this transform to the final composite transform.
        const typename TransformBaseType::ParametersType & fixedImage_paras =
          inputFixedTransform->GetNthTransform(i)->GetParameters();
        const typename TransformBaseType::ParametersType & fixedImage_fixed_paras =
          inputFixedTransform->GetNthTransform(i)->GetFixedParameters();
        subTransform->SetParameters( fixedImage_paras );
        subTransform->SetFixedParameters( fixedImage_fixed_paras );
        myFixedTransform->AddTransform( subTransform );
        }
      myFixedTransform->SetOnlyMostRecentTransformToOptimizeOn();
      fixedTransform = myFixedTransform;
      }
    else if( strcmp( inputMetric->GetFixedTransform()->GetNameOfClass(), "IdentityTransform" ) == 0 )
      {
      typedef typename itk::IdentityTransform<RealType, VImageDimension>    IdentityTransformType;
      typename IdentityTransformType::Pointer myFixedTransform = IdentityTransformType::New();
      fixedTransform = myFixedTransform;
      }
    else
      {
      itkExceptionMacro( "Fixed Transform should be either \"Composite\" or \"Identity\" transform." );
      }

    // Same procedure for the moving transform. Moving transform is always a Composite transform.
    typedef typename MetricType::MovingTransformType MovingTransformType;
    typename CompositeTransformType::Pointer movingTransform = CompositeTransformType::New();
    // TODO: Remove const_cast once ITKv4 is fixed to allow const Get Macro functions.
    typename CompositeTransformType::ConstPointer inputMovingTransform =
      dynamic_cast<CompositeTransformType *>( const_cast<MovingTransformType *>( inputMetric->GetMovingTransform() ) );
    const unsigned int N = inputMovingTransform->GetNumberOfTransforms();
    for( unsigned int i = 0; i < N; i++ )
      {
      typename TransformBaseType::Pointer subTransform(
        dynamic_cast<TransformBaseType *>( inputMovingTransform->GetNthTransform(i)->CreateAnother().GetPointer() ) );
      const typename TransformBaseType::ParametersType & moving_paras =
        inputMovingTransform->GetNthTransform(i)->GetParameters();
      const typename TransformBaseType::ParametersType & moving_fixed_paras =
        inputMovingTransform->GetNthTransform(i)->GetFixedParameters();
      subTransform->SetParameters( moving_paras );
      subTransform->SetFixedParameters( moving_fixed_paras );
      movingTransform->AddTransform( subTransform );
      }
    movingTransform->SetOnlyMostRecentTransformToOptimizeOn();

    metric->SetVirtualDomainFromImage( this->m_origFixedImage );
    metric->SetFixedImage( this->m_origFixedImage );
    metric->SetFixedTransform( fixedTransform );
    metric->SetMovingImage( this->m_origMovingImage );
    metric->SetMovingTransform( movingTransform );
    metric->Initialize();
    metricValue = metric->GetValue();
  }

  void WriteIntervalVolumes(itk::WeakPointer<OptimizerType> myOptimizer)
  {
    // Get the registration metric from the optimizer
    typename MetricType::ConstPointer inputMetric( dynamic_cast<MetricType const *>( myOptimizer->GetMetric() ) );

    // First, compute the moving transform
    typedef typename MetricType::MovingTransformType MovingTransformType;
    typename CompositeTransformType::Pointer movingTransform = CompositeTransformType::New();
    // TODO: Remove const_cast once ITKv4 is fixed to allow const Get Macro functions.
    typename CompositeTransformType::ConstPointer inputMovingTransform =
      dynamic_cast<CompositeTransformType *>( const_cast<MovingTransformType *>( inputMetric->GetMovingTransform() ) );
    const unsigned int N = inputMovingTransform->GetNumberOfTransforms();
    for( unsigned int i = 0; i < N; i++ )
      {
      typename TransformBaseType::Pointer subTransform(
        dynamic_cast<TransformBaseType *>( inputMovingTransform->GetNthTransform(i)->CreateAnother().GetPointer() ) );
      const typename TransformBaseType::ParametersType & moving_paras =
        inputMovingTransform->GetNthTransform(i)->GetParameters();
      const typename TransformBaseType::ParametersType & moving_fixed_paras =
        inputMovingTransform->GetNthTransform(i)->GetFixedParameters();
      subTransform->SetParameters( moving_paras );
      subTransform->SetFixedParameters( moving_fixed_paras );
      movingTransform->AddTransform( subTransform );
      }
    movingTransform->SetOnlyMostRecentTransformToOptimizeOn();

    // Now we apply this output transform to get warped image
    typedef itk::LinearInterpolateImageFunction<ImageType, RealType> LinearInterpolatorType;
    typename LinearInterpolatorType::Pointer linearInterpolator = LinearInterpolatorType::New();

    typedef itk::ResampleImageFilter<ImageType, ImageType, RealType> ResampleFilterType;
    typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    resampler->SetTransform( movingTransform );
    resampler->SetInput( this->m_origMovingImage );
    resampler->SetOutputParametersFromImage( this->m_origFixedImage );
    resampler->SetInterpolator( linearInterpolator );
    resampler->SetDefaultPixelValue( 0 );
    resampler->Update();

    // write the results to the disk
    // const unsigned int curLevel = this->m_Optimizer->GetCurrentLevel();
    const unsigned int curIter = this->m_Optimizer->GetCurrentIteration() + 1;
    if( curIter == 1 )
      {
      ++this->m_CurLevel;
      }
    std::stringstream currentFileName;
    currentFileName << "Stage" << this->m_CurrentStageNumber + 1 << "_level" << this->m_CurLevel;
    /*
    The name arrangement of written files are important to us.
    To prevent: "Iter1 Iter10 Iter2 Iter20" we use the following style.
    Then the order is: "Iter1 Iter2 ... Iters10 ... Itert20"
    */
    if( curIter > 9 )
      {
      currentFileName << "_Iters" << curIter << ".nii.gz";
      }
    else if( curIter > 19 )
      {
      currentFileName << "_Itert" << curIter << ".nii.gz";
      }
    else
      {
      currentFileName << "_Iter" << curIter << ".nii.gz";
      }
    std::cout << "*"; // The star befor each DIAGNOSTIC shows that its output is writtent out.

    typedef itk::ImageFileWriter<ImageType> WarpedImageWriterType;
    typename WarpedImageWriterType::Pointer writer = WarpedImageWriterType::New();
    writer->SetFileName( currentFileName.str().c_str() );
    writer->SetInput( resampler->GetOutput() );
    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      antscout << "Can't write warped image " << currentFileName.str().c_str() << std::endl;
      antscout << "Exception Object caught: " << std::endl;
      antscout << err << std::endl;
      }
  }

private:
  std::ostream & Logger() const
  {
    return *m_LogStream;
  }

  /**
   *  WeakPointer to the Optimizer
   */
  itk::WeakPointer<OptimizerType> m_Optimizer;

  std::vector<unsigned int>         m_NumberOfIterations;
  std::ostream *                    m_LogStream;
  itk::TimeProbe                    m_clock;
  itk::RealTimeClock::TimeStampType m_lastTotalTime;

  unsigned int m_ComputeFullScaleCCInterval;
  unsigned int m_WriteInterationsOutputsInIntervals;
  unsigned int m_CurrentStageNumber;
  unsigned int m_CurLevel;

  typename ImageType::Pointer       m_origFixedImage;
  typename ImageType::Pointer       m_origMovingImage;
};
}; // end namespace ants
#endif // antsRegistrationOptimizerCommandIterationUpdate__h_
