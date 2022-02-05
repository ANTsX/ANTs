#ifndef antsRegistrationCommandIterationUpdate__h_
#define antsRegistrationCommandIterationUpdate__h_

namespace ants
{
/** \class antsRegistrationCommandIterationUpdate
 *  \brief change parameters between iterations of registration
 */
template <typename TFilter>
class antsRegistrationCommandIterationUpdate final : public itk::Command
{
public:
  typedef antsRegistrationCommandIterationUpdate Self;
  typedef itk::Command                           Superclass;
  typedef itk::SmartPointer<Self>                Pointer;
  itkNewMacro(Self);

protected:
  antsRegistrationCommandIterationUpdate()
  {
    m_clock.Start();
    m_clock.Stop();
    const itk::RealTimeClock::TimeStampType now = m_clock.GetTotal();
    this->m_lastTotalTime = now;
    m_clock.Start();
    this->m_LogStream = &std::cout;
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

      typedef itk::GradientDescentOptimizerv4Template<typename TFilter::RealType> GradientDescentOptimizerType;
      GradientDescentOptimizerType *                                              optimizer =
        reinterpret_cast<GradientDescentOptimizerType *>(const_cast<TFilter *>(filter)->GetModifiableOptimizer());

      // TODO:  This looks very wrong.  There is a const_cast above, and then the change
      //       of the number of iterations here on what should be a const object.
      optimizer->SetNumberOfIterations(this->m_NumberOfIterations[currentLevel]);
    }
    else if (typeid(event) == typeid(itk::IterationEvent))
    {
      const unsigned int lCurrentIteration = filter->GetCurrentIteration();
      if (lCurrentIteration == 1)
      {
        // Print header line one time
        this->Logger() << "XDIAGNOSTIC,Iteration,metricValue,convergenceValue,ITERATION_TIME_INDEX,SINCE_LAST"
                       << std::flush << std::endl;
      }

      m_clock.Stop();
      const itk::RealTimeClock::TimeStampType now = m_clock.GetTotal();
      this->Logger() << "WDIAGNOSTIC, " << std::setw(5) << lCurrentIteration << ", " << std::scientific
                     << std::setprecision(12) << filter->GetCurrentMetricValue() << ", " << std::scientific
                     << std::setprecision(12) << filter->GetCurrentConvergenceValue() << ", " << std::setprecision(4)
                     << now << ", " << std::setprecision(4) << (now - this->m_lastTotalTime) << ", " << std::flush
                     << std::endl;
      this->m_lastTotalTime = now;
      m_clock.Start();
    }
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

private:
  std::ostream &
  Logger() const
  {
    return *m_LogStream;
  }

  std::vector<unsigned int>         m_NumberOfIterations;
  std::ostream *                    m_LogStream;
  itk::TimeProbe                    m_clock;
  itk::RealTimeClock::TimeStampType m_lastTotalTime;

  // typename ImageType::Pointer m_origFixedImage;
  // typename ImageType::Pointer m_origMovingImage;
};
} // end namespace ants
#endif // antsRegistrationCommandIterationUpdate__h_
