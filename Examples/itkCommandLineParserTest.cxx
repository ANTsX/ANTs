#include "itkPICSLAdvancedNormalizationToolKit.h"
#include "itkCommandLineParser.h"

int main( int argc, char *argv[] )
{
  typedef itk::PICSLAdvancedNormalizationToolKit<3> RegistrationType;
  RegistrationType::Pointer registration = RegistrationType::New();

  registration->ParseCommandLine( argc, argv );

  typedef itk::CommandLineParser ParserType;
  ParserType::Pointer parser = ParserType::New();

  typedef ParserType::OptionType OptionType;
  OptionType::Pointer fileOption = OptionType::New();
  fileOption->SetShortName( 'f' );
  fileOption->SetLongName( "file" );
  fileOption->SetDescription( "The fixed image file used in the image registration algorithm." );

  parser->AddOption( fileOption );

  typedef ParserType::OptionType OptionType;
  OptionType::Pointer metricOption = OptionType::New();
  metricOption->SetShortName( 'm' );
  metricOption->SetLongName( "metric" );
  metricOption->SetDescription( "The metric used by the image registration algorithm." );

  parser->AddOption( metricOption );

  typedef ParserType::OptionType OptionType;
  OptionType::Pointer verbosityOption = OptionType::New();
  verbosityOption->SetLongName( "verbosity" );
  verbosityOption->SetDescription( "Mundi vult decipi.  " );

  parser->AddOption( verbosityOption );

  typedef ParserType::OptionType OptionType;
  OptionType::Pointer dirOption = OptionType::New();
  dirOption->SetLongName( "directionality" );
  //  dirOption->SetShortName( 'd' );
  dirOption->SetDescription( "Mundi vult decipi.  " );

  parser->AddOption( dirOption );

  parser->Parse( argc, argv );

  if( !parser->GetOption( "directionality" ) )
    {
    std::cout << "N ULL" << std::endl;
    }
  else
    {
    for( unsigned int j = 0; j < 3; j++ )
      {
      std::cout << parser->ConvertVector<bool>(
        parser->GetOption( "directionality" )->GetValue( 0 ) )[j] << " x ";
      }
    }

  parser->PrintMenu( std::cout, 7 );

  //  std::cout << std::endl << std::endl << "--------------------------------------------"
  //            << std::endl << std::endl;

  //  parser->Print( std::cout << 7 );

  return 0;
};
