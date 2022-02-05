
#include "antsUtilities.h"
#include <algorithm>

#include "itkPICSLAdvancedNormalizationToolKit.h"
#include "itkCommandLineParser.h"

namespace ants
{
// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
itkCommandLineParserTest(std::vector<std::string> args, std::ostream * out_stream = nullptr)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "itkCommandLineParserTest");

  int     argc = args.size();
  char ** argv = new char *[args.size() + 1];
  for (unsigned int i = 0; i < args.size(); ++i)
  {
    // allocate space for the string plus a null character
    argv[i] = new char[args[i].length() + 1];
    std::strncpy(argv[i], args[i].c_str(), args[i].length());
    // place the null character in the end
    argv[i][args[i].length()] = '\0';
  }
  argv[argc] = 0;
  // class to automatically cleanup argv upon destruction
  class Cleanup_argv
  {
  public:
    Cleanup_argv(char ** argv_, int argc_plus_one_)
      : argv(argv_)
      , argc_plus_one(argc_plus_one_)
    {}

    ~Cleanup_argv()
    {
      for (unsigned int i = 0; i < argc_plus_one; ++i)
      {
        delete[] argv[i];
      }
      delete[] argv;
    }

  private:
    char **      argv;
    unsigned int argc_plus_one;
  };
  Cleanup_argv cleanup_argv(argv, argc + 1);

  // antscout->set_stream( out_stream );

  typedef itk::PICSLAdvancedNormalizationToolKit<3> RegistrationType;
  RegistrationType::Pointer                         registration = RegistrationType::New();

  registration->ParseCommandLine(argc, argv);

  typedef itk::CommandLineParser ParserType;
  ParserType::Pointer            parser = ParserType::New();

  typedef ParserType::OptionType OptionType;
  OptionType::Pointer            fileOption = OptionType::New();
  fileOption->SetShortName('f');
  fileOption->SetLongName("file");
  fileOption->SetDescription("The fixed image file used in the image registration algorithm.");

  parser->AddOption(fileOption);

  typedef ParserType::OptionType OptionType;
  OptionType::Pointer            metricOption = OptionType::New();
  metricOption->SetShortName('m');
  metricOption->SetLongName("metric");
  metricOption->SetDescription("The metric used by the image registration algorithm.");

  parser->AddOption(metricOption);

  typedef ParserType::OptionType OptionType;
  OptionType::Pointer            verbosityOption = OptionType::New();
  verbosityOption->SetLongName("verbosity");
  verbosityOption->SetDescription("Mundi vult decipi.  ");

  parser->AddOption(verbosityOption);

  typedef ParserType::OptionType OptionType;
  OptionType::Pointer            dirOption = OptionType::New();
  dirOption->SetLongName("directionality");
  //  dirOption->SetShortName( 'd' );
  dirOption->SetDescription("Mundi vult decipi.  ");

  parser->AddOption(dirOption);

  if (parser->Parse(argc, argv) == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  }

  if (!parser->GetOption("directionality"))
  {
    std::cout << "N ULL" << std::endl;
  }
  else
  {
    for (unsigned int j = 0; j < 3; j++)
    {
      std::cout << parser->ConvertVector<bool>(parser->GetOption("directionality")->GetValue(0))[j] << " x ";
    }
  }

  parser->PrintMenu(std::cout, 7);

  //  std::cout << std::endl << std::endl << "--------------------------------------------"
  //            << std::endl << std::endl;

  //  parser->Print( std::cout << 7 );

  return 0;
};
} // namespace ants
