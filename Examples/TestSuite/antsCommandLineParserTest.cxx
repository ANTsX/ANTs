/*=========================================================================
 *
 *  Regression test for bracket handling in itk::ants::CommandLineParser.
 *  An option value opened with '[' but never closed used to swallow every
 *  following argument and then be dropped without an error, so the option
 *  and all flags after it were silently ignored.
 *
 *=========================================================================*/
#include "antsCommandLineParser.h"
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>

namespace
{
using ParserType = itk::ants::CommandLineParser;

ParserType::Pointer
MakeParser()
{
  ParserType::Pointer parser = ParserType::New();
  {
    ParserType::OptionType::Pointer option = ParserType::OptionType::New();
    option->SetLongName("masks");
    option->SetShortName('x');
    parser->AddOption(option);
  }
  {
    ParserType::OptionType::Pointer option = ParserType::OptionType::New();
    option->SetLongName("float");
    parser->AddOption(option);
  }
  return parser;
}

int
Parse(ParserType::Pointer parser, std::vector<std::string> args)
{
  std::vector<char *> argv;
  for (auto & a : args)
  {
    argv.push_back(&a[0]);
  }
  return parser->Parse(static_cast<unsigned int>(argv.size()), argv.data());
}
} // namespace

int
main()
{
  int failures = 0;

  // Closed brackets, including a value split across arguments, parse normally.
  {
    ParserType::Pointer parser = MakeParser();
    if (Parse(parser, { "prog", "-x", "[fixed.nii,", "moving.nii]", "--float", "1" }) != EXIT_SUCCESS)
    {
      std::cerr << "FAIL: well-formed command line was rejected" << std::endl;
      ++failures;
    }
    else
    {
      ParserType::OptionType::Pointer masks = parser->GetOption("masks");
      ParserType::OptionType::Pointer useFloat = parser->GetOption("float");
      if (masks->GetNumberOfFunctions() != 1 || masks->GetFunction(0)->GetNumberOfParameters() != 2 ||
          masks->GetFunction(0)->GetParameter(0) != "fixed.nii" ||
          masks->GetFunction(0)->GetParameter(1) != "moving.nii")
      {
        std::cerr << "FAIL: --masks parameters not parsed as [fixed.nii,moving.nii]" << std::endl;
        ++failures;
      }
      if (useFloat->GetNumberOfFunctions() != 1 || useFloat->GetFunction(0)->GetName() != "1")
      {
        std::cerr << "FAIL: --float not parsed as 1" << std::endl;
        ++failures;
      }
    }
  }

  // A value missing its closing bracket at the end of the command line must be
  // reported, not silently dropped together with the arguments after it.
  {
    ParserType::Pointer parser = MakeParser();
    bool                thrown = false;
    try
    {
      Parse(parser, { "prog", "-x", "[fixed.nii,moving.nii", "--float", "1" });
    }
    catch (const itk::ExceptionObject &)
    {
      thrown = true;
    }
    if (!thrown)
    {
      ParserType::OptionType::Pointer masks = parser->GetOption("masks");
      std::cerr << "FAIL: unclosed '[' was accepted; --masks has " << masks->GetNumberOfFunctions()
                << " value(s), first is '" << (masks->GetNumberOfFunctions() ? masks->GetFunction(0)->GetName() : "")
                << "', --float has " << parser->GetOption("float")->GetNumberOfFunctions() << " value(s)" << std::endl;
      ++failures;
    }
  }

  if (failures)
  {
    std::cerr << failures << " failure(s)" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "antsCommandLineParserTest passed" << std::endl;
  return EXIT_SUCCESS;
}
