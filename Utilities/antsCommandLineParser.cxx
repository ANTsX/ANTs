/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "antsCommandLineParser.h"

#include <algorithm>

#include "itkMath.h"

namespace itk
{
namespace ants
{

std::string
ConvertToHumanReadable(const std::string & input)
{
  using TypeMapType = std::map<std::string, std::string>;
  TypeMapType cnvtMap;
  cnvtMap[typeid(signed char).name()] = "signed int";
  cnvtMap[typeid(unsigned char).name()] = "unsigned int";
  cnvtMap[typeid(signed short).name()] = "signed int";
  cnvtMap[typeid(unsigned short).name()] = "unsigned int";
  cnvtMap[typeid(signed int).name()] = "signed int";
  cnvtMap[typeid(unsigned int).name()] = "unsigned int";
  cnvtMap[typeid(signed long).name()] = "signed int";
  cnvtMap[typeid(unsigned long).name()] = "unsigned int";
  cnvtMap[typeid(float).name()] = "float";
  cnvtMap[typeid(double).name()] = "double";
  cnvtMap[typeid(std::string).name()] = "std::string";
  cnvtMap[typeid(char *).name()] = "char *";

  auto mi = cnvtMap.find(input);
  if (mi == cnvtMap.end())
  {
    return std::string("Unmapped Type");
  }
  return mi->second;
}
CommandLineParser ::CommandLineParser()
  : m_LeftDelimiter('[')
  , m_RightDelimiter(']')
{
  this->m_Options.clear();
  this->m_Command.clear();
  this->m_CommandDescription.clear();
  this->m_UnknownOptions.clear();
}

void
CommandLineParser ::AddOption(OptionType::Pointer option)
{
  if ((option->GetShortName() != '\0' || !this->GetOption(option->GetShortName())) ||
      (!option->GetLongName().empty() || !this->GetOption(option->GetLongName())))
  {
    this->m_Options.push_back(option);
  }
  else
  {
    if (option->GetShortName() != '\0' && this->GetOption(option->GetShortName()))
    {
      itkWarningMacro("Duplicate short option '-" << option->GetShortName() << "'");
    }
    if (!(option->GetLongName().empty()) && this->GetOption(option->GetLongName()))
    {
      itkWarningMacro("Duplicate long option '--" << option->GetLongName() << "'");
    }
  }
}

bool
CommandLineParser ::starts_with(const std::string & s1, const std::string & s2)
{
  return s2.length() <= s1.length() && s1.compare(0, s2.length(), s2) == 0;
}

int
CommandLineParser ::Parse(unsigned int argc, char ** argv)
{
  std::vector<std::string> arguments = this->RegroupCommandLineArguments(argc, argv);

  unsigned int n = 0;
  unsigned int order = 0;
  bool         allFlagsAreValid = true;

  if (arguments.size() > 1)
  {
    this->m_Command = arguments[n++];
  }

  while (n < arguments.size())
  {
    std::string argument = arguments[n++];
    std::string name;

    name.clear();
    if (starts_with(argument, "--"))
    {
      name = argument.substr(2, argument.length() - 1);
    }
    else if (starts_with(argument, "-") && !starts_with(argument, "--"))
    {
      name = argument.substr(1, 2);
    }
    if (!itk::Math::FloatAlmostEqual(atof(name.c_str()), 0.0))
    {
      continue;
    }
    if (!name.empty())
    {
      allFlagsAreValid &= this->ValidateFlag(name);
    }
    if ((!(name.empty())) && (!name.empty()))
    {
      OptionType::Pointer option = this->GetOption(name);
      if (!option)
      {
        OptionType::Pointer unknownOption = OptionType::New();
        if (name.length() > 1)
        {
          unknownOption->SetLongName(name);
        }
        else
        {
          unknownOption->SetShortName(name.at(0));
        }
        if (n == arguments.size())
        {
          unknownOption->AddFunction("1", this->m_LeftDelimiter, this->m_RightDelimiter, order++);
        }
        else
        {
          for (unsigned int m = n; m < arguments.size(); m++)
          {
            std::string function = arguments[m];
            if (!starts_with(function, "-"))
            {
              unknownOption->AddFunction(function, this->m_LeftDelimiter, this->m_RightDelimiter, order++);
            }
            else
            {
              if (m == n)
              {
                unknownOption->AddFunction("1", this->m_LeftDelimiter, this->m_RightDelimiter, order++);
              }
              n = m;
              break;
            }
          }
        }
        this->m_UnknownOptions.push_back(unknownOption);
      }
      else // the option exists
      {
        if (n == arguments.size())
        {
          option->AddFunction("1", this->m_LeftDelimiter, this->m_RightDelimiter, order++);
        }
        else
        {
          for (unsigned int m = n; m < arguments.size(); m++)
          {
            std::string function = arguments[m];
            if (!starts_with(function, "-") || !itk::Math::FloatAlmostEqual(atof(function.c_str()), 0.0))
            {
              option->AddFunction(function, this->m_LeftDelimiter, this->m_RightDelimiter, order++);
            }
            else
            {
              if (m == n)
              {
                option->AddFunction("1", this->m_LeftDelimiter, this->m_RightDelimiter, order++);
              }
              n = m;
              break;
            }
          }
        }
      }
    }
  }
  if (!allFlagsAreValid)
  {
    std::cerr << "ERROR:  Invalid command line flags found! Aborting execution." << std::endl;
    return EXIT_FAILURE;
  }
  this->AssignStages();
  return EXIT_SUCCESS;
}

std::vector<std::string>
CommandLineParser ::RegroupCommandLineArguments(unsigned int argc, char ** argv)
{
  /**
   * Inclusion of this function allows the user to use spaces inside
   * the left and right delimiters.  Also replace other left and right
   * delimiters.
   */
  std::vector<std::string> arguments;

  std::string currentArg("");
  bool        isArgOpen = false;

  for (unsigned int n = 0; n < argc; n++)
  {
    std::string a(argv[n]);

    if (n == 0)
    {
      arguments.push_back(a);
    }
    else
    {
      // replace left delimiters
      std::replace(a.begin(), a.begin() + 1, '{', '[');
      std::replace(a.begin(), a.begin() + 1, '(', '[');
      std::replace(a.begin(), a.begin() + 1, '<', '[');

      // replace right delimiters
      std::replace(a.end() - 1, a.end(), '}', ']');
      std::replace(a.end() - 1, a.end(), ')', ']');
      std::replace(a.end() - 1, a.end(), '>', ']');

      if (isArgOpen)
      {
        std::size_t leftDelimiterPosition = a.find(this->m_LeftDelimiter);
        if (leftDelimiterPosition != std::string::npos)
        {
          itkExceptionMacro("Incorrect command line specification. Missing leftDelimiterPosition? " << a);
        }

        std::size_t rightDelimiterPosition = a.find(this->m_RightDelimiter);
        if (rightDelimiterPosition != std::string::npos)
        {
          if (rightDelimiterPosition < a.length() - 1)
          {
            itkExceptionMacro("Incorrect command line specification. Missing rightDelimiterPosition? " << a);
          }
          else
          {
            currentArg += a;
            arguments.push_back(currentArg);
            currentArg.clear();
            isArgOpen = false;
          }
        }
        else
        {
          currentArg += a;
        }
      }
      else
      {
        std::size_t leftDelimiterPosition = a.find(this->m_LeftDelimiter);
        std::size_t rightDelimiterPosition = a.find(this->m_RightDelimiter);

        if (leftDelimiterPosition == std::string::npos)
        {
          if (rightDelimiterPosition == std::string::npos)
          {
            currentArg += a;
            arguments.push_back(currentArg);
            currentArg.clear();
          }
          else
          {
            itkExceptionMacro("Incorrect command line specification. " << a);
          }
        }
        else if (leftDelimiterPosition != std::string::npos && rightDelimiterPosition != std::string::npos &&
                 leftDelimiterPosition < rightDelimiterPosition)
        {
          if (rightDelimiterPosition < a.length() - 1)
          {
            itkExceptionMacro("Incorrect command line specification. " << a);
          }
          currentArg += a;
          arguments.push_back(currentArg);
          currentArg.clear();
          isArgOpen = false;
        }
        else if (rightDelimiterPosition == std::string::npos && leftDelimiterPosition != std::string::npos)
        {
          currentArg += a;
          isArgOpen = true;
        }
      }
    }
  }

  return arguments;
}

CommandLineParser::OptionType::Pointer
CommandLineParser ::GetOption(std::string name)
{
  if (name.length() == 1)
  {
    return this->GetOption(name.at(0));
  }

  OptionListType::iterator it;
  for (it = this->m_Options.begin(); it != this->m_Options.end(); ++it)
  {
    if (name.compare((*it)->GetLongName()) == 0)
    {
      return *it;
    }
  }
  return nullptr;
}

CommandLineParser::OptionType::Pointer
CommandLineParser ::GetOption(char name)
{
  OptionListType::iterator it;

  for (it = this->m_Options.begin(); it != this->m_Options.end(); ++it)
  {
    if (name == (*it)->GetShortName())
    {
      return *it;
    }
  }
  return nullptr;
}

bool
CommandLineParser::ValidateFlag(const std::string & currentFlag)
{
  bool validFlagFound = false;
  for (OptionListType::const_iterator it = this->m_Options.begin(); it != this->m_Options.end(); ++it)
  {
    const char        shortName = (*it)->GetShortName();
    const std::string longName = (*it)->GetLongName();
    if (((currentFlag.size() == 1) && (shortName == currentFlag[0])) || (longName == currentFlag))
    {
      validFlagFound = true;
    }
  }

  if ((!validFlagFound) && (!currentFlag.empty()) && (itk::Math::FloatAlmostEqual(atof(currentFlag.c_str()), 0.0)))
  {
    std::cout << "ERROR:  Invalid flag provided " << currentFlag << std::endl;
  }
  return validFlagFound;
}

void
CommandLineParser ::PrintMenu(std::ostream & os, Indent indent, bool printShortVersion) const
{
  os << std::endl;
  os << "COMMAND: " << std::endl;
  os << indent << this->m_Command << std::endl;
  if (!this->m_CommandDescription.empty() && !printShortVersion)
  {
    std::stringstream ss1;
    ss1 << indent << indent;

    std::stringstream ss2;
    ss2 << this->m_CommandDescription;

    std::string description = this->BreakUpStringIntoNewLines(ss2.str(), ss1.str(), 80);

    os << indent << indent << description << std::endl;
  }
  os << std::endl;
  os << "OPTIONS: " << std::endl;

  OptionListType::const_iterator it;
  for (it = this->m_Options.begin(); it != this->m_Options.end(); ++it)
  {
    os << indent;
    std::stringstream ss;
    ss << indent;

    if ((*it)->GetShortName() != '\0')
    {
      os << "-" << (*it)->GetShortName();
      ss << Indent(2);
      if (!((*it)->GetLongName()).empty())
      {
        os << ", "
           << "--" << (*it)->GetLongName() << " " << std::flush;
        ss << Indent(static_cast<std::size_t>(5) + ((*it)->GetLongName()).length());
      }
      else
      {
        os << " " << std::flush;
        ss << Indent(1);
      }
    }
    else
    {
      os << "--" << (*it)->GetLongName() << " " << std::flush;
      ss << Indent(static_cast<std::size_t>(3) + ((*it)->GetLongName()).length());
    }
    if ((*it)->GetNumberOfUsageOptions() > 0)
    {
      os << (*it)->GetUsageOption(0) << std::endl;
      for (unsigned int i = 1; i < (*it)->GetNumberOfUsageOptions(); i++)
      {
        os << ss.str() << (*it)->GetUsageOption(i) << std::endl;
      }
    }
    else
    {
      os << std::endl;
    }

    if (!((*it)->GetDescription().empty()) && !printShortVersion)
    {
      std::stringstream ss1;
      ss1 << indent << indent;

      std::stringstream ss2;
      ss2 << (*it)->GetDescription();

      std::string description = this->BreakUpStringIntoNewLines(ss2.str(), ss1.str(), 80);

      os << indent << indent << description << std::endl;
    }
    if (!printShortVersion)
    {
      if ((*it)->GetFunctions().size() == 1)
      {
        os << indent << indent << "<VALUES>: " << (*it)->GetFunction(0)->GetName();
        if (!(*it)->GetFunction(0)->GetParameters().empty())
        {
          os << "[";
          if ((*it)->GetFunction(0)->GetParameters().size() == 1)
          {
            os << (*it)->GetFunction(0)->GetParameter(0);
          }
          else
          {
            for (unsigned int i = 0; i < (*it)->GetFunction(0)->GetParameters().size() - 1; i++)
            {
              os << (*it)->GetFunction(0)->GetParameter(i) << ",";
            }
            os << (*it)->GetFunction(0)->GetParameter((*it)->GetFunction(0)->GetParameters().size() -
                                                      static_cast<std::size_t>(1));
          }
          os << "]";
        }
        os << std::endl;
      }
      else if ((*it)->GetFunctions().size() > 1)
      {
        os << indent << indent << "<VALUES>: ";
        for (unsigned int n = 0; n < (*it)->GetFunctions().size() - 1; n++)
        {
          os << (*it)->GetFunction(n)->GetName();
          if (!(*it)->GetFunction(n)->GetParameters().empty())
          {
            os << "[";
            if ((*it)->GetFunction(n)->GetParameters().size() == 1)
            {
              os << (*it)->GetFunction(n)->GetParameter(0) << "], ";
            }
            else
            {
              for (unsigned int i = 0; i < (*it)->GetFunction(n)->GetParameters().size() - 1; i++)
              {
                os << (*it)->GetFunction(n)->GetParameter(i) << ",";
              }
              os << (*it)->GetFunction(n)->GetParameter((*it)->GetFunction(n)->GetParameters().size() -
                                                        static_cast<std::size_t>(1))
                 << "], ";
            }
          }
          else
          {
            os << ", ";
          }
        }

        unsigned int nn = (*it)->GetFunctions().size() - static_cast<std::size_t>(1);

        os << (*it)->GetFunction(nn)->GetName();
        if (!(*it)->GetFunction(nn)->GetParameters().empty())
        {
          os << "[";
          if ((*it)->GetFunction(nn)->GetParameters().size() == 1)
          {
            os << (*it)->GetFunction(nn)->GetParameter(0) << "]";
          }
          else
          {
            for (unsigned int i = 0; i < (*it)->GetFunction(nn)->GetParameters().size() - 1; i++)
            {
              os << (*it)->GetFunction(nn)->GetParameter(i) << ",";
            }
            os << (*it)->GetFunction(nn)->GetParameter((*it)->GetFunction(nn)->GetParameters().size() -
                                                       static_cast<std::size_t>(1))
               << "]";
          }
        }
      }
      os << std::endl;
    }
  }
}

std::string
CommandLineParser ::BreakUpStringIntoNewLines(std::string  longString,
                                              std::string  indentString,
                                              unsigned int numberOfCharactersPerLine) const
{
  std::vector<std::string> tokens;

  this->TokenizeString(longString, tokens, " ");

  std::string  newString("");
  unsigned int currentTokenId = 0;
  unsigned int currentLineLength = 0;
  while (currentTokenId < tokens.size())
  {
    if (tokens[currentTokenId].length() >= numberOfCharactersPerLine)
    {
      newString += (std::string("\n") + tokens[currentTokenId] + std::string("\n"));
      currentTokenId++;
      currentLineLength = 0;
    }
    else if (currentTokenId < tokens.size() &&
             currentLineLength + tokens[currentTokenId].length() > numberOfCharactersPerLine)
    {
      newString += (std::string("\n") + indentString);
      currentLineLength = 0;
    }
    else
    {
      newString += (tokens[currentTokenId] + std::string(" "));
      currentLineLength += (tokens[currentTokenId].length() + 1);
      currentTokenId++;
    }
  }

  return newString;
}

void
CommandLineParser ::TokenizeString(std::string str, std::vector<std::string> & tokens, std::string delimiters) const
{
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  std::string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos)
  {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}

void
CommandLineParser ::AssignStages()
{
  OptionListType::const_iterator it;

  for (it = this->m_Options.begin(); it != this->m_Options.end(); ++it)
  {
    using OptionFunctionStackType = OptionType::FunctionStackType;
    OptionFunctionStackType functions = (*it)->GetFunctions();

    OptionFunctionStackType::const_iterator it2;

    unsigned int previousOrder = 0;
    unsigned int currentOrder = 0;
    for (it2 = functions.begin(); it2 != functions.end(); ++it2)
    {
      if (it2 == functions.begin())
      {
        previousOrder = (*it2)->GetArgOrder();
        (*it2)->SetStageID(0);
      }
      else
      {
        currentOrder = (*it2)->GetArgOrder();
        if (previousOrder == currentOrder + 1)
        {
          (*it2)->SetStageID(functions[it2 - functions.begin() - static_cast<std::size_t>(1)]->GetStageID());
        }
        else
        {
          (*it2)->SetStageID(functions[it2 - functions.begin() - static_cast<std::size_t>(1)]->GetStageID() + 1);
        }
        previousOrder = currentOrder;
      }
    }
  }
}

/**
 * Standard "PrintSelf" method
 */
void
CommandLineParser ::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Command: " << this->m_Command << std::endl;
  os << indent << "Options: " << std::endl;

  OptionListType::const_iterator it;
  for (it = this->m_Options.begin(); it != this->m_Options.end(); ++it)
  {
    (*it)->Print(os, indent);
  }

  if (!this->m_UnknownOptions.empty())
  {
    os << indent << "Unknown Options: " << std::endl;
    OptionListType::const_iterator its;
    for (its = this->m_UnknownOptions.begin(); its != this->m_UnknownOptions.end(); ++its)
    {
      (*its)->Print(os, indent);
    }
  }
}
} // end namespace ants
} // end namespace itk
