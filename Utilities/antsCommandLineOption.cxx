/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "antsCommandLineOption.h"

namespace itk
{
namespace ants
{
CommandLineOption
::CommandLineOption() : m_ShortName( '\0' ),
  m_LongName( "" ),
  m_Description( "" )
{
  this->m_OptionFunctions.clear();
  this->m_UsageOptions.clear();
}

void
CommandLineOption
::AddFunction( std::string functionString, char leftDelimiter, char rightDelimiter, unsigned int order )
{
  OptionFunctionType::Pointer optionFunction = OptionFunctionType::New();

  optionFunction->SetArgOrder( order );

  std::string::size_type leftDelimiterPos = functionString.find( leftDelimiter );
  std::string::size_type rightDelimiterPos = functionString.find( rightDelimiter );

  if( leftDelimiterPos == std::string::npos ||
      rightDelimiterPos == std::string::npos )
    {
    optionFunction->SetName( functionString );
    this->m_OptionFunctions.push_front( optionFunction );
    }
  else
    {
    OptionFunctionType::ParameterStackType parameters;

    optionFunction->SetName( functionString.substr( 0, leftDelimiterPos ) );

    std::string::size_type leftPos = leftDelimiterPos;
    std::string::size_type rightPos = functionString.find( ',', leftPos + 1 );
    while( rightPos != std::string::npos )
      {
      parameters.push_back( functionString.substr( leftPos + 1, rightPos - leftPos - 1 ) );
      leftPos = rightPos;
      rightPos = functionString.find( ',', leftPos + 1 );
      }

    rightPos = rightDelimiterPos;
    parameters.push_back( functionString.substr( leftPos + 1, rightPos - leftPos - 1 ) );

    optionFunction->SetParameters( parameters );

    this->m_OptionFunctions.push_front( optionFunction );
    }

  this->Modified();
}

void
CommandLineOption
::SetUsageOption( unsigned int i, std::string usage )
{
  if( i >= this->m_UsageOptions.size() )
    {
    this->m_UsageOptions.resize( i + 1 );
    }
  this->m_UsageOptions[i] = usage;
  this->Modified();
}
} // end namespace ants
} // end namespace itk
