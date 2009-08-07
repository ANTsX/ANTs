/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkWASPCommandLineParser.cxx,v $
  Language:  C++
  Date:      $Date: 2009/01/22 22:43:11 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkWASPCommandLineParser.h"

namespace itk
{
WASPCommandLineParser
::WASPCommandLineParser()
{
  this->m_Options.clear();
  this->m_Command.clear();
  this->m_CommandDescription.clear();
  this->m_UnknownOptions.clear();
}

void
WASPCommandLineParser
::AddOption( OptionType::Pointer option )
{
  if( ( option->GetShortName() != '\0' ||
        !this->GetOption( option->GetShortName() ) )
      || ( !option->GetLongName().empty() ||
           !this->GetOption( option->GetLongName() ) ) )
    {
    this->m_Options.push_back( option );
    }
  else
    {
    if( option->GetShortName() != '\0' &&
        this->GetOption( option->GetShortName() ) )
      {
      itkWarningMacro( "Duplicate short option '-"
                       << option->GetShortName() << "'" );
      }
    if( !( option->GetLongName().empty() ) &&
        this->GetOption( option->GetLongName() ) )
      {
      itkWarningMacro( "Duplicate long option '--"
                       << option->GetLongName() << "'" );
      }
    }
}

void
WASPCommandLineParser
::Parse( unsigned int argc, char * *argv )
{
  unsigned int n = 0;

  this->m_Command = std::string( argv[n++] );

  while( n < argc )
    {
    std::string argument = std::string( argv[n++] );
    std::string name;

    name.clear();
    if( argument.find( "--" ) == 0 )
      {
      name = argument.substr( 2, argument.length() - 1 );
      }
    else if( argument.find( "-" ) == 0 && argument.find( "--" ) > 0 )
      {
      name = argument.substr( 1, 2 );
      }

    if( !( name.empty() ) )
      {
      OptionType::Pointer option = this->GetOption( name );
      if( !option )
        {
        OptionType::Pointer unknownOption = OptionType::New();
        if( name.length() > 1 )
          {
          unknownOption->SetLongName( name );
          }
        else
          {
          unknownOption->SetShortName( name.at( 0 ) );
          }
        if( n == argc )
          {
          unknownOption->AddValue( "1" );
          }
        else
          {
          for( unsigned int m = n; m < argc; m++ )
            {
            std::string value = std::string( argv[m] );
            if( value.find( "-" ) != 0 )
              {
              unknownOption->AddValue( value );
              }
            else
              {
              if( m == n )
                {
                unknownOption->AddValue( "1" );
                }
              n = m;
              break;
              }
            }
          }
        this->m_UnknownOptions.push_back( unknownOption );
        }
      else  // the option exists
        {
        if( n == argc )
          {
          option->AddValue( "1" );
          }
        else
          {
          for( unsigned int m = n; m < argc; m++ )
            {
            std::string value = std::string( argv[m] );
            if( value.find( "-" ) != 0 )
              {
              option->AddValue( value );
              }
            else
              {
              if( m == n )
                {
                option->AddValue( "1" );
                }
              n = m;
              break;
              }
            }
          }
        }
      }
    }
}

WASPCommandLineParser::OptionType::Pointer
WASPCommandLineParser
::GetOption( std::string name )
{
  if( name.length() == 1 )
    {
    return this->GetOption( name.at( 0 ) );
    }

  OptionListType::iterator it;
  for( it = this->m_Options.begin(); it != this->m_Options.end(); it++ )
    {
    if( name.compare( (*it)->GetLongName() ) == 0 )
      {
      return *it;
      }
    }
  return NULL;
}

WASPCommandLineParser::OptionType::Pointer
WASPCommandLineParser
::GetOption( char name )
{
  OptionListType::iterator it;

  for( it = this->m_Options.begin(); it != this->m_Options.end(); it++ )
    {
    if( name == (*it)->GetShortName() )
      {
      return *it;
      }
    }
  return NULL;
}

void
WASPCommandLineParser
::PrintMenu( std::ostream& os, Indent indent ) const
{
  os << std::endl;
  os << "COMMAND: " << std::endl;
  os << indent << this->m_Command << std::endl;
  if( !this->m_CommandDescription.empty() )
    {
    os << indent << indent << this->m_CommandDescription << std::endl;
    }
  os << std::endl;
  os << "OPTIONS: " << std::endl;

  OptionListType::const_iterator it;
  for( it = this->m_Options.begin(); it != this->m_Options.end(); it++ )
    {
    os << indent;
    if( (*it)->GetShortName() != '\0' )
      {
      os << "-" << (*it)->GetShortName();
      if( !( (*it)->GetLongName() ).empty() )
        {
        os << ", " << "--" << (*it)->GetLongName() << ": " << std::endl;
        }
      else
        {
        os << ": " << std::endl;
        }
      }
    else
      {
      os << "--" << (*it)->GetLongName() << ": " << std::endl;
      }

    if( !( (*it)->GetDescription().empty() ) )
      {
      os << indent << indent << (*it)->GetDescription() << std::endl;
      }
    if( (*it)->GetValues().size() == 1 )
      {
      os << indent << indent << "<VALUES>: " << (*it)->GetValue( 0 );
      if( (*it)->GetParameters( 0 ).size() > 0 )
        {
        os << "[";
        if( (*it)->GetParameters( 0 ).size() == 1 )
          {
          os << (*it)->GetParameter( 0, 0 );
          }
        else
          {
          for( unsigned int i = 0;
               i < (*it)->GetParameters( 0 ).size() - 1; i++ )
            {
            os << (*it)->GetParameter( 0, i ) << ",";
            }
          os << (*it)->GetParameter( 0, (*it)->GetParameters( 0 ).size() - 1 );
          }
        os << "]";
        }
      os << std::endl;
      }
    else if( (*it)->GetValues().size() > 1 )
      {
      os << indent << indent << "<VALUES>: ";
      for( unsigned int n = 0; n < (*it)->GetValues().size() - 1; n++ )
        {
        os << (*it)->GetValue( n );
        if( (*it)->GetParameters( n ).size() > 0 )
          {
          os << "[";
          if( (*it)->GetParameters( n ).size() == 1 )
            {
            os << (*it)->GetParameter( n, 0 ) << "], ";
            }
          else
            {
            for( unsigned int i = 0;
                 i < (*it)->GetParameters( n ).size() - 1; i++ )
              {
              os << (*it)->GetParameter( n, i ) << ",";
              }
            os << (*it)->GetParameter( n,
                                       (*it)->GetParameters( n ).size() - 1 ) << "], ";
            }
          }
        else
          {
          os << ", ";
          }
        }

      unsigned int n = (*it)->GetValues().size() - 1;

      os << (*it)->GetValue( n );
      if( (*it)->GetParameters( n ).size() > 0 )
        {
        os << "[";
        if( (*it)->GetParameters( n ).size() == 1 )
          {
          os << (*it)->GetParameter( n, 0 ) << "]";
          }
        else
          {
          for( unsigned int i = 0;
               i < (*it)->GetParameters( n ).size() - 1; i++ )
            {
            os << (*it)->GetParameter( n, i ) << ",";
            }
          os << (*it)->GetParameter( n,
                                     (*it)->GetParameters( n ).size() - 1 ) << "]";
          }
        }
      }
    os << std::endl;
    }
}

/**
 * Standard "PrintSelf" method
 */
void
WASPCommandLineParser
::PrintSelf( std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Command: " << this->m_Command << std::endl;
  os << indent << "Options: " << std::endl;

  OptionListType::const_iterator it;
  for( it = this->m_Options.begin(); it != this->m_Options.end(); it++ )
    {
    (*it)->Print( os, indent );
    }

  if( this->m_UnknownOptions.size() )
    {
    os << indent << "Unknown Options: " << std::endl;
    OptionListType::const_iterator its;
    for( its = this->m_UnknownOptions.begin();
         its != this->m_UnknownOptions.end(); its++ )
      {
      (*its)->Print( os, indent );
      }
    }
}
} // end namespace itk
