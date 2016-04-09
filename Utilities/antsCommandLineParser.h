/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __antsCommandLineParser_h
#define __antsCommandLineParser_h

#include "antsCommandLineOption.h"

#include "itkDataObject.h"
#include "itkObjectFactory.h"
#include "itkMacro.h"
#include "itkNumericTraits.h"

#include <list>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>

#include <typeinfo>

namespace itk
{
namespace ants
{
  /**
   * A untilty function to convert internal typeid(???).name() to
   * the human readable equivalent format.
   */
  extern std::string ConvertToHumanReadable(const std::string & input);

/** \class CommandLineParser
    \brief Simple command line parser.
    \par
    Parses the standard ( argc, argv ) variables which are stored
    as options in the helper class antsCommandLineOption.  Also contains
    routines for converting types including std::vectors using 'x' as a
    delimiter.  For example, I can specify the 3-element std::vector
    {10, 20, 30} as "10x20x30".
*/

class CommandLineParser
  : public       DataObject
{
public:
  /** Standard class typedefs. */
  typedef CommandLineParser        Self;
  typedef DataObject               Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( CommandLineParser, DataObject );

  typedef CommandLineOption              OptionType;
  typedef std::list<OptionType::Pointer> OptionListType;
  typedef std::list<std::string>         StringListType;

  /**
   * Interface routines
   */

  OptionType::Pointer GetOption( char );

  OptionType::Pointer GetOption( std::string );

  bool starts_with( const std::string &, const std::string & );
  int Parse( unsigned int, char * * );
  bool ValidateFlag(const std::string & currentFlag);

  void AddOption( OptionType::Pointer );

  void PrintMenu( std::ostream& os, Indent indent, bool printShortVersion = false ) const;

  itkSetStringMacro( Command );
  itkGetStringMacro( Command );

  itkSetStringMacro( CommandDescription );
  itkGetStringMacro( CommandDescription );

  OptionListType GetOptions() const
  {
    return this->m_Options;
  }

  OptionListType GetUnknownOptions() const
  {
    return this->m_UnknownOptions;
  }

  /**
   * This feature is designed for a more advanced command line usage
   * where multiple option values are used per stage (e.g.
   * antsRegistration).  Multiple option value are considered to be of
   * the same stage if they are situated adjacently on the command line.
   */
  void AssignStages();


  template <class TValue>
  TValue Convert( std::string optionString ) const
    {
    //Strip whitespace at end
    optionString.erase(optionString.find_last_not_of(" \n\r\t")+1);
    TValue             value;
    std::istringstream iss( optionString );
    if (!(iss >> value)  //Conversion did not fail
      || !( iss.peek() == EOF ) // All content parsed
    )
      {
      std::string internalTypeName( typeid(value).name() );
      itkExceptionMacro( "ERROR: Parse error occured during command line argument processing\n"
        << "ERROR: Unable to convert '" << optionString
        << "' to type '" << internalTypeName << "' as " << ConvertToHumanReadable(internalTypeName) << std::endl);
      }
    return value;
    }

  template <class TValue>
  std::vector<TValue> ConvertVector( std::string optionString ) const
  {
    //Strip whitespace at end
    optionString.erase(optionString.find_last_not_of(" \n\r\t")+1);

    std::vector<std::string> optionElementString;
    std::istringstream f(optionString);
    std::string s;
    while( std::getline(f, s, 'x'))
      {
      optionElementString.push_back(s);
      }

    std::vector< TValue > values;
    for ( std::vector< std::string >::const_iterator oESit = optionElementString.begin();
      oESit != optionElementString.end(); ++oESit)
      {
      const TValue & value = this->Convert<TValue>( *oESit );
      values.push_back ( value );
      }
    return values;
  }

protected:
  CommandLineParser();
  virtual ~CommandLineParser()
  {
  }

  void PrintSelf( std::ostream& os, Indent indent ) const ITK_OVERRIDE;

private:
  CommandLineParser( const Self & ); // purposely not implemented
  void operator=( const Self & );    // purposely not implemented

  std::vector<std::string> RegroupCommandLineArguments( unsigned int, char * * );

  std::string BreakUpStringIntoNewLines( std::string, const std::string, unsigned int ) const;

  void TokenizeString( std::string, std::vector<std::string> &, std::string ) const;

  OptionListType m_Options;
  std::string    m_Command;
  std::string    m_CommandDescription;
  OptionListType m_UnknownOptions;

  char m_LeftDelimiter;
  char m_RightDelimiter;
};
} // end namespace ants
} // end namespace itk

#endif
