/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkCommandLineParser.h,v $
  Language:  C++
  Date:      $Date: 2008/11/15 23:46:06 $
  Version:   $Revision: 1.16 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkCommandLineParser_h
#define __itkCommandLineParser_h

#include "itkDataObject.h"
#include "itkObjectFactory.h"

#include "itkCommandLineOption.h"
#include "itkMacro.h"
#include "itkNumericTraits.h"

#include <list>
#include <string>
#include <vector>

namespace itk
{
/**
 * Definition for the command line parsers
 */
class ITK_EXPORT CommandLineParser
  : public DataObject
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

  void Parse( unsigned int, char * * );

  void AddOption( OptionType::Pointer );

  void PrintMenu( std::ostream& os, Indent indent ) const;

  itkSetMacro( CommandDescription, std::string );
  itkGetMacro( CommandDescription, std::string );

  template <class TValue>
  TValue Convert( std::string optionString )
  {
    TValue             value;
    std::istringstream iss( optionString );

    iss >> value;
    return value;
  }

  template <class TValue>
  std::vector<TValue> ConvertVector( std::string optionString )
  {
    std::vector<TValue>    values;
    std::string::size_type crosspos = optionString.find( 'x', 0 );

    if( crosspos == std::string::npos )
      {
      values.push_back( this->Convert<TValue>( optionString ) );
      }
    else
      {
      std::string        element = optionString.substr( 0, crosspos );
      TValue             value;
      std::istringstream iss( element );
      iss >> value;
      values.push_back( value );
      while( crosspos != std::string::npos )
        {
        std::string::size_type crossposfrom = crosspos;
        crosspos = optionString.find( 'x', crossposfrom + 1 );
        if( crosspos == std::string::npos )
          {
          element = optionString.substr( crossposfrom + 1, optionString.length() );
          }
        else
          {
          element = optionString.substr( crossposfrom + 1, crosspos );
          }
        std::istringstream iss( element );
        iss >> value;
        values.push_back( value );
        }
      }
    return values;
  }

protected:
  CommandLineParser();
  virtual ~CommandLineParser()
  {
  }

  void PrintSelf( std::ostream& os, Indent indent ) const;

private:
  CommandLineParser( const Self & ); // purposely not implemented
  void operator=( const Self & );    // purposely not implemented

  OptionListType m_Options;
  std::string    m_Command;
  std::string    m_CommandDescription;
  OptionListType m_UnknownOptions;
};
} // end namespace itk

#endif
