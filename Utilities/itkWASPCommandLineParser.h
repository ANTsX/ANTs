/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkWASPCommandLineParser.h,v $
  Language:  C++
  Date:      $Date: 2009/01/22 22:43:11 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkWASPCommandLineParser_h
#define __itkWASPCommandLineParser_h

#include "itkDataObject.h"
#include "itkObjectFactory.h"

#include "itkWASPCommandLineOption.h"
#include "itkMacro.h"
#include "itkNumericTraits.h"

#include <list>
#include <string>
#include <vector>

namespace itk
{
/** \class WASPCommandLineParser
    \brief Simple command line parser.
    \par
    Parses the standard ( argc, argv ) variables which are stored
    as options in the helper class itkWASPCommandLineOption.  Also contains
    routines for converting types including std::vectors using 'x' as a
    delimiter.  For example, I can specify the 3-element std::vector
    {10, 20, 30} as "10x20x30".
*/

class ITK_EXPORT WASPCommandLineParser
  : public DataObject
{
public:
  /** Standard class typedefs. */
  typedef WASPCommandLineParser    Self;
  typedef DataObject               Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( WASPCommandLineParser, DataObject );

  typedef WASPCommandLineOption          OptionType;
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

  itkSetStringMacro( Command );
  itkGetStringMacro( Command );

  itkSetStringMacro( CommandDescription );
  itkGetStringMacro( CommandDescription );

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
          element = optionString.substr(
              crossposfrom + 1, optionString.length() );
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
  WASPCommandLineParser();
  virtual ~WASPCommandLineParser()
  {
  }

  void PrintSelf( std::ostream& os, Indent indent ) const;

private:
  WASPCommandLineParser( const Self & ); // purposely not implemented
  void operator=( const Self & );        // purposely not implemented

  OptionListType m_Options;
  std::string    m_Command;
  std::string    m_CommandDescription;
  OptionListType m_UnknownOptions;
};
} // end namespace itk

#endif
