/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkCommandLineOption.h,v $
  Language:  C++
  Date:      $Date: 2008/11/15 23:46:06 $
  Version:   $Revision: 1.16 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkCommandLineOption_h
#define __itkCommandLineOption_h

#include "itkDataObject.h"
#include "itkObjectFactory.h"

#include "itkMacro.h"
#include "itkNumericTraits.h"

#include <list>
#include <string>
#include <deque>

namespace itk
{
class ITK_EXPORT CommandLineOption
  : public DataObject
{
public:
  typedef CommandLineOption  Self;
  typedef DataObject         Superclass;
  typedef SmartPointer<Self> Pointer;

  itkNewMacro( Self );

  itkTypeMacro( Option, DataObject );

  typedef std::string                ValueType;
  typedef std::deque<ValueType>      ValueStackType;
  typedef std::deque<ValueStackType> ParameterStackType;

  ValueStackType GetValues()
  {
    return this->m_Values;
  }

  unsigned int GetNumberOfValues()
  {
    return this->m_Values.size();
  }

  std::string GetValue( unsigned int i = 0 )
  {
    if( i < this->m_Values.size() )
      {
      return this->m_Values[i];
      }
    else
      {
      return std::string( "" );
      }
  }

  ValueStackType GetParameters( unsigned int i = 0 )
  {
    if( i < this->m_Parameters.size() )
      {
      return this->m_Parameters[i];
      }
    else
      {
      ValueStackType empty;
      return empty;
      }
  }

  std::string GetParameter( unsigned int i = 0, unsigned int j = 0 )
  {
    if( i < this->m_Parameters.size() && j < this->m_Parameters[i].size() )
      {
      return this->m_Parameters[i][j];
      }
    return std::string( "" );
  }

  unsigned int GetNumberOfParameters( unsigned int i = 0 )
  {
    if( i < this->m_Parameters.size() )
      {
      return this->m_Parameters[i].size();
      }
    return 0;
  }

  itkSetMacro( ShortName, char );
  itkGetMacro( ShortName, char );

  itkSetMacro( LongName, std::string );
  itkGetMacro( LongName, std::string );

  itkSetMacro( Description, std::string );
  itkGetMacro( Description, std::string );

  void AddValue( std::string );

  void SetValue( unsigned int, std::string );

protected:
  CommandLineOption();
  virtual ~CommandLineOption()
  {
  };
private:

  char               m_ShortName;
  std::string        m_LongName;
  std::string        m_Description;
  ValueStackType     m_Values;
  ParameterStackType m_Parameters;
};
} // end namespace itk

#endif
