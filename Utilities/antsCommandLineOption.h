/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __antsCommandLineOption_h
#define __antsCommandLineOption_h

#include "itkDataObject.h"
#include "itkObjectFactory.h"
#include "itkMacro.h"
#include "itkNumericTraits.h"

#include <list>
#include <string>
#include <deque>

namespace itk
{
namespace ants
{
/** \class CommandLineOption
    \brief Simple data structure for holding command line options.
    An option can have multiple values with each value holding 0 or more
    parameters.  E.g. suppose we were creating an image registration program
    which has several transformation model options such as 'rigid', 'affine',
    and 'deformable'.  A instance of the command line option could have a
    long name of "transformation", short name 't', and description
    "Transformation model---rigid, affine, or deformable".  The values for
    this option would be "rigid", "affine", and "deformable".  Each value
    would then hold parameters that relate to that value.  For example, a
    possible subsection of the command line would be

    " --transformation rigid[parameter1,parameter2,etc.]
      -m mutualinformation[parameter1] --optimization gradientdescent"
*/

class OptionFunction
  : public       DataObject
{
public:
  OptionFunction() :
    m_Name( "" ),
    m_ArgOrder( 0 ),
    m_StageID( 0 )
  {
  };
  ~OptionFunction()
  {
  };

  typedef OptionFunction     Self;
  typedef DataObject         Superclass;
  typedef SmartPointer<Self> Pointer;

  itkNewMacro( Self );

  itkTypeMacro( Option, DataObject );

  typedef std::deque<std::string> ParameterStackType;

  itkSetStringMacro( Name );
  itkGetStringMacro( Name );
  itkGetMacro( Name, std::string );

  itkSetMacro( ArgOrder, unsigned int );
  itkGetConstMacro( ArgOrder, unsigned int );

  itkSetMacro( StageID, unsigned int );
  itkGetConstMacro( StageID, unsigned int );

  ParameterStackType GetParameters()
  {
    return this->m_Parameters;
  }

  void SetParameters( ParameterStackType parameters )
  {
    this->m_Parameters = parameters;
    this->Modified();
  }

  std::string GetParameter( unsigned int i = 0 )
  {
    if( i < this->m_Parameters.size() )
      {
      return this->m_Parameters[i];
      }
    else
      {
      std::string empty( "" );
      return empty;
      }
  }

  unsigned int GetNumberOfParameters()
  {
    return this->m_Parameters.size();
  }

private:
  std::string        m_Name;
  unsigned int       m_ArgOrder;
  unsigned int       m_StageID;
  ParameterStackType m_Parameters;
};

class CommandLineOption
  : public       DataObject
{
public:
  typedef CommandLineOption  Self;
  typedef DataObject         Superclass;
  typedef SmartPointer<Self> Pointer;

  itkNewMacro( Self );

  itkTypeMacro( Option, DataObject );

  typedef OptionFunction OptionFunctionType;

  typedef std::deque<OptionFunctionType::Pointer> FunctionStackType;
  typedef std::deque<std::string>                 UsageOptionStackType;

  FunctionStackType GetFunctions()
  {
    return this->m_OptionFunctions;
  }

  unsigned int GetNumberOfFunctions()
  {
    return this->m_OptionFunctions.size();
  }

  OptionFunction::Pointer GetFunction( unsigned int i = 0 )
  {
    if( i < this->m_OptionFunctions.size() )
      {
      return this->m_OptionFunctions[i];
      }
    else
      {
      return ITK_NULLPTR;
      }
  }

  UsageOptionStackType GetUsageOptions()
  {
    return this->m_UsageOptions;
  }

  unsigned int GetNumberOfUsageOptions()
  {
    return this->m_UsageOptions.size();
  }

  std::string GetUsageOption( unsigned int i = 0 )
  {
    if( i < this->m_UsageOptions.size() )
      {
      return this->m_UsageOptions[i];
      }
    else
      {
      return std::string( "" );
      }
  }

  itkSetMacro( ShortName, char );
  itkGetConstMacro( ShortName, char );

  itkSetStringMacro( LongName );
  itkGetConstMacro( LongName, std::string );

  itkSetStringMacro( Description );
  itkGetMacro( Description, std::string );

  void AddFunction( std::string, char, char, unsigned int order = 0 );

  void AddFunction( std::string s )
  {
    this->AddFunction( s, '[', ']' );
  }

  void SetUsageOption( unsigned int, std::string );

protected:
  CommandLineOption();
  virtual ~CommandLineOption()
  {
  };
private:
  char                 m_ShortName;
  std::string          m_LongName;
  std::string          m_Description;
  UsageOptionStackType m_UsageOptions;
  FunctionStackType    m_OptionFunctions;
};
} // end namespace ants
} // end namespace itk

#endif
