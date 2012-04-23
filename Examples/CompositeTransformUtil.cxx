#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "antsUtilities.h"
#include <algorithm>

#include "itkTransform.h"
#include "itkCompositeTransform.h"
#include "itkDisplacementFieldTransform.h"
#include "antsCommandLineParser.h"
#include "itkantsReadWriteTransform.h"
#include "itkTransformFactory.h"

namespace ants
{
bool MatOffRegistered[2] = { false, false };

template <unsigned int VImageDimension>
void RegisterMatOff()
{
  if( !MatOffRegistered[VImageDimension - 2] )
    {
    MatOffRegistered[VImageDimension - 2] = true;
    // Register the matrix offset transform base class to the
    // transform factory for compatibility with the current ANTs.
    typedef itk::MatrixOffsetTransformBase<double, VImageDimension,
                                           VImageDimension> MatrixOffsetTransformType;
    itk::TransformFactory<MatrixOffsetTransformType>::RegisterTransform();
    }
}

/**
 * print the usage and exit
 */
void UsageAndExit()
{
  antscout << "Usage: CompositeTransformUtil --disassemble <CompositeTransform FileName>"
           << " <transform name prefix>" << std::endl
           << "or" << std::endl
           << "CompositeTransformUtil  --assemble <CompositeTransform> "
           << "<transform 1> <transform 2 > ... <transform N>" << std::endl;
  throw std::exception();
}

/**
 * given a composite transform, write out all its component
 * transforms.
 */
template <unsigned int VImageDimension>
int
Disassemble(itk::TransformBase *transform, const std::string & transformName, const std::string & prefix)
{
  typedef itk::CompositeTransform<double, VImageDimension>                  CompositeTransformType;
  typedef typename CompositeTransformType::TransformTypePointer             TransformPointer;
  typedef typename itk::DisplacementFieldTransform<double, VImageDimension> DisplacementFieldTransformType;
  typedef typename DisplacementFieldTransformType::DisplacementFieldType    DisplacementFieldType;

  CompositeTransformType *composite =
    dynamic_cast<CompositeTransformType *>(transform);
  if( composite == 0 )
    {
    antscout << "Transform File " << transformName << " is a "
             << transform->GetNameOfClass() << " not a Composite Transform."
             << std::endl;
    return EXIT_FAILURE;
    }

  unsigned int numTransforms = composite->GetNumberOfTransforms();
  for( unsigned int i = 0; i < numTransforms; ++i )
    {
    TransformPointer                curXfrm = composite->GetNthTransform(i);
    DisplacementFieldTransformType *dispXfrm =
      dynamic_cast<DisplacementFieldTransformType *>(curXfrm.GetPointer() );
    std::stringstream fname;
    fname << std::setfill('0') << std::setw(2) << i
          << "_" << prefix << "_"
          << curXfrm->GetNameOfClass();
    if( dispXfrm != 0 )
      {
      fname << ".nii.gz";    // if it's a displacement field transform
      }
    else
      {
      fname << ".txt";
      }
    itk::ants::WriteTransform<VImageDimension>(curXfrm, fname.str() );
    }
  return EXIT_SUCCESS;
}

int Disassemble(const std::string & CompositeName,
                const std::string & Prefix)
{
  itk::TransformBase::Pointer transform = itk::ants::ReadTransform<2>(CompositeName).GetPointer();

  if( transform.IsNull() )
    {
    transform = itk::ants::ReadTransform<3>(CompositeName).GetPointer();
    if( transform.IsNull() )
      {
      return EXIT_FAILURE; // ReadTransform prints error messages on
                           // failure.
      }
    }
  unsigned int inDim(transform->GetInputSpaceDimension() ),
  outDim(transform->GetOutputSpaceDimension() );
  if( inDim != outDim )
    {
    antscout << "Can't handle mixed input & output dimension: input("
             << inDim << ") output (" << outDim << ")" << std::endl;
    return EXIT_FAILURE;
    }

  switch( inDim )
    {
    case 2:
      return Disassemble<2>(transform, CompositeName, Prefix);
    case 3:
      return Disassemble<3>(transform, CompositeName, Prefix);
    default:
      antscout << "Unknown dimension " << inDim << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}

template <unsigned int VImageDimension>
int
Assemble(const std::string & CompositeName,
         const std::vector<std::string> & transformNames,
         typename itk::CompositeTransform<double, VImageDimension>::TransformType::Pointer & firstTransform)
{
  typedef itk::CompositeTransform<double, VImageDimension> CompositeTransformType;
  typedef typename CompositeTransformType::TransformType   TransformType;
  typename CompositeTransformType::Pointer composite =
    CompositeTransformType::New();
  composite->AddTransform(firstTransform);
  for( unsigned i = 1; i < transformNames.size(); ++i )
    {
    typename TransformType::Pointer curXfrm = itk::ants::ReadTransform<VImageDimension>(transformNames[i]);
    if( curXfrm.IsNull() )
      {
      return EXIT_FAILURE; // ReadTransform will complain if anything goes wrong.
      }
    composite->AddTransform(curXfrm);
    }
  typename TransformType::Pointer genericXfrmPtr = composite.GetPointer();
  return itk::ants::WriteTransform<VImageDimension>(genericXfrmPtr, CompositeName);
}

int Assemble(const std::string & CompositeName,
             const std::vector<std::string> & transformNames)
{
  typedef itk::CompositeTransform<double, 2>::TransformType Transform2DType;
  typedef itk::CompositeTransform<double, 3>::TransformType Transform3DType;

  Transform2DType::Pointer transform2D = itk::ants::ReadTransform<2>(transformNames[0]);
  if( transform2D.IsNotNull() )
    {
    return Assemble<2>(CompositeName, transformNames, transform2D);
    }
  Transform3DType::Pointer transform3D = itk::ants::ReadTransform<3>(transformNames[0]);
  if( transform3D.IsNotNull() )
    {
    return Assemble<3>(CompositeName, transformNames, transform3D);
    }
  return EXIT_FAILURE;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int CompositeTransformUtil( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "CompositeTransformUtil" );
  std::remove( args.begin(), args.end(), std::string( "" ) );
  std::remove( args.begin(), args.end(), std::string( "" ) );
  int     argc = args.size();
  char* * argv = new char *[args.size() + 1];
  for( unsigned int i = 0; i < args.size(); ++i )
    {
    // allocate space for the string plus a null character
    argv[i] = new char[args[i].length() + 1];
    std::strncpy( argv[i], args[i].c_str(), args[i].length() );
    // place the null character in the end
    argv[i][args[i].length()] = '\0';
    }
  argv[argc] = 0;
  // class to automatically cleanup argv upon destruction
  class Cleanup_argv
  {
public:
    Cleanup_argv( char* * argv_, int argc_plus_one_ ) : argv( argv_ ), argc_plus_one( argc_plus_one_ )
    {
    }

    ~Cleanup_argv()
    {
      for( unsigned int i = 0; i < argc_plus_one; ++i )
        {
        delete[] argv[i];
        }
      delete[] argv;
    }

private:
    char* *      argv;
    unsigned int argc_plus_one;
  };
  Cleanup_argv cleanup_argv( argv, argc + 1 );

  antscout->set_stream( out_stream );

  if( argc < 2 )
    {
    UsageAndExit();
    }
  ++argv; --argc;
  std::string action(*argv);
  ++argv; --argc;
  if( argc == 0 )
    {
    antscout << "Missing CompositeTransformName" << std::endl;
    UsageAndExit();
    }

  RegisterMatOff<2>();
  RegisterMatOff<3>();

  std::string CompositeName(*argv);
  ++argv; --argc;
  if( action == "--disassemble" )
    {
    if( argc == 0 )
      {
      antscout << "Missing output transforms prefix" << std::endl;
      UsageAndExit();
      }
    std::string Prefix(*argv);
    return Disassemble(CompositeName, Prefix);
    }
  else if( action != "--assemble" )
    {
    antscout << "Unknown action " << action << std::endl;
    UsageAndExit();
    }
  std::vector<std::string> transformNames;
  while( argc > 0 )
    {
    std::string transformName(*argv);
    ++argv; --argc;
    transformNames.push_back(transformName);
    }

  if( transformNames.size() < 1 )
    {
    antscout << "Missing transform names to "
             << "assemble into a composite transform"
             << std::endl;
    UsageAndExit();
    }
  return Assemble(CompositeName, transformNames);
}
} // namespace ants
