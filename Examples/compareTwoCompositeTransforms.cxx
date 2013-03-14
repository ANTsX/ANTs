#include "antsUtilities.h"
#include "itkTransform.h"
#include "itkCompositeTransform.h"
#include "itkantsReadWriteTransform.h"

namespace ants
{
template <unsigned int VImageDimension>
int compareComposites( const typename itk::Transform<double, VImageDimension, VImageDimension>::Pointer & firstTransform,
                       const typename itk::Transform<double, VImageDimension, VImageDimension>::Pointer & secondTransform )
{
  static const std::string CompositeTransformID("CompositeTransform");
  if( firstTransform->GetNameOfClass() == CompositeTransformID && secondTransform->GetNameOfClass() == CompositeTransformID )
    {
    const typename itk::CompositeTransform<double, VImageDimension>::ConstPointer Comp1 =
                                      dynamic_cast<const itk::CompositeTransform<double, VImageDimension> *>( firstTransform.GetPointer() );
    const typename itk::CompositeTransform<double, VImageDimension>::ConstPointer Comp2 =
                                      dynamic_cast<const itk::CompositeTransform<double, VImageDimension> *>( secondTransform.GetPointer() );
    if( Comp1 != Comp2 )
      {
      std::cerr << "The input composite transforms are not equal!" << std::endl;
      return EXIT_FAILURE;
      }
    }
  else
    {
    std::cerr << "The input transforms MUST be composite transforms" << std::endl;
    return EXIT_FAILURE;
    }
  antscout << "Two input composite transforms are the same!" << std::endl;
  return EXIT_SUCCESS;
}

int compareTwoCompositeTransforms( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // the arguments coming in as 'args' is a replacement for the standard (argc,argv) format
  // Just notice that the argv[i] equals to args[i-1]
  // and the argc equals:
  int argc = args.size() + 1;
  
  if( argc != 3 )
    {
    std::cerr << "Usage: compareTwoCompositeTransform\n"
    <<
    "<First Composite Transform> , <Second Composite Transform>"
    << std::endl;
    return EXIT_FAILURE;
    }
  
  antscout->set_stream( out_stream );
  {
  itk::Transform<double, 2, 2>::Pointer firstTransform = itk::ants::ReadTransform<2>(args[0]);
  itk::Transform<double, 2, 2>::Pointer secondTransform = itk::ants::ReadTransform<2>(args[1]);
  if( firstTransform.IsNotNull() && secondTransform.IsNotNull() )
    {
    return compareComposites<2>(firstTransform, secondTransform);
    }
  }
  {
  itk::Transform<double, 3, 3>::Pointer firstTransform = itk::ants::ReadTransform<3>(args[0]);
  itk::Transform<double, 3, 3>::Pointer secondTransform = itk::ants::ReadTransform<3>(args[1]);
  if( firstTransform.IsNotNull() && secondTransform.IsNotNull() )
    {
    return compareComposites<3>(firstTransform, secondTransform);
    }
  }
  std::cerr << "Can't read input transforms" << std::endl;
  return EXIT_FAILURE;
}
} // namespace ants





