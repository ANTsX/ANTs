#include "itkImage.h"
#include "itkTransform.h"
#include "itkCompositeTransform.h"
#include "itkantsReadWriteTransform.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "antsUtilities.h"

namespace ants
{

class CompareFloat{
    //  Float           Memory          Bias (unsigned)
    //  -----           ------          ---------------
    //   NaN            0xFFFFFFFF      0xFF800001
    //   NaN            0xFF800001      0xFFFFFFFF
    //  -Infinity       0xFF800000      0x00000000 ---
    //  -3.40282e+038   0xFF7FFFFF      0x00000001    |
    //  -1.40130e-045   0x80000001      0x7F7FFFFF    |
    //  -0.0            0x80000000      0x7F800000    |--- Valid <= 0xFF000000.
    //   0.0            0x00000000      0x7F800000    |    NaN > 0xFF000000
    //   1.40130e-045   0x00000001      0x7F800001    |
    //   3.40282e+038   0x7F7FFFFF      0xFEFFFFFF    |
    //   Infinity       0x7F800000      0xFF000000 ---
    //   NaN            0x7F800001      0xFF000001
    //   NaN            0x7FFFFFFF      0xFF7FFFFF
    //
    //   Either value of NaN returns false.
    //   -Infinity and +Infinity are not "close".
    //   -0 and +0 are equal.

public:
  union{
    float          m_f32;
    unsigned int   m_u32;
  };

  static bool Is_Close( float A, float B, unsigned int unitsDelta = 4 )
  {
  unsigned int a = GetBiased( A );
  unsigned int b = GetBiased( B );

  if( (a > 0xFF000000) || (b > 0xFF000000) )
    {
    return false;
    }
  return( (static_cast<unsigned int>(abs( long(a - b) ))) < unitsDelta );
  }

protected:

  static unsigned int GetBiased( float f )
  {
  unsigned int r = ((CompareFloat*)&f)->m_u32;

  if( r & 0x80000000 )
    {
    return( ~r - 0x007FFFFF );
    }
  return( r + 0x7F800000 );
  }

};

template <unsigned int VImageDimension>
int compareComposites( const typename itk::Transform<double, VImageDimension, VImageDimension>::Pointer & firstTransform,
                       const typename itk::Transform<double, VImageDimension, VImageDimension>::Pointer & secondTransform )
{
  typedef typename itk::CompositeTransform<double, VImageDimension>            CompositeTransformType;
  typedef typename CompositeTransformType::ScalarType                          RealType;

  typedef typename itk::DisplacementFieldTransform<RealType, VImageDimension>  DisplacementFieldTransformType;
  typedef typename DisplacementFieldTransformType::DisplacementFieldType       DisplacementFieldType;

  const std::string CompositeTransformID ("CompositeTransform");

  if( firstTransform->GetNameOfClass() == CompositeTransformID && secondTransform->GetNameOfClass() == CompositeTransformID )
    {
    const typename CompositeTransformType::ConstPointer Comp1 =
                                      dynamic_cast<const itk::CompositeTransform<double, VImageDimension> *>( firstTransform.GetPointer() );
    const typename CompositeTransformType::ConstPointer Comp2 =
                                      dynamic_cast<const itk::CompositeTransform<double, VImageDimension> *>( secondTransform.GetPointer() );

    if( Comp1->GetNumberOfTransforms() == Comp2->GetNumberOfTransforms() )
      {
      const unsigned int N = Comp1->GetNumberOfTransforms();
      for(unsigned int i = 0; i < N; ++i)
         {
         if( Comp1->GetNthTransform(i)->GetNameOfClass() == Comp2->GetNthTransform(i)->GetNameOfClass() )
           {
           if( strcmp( Comp1->GetNthTransform(i)->GetNameOfClass(), "DisplacementFieldTransform" ) )
             {
             if( Comp1->GetNthTransform(i)->GetFixedParameters() != Comp2->GetNthTransform(i)->GetFixedParameters() ||
                 Comp1->GetNthTransform(i)->GetParameters() != Comp2->GetNthTransform(i)->GetParameters() )
               {
               std::cerr << "The input composite transforms are not equal! The transform parameters are different!" << std::endl;
               return EXIT_FAILURE;
               }
             }
           else
             {
             // compare two displacement field transforms by considering a tolerance
             const typename DisplacementFieldTransformType::ConstPointer DispTrans1 = dynamic_cast<DisplacementFieldTransformType *>( Comp1->GetNthTransform(i).GetPointer() );
             const typename DisplacementFieldTransformType::ConstPointer DispTrans2 = dynamic_cast<DisplacementFieldTransformType *>( Comp2->GetNthTransform(i).GetPointer() );

             const typename DisplacementFieldType::ConstPointer DispField1 = DispTrans1->GetDisplacementField();
             const typename DisplacementFieldType::ConstPointer DispField2 = DispTrans2->GetDisplacementField();

             typedef itk::ImageRegionConstIteratorWithIndex<DisplacementFieldType> DispIteratorType;
             DispIteratorType dit1( DispField1, DispField1->GetLargestPossibleRegion() );
             DispIteratorType dit2( DispField2, DispField2->GetLargestPossibleRegion() );

             dit1.GoToBegin();
             dit2.GoToBegin();

             while( !dit1.IsAtEnd() && !dit2.IsAtEnd() )
                  {
                  typename DisplacementFieldType::PixelType v1 = dit1.Get();
                  typename DisplacementFieldType::PixelType v2 = dit2.Get();
                  for (unsigned int index=0; index<VImageDimension; ++index) {
                    if( !CompareFloat::Is_Close( v1[index], v2[index]) ) // Compares two float numbers.
                      {
                      std::cerr << "The input composite transforms are not equal! The diffeomorphic transform parameters are different!" << std::endl;
                      return EXIT_FAILURE;
                      }
                  }
                  ++dit1; ++dit2;
                  }
             }
           }
         else
           {
           std::cerr << "The input composite transforms are not equal! They contain different types of transform!" << std::endl;
           return EXIT_FAILURE;
           }
         }
      }
    else
      {
      std::cerr << "The input composite transforms are not equal! They do not have the same number of transforms!" << std::endl;
      return EXIT_FAILURE;
      }
    }
  else
    {
    std::cerr << "The input transforms MUST be composite transforms." << std::endl;
    return EXIT_FAILURE;
    }
  std::cout << "Two input composite transforms are the same!" << std::endl;
  return EXIT_SUCCESS;
}

int compareTwoCompositeTransforms( std::vector<std::string> args, std::ostream* /*out_stream = NULL */ )
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

  // antscout->set_stream( out_stream );
  {
  itk::Transform<double, 2, 2>::Pointer firstTransform = itk::ants::ReadTransform<double, 2>(args[0]);
  itk::Transform<double, 2, 2>::Pointer secondTransform = itk::ants::ReadTransform<double, 2>(args[1]);
  if( firstTransform.IsNotNull() && secondTransform.IsNotNull() )
    {
    return compareComposites<2>(firstTransform, secondTransform);
    }
  }
  {
  itk::Transform<double, 3, 3>::Pointer firstTransform = itk::ants::ReadTransform<double, 3>(args[0]);
  itk::Transform<double, 3, 3>::Pointer secondTransform = itk::ants::ReadTransform<double, 3>(args[1]);
  if( firstTransform.IsNotNull() && secondTransform.IsNotNull() )
    {
    return compareComposites<3>(firstTransform, secondTransform);
    }
  }
  std::cerr << "Can't read input transforms" << std::endl;
  return EXIT_FAILURE;
}
} // namespace ants
