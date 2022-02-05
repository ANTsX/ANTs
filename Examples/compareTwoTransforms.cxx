#include "itkImage.h"
#include "itkTransform.h"
#include "itkCompositeTransform.h"
#include "itkantsReadWriteTransform.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkBSplineTransform.h"
#include <vcl_compiler.h>
#include <iostream>
#include <algorithm>

#include "antsUtilities.h"

namespace ants
{

template <unsigned int VImageDimension>
int
compareTransforms(const typename itk::Transform<double, VImageDimension, VImageDimension>::Pointer & firstTransform,
                  const typename itk::Transform<double, VImageDimension, VImageDimension>::Pointer & secondTransform)
{
  using CompositeTransformType = typename itk::CompositeTransform<double, VImageDimension>;
  using RealType = typename CompositeTransformType::ScalarType;

  using DisplacementFieldTransformType = typename itk::DisplacementFieldTransform<RealType, VImageDimension>;
  using DisplacementFieldType = typename DisplacementFieldTransformType::DisplacementFieldType;

  const std::string CompositeTransformID("CompositeTransform");

  if (firstTransform->GetNameOfClass() == CompositeTransformID &&
      secondTransform->GetNameOfClass() == CompositeTransformID)
  {
    std::cout << "The input transforms are composite transform types." << std::endl;

    const typename CompositeTransformType::ConstPointer Comp1 =
      dynamic_cast<const itk::CompositeTransform<double, VImageDimension> *>(firstTransform.GetPointer());
    const typename CompositeTransformType::ConstPointer Comp2 =
      dynamic_cast<const itk::CompositeTransform<double, VImageDimension> *>(secondTransform.GetPointer());

    if (Comp1->GetNumberOfTransforms() == Comp2->GetNumberOfTransforms())
    {
      const float        coordinateTolerance = 1e-0;
      const unsigned int N = Comp1->GetNumberOfTransforms();
      for (unsigned int i = 0; i < N; ++i)
      {
        if (Comp1->GetNthTransform(i)->GetNameOfClass() == Comp2->GetNthTransform(i)->GetNameOfClass())
        {
          if (strcmp(Comp1->GetNthTransform(i)->GetNameOfClass(), "DisplacementFieldTransform"))
          {
            if (!Comp1->GetNthTransform(i)->GetFixedParameters().is_equal(
                  Comp2->GetNthTransform(i)->GetFixedParameters(), coordinateTolerance) ||
                !Comp1->GetNthTransform(i)->GetParameters().is_equal(Comp2->GetNthTransform(i)->GetParameters(),
                                                                     coordinateTolerance))
            {
              std::cerr << i << ": " << std::endl;
              std::cerr << Comp1->GetNthTransform(i)->GetParameters() << std::endl;
              std::cerr << Comp2->GetNthTransform(i)->GetParameters() << std::endl;
              std::cerr << "The input composite transforms are not equal! The transform parameters are different!"
                        << std::endl;
              return EXIT_FAILURE;
            }
          }
          else
          {
            // compare two displacement field transforms by considering a tolerance
            const typename DisplacementFieldTransformType::ConstPointer DispTrans1 =
              dynamic_cast<DisplacementFieldTransformType *>(Comp1->GetNthTransform(i).GetPointer());
            const typename DisplacementFieldTransformType::ConstPointer DispTrans2 =
              dynamic_cast<DisplacementFieldTransformType *>(Comp2->GetNthTransform(i).GetPointer());

            const typename DisplacementFieldType::ConstPointer DispField1 = DispTrans1->GetDisplacementField();
            const typename DisplacementFieldType::ConstPointer DispField2 = DispTrans2->GetDisplacementField();

            using DispIteratorType = itk::ImageRegionConstIteratorWithIndex<DisplacementFieldType>;
            DispIteratorType dit1(DispField1, DispField1->GetLargestPossibleRegion());
            DispIteratorType dit2(DispField2, DispField2->GetLargestPossibleRegion());

            dit1.GoToBegin();
            dit2.GoToBegin();

            while (!dit1.IsAtEnd() && !dit2.IsAtEnd())
            {
              typename DisplacementFieldType::PixelType v1 = dit1.Get();
              typename DisplacementFieldType::PixelType v2 = dit2.Get();
              for (unsigned int index = 0; index < VImageDimension; ++index)
              {
                if (!itk::Math::FloatAlmostEqual(v1[index], v2[index], 4, 20)) // Compares two float numbers.
                {
                  std::cerr << v1[index] << " != " << v2[index] << std::endl;
                  std::cerr << "The input composite transforms are not equal! The diffeomorphic transform parameters "
                               "are different!"
                            << std::endl;
                  return EXIT_FAILURE;
                }
              }
              ++dit1;
              ++dit2;
            }
          }
        }
        else
        {
          std::cerr << "The input composite transforms are not equal! They contain different types of transform!"
                    << std::endl;
          return EXIT_FAILURE;
        }
      }
    }
    else
    {
      std::cerr << "The input composite transforms are not equal! They do not have the same number of transforms!"
                << std::endl;
      return EXIT_FAILURE;
    }
  }
  else
  {
    std::cout << "The input transforms are not composite transforms." << std::endl;
    const std::string DisplacementFieldTransformID("DisplacementFieldTransform");

    if (firstTransform->GetNameOfClass() == DisplacementFieldTransformID &&
        secondTransform->GetNameOfClass() == DisplacementFieldTransformID)
    {
      std::cout << "The input transforms are displacement field transform type." << std::endl;

      // compare two displacement field transforms by considering a tolerance
      const typename DisplacementFieldTransformType::ConstPointer DispTrans1 =
        dynamic_cast<DisplacementFieldTransformType *>(firstTransform.GetPointer());
      const typename DisplacementFieldTransformType::ConstPointer DispTrans2 =
        dynamic_cast<DisplacementFieldTransformType *>(secondTransform.GetPointer());

      const typename DisplacementFieldType::ConstPointer DispField1 = DispTrans1->GetDisplacementField();
      const typename DisplacementFieldType::ConstPointer DispField2 = DispTrans2->GetDisplacementField();

      using DispIteratorType = itk::ImageRegionConstIteratorWithIndex<DisplacementFieldType>;
      DispIteratorType dit1(DispField1, DispField1->GetLargestPossibleRegion());
      DispIteratorType dit2(DispField2, DispField2->GetLargestPossibleRegion());

      dit1.GoToBegin();
      dit2.GoToBegin();

      while (!dit1.IsAtEnd() && !dit2.IsAtEnd())
      {
        typename DisplacementFieldType::PixelType v1 = dit1.Get();
        typename DisplacementFieldType::PixelType v2 = dit2.Get();
        for (unsigned int index = 0; index < VImageDimension; ++index)
        {
          if (!itk::Math::FloatAlmostEqual(v1[index], v2[index], 4, 20)) // Compares two float numbers.
          {
            std::cerr << v1[index] << " != " << v2[index] << std::endl;
            std::cerr << "The input displacement field transforms are not equal! The diffeomorphic transform "
                         "parameters are different!"
                      << std::endl;
            return EXIT_FAILURE;
          }
        }
        ++dit1;
        ++dit2;
      }
    }
    else
    {
      std::cout << "The input transforms are neither composite transform nor displacement field transform."
                << std::endl;
      std::cout << "First Transform Type: " << firstTransform->GetNameOfClass() << std::endl;
      std::cout << "Second Transform Type: " << secondTransform->GetNameOfClass() << std::endl;
      if (firstTransform->GetFixedParameters() != secondTransform->GetFixedParameters() ||
          firstTransform->GetParameters() != secondTransform->GetParameters())
      {
        std::cerr << "The input transforms are not equal! The transform parameters are different!" << std::endl;
        return EXIT_FAILURE;
      }
    }
  }
  std::cout << "Two input transforms are the same!" << std::endl;
  return EXIT_SUCCESS;
}

int
compareTwoTransforms(std::vector<std::string> args, std::ostream * /* out_stream = nullptr */)
{
  // the arguments coming in as 'args' is a replacement for the standard (argc,argv) format
  // Just notice that the argv[i] equals to args[i-1]
  // and the argc equals:
  int argc = args.size() + 1;

  if (argc != 3)
  {
    std::cerr << "Usage: compareTwoTransforms\n"
              << "<First Transform> , <Second Transform>" << std::endl;
    return EXIT_FAILURE;
  }

  {
    itk::Transform<double, 2, 2>::Pointer firstTransform = itk::ants::ReadTransform<double, 2>(args[0]);
    itk::Transform<double, 2, 2>::Pointer secondTransform = itk::ants::ReadTransform<double, 2>(args[1]);
    if (firstTransform.IsNotNull() && secondTransform.IsNotNull())
    {
      using BSplineTransformType = itk::BSplineTransform<double, 2, 2>;
      BSplineTransformType::Pointer bsplineInput1 = dynamic_cast<BSplineTransformType *>(firstTransform.GetPointer());
      BSplineTransformType::Pointer bsplineInput2 = dynamic_cast<BSplineTransformType *>(secondTransform.GetPointer());
      if (bsplineInput1.IsNull() && bsplineInput2.IsNull())
      {
        return compareTransforms<2>(firstTransform, secondTransform);
      }
      else
      {
        std::cerr << "BSpline transform type is not supported." << std::endl;
        return EXIT_FAILURE;
      }
    }
  }
  {
    itk::Transform<double, 3, 3>::Pointer firstTransform = itk::ants::ReadTransform<double, 3>(args[0]);
    itk::Transform<double, 3, 3>::Pointer secondTransform = itk::ants::ReadTransform<double, 3>(args[1]);
    if (firstTransform.IsNotNull() && secondTransform.IsNotNull())
    {
      using BSplineTransformType = itk::BSplineTransform<double, 3, 3>;
      BSplineTransformType::Pointer bsplineInput1 = dynamic_cast<BSplineTransformType *>(firstTransform.GetPointer());
      BSplineTransformType::Pointer bsplineInput2 = dynamic_cast<BSplineTransformType *>(secondTransform.GetPointer());
      if (bsplineInput1.IsNull() && bsplineInput2.IsNull())
      {
        return compareTransforms<3>(firstTransform, secondTransform);
      }
      else
      {
        std::cerr << "BSpline transform type is not supported." << std::endl;
        return EXIT_FAILURE;
      }
    }
  }
  std::cerr << "Can't read input transforms" << std::endl;
  return EXIT_FAILURE;
}
} // namespace ants
