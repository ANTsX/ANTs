/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include "antsUtilities.h"

#include "antsCommandLineParser.h"
#include "itkTransformFileReader.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkCompositeTransform.h"


#include <sstream>

namespace ants
{


// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
antsTransformInfo(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "antsTransformInfo");

  int     argc = args.size();
  char ** argv = new char *[args.size() + 1];
  for (unsigned int i = 0; i < args.size(); ++i)
  {
    // allocate space for the string plus a null character
    argv[i] = new char[args[i].length() + 1];
    std::strncpy(argv[i], args[i].c_str(), args[i].length());
    // place the null character in the end
    argv[i][args[i].length()] = '\0';
  }
  argv[argc] = nullptr;
  // class to automatically cleanup argv upon destruction
  class Cleanup_argv
  {
  public:
    Cleanup_argv(char ** argv_, int argc_plus_one_)
      : argv(argv_)
      , argc_plus_one(argc_plus_one_)
    {}

    ~Cleanup_argv()
    {
      for (unsigned int i = 0; i < argc_plus_one; ++i)
      {
        delete[] argv[i];
      }
      delete[] argv;
    }

  private:
    char **      argv;
    unsigned int argc_plus_one;
  };
  Cleanup_argv cleanup_argv(argv, argc + 1);

  using ReadScalarType = double;

  for (int i = 1; i < argc; i++)
  {

    std::cout << "Transform file: " << argv[i] << std::endl;

    using TransformReaderType = itk::TransformFileReaderTemplate<ReadScalarType>;
    auto reader = TransformReaderType::New();

    reader->SetFileName(argv[i]);
    try
    {
      reader->Update();
    }
    catch (const itk::ExceptionObject & excp)
    {
      std::cerr << "Error while reading the transform file" << std::endl;
      std::cerr << excp << std::endl;
      std::cerr << "[FAILED]" << std::endl;
      return EXIT_FAILURE;
    }

    const TransformReaderType::TransformListType * transforms = reader->GetTransformList();
    std::cout << "Number of transforms = " << transforms->size() << std::endl;

    unsigned int Dimension = transforms->front()->GetInputSpaceDimension();


    for (auto it = transforms->begin(); it != transforms->end(); ++it)
    {

      if (Dimension == 3)
      {
        using ReadCompositeTransformType3D = itk::CompositeTransform<ReadScalarType, 3>;
        if (!strcmp((*it)->GetNameOfClass(), "CompositeTransform"))
        {
          ReadCompositeTransformType3D::Pointer compositeRead =
            static_cast<ReadCompositeTransformType3D *>((*it).GetPointer());
          const ReadCompositeTransformType3D::TransformQueueType & compositeTransformQueue =
            compositeRead->GetTransformQueue();
          std::cout << "Number of transforms in composite queue = " << compositeTransformQueue.size() << std::endl;

          for (auto it2 = compositeTransformQueue.begin(); it2 != compositeTransformQueue.end(); ++it2)
          {

            using TransformType3D = itk::MatrixOffsetTransformBase<double, 3, 3>;
            TransformType3D * itktx3d = dynamic_cast<TransformType3D *>((*it2).GetPointer());

            if (itktx3d)
            {
              itktx3d->Print(std::cout);
              std::cout << "Determinant: " << vnl_determinant(itktx3d->GetMatrix().GetVnlMatrix()) << std::endl;
            }
            else
            {
              it2->Print(std::cout);
            }
          }
        }
        else
        {
          using TransformType3D = itk::MatrixOffsetTransformBase<double, 3, 3>;
          TransformType3D * itktx3d = dynamic_cast<TransformType3D *>((*it).GetPointer());
          if (itktx3d)
          {
              itktx3d->Print(std::cout);
              std::cout << "Determinant: " << vnl_determinant(itktx3d->GetMatrix().GetVnlMatrix()) << std::endl;
          }
          else
          {
              (*it)->Print(std::cout);
          }
        }
      }
      else if (Dimension == 2)
      {
        using ReadCompositeTransformType2D = itk::CompositeTransform<ReadScalarType, 2>;
        if (!strcmp((*it)->GetNameOfClass(), "CompositeTransform"))
        {
          ReadCompositeTransformType2D::Pointer compositeRead =
            static_cast<ReadCompositeTransformType2D *>((*it).GetPointer());
          const ReadCompositeTransformType2D::TransformQueueType & compositeTransformQueue =
            compositeRead->GetTransformQueue();
          std::cout << "Number of transforms in composite queue = " << compositeTransformQueue.size() << std::endl;

          for (auto it2 = compositeTransformQueue.begin(); it2 != compositeTransformQueue.end(); ++it2)
          {

            using TransformType2D = itk::MatrixOffsetTransformBase<double, 2, 2>;
            TransformType2D * itktx2d = dynamic_cast<TransformType2D *>((*it2).GetPointer());

            if (itktx2d)
            {
              itktx2d->Print(std::cout);
              std::cout << "Determinant: " << vnl_determinant(itktx2d->GetMatrix().GetVnlMatrix()) << std::endl;
            }
            else
            {
              it2->Print(std::cout);
            }
          }
        }
        else
        {
          using TransformType2D = itk::MatrixOffsetTransformBase<double, 2, 2>;
          TransformType2D * itktx2d = dynamic_cast<TransformType2D *>((*it).GetPointer());
          if (itktx2d)
          {
            itktx2d->Print(std::cout);
            std::cout << "Determinant: " << vnl_determinant(itktx2d->GetMatrix().GetVnlMatrix()) << std::endl;
          }
          else
          {
            (*it)->Print(std::cout);
          }
        }
      }
    }
  }

  return EXIT_SUCCESS;
}

} // namespace ants
