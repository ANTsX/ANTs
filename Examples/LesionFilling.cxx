#include <iostream>
#include <ostream>
#include <sstream>
#include <vector>
#include <iomanip>

#include "antsUtilities.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkMaskImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageRegionIterator.h"
// LesionFilling dimension t1.nii.gz lesionmask output.nii.gz
namespace ants
{
template <unsigned int ImageDimension>
int
LesionFilling(int argc, char * argv[])
{
  const char * T1FileName = argv[2];
  const char * LesionMapFileName = argv[3];
  const char * OutputFileName = argv[4];
  if (argc < 3)
  {
    std::cout << "Missing arguments, see usage." << std::endl;
    throw;
  }
  if (argc < 4)
  {
    std::cout << "2 more arguments are necessary, see usage." << std::endl;
    throw;
  }
  if (argc < 5)
  {
    std::cout << "Missing output filename." << std::endl;
    throw;
  }
  using T1ImageType = itk::Image<double, ImageDimension>;
  using LesionImageType = itk::Image<unsigned short, ImageDimension>;
  using T1ImageReaderType = itk::ImageFileReader<T1ImageType>;
  using LesionImageReaderType = itk::ImageFileReader<LesionImageType>;
  typename LesionImageReaderType::Pointer LesionReader = LesionImageReaderType::New();
  LesionReader->SetFileName(LesionMapFileName);
  try
  {
    LesionReader->Update();
  }
  catch (const itk::ExceptionObject & itkNotUsed(excp))
  {
    std::cout << "no lesion mask that can be read" << std::endl;
    return 0;
  }
  typename T1ImageReaderType::Pointer T1Reader = T1ImageReaderType::New();
  T1Reader->SetFileName(T1FileName);
  try
  {
    T1Reader->Update();
  }
  catch (const itk::ExceptionObject & itkNotUsed(excp))
  {
    std::cout << "no T1 image that can be read" << std::endl;
    return 0;
  }
  typename T1ImageType::Pointer outImage = nullptr;
  outImage = T1Reader->GetOutput();
  using IteratorType = itk::ImageRegionIterator<T1ImageType>;
  using BinaryThresholdImageFilterType = itk::BinaryThresholdImageFilter<T1ImageType, T1ImageType>;
  using StructuringElementType = itk::BinaryBallStructuringElement<double, ImageDimension>;
  using DilateFilterType = itk::BinaryDilateImageFilter<T1ImageType, T1ImageType, StructuringElementType>;
  using RelabelComponentFilterType = itk::RelabelComponentImageFilter<LesionImageType, LesionImageType>;
  using ConnectedComponentFilterType = itk::ConnectedComponentImageFilter<LesionImageType, LesionImageType>;

  typename ConnectedComponentFilterType::Pointer connected = ConnectedComponentFilterType::New();
  typename RelabelComponentFilterType::Pointer   relabeled = RelabelComponentFilterType::New();
  connected->SetInput(LesionReader->GetOutput());
  connected->Update();
  relabeled->SetInput(connected->GetOutput());
  relabeled->Update();

  const int LesionNumber = relabeled->GetNumberOfObjects();
  std::cout << "Number of lesions: " << LesionNumber << std::endl;
  for (int i = 1; i < LesionNumber + 1; i++)
  {
    using FilterType = itk::CastImageFilter<LesionImageType, T1ImageType>;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput(relabeled->GetOutput());
    filter->Update();
    typename BinaryThresholdImageFilterType::Pointer thresholdFilter = BinaryThresholdImageFilterType::New();
    thresholdFilter->SetInput(filter->GetOutput());
    thresholdFilter->SetLowerThreshold((double)i - 0.1);
    thresholdFilter->SetUpperThreshold((double)i + 0.1);
    thresholdFilter->SetInsideValue(1);
    thresholdFilter->SetOutsideValue(0);
    thresholdFilter->Update();
    // Neighbouring voxel
    // filling lesions with the voxels surrounding them
    // first finding the edges of lesions
    // by subtracting dilated lesion map from lesion map itself
    typename DilateFilterType::Pointer binaryDilate = DilateFilterType::New();

    int elementRadius = 1;

    // dilate more around very small lesions to improve the chance of finding
    // suitable replacement voxels
    // index is offset because background is not counted
    if (relabeled->GetSizeOfObjectsInPixels()[i - 1] < 5)
    {
      elementRadius = 2;
    }

    StructuringElementType structuringElement;
    structuringElement.SetRadius(elementRadius); // 3x3 structuring element
    structuringElement.CreateStructuringElement();
    binaryDilate->SetKernel(structuringElement);
    binaryDilate->SetInput(thresholdFilter->GetOutput());
    binaryDilate->SetDilateValue(1);
    binaryDilate->Update();
    // subtract dilated image form non-dilated one
    using SubtractImageFilterType = itk::SubtractImageFilter<T1ImageType, T1ImageType>;
    typename SubtractImageFilterType::Pointer subtractFilter = SubtractImageFilterType::New();
    // output = image1 - image2
    subtractFilter->SetInput1(binaryDilate->GetOutput());
    subtractFilter->SetInput2(thresholdFilter->GetOutput());
    subtractFilter->Update();
    // multiply the outer lesion mask with T1 to get only the neighbouring voxels
    using MaskFilterType = itk::MaskImageFilter<T1ImageType, T1ImageType>;
    typename MaskFilterType::Pointer maskFilter = MaskFilterType::New();
    maskFilter->SetInput(outImage);
    maskFilter->SetMaskImage(subtractFilter->GetOutput());
    maskFilter->Update();
    typename T1ImageType::Pointer LesionEdge = nullptr;
    LesionEdge = maskFilter->GetOutput();

    // calculating mean lesion intesity
    // Note: lesions should not be filled with values
    // less than their originial values, this is a
    // trick to exclude any CSF voxels in surronding voxels (if any)
    typename MaskFilterType::Pointer maskFilterLesion = MaskFilterType::New();
    maskFilterLesion->SetInput(T1Reader->GetOutput());
    // do we need to pass a double lesion
    maskFilterLesion->SetMaskImage(thresholdFilter->GetOutput());
    maskFilterLesion->Update();
    IteratorType it(maskFilterLesion->GetOutput(), maskFilterLesion->GetOutput()->GetLargestPossibleRegion());
    it.GoToBegin();
    /** Walk over the image. */
    int    counter = 0;
    double meanInsideLesion = 0;
    while (!it.IsAtEnd())
    {
      if (!itk::Math::FloatAlmostEqual(it.Value(), 0.0))
      {
        // counting number of voxels inside lesion
        counter++;
        meanInsideLesion += it.Get();
      }
      ++it;
    }
    if (counter > 0)
    {
      meanInsideLesion /= (double)counter;
    }
    else
    {
      std::cerr << "Intensity in lesion " << i << " of " << LesionNumber << " is zero, will not fill" << std::endl;
      continue;
    }

    // check that all outer voxels are more than the mean
    // intensity of the lesion, i.e. not including CSF voxels
    IteratorType itNoCSF(maskFilter->GetOutput(), maskFilter->GetOutput()->GetLargestPossibleRegion());
    itNoCSF.GoToBegin();
    std::vector<double> outerWMVoxels;
    while (!itNoCSF.IsAtEnd())
    {
      if (itNoCSF.Get() >= meanInsideLesion)
      {
        outerWMVoxels.push_back(itNoCSF.Get());
      } // end if
      ++itNoCSF;
    } // end while
    // walk through original T1
    // and change inside the lesion with a random pick from
    // collected normal appearing WM voxels (outerWMVoxels)
    IteratorType it4(outImage, outImage->GetLargestPossibleRegion());
    IteratorType itL(thresholdFilter->GetOutput(), thresholdFilter->GetOutput()->GetLargestPossibleRegion());
    int          max = outerWMVoxels.size();

    if (max == 0)
    {
      std::cerr << "Intensity surrounding lesion " << i << " of " << LesionNumber
                << " is less than mean lesion intensity, will not fill" << std::endl;
      continue;
    }

    int min = 0;
    it4.GoToBegin();
    itL.GoToBegin();
    while (!it4.IsAtEnd())
    {
      if (itk::Math::FloatAlmostEqual(itL.Get(), 1.0))
      {
        int index = min + (rand() % (int)(max - min));
        it4.Set(outerWMVoxels[index]);
      } // end if
      ++it4;
      ++itL;
    } // end while
  }   // end of loop for lesions
  using WriterType = itk::ImageFileWriter<T1ImageType>;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput(outImage);
  writer->SetFileName(OutputFileName);
  try
  {
    writer->Update();
  }
  catch (const itk::ExceptionObject & err)
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
} // main int

int
LesionFilling(std::vector<std::string> args, std::ostream * itkNotUsed(out_stream))
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "LesionFilling");

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

  // antscout->set_stream( out_stream );

  // LesionFilling dimension t1.nii.gz lesionmask output.nii.gz
  if (argc < 3)
  {
    std::cout << "Example usage: " << argv[0]
              << " imageDimension T1_image.nii.gz lesion_mask.nii.gz output_lesion_filled.nii.gz" << std::endl;

    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  switch (std::stoi(argv[1]))
  {
    case 2:
    {
      return LesionFilling<2>(argc, argv);
    }
    break;
    case 3:
    {
      return LesionFilling<3>(argc, argv);
    }
    break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
} // int LesionFilling std::vector
} // namespace ants
