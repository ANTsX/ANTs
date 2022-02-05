/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "antsUtilities.h"
#include <algorithm>

#include "itkImageFileReader.h"

#include "itkImageRegionIteratorWithIndex.h"

#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkFloatArray.h"
#include "vtkPoints.h"
#include "vtkPointData.h"

namespace ants
{
// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
ConvertVectorFieldToVTK(std::vector<std::string> args, std::ostream * out_stream = nullptr)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "ConvertVectorFieldToVTK");

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
  argv[argc] = 0;
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

  if (argc < 3)
  {
    std::cout << "Usage: " << argv[0]
              << " inputDisplacementField outputVTKFile maskImage(optional) slice(optional) whichAxis(optional)"
              << std::endl;
    return EXIT_FAILURE;
  }

  typedef float          PixelType;
  constexpr unsigned int ImageDimension = 3;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<int, ImageDimension>       MaskImageType;

  typedef double                                 RealType;
  typedef itk::Vector<RealType, ImageDimension>  VectorType;
  typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;

  typedef itk::ImageFileReader<DisplacementFieldType> ReaderType;
  ReaderType::Pointer                                 reader = ReaderType::New();
  reader->SetFileName(argv[1]);
  //  reader->SetUseAvantsNamingConvention( true );
  reader->Update();

  MaskImageType::Pointer mask;
  if (argc >= 4)
  {
    typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
    MaskReaderType::Pointer                     maskreader = MaskReaderType::New();
    maskreader->SetFileName(argv[3]);
    maskreader->Update();
    mask = maskreader->GetOutput();
  }
  else
  {
    // ORIENTATION ALERT  -- the original code here
    // set the region, spacing, and origin without setting directions.
    mask = AllocImage<MaskImageType>(reader->GetOutput(), 1);
  }

  int size[ImageDimension];
  int totalsize = 1;
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    size[i] = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[i];
    if (argc > 4 && std::stoi(argv[5]) == (int)i)
    {
      size[i] = 1;
    }
    totalsize *= size[i];
  }

  int totalPoints = totalsize;

  vtkUnstructuredGrid * field = vtkUnstructuredGrid::New();
  vtkPoints *           points = vtkPoints::New();
  points->Allocate(totalPoints);
  vtkFloatArray * vectors = vtkFloatArray::New();
  vectors->SetNumberOfComponents(3);
  vectors->SetNumberOfTuples(totalPoints);

  float x[3], v[3];
  int   offset = 0;

  itk::ImageRegionIteratorWithIndex<MaskImageType> It(mask, mask->GetLargestPossibleRegion());
  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    DisplacementFieldType::IndexType idx = It.GetIndex();

    if ((argc > 4 && idx[atoi(argv[5])] != std::stoi(argv[4])) || It.Get() == 0)
    {
      continue;
    }
    DisplacementFieldType::PointType point;
    reader->GetOutput()->TransformIndexToPhysicalPoint(idx, point);

    VectorType V = reader->GetOutput()->GetPixel(idx);
    for (unsigned int i = 0; i < ImageDimension; i++)
    {
      x[i] = point[i];
      v[i] = V[i];
    }
    //    offset = idx[0] + idx[1]*(size[1]+1) + idx[2]*(size[1]+1)*(size[3]+1);
    points->InsertPoint(offset, x);
    vectors->InsertTuple(offset++, v);
  }

  field->SetPoints(points);
  field->GetPointData()->SetVectors(vectors);

  points->Delete();
  vectors->Delete();

  vtkUnstructuredGridWriter * writer = vtkUnstructuredGridWriter::New();
  writer->SetInput(field);
  //  writer->SetFileTypeToBinary();
  writer->SetFileName(argv[2]);
  writer->Write();
  return 0;
}
} // namespace ants
