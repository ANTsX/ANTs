/*=========================================================================

  Program:   ITK General
  Module:    $RCSfile: StackSlices.cxx,v $
  Language:  C++
  Date:      $Date: 2009/04/02 19:55:26 $
  Version:   $Revision: 1.22 $

  Author: Jeffrey T. Duda (jtduda@seas.upenn.edu)
  Institution: PICSL

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

=========================================================================*/
#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "ReadWriteImage.h"
/* FlipScalarVolume
 * This program takes a volume and flips it along the
 * indicated axes
 */
int main( int argc, char *argv[] )
{
  // Pixel and Image typedefs
  typedef float PixelType;

  typedef itk::Image<PixelType, 3> ImageType;
  typedef itk::Image<PixelType, 2> SliceType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::ImageFileWriter<ImageType> WriterType;

  typedef itk::ExtractImageFilter<ImageType, SliceType> ExtractFilterType;

  typedef itk::ImageRegionIteratorWithIndex<ImageType> ImageIt;
  typedef itk::ImageRegionIteratorWithIndex<SliceType> SliceIt;

  // Check for valid input paramters
  if( argc < 5 )
    {
    std::cout << "Usage: " << argv[0] << " outputvolume x y z inputvolumes" << std::endl;
    return 1;
    }

  char * stackName = argv[1];
  int    dimVars[3];
  dimVars[0] = atoi(argv[2]);
  dimVars[1] = atoi(argv[3]);
  dimVars[2] = atoi(argv[4]);

  int dim = -1;
  int slice = -1;
  // which dim to extract slice from
  for( unsigned int i = 0; i < 3; i++ )
    {
    if( dimVars[i] > -1 )
      {
      if( (dim > -1) || (slice > -1) )
        {
        std::cerr << "Can only choose slice from 1 dimension" << std::endl;
        return EXIT_FAILURE;
        }
      dim = i;
      slice = dimVars[i];
      }
    }

  unsigned long nSlices = argc - 5;
  // std::cout << nSlices << std::endl;

  ReaderType::Pointer firstReader = ReaderType::New();
  firstReader->SetFileName( argv[5] );
  firstReader->Update();
  std::cout << " Slice 0 :: " << std::string(argv[5]) << std::endl;

  ImageType::Pointer    stack = ImageType::New();
  ImageType::RegionType region = firstReader->GetOutput()->GetLargestPossibleRegion();
  region.SetSize(dim, nSlices);
  stack->SetRegions(region);
  stack->Allocate();

  SliceType::SizeType size;
  size.Fill(0);
  SliceType::IndexType index2d;
  if( dim == 0 )
    {
    size[0] = firstReader->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
    size[1] = firstReader->GetOutput()->GetLargestPossibleRegion().GetSize()[2];
    }
  if( dim == 1 )
    {
    size[0] = firstReader->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
    size[1] = firstReader->GetOutput()->GetLargestPossibleRegion().GetSize()[2];
    }
  if( dim == 2 )
    {
    size[0] = firstReader->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
    size[1] = firstReader->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
    }
  SliceType::Pointer    stack2 = SliceType::New();
  SliceType::RegionType region2;
  region2.SetSize(size);
  stack2->SetRegions(region2);
  stack2->Allocate();
  //  std::cout << region << std::endl;

  ImageType::RegionType extractRegion = stack->GetLargestPossibleRegion();
  extractRegion.SetSize(dim, 0);
  extractRegion.SetIndex(dim, slice);

  ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
  extractFilter->SetInput( firstReader->GetOutput() );
  extractFilter->SetDirectionCollapseToIdentity();
  extractFilter->SetExtractionRegion( extractRegion );
  extractFilter->Update();

  SliceIt it1(extractFilter->GetOutput(), extractFilter->GetOutput()->GetLargestPossibleRegion() );
  float   mean = 0.0, ct = 1;
  while( !it1.IsAtEnd() )
    {
    float val = it1.Value();
    if( val > 0 )
      {
      mean += it1.Value();
      ct++;
      }
    ++it1;
    }

  mean /= ct;    std::cout << " Mean " << mean << std::endl;
  mean = 1;
  it1.GoToBegin();
  while( !it1.IsAtEnd() )
    {
    ImageType::IndexType index;
    index.Fill(0);
    if( dim == 0 )
      {
      index[0] = 0;
      index[1] = it1.GetIndex()[0];
      index[2] = it1.GetIndex()[1];
      }
    if( dim == 1 )
      {
      index[0] = it1.GetIndex()[0];
      index[1] = 0;
      index[2] = it1.GetIndex()[1];
      }
    if( dim == 2 )
      {
      index[0] = it1.GetIndex()[0];
      index[1] = it1.GetIndex()[1];
      index[2] = 0;
      }
    index2d[0] = it1.GetIndex()[0];
    index2d[1] = it1.GetIndex()[1];
    stack->SetPixel(index, it1.Value() / mean);
    stack2->SetPixel(index2d, it1.Value() / mean);
    ++it1;
    }

  if( nSlices == 1 )
    {
    std::cout << " write slice " << std::endl;
    WriteImage<SliceType>(stack2, stackName);
    return 0;
    }
  // std::cout << "Stacking." << std::flush;
  for( unsigned int i = 1; i < nSlices; i++ )
    {
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[5 + i] );
    reader->Update();

    std::cout << " Slice " << i << " :: " << std::string(argv[5 + i]) << std::endl;
    //    std::cout << reader->GetOutput()->GetLargestPossibleRegion().GetSize() << std::endl;
    ExtractFilterType::Pointer extract = ExtractFilterType::New();
    extract->SetInput( reader->GetOutput() );
    extract->SetDirectionCollapseToIdentity();
    extract->SetExtractionRegion( extractRegion );
    extract->Update();

    SliceIt it(extract->GetOutput(), extract->GetOutput()->GetLargestPossibleRegion() );
    float   mean = 0.0, ct = 1;
    while( !it.IsAtEnd() )
      {
      mean += it.Value();
      ct++;
      ++it;
      }

    mean /= ct;
    std::cout << " Mean " << mean << std::endl;
    mean = 1;
    it.GoToBegin();
    while( !it.IsAtEnd() )
      {
      ImageType::IndexType index;
      index.Fill(0);
      if( dim == 0 )
        {
        index[0] = i;
        index[1] = it.GetIndex()[0];
        index[2] = it.GetIndex()[1];
        }
      if( dim == 1 )
        {
        index[0] = it.GetIndex()[0];
        index[1] = i;
        index[2] = it.GetIndex()[1];
        }
      if( dim == 2 )
        {
        index[0] = it.GetIndex()[0];
        index[1] = it.GetIndex()[1];
        index[2] = i;
        }
      stack->SetPixel(index, it.Value() / mean);
      ++it;
      }

    // std::cout << i << "." << std::flush;
    }
  // std::cout << "Done" << std::endl;

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( stackName );
  writer->SetInput( stack );
  writer->Update();

  // Input parameters
//   char * inputName  = argv[1];
//   unsigned int flip_x = atoi( argv[2] );
//   unsigned int flip_y = atoi( argv[3] );
//   unsigned int flip_z = atoi( argv[4] );
//   char * outputName = argv[5];

//   ReaderType::Pointer reader = ReaderType::New();
//   reader->SetFileName( inputName );
//   reader->Update();

//   // flip in desired directions for correct display
//   FlipFilterType::Pointer flip = FlipFilterType::New();
//   FlipFilterType::FlipAxesArrayType flipOver;
//   flipOver[0] = flip_x;
//   flipOver[1] = flip_y;
//   flipOver[2] = flip_z;
//   flip->SetInput( reader->GetOutput() );
//   flip->SetFlipAxes( flipOver );
//   flip->Update();

//   // write output
//   WriterType::Pointer writer = WriterType::New();
//   writer->SetInput( flip->GetOutput() );
//   writer->SetFileName( outputName );
//   writer->Update();

  return 0;
}
