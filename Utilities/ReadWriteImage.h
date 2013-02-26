#ifndef __ReadWriteImage_h_
#define __ReadWriteImage_h_
#include <antsAllocImage.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "itkVector.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkWarpImageFilter.h"
// #include "itkInverseWarpImageFilter.h"
#include "itkAffineTransform.h"
#include "itkImageRegionIterator.h"
#include "itkResampleImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkLogTensorImageFilter.h"
#include "itkExpTensorImageFilter.h"
#include "itkCastImageFilter.h"
#include <sys/stat.h>

extern bool ANTSFileExists(const std::string & strFilename);

// Nifti stores DTI values in lower tri format but itk uses upper tri
// currently, nifti io does nothing to deal with this. if this changes
// the function below should be modified/eliminated.

#if 1                           // currrently unimplemented
template <class TImageType>
void NiftiDTICheck(itk::SmartPointer<TImageType> &, const char *, bool )
{
}

#else
template <class TImageType>
void NiftiDTICheck(itk::SmartPointer<TImageType> & target, const char *file, bool makeLower)
{
  typedef typename TImageType::PixelType PixType;

  return;
  // typedef itk::ImageFileWriter<TImageType> Writer;
  // typename Writer::Pointer writer = Writer::New();
  // writer->SetInput( target );
  // writer->SetFileName( "testdt.nii" );
  // writer->Update();

  // Check for nifti file
  std::string            filename = file;
  std::string::size_type pos1 = filename.find( ".nii" );
  std::string::size_type pos2 = filename.find( ".nia" );
  if( (pos1 == std::string::npos) && (pos2 == std::string::npos) )
    {
    return;
    }

  if( PixType::Length != 6 )
    {
    return;
    }

  ::ants::antscout << "Performing lower/upper triangular format check for Nifti DTI" << std::endl;

  // swap elements 2 and 3 for lower<->upper conversion
  itk::ImageRegionIteratorWithIndex<TImageType>
  iter(target, target->GetLargestPossibleRegion() );

  unsigned int looksLikeLower = 0;
  unsigned int looksLikeUpper = 0;
  unsigned int nBadVoxels = 0;
  unsigned int count = 0;

  unsigned int el2Neg = 0;
  unsigned int el3Neg = 0;

  while( !iter.IsAtEnd() )
    {
    bool isValid = true;
    for( unsigned int i = 0; i < 6; i++ )
      {
      if( iter.Get()[i] != iter.Get()[i] )
        {
        ++nBadVoxels;
        isValid = false;
        }
      }

    double el2 = iter.Get()[2];
    double el3 = iter.Get()[3];

    if( el2 < 0 )
      {
      ++el2Neg;
      }
    if( el3 < 0 )
      {
      ++el3Neg;
      }

    if( isValid )
      {
      if( el2 > el3 )
        {
        ++looksLikeLower;
        }
      else
        {
        ++looksLikeUpper;
        }
      }

    ++count;
    ++iter;
    }

  // ::ants::antscout << "Invalid: " << nBadVoxels << "/" << count << std::endl;
  // ::ants::antscout << "Lower: " << looksLikeLower << ", Upper: " << looksLikeUpper << std::endl;
  // ::ants::antscout << "el2Neg: " << el2Neg << ", el3Neg: " << el3Neg << std::endl;

  if( ( (looksLikeUpper > looksLikeLower) && makeLower) || ( (looksLikeLower > looksLikeUpper) && !makeLower) )
    {
    ::ants::antscout << "Performing lower/upper triangular format swap for Nifti DTI" << std::endl;

    iter.GoToBegin();
    while( !iter.IsAtEnd() )
      {
      PixType pix = iter.Get();
      typename PixType::ValueType temp;
      temp = pix[2];
      pix[2] = pix[3];
      pix[3] = temp;
      iter.Set(pix);
      ++iter;
      }
    }
}

#endif

template <class TImageType>
void ReadTensorImage(itk::SmartPointer<TImageType> & target, const char *file, bool takelog = true)
{
  if( !ANTSFileExists(std::string(file) ) )
    {
    ::ants::antscout << " file " << std::string(file) << " does not exist . " << std::endl;  return;
    }

  typedef TImageType                      ImageType;
  typedef itk::ImageFileReader<ImageType> FileSourceType;
  typedef typename ImageType::PixelType   PixType;

  typedef itk::LogTensorImageFilter<ImageType, ImageType> LogFilterType;
  typename FileSourceType::Pointer reffilter = NULL;
  if( file[0] == '0' && file[1] == 'x' )
    {
    std::stringstream strstream;
    strstream << file;
    void* ptr;
    strstream >> ptr;
    target = *( static_cast<typename TImageType::Pointer *>( ptr ) );
    }
  else
    {
    // Read the image files begin

    reffilter = FileSourceType::New();
    reffilter->SetFileName( file );
    try
      {
      reffilter->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      ::ants::antscout << "Exception caught during reference file reading " << std::endl;
      ::ants::antscout << e << " file " << file << std::endl;
      target = NULL;
      return;
      }

    target = reffilter->GetOutput();
    }
  NiftiDTICheck<ImageType>(target, file, false);

  if( takelog )
    {
    typename LogFilterType::Pointer logFilter = LogFilterType::New();
    logFilter->SetInput( reffilter->GetOutput() );
    logFilter->Update();
    target = logFilter->GetOutput();
    ::ants::antscout << "Returning Log(D) for log-euclidean math ops" << std::endl;
    }
}

template <class TImageType>
// void ReadImage(typename TImageType::Pointer target, const char *file)
bool ReadImage(itk::SmartPointer<TImageType> & target, const char *file)
{
  typedef typename TImageType::PixelType PixelType;
  enum { ImageDimension = TImageType::ImageDimension };
  if( std::string(file).length() < 3 )
    {
    ::ants::antscout << " bad file name " << std::string(file) << std::endl;    target = NULL;  return false;
    }

  std::string comparetype1 = std::string( "0x" );
  std::string comparetype2 = std::string( file );
  comparetype2 = comparetype2.substr( 0, 2 );
  // Read the image files begin
  if(  comparetype1 == comparetype2  )
    {
    typedef TImageType RImageType;
    std::stringstream strstream;
    strstream << file;
    void* ptr;
    strstream >> ptr;
    typename RImageType::Pointer Rimage = *( static_cast<typename RImageType::Pointer *>( ptr ) );
    /** more needs to be done here to cast the pointer to an image type --- this is a work-around */
    typedef itk::CastImageFilter<RImageType, TImageType> CastFilterType;
    typename CastFilterType::Pointer caster = CastFilterType::New();
    caster->SetInput( Rimage );
    caster->UpdateLargestPossibleRegion();
    target = caster->GetOutput();
    }
  else
    {
    if( !ANTSFileExists(std::string(file) ) )
      {
      ::ants::antscout << " file " << std::string(file) << " does not exist . " << std::endl; target = NULL;
      return false;
      }
    typedef TImageType                      ImageType;
    typedef itk::ImageFileReader<ImageType> FileSourceType;
    typedef typename ImageType::PixelType   PixType;
    //    const unsigned int ImageDimension=ImageType::ImageDimension;
    typename FileSourceType::Pointer reffilter = FileSourceType::New();
    reffilter->SetFileName( file );
    try
      {
      reffilter->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      ::ants::antscout << "Exception caught during reference file reading " << std::endl;
      ::ants::antscout << e << " file " << file << std::endl;
      target = NULL;
      std::exception();
      return false;
      }

    // typename ImageType::DirectionType dir;
    // dir.SetIdentity();
    //  reffilter->GetOutput()->SetDirection(dir);

    // ::ants::antscout << " setting pointer " << std::endl;
    target = reffilter->GetOutput();
    }
  return true;
}

template <class ImageType>
typename ImageType::Pointer ReadImage(char* fn )
{
  // Read the image files begin
  typedef itk::ImageFileReader<ImageType> FileSourceType;
  typedef typename ImageType::PixelType   PixType;
//   const unsigned int ImageDimension=ImageType::ImageDimension;
  typename FileSourceType::Pointer reffilter = FileSourceType::New();
  reffilter->SetFileName( fn );
  try
    {
    reffilter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    ::ants::antscout << "Exception caught during reference file reading " << std::endl;
    ::ants::antscout << e << std::endl;
    return NULL;
    }

  typename ImageType::DirectionType dir;
  dir.SetIdentity();
  //  reffilter->GetOutput()->SetDirection(dir);

  typename ImageType::Pointer target = reffilter->GetOutput();
  // if (reffilter->GetImageIO->GetNumberOfComponents() == 6)
  // NiftiDTICheck<ImageType>(target,fn);

  return target;
}

template <class ImageType>
typename ImageType::Pointer ReadTensorImage(char* fn, bool takelog = true )
{
  // Read the image files begin
  typedef itk::ImageFileReader<ImageType>                 FileSourceType;
  typedef typename ImageType::PixelType                   PixType;
  typedef itk::LogTensorImageFilter<ImageType, ImageType> LogFilterType;
//   const unsigned int ImageDimension=ImageType::ImageDimension;
  typename FileSourceType::Pointer reffilter = FileSourceType::New();
  reffilter->SetFileName( fn );
  try
    {
    reffilter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    ::ants::antscout << "Exception caught during reference file reading " << std::endl;
    ::ants::antscout << e << std::endl;
    return NULL;
    }

  typename ImageType::Pointer target = reffilter->GetOutput();

  NiftiDTICheck<ImageType>(target, fn, false);

  if( takelog )
    {
    typename LogFilterType::Pointer logFilter = LogFilterType::New();
    logFilter->SetInput( target->GetOutput() );
    logFilter->Update();
    target = logFilter->GetOutput();
    }

  return target;
}

template <class TImageType>
bool WriteImage(itk::SmartPointer<TImageType> image, const char *file)
{
  if( std::string(file).length() < 3 )
    {
    return false;
    }

  //  typename TImageType::DirectionType dir;
  // dir.SetIdentity();
  // image->SetDirection(dir);
  //  ::ants::antscout << " now Write direction " << image->GetOrigin() << std::endl;

  // if (writer->GetImageIO->GetNumberOfComponents() == 6)
  // NiftiDTICheck<TImageType>(image,file);

  if( file[0] == '0' && file[1] == 'x' )
    {
    std::stringstream strstream;
    strstream << file;
    void* ptr;
    strstream >> ptr;
    *( static_cast<typename TImageType::Pointer *>( ptr ) ) = image;
    }
  else
    {
    typename itk::ImageFileWriter<TImageType>::Pointer writer =
      itk::ImageFileWriter<TImageType>::New();
    writer->SetFileName(file);
    if( !image )
      {
      ::ants::antscout << " file " << file << " does not exist " << std::endl;
      std::exception();
      }
    writer->SetInput(image);
    writer->Update();
    }
  return true;
}

template <class TImageType>
void WriteTensorImage(itk::SmartPointer<TImageType> image, const char *file, bool takeexp = true)
{
  typedef typename TImageType::PixelType                    PixType;
  typedef itk::ExpTensorImageFilter<TImageType, TImageType> ExpFilterType;
  typename itk::ImageFileWriter<TImageType>::Pointer writer =
    itk::ImageFileWriter<TImageType>::New();
  writer->SetFileName(file);

  typename TImageType::Pointer writeImage = image;

  if( takeexp )
    {
    typename ExpFilterType::Pointer expFilter = ExpFilterType::New();
    expFilter->SetInput( image );
    expFilter->Update();
    writeImage = expFilter->GetOutput();
    ::ants::antscout << "Taking Exp(D) before writing" << std::endl;
    }

  // convert from upper tri to lower tri
  NiftiDTICheck<TImageType>(writeImage, file, true); // BA May 30 2009 -- remove b/c ITK fixed NIFTI reader

  if( file[0] == '0' && file[1] == 'x' )
    {
    std::stringstream strstream;
    strstream << file;
    void* ptr;
    strstream >> ptr;
    *( static_cast<typename TImageType::Pointer *>( ptr ) ) = writeImage;
    }
  else
    {
    writer->SetInput(writeImage);
    writer->Update();
    }
}

template <class TImage, class TField>
typename TField::Pointer
ReadWarpFromFile( std::string warpfn, std::string ext)
{
  typedef TField                        FieldType;
  typedef typename FieldType::PixelType VectorType;
  enum { ImageDimension = FieldType::ImageDimension };

  typedef itk::Image<float, ImageDimension> RealImageType;
  typedef RealImageType                     ImageType;

//  typedef itk::Vector<float,itkGetStaticConstMacro(ImageDimension)>         VectorType;
//  typedef itk::Image<VectorType,itkGetStaticConstMacro(ImageDimension)>     FieldType;
// ::ants::antscout << " warp file name " << warpfn + ext << std::endl;

// First - read the vector fields
// NOTE : THE FIELD SHOULD WARP INPUT1 TO INPUT2, THUS SHOULD POINT
// FROM INPUT2 TO INPUT1
  std::string fn = warpfn + "x" + ext;
  typename RealImageType::Pointer xvec = ReadImage<ImageType>( (char *)fn.c_str() );
  //  ::ants::antscout << " done reading " << fn << std::endl;
  fn = warpfn + "y" + ext;
  typename RealImageType::Pointer yvec = ReadImage<ImageType>( (char *)fn.c_str() );
  // ::ants::antscout << " done reading " << fn << std::endl;
  fn = warpfn + "z" + ext;
  typename RealImageType::Pointer zvec = NULL;
  // ::ants::antscout << " done reading " << fn << std::endl;
  if( ImageDimension == 3 )
    {
    zvec = ReadImage<ImageType>( (char *)fn.c_str() );
    }

  typename FieldType::Pointer field = AllocImage<FieldType>(xvec);

  itk::ImageRegionIteratorWithIndex<RealImageType>
  it( xvec, xvec->GetLargestPossibleRegion() );

  //  ::ants::antscout << " spacing xv " << xvec->GetSpacing()[0]
  // << " field " << field->GetSpacing()[0] << std::endl;

  unsigned int ct = 0;
  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    ct++;
    typename ImageType::IndexType index = it.GetIndex();

    VectorType disp;
    disp[0] = xvec->GetPixel(index);
    disp[1] = yvec->GetPixel(index);
    if( ImageDimension == 3 )
      {
      disp[2] = zvec->GetPixel(index);
      }

    field->SetPixel(index, disp);

    //    if (ct == 10000) ::ants::antscout << " 10000th pix " << disp << std::endl;
    }

  return field;
}

template <class TImage>
typename TImage::Pointer
MakeNewImage(typename TImage::Pointer image1, typename TImage::PixelType initval )
{
  typedef itk::ImageRegionIteratorWithIndex<TImage> Iterator;
  typename TImage::Pointer varimage = AllocImage<TImage>(image1);
  Iterator vfIter2( varimage,  varimage->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    if( initval >= 0 )
      {
      vfIter2.Set(initval);
      }
    else
      {
      vfIter2.Set(image1->GetPixel(vfIter2.GetIndex() ) );
      }
    }

  return varimage;
}

template <class TField>
void
WriteDisplacementField(TField* field, std::string filename)
{
  typedef TField                        FieldType;
  typedef typename FieldType::PixelType VectorType;
  enum { ImageDimension = FieldType::ImageDimension };

  typedef itk::Image<float, ImageDimension> RealImageType;

  // Initialize the caster to the displacement field
  typedef itk::VectorIndexSelectionCastImageFilter<FieldType, RealImageType> IndexSelectCasterType;
  for( unsigned int dim = 0; dim < ImageDimension; dim++ )
    {
    typename IndexSelectCasterType::Pointer fieldCaster = IndexSelectCasterType::New();
    fieldCaster->SetInput( field );

    fieldCaster->SetIndex( dim );
    fieldCaster->Update();

    // Set up the output filename
    std::string outfile = filename + static_cast<char>('x' + dim) + std::string("vec.nii.gz");
    ::ants::antscout << "Writing displacements to " << outfile << " spacing "
                     << field->GetSpacing()[0] << std::endl;
    typename RealImageType::Pointer fieldcomponent = fieldCaster->GetOutput();
    fieldcomponent->SetSpacing(field->GetSpacing() );
    fieldcomponent->SetOrigin(field->GetOrigin() );
    fieldcomponent->SetDirection(field->GetDirection() );

    WriteImage<RealImageType>(fieldcomponent, outfile.c_str() );
    }
  ::ants::antscout << "...done" << std::endl;
  return;
}

template <class TField>
void
WriteDisplacementField2(TField* field, std::string filename, std::string app)
{
  typedef TField                        FieldType;
  typedef typename FieldType::PixelType VectorType;
  enum { ImageDimension = FieldType::ImageDimension };

  typedef itk::Image<float, ImageDimension> RealImageType;

  // Initialize the caster to the displacement field
  typedef itk::VectorIndexSelectionCastImageFilter<FieldType, RealImageType> IndexSelectCasterType;
  for( unsigned int dim = 0; dim < ImageDimension; dim++ )
    {
    typename IndexSelectCasterType::Pointer fieldCaster = IndexSelectCasterType::New();
    fieldCaster->SetInput( field );

    fieldCaster->SetIndex( dim );
    fieldCaster->Update();

    // Set up the output filename
    std::string outfile = filename + static_cast<char>('x' + dim) + std::string(app);
    ::ants::antscout << "Writing displacements to " << outfile << " spacing "
                     << field->GetSpacing()[0] << std::endl;
    typename RealImageType::Pointer fieldcomponent = fieldCaster->GetOutput();
    fieldcomponent->SetSpacing(field->GetSpacing() );
    fieldcomponent->SetOrigin(field->GetOrigin() );

    WriteImage<RealImageType>(fieldcomponent, outfile.c_str() );
    }
  ::ants::antscout << "...done" << std::endl;
  return;
}

#endif
