#ifndef ReadWriteData_h_
#define ReadWriteData_h_

#include "antsAllocImage.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include "itkVector.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLabeledPointSetFileReader.h"
#include "itkLabeledPointSetFileWriter.h"
#include "itkImageIntensityAndGradientToPointSetFilter.h"
#include "itkWarpImageFilter.h"
// #include "itkInverseWarpImageFilter.h"
#include "itkAffineTransform.h"
#include "itkImageRegionIterator.h"
#include "itkResampleImageFilter.h"
#include "itkVariableLengthVector.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkLogTensorImageFilter.h"
#include "itkExpTensorImageFilter.h"
#include "itkCastImageFilter.h"
#include <sys/stat.h>

extern bool
ANTSFileExists(const std::string & strFilename);

extern bool
ANTSFileIsImage(const std::string & filename);

// Nifti stores DTI values in lower tri format but itk uses upper tri
// currently, nifti io does nothing to deal with this. if this changes
// the function below should be modified/eliminated.

#if 1 // currrently unimplemented
template <typename TImageType>
void
NiftiDTICheck(itk::SmartPointer<TImageType> &, const char *, bool)
{}

#else
template <typename TImageType>
void
NiftiDTICheck(itk::SmartPointer<TImageType> & target, const char * file, bool makeLower)
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
  std::string::size_type pos1 = filename.find(".nii");
  std::string::size_type pos2 = filename.find(".nia");
  if ((pos1 == std::string::npos) && (pos2 == std::string::npos))
  {
    return;
  }

  if (PixType::Length != 6)
  {
    return;
  }

  std::cout << "Performing lower/upper triangular format check for Nifti DTI" << std::endl;

  // swap elements 2 and 3 for lower<->upper conversion
  itk::ImageRegionIteratorWithIndex<TImageType> iter(target, target->GetLargestPossibleRegion());

  unsigned int looksLikeLower = 0;
  unsigned int looksLikeUpper = 0;
  unsigned int nBadVoxels = 0;
  unsigned int count = 0;

  unsigned int el2Neg = 0;
  unsigned int el3Neg = 0;

  while (!iter.IsAtEnd())
  {
    bool isValid = true;
    for (unsigned int i = 0; i < 6; i++)
    {
      if (iter.Get()[i] != iter.Get()[i])
      {
        ++nBadVoxels;
        isValid = false;
      }
    }

    double el2 = iter.Get()[2];
    double el3 = iter.Get()[3];

    if (el2 < 0)
    {
      ++el2Neg;
    }
    if (el3 < 0)
    {
      ++el3Neg;
    }

    if (isValid)
    {
      if (el2 > el3)
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

  // std::cout << "Invalid: " << nBadVoxels << "/" << count << std::endl;
  // std::cout << "Lower: " << looksLikeLower << ", Upper: " << looksLikeUpper << std::endl;
  // std::cout << "el2Neg: " << el2Neg << ", el3Neg: " << el3Neg << std::endl;

  if (((looksLikeUpper > looksLikeLower) && makeLower) || ((looksLikeLower > looksLikeUpper) && !makeLower))
  {
    std::cout << "Performing lower/upper triangular format swap for Nifti DTI" << std::endl;

    iter.GoToBegin();
    while (!iter.IsAtEnd())
    {
      PixType                     pix = iter.Get();
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

template <typename TImageType>
void
ReadTensorImage(itk::SmartPointer<TImageType> & target, const char * file, bool takelog = true)
{
  if (!ANTSFileExists(std::string(file)))
  {
    std::cerr << " file " << std::string(file) << " does not exist . " << std::endl;
    target = nullptr;
    return;
  }

  if (!ANTSFileIsImage(file))
  {
    std::cerr << " file " << std::string(file) << " is not recognized as a supported image format . " << std::endl;
    target = nullptr;
    return;
  }

  typedef TImageType                      ImageType;
  typedef itk::ImageFileReader<ImageType> FileSourceType;

  typedef itk::LogTensorImageFilter<ImageType, ImageType> LogFilterType;
  typename FileSourceType::Pointer                        reffilter = nullptr;
  if (file[0] == '0' && file[1] == 'x')
  {
    void * ptr;
    sscanf(file, "%p", (void **)&ptr);
    target = *(static_cast<typename TImageType::Pointer *>(ptr));
  }
  else
  {
    // Read the image files begin

    reffilter = FileSourceType::New();
    reffilter->SetFileName(file);
    try
    {
      reffilter->Update();
    }
    catch (const itk::ExceptionObject & e)
    {
      std::cerr << "Exception caught during reference file reading " << std::endl;
      std::cerr << e << " file " << file << std::endl;
      target = nullptr;
      return;
    }

    target = reffilter->GetOutput();
  }

  // NiftiDTICheck<ImageType>(target, file, false);

  if (takelog)
  {
    typename LogFilterType::Pointer logFilter = LogFilterType::New();
    logFilter->SetInput(reffilter->GetOutput());
    try
    {
      logFilter->Update();
    }
    catch (const itk::ExceptionObject & e)
    {
      std::cerr << "Exception caught during log tensor filter " << std::endl;
      std::cerr << e << " file " << file << std::endl;
      target = nullptr;
      return;
    }
    target = logFilter->GetOutput();
    std::cout << "Returning Log(D) for log-euclidean math ops" << std::endl;
  }
}

// function to determine if a file name is a memory address rather than a file
inline bool
FileIsPointer(const char *file)
{
  bool fileIsPointer = false;

  std::string comparetype1 = std::string("0x");
  std::string comparetype2 = std::string(file);
  comparetype2 = comparetype2.substr(0, 2);

  if (comparetype1 == comparetype2)
  {
    fileIsPointer = true;
  }
  else
  {
    std::string fileString = std::string(file);
    // Also treat file as a pointer if it is composed entirely of hex characters, and has
    // appropriate length
    // Windows pointers are just hex numbers without the leading 0x
    bool fileIsHexChars = true;

    std::string::iterator fileIt;

    for (fileIt = fileString.begin(); fileIt < fileString.end(); ++fileIt)
    {
      if (!std::isxdigit(*fileIt))
      {
        fileIsHexChars = false;
      }
    }

    if (fileIsHexChars && fileString.length() == (sizeof(int*) * 2))
    {
      fileIsPointer = true;
    }
  }
  return fileIsPointer;
}

template <typename TImageType>
// void ReadImage(typename TImageType::Pointer target, const char *file)
bool
ReadImage(itk::SmartPointer<TImageType> & target, const char * file)
{
  enum
  {
    ImageDimension = TImageType::ImageDimension
  };
  if (std::string(file).length() < 3)
  {
    target = nullptr;
    return false;
  }

  if (FileIsPointer(file))
  {
    typedef TImageType RImageType;
    void *             ptr;
    sscanf(file, "%p", (void **)&ptr);
    typename RImageType::Pointer Rimage = *(static_cast<typename RImageType::Pointer *>(ptr));
    /** more needs to be done here to cast the pointer to an image type --- this is a work-around */
    typedef itk::CastImageFilter<RImageType, TImageType> CastFilterType;
    typename CastFilterType::Pointer                     caster = CastFilterType::New();
    caster->SetInput(Rimage);
    caster->UpdateLargestPossibleRegion();
    target = caster->GetOutput();
  }
  else
  {
    if (!ANTSFileExists(std::string(file)))
    {
      std::cerr << " file " << std::string(file) << " does not exist . " << std::endl;
      target = nullptr;
      return false;
    }

    if (!ANTSFileIsImage(file))
    {
      std::cerr << " file " << std::string(file) << " is not recognized as a supported image format . " << std::endl;
      target = nullptr;
      return false;
    }

    typedef TImageType                      ImageType;
    typedef itk::ImageFileReader<ImageType> FileSourceType;

    typename FileSourceType::Pointer reffilter = FileSourceType::New();
    reffilter->SetFileName(file);
    try
    {
      reffilter->Update();
    }
    catch (const itk::ExceptionObject & e)
    {
      std::cerr << "Exception caught during reference file reading " << std::endl;
      std::cerr << e << " file " << file << std::endl;
      target = nullptr;
      std::exception();
      return false;
    }

    // typename ImageType::DirectionType dir;
    // dir.SetIdentity();
    //  reffilter->GetOutput()->SetDirection(dir);

    // std::cout << " setting pointer " << std::endl;
    target = reffilter->GetOutput();
  }
  return true;
}

template <typename ImageType>
typename ImageType::Pointer
ReadImage(char * fn)
{
  // Read the image files begin
  typedef itk::ImageFileReader<ImageType> FileSourceType;

  typename FileSourceType::Pointer reffilter = FileSourceType::New();
  reffilter->SetFileName(fn);
  try
  {
    reffilter->Update();
  }
  catch (const itk::ExceptionObject & e)
  {
    std::cerr << "Exception caught during image reference file reading " << std::endl;
    std::cerr << e << std::endl;
    return nullptr;
  }

  // typename ImageType::DirectionType dir;
  // dir.SetIdentity();
  //  reffilter->GetOutput()->SetDirection(dir);

  typename ImageType::Pointer target = reffilter->GetOutput();
  // if (reffilter->GetImageIO->GetNumberOfComponents() == 6)
  // NiftiDTICheck<ImageType>(target,fn);

  return target;
}

template <typename ImageType>
typename ImageType::Pointer
ReadTensorImage(char * fn, bool takelog = true)
{
  // Read the image files begin
  typedef itk::ImageFileReader<ImageType>                 FileSourceType;
  typedef itk::LogTensorImageFilter<ImageType, ImageType> LogFilterType;

  typename FileSourceType::Pointer reffilter = FileSourceType::New();
  reffilter->SetFileName(fn);
  try
  {
    reffilter->Update();
  }
  catch (const itk::ExceptionObject & e)
  {
    std::cerr << "Exception caught during tensor image reference file reading " << std::endl;
    std::cerr << e << std::endl;
    return nullptr;
  }

  typename ImageType::Pointer target = reffilter->GetOutput();

  NiftiDTICheck<ImageType>(target, fn, false);

  if (takelog)
  {
    typename LogFilterType::Pointer logFilter = LogFilterType::New();
    logFilter->SetInput(target->GetOutput());
    logFilter->Update();
    target = logFilter->GetOutput();
  }

  return target;
}

template <typename TPointSet>
// void ReadImage(typename TPointSet::Pointer target, const char *file)
bool
ReadLabeledPointSet(itk::SmartPointer<TPointSet> & target,
                    const char *                   file,
                    bool                           boundaryPointsOnly = false,
                    float                          samplingPercentage = 1.0)
{
  if (std::string(file).length() < 3)
  {
    std::cerr << " bad file name " << std::string(file) << std::endl;
    target = nullptr;
    return false;
  }

  if (!ANTSFileExists(std::string(file)))
  {
    std::cerr << " file " << std::string(file) << " does not exist . " << std::endl;
    target = nullptr;
    return false;
  }

  // Read the image files begin
  typedef itk::LabeledPointSetFileReader<TPointSet> FileSourceType;
  typename FileSourceType::Pointer                  reffilter = FileSourceType::New();
  reffilter->SetFileName(file);
  reffilter->SetExtractBoundaryPoints(boundaryPointsOnly);
  reffilter->SetRandomPercentage(samplingPercentage);
  try
  {
    reffilter->Update();
  }
  catch (const itk::ExceptionObject & e)
  {
    std::cerr << "Exception caught during point set reference file reading " << std::endl;
    std::cerr << e << std::endl;
    target = nullptr;
    return false;
  }

  target = reffilter->GetOutput();

  return true;
}

template <typename TImage, typename TMask, typename TPointSet>
bool
ReadImageIntensityPointSet(itk::SmartPointer<TPointSet> & target,
                           const char *                   imageFile,
                           const char *                   maskFile,
                           std::vector<unsigned int>      neighborhoodRadius,
                           double                         sigma)
{
  if (std::string(imageFile).length() < 3)
  {
    std::cerr << " bad image file name " << std::string(imageFile) << std::endl;
    target = nullptr;
    return false;
  }

  if (!ANTSFileExists(std::string(imageFile)))
  {
    std::cerr << " image file " << std::string(imageFile) << " does not exist . " << std::endl;
    target = nullptr;
    return false;
  }

  if (!ANTSFileIsImage(imageFile))
  {
    std::cerr << " file " << std::string(imageFile) << " is not recognized as a supported image format . " << std::endl;
    target = nullptr;
    return false;
  }

  if (std::string(maskFile).length() < 3)
  {
    std::cerr << " bad mask file name " << std::string(maskFile) << std::endl;
    target = nullptr;
    return false;
  }

  if (!ANTSFileExists(std::string(maskFile)))
  {
    std::cerr << " mask file " << std::string(maskFile) << " does not exist . " << std::endl;
    target = nullptr;
    return false;
  }

  if (!ANTSFileIsImage(maskFile))
  {
    std::cerr << " file " << std::string(maskFile) << " is not recognized as a supported image format . " << std::endl;
    target = nullptr;
    return false;
  }

  if (neighborhoodRadius.size() != TImage::ImageDimension)
  {
    std::cerr << " size of the neighborhood radius is not equal to the image dimension." << std::endl;
    target = nullptr;
    return false;
  }

  typename TImage::Pointer intensityImage = ReadImage<TImage>(const_cast<char *>(imageFile));
  typename TMask::Pointer  maskImage = ReadImage<TMask>(const_cast<char *>(maskFile));

  typedef itk::ImageIntensityAndGradientToPointSetFilter<TImage, TMask, TPointSet> FilterType;

  typename FilterType::NeighborhoodRadiusType radius;
  for (unsigned int d = 0; d < TImage::ImageDimension; d++)
  {
    radius[d] = neighborhoodRadius[d];
  }

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1(intensityImage);
  filter->SetInput2(maskImage);
  filter->SetSigma(sigma);
  filter->SetNeighborhoodRadius(radius);
  filter->Update();

  target = filter->GetOutput();

  return true;
}

template <typename TPointSet>
typename TPointSet::Pointer
ReadLabeledPointSet(char * fn)
{
  if (!ANTSFileExists(std::string(fn)))
  {
    std::cerr << " file " << std::string(fn) << " does not exist . " << std::endl;
    return;
  }

  // Read the image files begin
  typedef itk::LabeledPointSetFileReader<TPointSet> FileSourceType;
  typename FileSourceType::Pointer                  reffilter = FileSourceType::New();
  reffilter->SetFileName(fn);
  try
  {
    reffilter->Update();
  }
  catch (const itk::ExceptionObject & e)
  {
    std::cerr << "Exception caught during point set reference file reading " << std::endl;
    std::cerr << e << std::endl;
    return nullptr;
  }

  typename TPointSet::Pointer target = reffilter->GetOutput();

  return target;
}

template <typename TPointSet>
bool
WritePointSet(itk::SmartPointer<TPointSet> pointSet, const char * file)
{
  if (std::string(file).length() < 3)
  {
    return false;
  }

  typename itk::LabeledPointSetFileWriter<TPointSet>::Pointer writer = itk::LabeledPointSetFileWriter<TPointSet>::New();
  writer->SetFileName(file);
  if (!pointSet)
  {
    std::cerr << " Point set is nullptr." << std::endl;
    std::exception();
  }
  writer->SetInput(pointSet);
  writer->Update();

  return true;
}

namespace ANTs
{
template <typename TImageType>
bool
WriteImage(const itk::SmartPointer<TImageType> image, const char * file)
{
  if (std::string(file).length() < 3)
  {
    return false;
  }

  //  typename TImageType::DirectionType dir;
  // dir.SetIdentity();
  // image->SetDirection(dir);
  //  std::cout << " now Write direction " << image->GetOrigin() << std::endl;

  // if (writer->GetImageIO->GetNumberOfComponents() == 6)
  // NiftiDTICheck<TImageType>(image,file);

  if (FileIsPointer(file))
  {
    void * ptr;
    sscanf(file, "%p", (void **)&ptr);
    *(static_cast<typename TImageType::Pointer *>(ptr)) = image;
  }
  else
  {
    typename itk::ImageFileWriter<TImageType>::Pointer writer = itk::ImageFileWriter<TImageType>::New();
    writer->SetFileName(file);
    if (!image)
    {
      std::cerr << "Image is nullptr." << std::endl;
      std::exception();
    }
    writer->SetInput(image);
    writer->SetUseCompression(true);
    writer->Update();
  }
  return true;
}
} // namespace ANTs

template <typename TImageType>
void
WriteTensorImage(itk::SmartPointer<TImageType> image, const char * file, bool takeexp = true)
{
  typedef itk::ExpTensorImageFilter<TImageType, TImageType> ExpFilterType;
  typename itk::ImageFileWriter<TImageType>::Pointer        writer = itk::ImageFileWriter<TImageType>::New();
  writer->SetFileName(file);

  typename TImageType::Pointer writeImage = image;

  if (takeexp)
  {
    typename ExpFilterType::Pointer expFilter = ExpFilterType::New();
    expFilter->SetInput(image);
    expFilter->Update();
    writeImage = expFilter->GetOutput();
    std::cout << "Taking Exp(D) before writing" << std::endl;
  }

  // convert from upper tri to lower tri
  NiftiDTICheck<TImageType>(writeImage, file, true); // BA May 30 2009 -- remove b/c ITK fixed NIFTI reader

  if (file[0] == '0' && file[1] == 'x')
  {
    void * ptr;
    sscanf(file, "%p", (void **)&ptr);
    *(static_cast<typename TImageType::Pointer *>(ptr)) = writeImage;
  }
  else
  {
    writer->SetInput(writeImage);
    writer->SetUseCompression(true);
    writer->Update();
  }
}

template <typename TImage, typename TField>
typename TField::Pointer
ReadWarpFromFile(std::string warpfn, std::string ext)
{
  typedef TField                        FieldType;
  typedef typename FieldType::PixelType VectorType;
  enum
  {
    ImageDimension = FieldType::ImageDimension
  };

  typedef itk::Image<float, ImageDimension> RealImageType;
  typedef RealImageType                     ImageType;

  //  typedef itk::Vector<float,Self::ImageDimension>         VectorType;
  //  typedef itk::Image<VectorType,Self::ImageDimension>     FieldType;
  // std::cout << " warp file name " << warpfn + ext << std::endl;

  // First - read the vector fields
  // NOTE : THE FIELD SHOULD WARP INPUT1 TO INPUT2, THUS SHOULD POINT
  // FROM INPUT2 TO INPUT1
  std::string                     fn = warpfn + "x" + ext;
  typename RealImageType::Pointer xvec = ReadImage<ImageType>((char *)fn.c_str());
  //  std::cout << " done reading " << fn << std::endl;
  fn = warpfn + "y" + ext;
  typename RealImageType::Pointer yvec = ReadImage<ImageType>((char *)fn.c_str());
  // std::cout << " done reading " << fn << std::endl;
  fn = warpfn + "z" + ext;
  typename RealImageType::Pointer zvec = nullptr;
  // std::cout << " done reading " << fn << std::endl;
  if (ImageDimension == 3)
  {
    zvec = ReadImage<ImageType>((char *)fn.c_str());
  }

  typename FieldType::Pointer field = AllocImage<FieldType>(xvec);

  itk::ImageRegionIteratorWithIndex<RealImageType> it(xvec, xvec->GetLargestPossibleRegion());

  //  std::cout << " spacing xv " << xvec->GetSpacing()[0]
  // << " field " << field->GetSpacing()[0] << std::endl;

  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    typename ImageType::IndexType index = it.GetIndex();

    VectorType disp;
    disp[0] = xvec->GetPixel(index);
    disp[1] = yvec->GetPixel(index);
    if (ImageDimension == 3)
    {
      disp[2] = zvec->GetPixel(index);
    }

    field->SetPixel(index, disp);

    //    if (ct == 10000) std::cout << " 10000th pix " << disp << std::endl;
  }

  return field;
}


template <typename TImage>
typename TImage::Pointer
MakeNewImage(typename TImage::Pointer image1, typename TImage::PixelType initval)
{
  typedef itk::ImageRegionIteratorWithIndex<TImage> Iterator;
  typename TImage::Pointer                          varimage = AllocImage<TImage>(image1);
  Iterator                                          vfIter2(varimage, varimage->GetLargestPossibleRegion());
  for (vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2)
  {
    if (initval >= 0)
    {
      vfIter2.Set(initval);
    }
    else
    {
      vfIter2.Set(image1->GetPixel(vfIter2.GetIndex()));
    }
  }

  return varimage;
}

template <typename TField>
void
WriteDisplacementField(TField * field, std::string filename)
{
  typedef TField FieldType;
  enum
  {
    ImageDimension = FieldType::ImageDimension
  };

  typedef itk::Image<float, ImageDimension> RealImageType;

  // Initialize the caster to the displacement field
  typedef itk::VectorIndexSelectionCastImageFilter<FieldType, RealImageType> IndexSelectCasterType;
  for (unsigned int dim = 0; dim < ImageDimension; dim++)
  {
    typename IndexSelectCasterType::Pointer fieldCaster = IndexSelectCasterType::New();
    fieldCaster->SetInput(field);

    fieldCaster->SetIndex(dim);
    fieldCaster->Update();

    // Set up the output filename
    std::string outfile = filename + static_cast<char>('x' + dim) + std::string("vec.nii.gz");
    std::cout << "Writing displacements to " << outfile << " spacing " << field->GetSpacing()[0] << std::endl;
    typename RealImageType::Pointer fieldcomponent = fieldCaster->GetOutput();
    fieldcomponent->SetSpacing(field->GetSpacing());
    fieldcomponent->SetOrigin(field->GetOrigin());
    fieldcomponent->SetDirection(field->GetDirection());

    ANTs::WriteImage<RealImageType>(fieldcomponent, outfile.c_str());
  }
  std::cout << "...done" << std::endl;
  return;
}

template <typename TField>
void
WriteDisplacementField2(TField * field, std::string filename, std::string app)
{
  typedef TField FieldType;
  enum
  {
    ImageDimension = FieldType::ImageDimension
  };

  typedef itk::Image<float, ImageDimension> RealImageType;

  // Initialize the caster to the displacement field
  typedef itk::VectorIndexSelectionCastImageFilter<FieldType, RealImageType> IndexSelectCasterType;
  for (unsigned int dim = 0; dim < ImageDimension; dim++)
  {
    typename IndexSelectCasterType::Pointer fieldCaster = IndexSelectCasterType::New();
    fieldCaster->SetInput(field);

    fieldCaster->SetIndex(dim);
    fieldCaster->Update();

    // Set up the output filename
    std::string outfile = filename + static_cast<char>('x' + dim) + std::string(app);
    std::cout << "Writing displacements to " << outfile << " spacing " << field->GetSpacing()[0] << std::endl;
    typename RealImageType::Pointer fieldcomponent = fieldCaster->GetOutput();
    fieldcomponent->SetSpacing(field->GetSpacing());
    fieldcomponent->SetOrigin(field->GetOrigin());

    ANTs::WriteImage<RealImageType>(fieldcomponent, outfile.c_str());
  }
  std::cout << "...done" << std::endl;
  return;
}

template <typename TTimeSeriesImageType, typename MultiChannelImageType>
typename MultiChannelImageType::Pointer
ConvertTimeSeriesImageToMultiChannelImage(TTimeSeriesImageType * timeSeriesImage)
{
  enum
  {
    ImageDimension = MultiChannelImageType::ImageDimension
  };

  typename MultiChannelImageType::SpacingType          spacing;
  typename MultiChannelImageType::PointType            origin;
  typename MultiChannelImageType::RegionType::SizeType size;
  typename MultiChannelImageType::DirectionType        direction;

  typename TTimeSeriesImageType::SpacingType          timeSeriesSpacing = timeSeriesImage->GetSpacing();
  typename TTimeSeriesImageType::PointType            timeSeriesOrigin = timeSeriesImage->GetOrigin();
  typename TTimeSeriesImageType::RegionType::SizeType timeSeriesSize = timeSeriesImage->GetRequestedRegion().GetSize();
  typename TTimeSeriesImageType::DirectionType        timeSeriesDirection = timeSeriesImage->GetDirection();

  for (itk::SizeValueType d = 0; d < ImageDimension; d++)
  {
    spacing[d] = timeSeriesSpacing[d];
    origin[d] = timeSeriesOrigin[d];
    size[d] = timeSeriesSize[d];
    for (itk::SizeValueType e = 0; e < ImageDimension; e++)
    {
      direction(d, e) = timeSeriesDirection(d, e);
    }
  }

  typename MultiChannelImageType::Pointer multiChannelImage = MultiChannelImageType::New();
  multiChannelImage->SetRegions(size);
  multiChannelImage->SetSpacing(spacing);
  multiChannelImage->SetOrigin(origin);
  multiChannelImage->SetDirection(direction);
  multiChannelImage->SetVectorLength(timeSeriesSize[ImageDimension]);
  multiChannelImage->AllocateInitialized();

  itk::ImageRegionIteratorWithIndex<MultiChannelImageType> It(multiChannelImage,
                                                              multiChannelImage->GetRequestedRegion());
  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    typename MultiChannelImageType::IndexType index = It.GetIndex();

    typename MultiChannelImageType::PixelType multiChannelVoxel;
    multiChannelVoxel.SetSize(timeSeriesSize[ImageDimension]);

    for (itk::SizeValueType n = 0; n < timeSeriesSize[ImageDimension]; n++)
    {
      typename TTimeSeriesImageType::IndexType timeSeriesIndex;
      for (itk::SizeValueType d = 0; d < ImageDimension; d++)
      {
        timeSeriesIndex[d] = index[d];
      }
      timeSeriesIndex[ImageDimension] = n;
      multiChannelVoxel[n] = timeSeriesImage->GetPixel(timeSeriesIndex);
    }
    It.Set(multiChannelVoxel);
  }

  return multiChannelImage;
}

template <typename MultiChannelImageType, typename TimeSeriesImageType>
typename TimeSeriesImageType::Pointer
ConvertMultiChannelImageToTimeSeriesImage(MultiChannelImageType * multiChannelImage)
{
  enum
  {
    ImageDimension = MultiChannelImageType::ImageDimension
  };

  typename MultiChannelImageType::SpacingType          spacing = multiChannelImage->GetSpacing();
  typename MultiChannelImageType::PointType            origin = multiChannelImage->GetOrigin();
  typename MultiChannelImageType::RegionType::SizeType size = multiChannelImage->GetRequestedRegion().GetSize();
  typename MultiChannelImageType::DirectionType        direction = multiChannelImage->GetDirection();

  typename TimeSeriesImageType::SpacingType          timeSeriesSpacing;
  typename TimeSeriesImageType::PointType            timeSeriesOrigin;
  typename TimeSeriesImageType::RegionType::SizeType timeSeriesSize;
  typename TimeSeriesImageType::DirectionType        timeSeriesDirection;
  timeSeriesDirection.SetIdentity();

  typename MultiChannelImageType::IndexType index;
  index.Fill(0);
  typename MultiChannelImageType::PixelType multiChannelVoxel = multiChannelImage->GetPixel(index);

  for (itk::SizeValueType d = 0; d < ImageDimension; d++)
  {
    timeSeriesSpacing[d] = spacing[d];
    timeSeriesOrigin[d] = origin[d];
    timeSeriesSize[d] = size[d];
    for (itk::SizeValueType e = 0; e < ImageDimension; e++)
    {
      timeSeriesDirection(d, e) = direction(d, e);
    }
  }
  timeSeriesSpacing[ImageDimension] = 1;
  timeSeriesOrigin[ImageDimension] = 0;
  timeSeriesSize[ImageDimension] = multiChannelVoxel.GetSize();

  typename TimeSeriesImageType::Pointer timeSeriesImage =
    AllocImage<TimeSeriesImageType>(timeSeriesSize, timeSeriesSpacing, timeSeriesOrigin, timeSeriesDirection);

  itk::ImageRegionIteratorWithIndex<MultiChannelImageType> It(multiChannelImage,
                                                              multiChannelImage->GetRequestedRegion());
  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    index = It.GetIndex();
    multiChannelVoxel = It.Get();

    for (itk::SizeValueType n = 0; n < timeSeriesSize[ImageDimension]; n++)
    {
      typename TimeSeriesImageType::IndexType timeSeriesIndex;
      for (itk::SizeValueType d = 0; d < ImageDimension; d++)
      {
        timeSeriesIndex[d] = index[d];
      }
      timeSeriesIndex[ImageDimension] = n;
      timeSeriesImage->SetPixel(timeSeriesIndex, multiChannelVoxel[n]);
    }
  }

  return timeSeriesImage;
}

template <typename TTimeSeriesImageType, typename FiveDimensionalImageType>
typename FiveDimensionalImageType::Pointer
ConvertTimeSeriesImageToFiveDimensionalImage(TTimeSeriesImageType * timeSeriesImage)
{
  enum
  {
    ImageDimension = TTimeSeriesImageType::ImageDimension - 1
  };

  typename FiveDimensionalImageType::SpacingType          spacing;
  typename FiveDimensionalImageType::PointType            origin;
  typename FiveDimensionalImageType::RegionType::SizeType size;
  typename FiveDimensionalImageType::DirectionType        direction;

  origin.Fill(0.0);
  spacing.Fill(1.0);
  direction.SetIdentity();

  typename TTimeSeriesImageType::SpacingType          timeSeriesSpacing = timeSeriesImage->GetSpacing();
  typename TTimeSeriesImageType::PointType            timeSeriesOrigin = timeSeriesImage->GetOrigin();
  typename TTimeSeriesImageType::RegionType::SizeType timeSeriesSize = timeSeriesImage->GetRequestedRegion().GetSize();
  typename TTimeSeriesImageType::DirectionType        timeSeriesDirection = timeSeriesImage->GetDirection();

  for (itk::SizeValueType d = 0; d < ImageDimension; d++)
  {
    spacing[d] = timeSeriesSpacing[d];
    origin[d] = timeSeriesOrigin[d];
    size[d] = timeSeriesSize[d];
    for (itk::SizeValueType e = 0; e < ImageDimension; e++)
    {
      direction(d, e) = timeSeriesDirection(d, e);
    }
  }
  size[3] = itk::NumericTraits<itk::SizeValueType>::OneValue();
  size[4] = timeSeriesSize[ImageDimension];

  typename FiveDimensionalImageType::Pointer FiveDimensionalImage = FiveDimensionalImageType::New();
  FiveDimensionalImage->SetRegions(size);
  FiveDimensionalImage->SetSpacing(spacing);
  FiveDimensionalImage->SetOrigin(origin);
  FiveDimensionalImage->SetDirection(direction);
  FiveDimensionalImage->AllocateInitialized();

  itk::ImageRegionIteratorWithIndex<FiveDimensionalImageType> It(FiveDimensionalImage,
                                                                 FiveDimensionalImage->GetRequestedRegion());
  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    typename FiveDimensionalImageType::IndexType index = It.GetIndex();
    typename TTimeSeriesImageType::IndexType     timeSeriesIndex;

    timeSeriesIndex[0] = index[0];
    timeSeriesIndex[1] = index[1];
    timeSeriesIndex[2] = index[2];
    timeSeriesIndex[3] = index[4];

    It.Set(timeSeriesImage->GetPixel(timeSeriesIndex));
  }

  return FiveDimensionalImage;
}

template <typename FiveDimensionalImageType, typename TimeSeriesImageType>
typename TimeSeriesImageType::Pointer
ConvertFiveDimensionalImageToTimeSeriesImage(FiveDimensionalImageType * FiveDimensionalImage)
{
  enum
  {
    ImageDimension = TimeSeriesImageType::ImageDimension - 1
  };

  typename FiveDimensionalImageType::SpacingType          spacing = FiveDimensionalImage->GetSpacing();
  typename FiveDimensionalImageType::PointType            origin = FiveDimensionalImage->GetOrigin();
  typename FiveDimensionalImageType::RegionType::SizeType size = FiveDimensionalImage->GetRequestedRegion().GetSize();
  typename FiveDimensionalImageType::DirectionType        direction = FiveDimensionalImage->GetDirection();

  typename TimeSeriesImageType::SpacingType          timeSeriesSpacing;
  typename TimeSeriesImageType::PointType            timeSeriesOrigin;
  typename TimeSeriesImageType::RegionType::SizeType timeSeriesSize;
  typename TimeSeriesImageType::DirectionType        timeSeriesDirection;

  timeSeriesOrigin.Fill(0.0);
  timeSeriesSpacing.Fill(1.0);
  timeSeriesDirection.SetIdentity();

  for (itk::SizeValueType d = 0; d < ImageDimension; d++)
  {
    timeSeriesSpacing[d] = spacing[d];
    timeSeriesOrigin[d] = origin[d];
    timeSeriesSize[d] = size[d];
    for (itk::SizeValueType e = 0; e < ImageDimension; e++)
    {
      timeSeriesDirection(d, e) = direction(d, e);
    }
  }
  timeSeriesSize[ImageDimension] = size[4];

  typename TimeSeriesImageType::Pointer timeSeriesImage =
    AllocImage<TimeSeriesImageType>(timeSeriesSize, timeSeriesSpacing, timeSeriesOrigin, timeSeriesDirection);

  itk::ImageRegionIteratorWithIndex<FiveDimensionalImageType> It(FiveDimensionalImage,
                                                                 FiveDimensionalImage->GetRequestedRegion());
  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    typename FiveDimensionalImageType::IndexType index = It.GetIndex();
    typename TimeSeriesImageType::IndexType      timeSeriesIndex;

    timeSeriesIndex[0] = index[0];
    timeSeriesIndex[1] = index[1];
    timeSeriesIndex[2] = index[2];
    timeSeriesIndex[3] = index[4];

    timeSeriesImage->SetPixel(timeSeriesIndex, It.Get());
  }

  return timeSeriesImage;
}

class nullBuf : public std::streambuf
{
public:
  std::streamsize
  xsputn(const char * itkNotUsed(s), std::streamsize n) override
  {
    return n;
  }

  int
  overflow(int itkNotUsed(c)) override
  {
    return 1;
  }
};

class nullStream : public std::ostream
{
public:
  nullStream()
    : std::ostream(&buf)
  {}

private:
  nullBuf buf;
};

#endif
