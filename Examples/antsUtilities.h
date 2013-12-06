/**
 * There were several functions that had been copied and
 * pasted over and over again in this library.  This
 * header files is contains a common definitoin for
 * those file.
 * \author Hans J. Johnson
 */
#ifndef __antsUtilities_h__
#define __antsUtilities_h__

// #include "antscout.hxx"
#include "antsAllocImage.h"
#include <string>
#include <vector>
#include <algorithm>

#include <stdlib.h>
#include <errno.h>
#include <cmath>
#include <iostream>

#include "itkVector.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"

// We need to ensure that only one of these exists!
namespace ants
{
// extern boost::iostreams::stream<ants_Sink> std::cout;

template <class TImage>
bool IsInside( typename TImage::Pointer input, typename TImage::IndexType index )
{
  /** FIXME - should use StartIndex - */
  typedef TImage ImageType;
  enum { ImageDimension = ImageType::ImageDimension };
  bool isinside = true;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    float shifted = index[i];
    if( shifted < 0 || shifted >  input->GetLargestPossibleRegion().GetSize()[i] - 1  )
      {
      isinside = false;
      }
    }
  return isinside;
}

// ##########################################################################
// ##########################################################################
// ##########################################################################
// ##########################################################################
// ##########################################################################
// ##########################################################################
// ##########################################################################
// Templates
template <class TImage>
typename TImage::Pointer  Morphological( typename TImage::Pointer input, float rad, unsigned int option,
                                         float dilateval)
{
  typedef TImage ImageType;
  enum { ImageDimension = TImage::ImageDimension };
  typedef typename TImage::PixelType PixelType;

  if( option == 0 )
    {
    std::cout << " binary eroding the image " << std::endl;
    }
  else if( option == 1 )
    {
    std::cout << " binary dilating the image " << std::endl;
    }
  else if( option == 2 )
    {
    std::cout << " binary opening the image " << std::endl;
    }
  else if( option == 3 )
    {
    std::cout << " binary closing the image " << std::endl;
    }
  else if( option == 4 )
    {
    std::cout << " grayscale eroding the image " << std::endl;
    }
  else if( option == 5 )
    {
    std::cout << " grayscale dilating the image " << std::endl;
    }
  else if( option == 6 )
    {
    std::cout << " grayscale opening the image " << std::endl;
    }
  else if( option == 7 )
    {
    std::cout << " grayscale closing the image " << std::endl;
    }

  typedef itk::BinaryBallStructuringElement<
      PixelType,
      ImageDimension>             StructuringElementType;

  typedef itk::BinaryErodeImageFilter<
      TImage,
      TImage,
      StructuringElementType>  ErodeFilterType;

  typedef itk::BinaryDilateImageFilter<
      TImage,
      TImage,
      StructuringElementType>  DilateFilterType;
  typename ErodeFilterType::Pointer  binaryErode  = ErodeFilterType::New();
  typename DilateFilterType::Pointer binaryDilate = DilateFilterType::New();

  StructuringElementType structuringElement;

  structuringElement.SetRadius( (unsigned long) rad );  // 3x3x3 structuring element

  structuringElement.CreateStructuringElement();

  binaryErode->SetKernel(  structuringElement );
  binaryDilate->SetKernel( structuringElement );
  // binaryOpen->SetKernal( structuringElement );
  // binaryClose->SetKernel( structuringElement );
  //
  // typename OpeningFilterType::Pointer  binaryOpen  = OpeningFilterType::New();
  // typename ClosingFilterType::Pointer binaryClose = ClosingFilterType::New();
  typedef itk::GrayscaleErodeImageFilter<
      TImage,
      TImage,
      StructuringElementType> GrayscaleErodeFilterType;

  typedef itk::GrayscaleDilateImageFilter<
      TImage,
      TImage,
      StructuringElementType> GrayscaleDilateFilterType;

  typename GrayscaleErodeFilterType::Pointer grayscaleErode = GrayscaleErodeFilterType::New();
  typename GrayscaleDilateFilterType::Pointer grayscaleDilate = GrayscaleDilateFilterType::New();
  grayscaleErode->SetKernel( structuringElement );
  grayscaleDilate->SetKernel( structuringElement );

  //  It is necessary to define what could be considered objects on the binary
  //  images. This is specified with the methods \code{SetErodeValue()} and
  //  \code{SetDilateValue()}. The value passed to these methods will be
  //  considered the value over which the dilation and erosion rules will apply
  binaryErode->SetErodeValue( (unsigned int ) dilateval );
  binaryDilate->SetDilateValue(  (unsigned int ) dilateval );

  typename TImage::Pointer temp;
  if( option == 1 )
    {
    std::cout << " Dilate " << rad << std::endl;
    binaryDilate->SetInput( input );
    binaryDilate->Update();
    temp = binaryDilate->GetOutput();
    }
  else if( option == 0 )
    {
    std::cout << " Erode " << rad << std::endl;
    binaryErode->SetInput( input );  // binaryDilate->GetOutput() );
    binaryErode->Update();
    temp = binaryErode->GetOutput();
    }
  else if( option == 2 )
    {
    // dilate(erode(img))
    std::cout << " Binary Open " << rad << std::endl;
    // binaryOpen->SetInput( input );//binaryDilate->GetOutput() );
    // binaryOpen->Update();
    binaryErode->SetInput( input );
    binaryDilate->SetInput( binaryErode->GetOutput() );
    binaryDilate->Update();
    temp = binaryDilate->GetOutput();
    }
  else if( option == 3 )
    {
    std::cout << " Binary Close " << rad << std::endl;
    // binaryClose->SetInput( input );//binaryDilate->GetOutput() );
    // binaryClose->Update();
    binaryDilate->SetInput( input );
    binaryErode->SetInput( binaryDilate->GetOutput() );
    binaryErode->Update();
    temp = binaryErode->GetOutput();
    }
  else if( option == 4 )
    {
    std::cout << " Grayscale Erode " << rad << std::endl;
    grayscaleErode->SetInput( input ); // binaryDilate->GetOutput() );
    grayscaleErode->Update();
    temp = grayscaleErode->GetOutput();
    }
  else if( option == 5 )
    {
    std::cout << " Grayscale Dilate " << rad << std::endl;
    grayscaleDilate->SetInput( input ); // binaryDilate->GetOutput() );
    grayscaleDilate->Update();
    temp = grayscaleDilate->GetOutput();
    std::cout << " Grayscale Dilate Done " << temp << std::endl;
    }
  else if( option == 6 )
    {
    std::cout << " Grayscale Open " << rad << std::endl;
    grayscaleErode->SetInput( input ); // binaryDilate->GetOutput() );
    grayscaleErode->Update();
    grayscaleDilate->SetInput( grayscaleErode->GetOutput() );
    grayscaleDilate->Update();
    temp = grayscaleDilate->GetOutput();
    }
  else if( option == 7 )
    {
    std::cout << " Grayscale Close " << rad << std::endl;
    grayscaleDilate->SetInput( input ); // binaryDilate->GetOutput() );
    grayscaleDilate->Update();
    grayscaleErode->SetInput( grayscaleDilate->GetOutput() );
    grayscaleErode->Update();
    temp = grayscaleErode->GetOutput();
    }

  if( option == 0 )
    {
    // FIXME - replace with threshold filter?
    typedef itk::ImageRegionIteratorWithIndex<ImageType> ImageIteratorType;
    ImageIteratorType o_iter( temp, temp->GetLargestPossibleRegion() );
    o_iter.GoToBegin();
    while( !o_iter.IsAtEnd() )
      {
      if( o_iter.Get() > 0.5 && input->GetPixel(o_iter.GetIndex() ) > 0.5 )
        {
        o_iter.Set(1);
        }
      else
        {
        o_iter.Set(0);
        }
      ++o_iter;
      }
    }

  return temp;
}

#if 0
// TODO:  I am pretty sure that this can be completely
// replaced by the Morphological template above
// with option = true, flase, and
template <class TImage>
typename TImage::Pointer  MorphologicalBinary( typename TImage::Pointer input, float rad, bool option)
{
  typedef TImage ImageType;
  enum { ImageDimension = TImage::ImageDimension };
  typedef typename TImage::PixelType PixelType;

  if( !option )
    {
    std::cout << " eroding the image " << std::endl;
    }
  else
    {
    std::cout << " dilating the image " << std::endl;
    }
  typedef itk::BinaryBallStructuringElement<
      PixelType,
      ImageDimension>             StructuringElementType;

  typedef itk::BinaryErodeImageFilter<
      TImage,
      TImage,
      StructuringElementType>  ErodeFilterType;

  typedef itk::BinaryDilateImageFilter<
      TImage,
      TImage,
      StructuringElementType>  DilateFilterType;

  typename ErodeFilterType::Pointer  binaryErode  = ErodeFilterType::New();
  typename DilateFilterType::Pointer binaryDilate = DilateFilterType::New();

  StructuringElementType structuringElement;

  structuringElement.SetRadius( (unsigned long) rad );  // 3x3x3 structuring element

  structuringElement.CreateStructuringElement();

  binaryErode->SetKernel(  structuringElement );
  binaryDilate->SetKernel( structuringElement );

  //  It is necessary to define what could be considered objects on the binary
  //  images. This is specified with the methods \code{SetErodeValue()} and
  //  \code{SetDilateValue()}. The value passed to these methods will be
  //  considered the value over which the dilation and erosion rules will apply
  binaryErode->SetErodeValue( 1 );
  binaryDilate->SetDilateValue( 1 );

  typename TImage::Pointer temp;
  if( option )
    {
    binaryDilate->SetInput( input );
    binaryDilate->Update();
    temp = binaryDilate->GetOutput();
    }
  else
    {
    binaryErode->SetInput( input );  // binaryDilate->GetOutput() );
    binaryErode->Update();
    temp = binaryErode->GetOutput();

    typedef itk::ImageRegionIteratorWithIndex<ImageType> ImageIteratorType;
    ImageIteratorType o_iter( temp, temp->GetLargestPossibleRegion() );
    o_iter.GoToBegin();
    while( !o_iter.IsAtEnd() )
      {
      if( o_iter.Get() > 0.5 && input->GetPixel(o_iter.GetIndex() ) > 0.5 )
        {
        o_iter.Set(1);
        }
      else
        {
        o_iter.Set(0);
        }
      ++o_iter;
      }
    }

  return temp;
}

#endif

template <class TImage>
typename TImage::Pointer BinaryThreshold(
  typename TImage::PixelType low,
  typename TImage::PixelType high,
  typename TImage::PixelType replaceval, typename TImage::Pointer input)
{
  typedef typename TImage::PixelType PixelType;
  // Begin Threshold Image
  typedef itk::BinaryThresholdImageFilter<TImage, TImage> InputThresholderType;
  typename InputThresholderType::Pointer inputThresholder =
    InputThresholderType::New();

  inputThresholder->SetInput( input );
  inputThresholder->SetInsideValue(  replaceval );
  int outval = 0;
  if( (float) replaceval == (float) -1 )
    {
    outval = 1;
    }
  inputThresholder->SetOutsideValue( outval );

  if( high < low )
    {
    high = 255;
    }
  inputThresholder->SetLowerThreshold( (PixelType) low );
  inputThresholder->SetUpperThreshold( (PixelType) high );
  inputThresholder->Update();

  return inputThresholder->GetOutput();
}

template <class TPixel, unsigned int VDim>
class VectorPixelCompare
{
public:
  bool operator()( const itk::Vector<TPixel, VDim> & v1,
                   const itk::Vector<TPixel, VDim> & v2 ) const
  {
    // Ordering of vectors based on 1st component, then second, etc.
    for( size_t i = 0; i < VDim; i++ )
      {
      if( v1[i] < v2[i] )
        {
        return true;
        }
      else if( v1[i] > v2[i] )
        {
        return false;
        }
      }
    return false;
  }
};

template <class ImageType, class AffineTransform>
void GetAffineTransformFromImage(const typename ImageType::Pointer& img,
                                 typename AffineTransform::Pointer & aff)
{
  typedef typename ImageType::DirectionType         DirectionType;
  typedef typename ImageType::PointType             PointType;
  typedef typename ImageType::SpacingType           SpacingType;
  typedef typename AffineTransform::TranslationType VectorType;

  DirectionType direction = img->GetDirection();

  VectorType translation;
  // translation.Fill(0);
  for( unsigned int i = 0; i < ImageType::GetImageDimension(); i++ )
    {
    translation[i] = img->GetOrigin()[i];
    }

  aff->SetMatrix(direction);
  // aff->SetCenter(pt);
  PointType pt; pt.Fill(0);
  aff->SetOffset(translation);
  aff->SetCenter(pt);

  std::cout << "aff from image:" << aff << std::endl;
}

template <class WarperType, class ImageType>
void GetLargestSizeAfterWarp(typename WarperType::Pointer & warper,
                             typename ImageType::Pointer & img,
                             typename ImageType::SizeType & largest_size,
                             typename ImageType::PointType & origin_warped)
{
  typedef typename ImageType::SizeType  SizeType;
  typedef typename ImageType::PointType PointType;

  const int ImageDimension = ImageType::GetImageDimension();
  // typedef typename ImageType::PointType PointType;
  typedef typename std::vector<PointType> PointList;

  typedef typename ImageType::IndexType IndexType;

  // PointList pts_orig;
  PointList pts_warped;

  typename ImageType::SizeType imgsz;
  imgsz = img->GetLargestPossibleRegion().GetSize();

  typename ImageType::SpacingType spacing;
  spacing = img->GetSpacing();

  pts_warped.clear();
  if( ImageDimension == 3 )
    {
    for( int i = 0; i < 8; i++ )
      {
      IndexType ind;

      switch( i )
        {
        case 0:
          {
          ind[0] = 0; ind[1] = 0; ind[2] = 0;
          }
          break;
        case 1:
          {
          ind[0] = imgsz[0] - 1; ind[1] = 0; ind[2] = 0;
          }
          break;
        case 2:
          {
          ind[0] = 0; ind[1] = imgsz[1] - 1; ind[2] = 0;
          }
          break;
        case 3:
          {
          ind[0] = imgsz[0] - 1; ind[1] = imgsz[1] - 1; ind[2] = 0;
          }
          break;
        case 4:
          {
          ind[0] = 0; ind[1] = 0; ind[2] = imgsz[2] - 1;
          }
          break;
        case 5:
          {
          ind[0] = imgsz[0] - 1; ind[1] = 0; ind[2] = imgsz[2] - 1;
          }
          break;
        case 6:
          {
          ind[0] = 0; ind[1] = imgsz[1] - 1; ind[2] = imgsz[2] - 1;
          }
          break;
        case 7:
          {
          ind[0] = imgsz[0] - 1; ind[1] = imgsz[1] - 1; ind[2] = imgsz[2] - 1;
          }
          break;
        }
      PointType pt_orig, pt_warped;
      img->TransformIndexToPhysicalPoint(ind, pt_orig);
      if( warper->MultiInverseAffineOnlySinglePoint(pt_orig, pt_warped) == false )
        {
        std::cout << "ERROR: outside of numeric boundary with affine transform." << std::endl;
        throw std::exception();
        }
      pts_warped.push_back(pt_warped);
      std::cout << '[' << i << ']' << ind << ',' << pt_orig << "->" << pt_warped << std::endl;
      }
    }
  else if( ImageDimension == 2 )
    {
    for( int i = 0; i < 4; i++ )
      {
      IndexType ind;

      switch( i )
        {
        case 0:
          {
          ind[0] = 0; ind[1] = 0;
          }
          break;
        case 1:
          {
          ind[0] = imgsz[0] - 1; ind[1] = 0;
          }
          break;
        case 2:
          {
          ind[0] = 0; ind[1] = imgsz[1] - 1;
          }
          break;
        case 3:
          {
          ind[0] = imgsz[0] - 1; ind[1] = imgsz[1] - 1;
          }
          break;
        }
      PointType pt_orig, pt_warped;
      img->TransformIndexToPhysicalPoint(ind, pt_orig);
      if( warper->MultiInverseAffineOnlySinglePoint(pt_orig, pt_warped) == false )
        {
        std::cout << "ERROR: outside of numeric boundary with affine transform." << std::endl;
        throw std::exception();
        }
      pts_warped.push_back(pt_warped);
      std::cout << '[' << i << ']' << ind << ',' << pt_orig << "->" << pt_warped << std::endl;
      }
    }
  else
    {
    std::cout << "could not determine the dimension after warping for non 2D/3D volumes" << std::endl;
    throw std::exception();
    }

  PointType pt_min, pt_max;
  pt_min = pts_warped[0];
  pt_max = pts_warped[0];
  for( unsigned int k = 0; k < pts_warped.size(); k++ )
    {
    for( int i = 0; i < ImageDimension; i++ )
      {
      pt_min[i] = (pt_min[i] < pts_warped[k][i]) ? (pt_min[i]) : (pts_warped[k][i]);
      pt_max[i] = (pt_max[i] > pts_warped[k][i]) ? (pt_max[i]) : (pts_warped[k][i]);
      }
    }
  for( int i = 0; i < ImageDimension; i++ )
    {
    largest_size[i] = (int) (ceil( (pt_max[i] - pt_min[i]) / spacing[i]) + 1);
    }

  origin_warped = pt_min;
  std::cout << "origin_warped: " << origin_warped << std::endl;
  std::cout << "pt_min: " << pt_min << " pt_max:" << pt_max << " largest_size:" << largest_size << std::endl;
}

template <class TImageIn, class TImageOut>
typename TImageOut::Pointer
arCastImage( typename TImageIn::Pointer Rimage )
{
  typedef itk::CastImageFilter<TImageIn, TImageOut> CastFilterType;
  typename CastFilterType::Pointer caster = CastFilterType::New();
  caster->SetInput( Rimage );
  caster->Update();
  return caster->GetOutput();
}

} // end namespace

// ##########################################################################
// TODO: KENT:  This block feels like it could be better encapsulated as a c++ class
//
typedef enum
  {
  INVALID_FILE = 1,
  AFFINE_FILE,
  DEFORMATION_FILE,
  IMAGE_AFFINE_HEADER,
  IDENTITY_TRANSFORM
  } TRAN_FILE_TYPE;

// TODO: This should be a class.
typedef struct
  {
  //    char *filename;
  std::string filename;
  TRAN_FILE_TYPE file_type;
  bool do_affine_inv;
  //    void SetValue(char *filename, TRAN_FILE_TYPE file_type, bool do_affine_inv){
  //        this.filename = filename;
  //        this.file_type = file_type;
  //        this.do_affine_inv = do_affine_inv;
  //    };
  double weight;   // for average
  } TRAN_OPT;

typedef std::vector<TRAN_OPT> TRAN_OPT_QUEUE;

typedef struct
  {
  bool physical_units;
  std::vector<double> sigma;
  } MLINTERP_OPT;

typedef struct
  {
  bool use_NN_interpolator;
  bool use_MultiLabel_interpolator;
  bool use_BSpline_interpolator;
  bool use_TightestBoundingBox;
  char * reference_image_filename;
  bool use_RotationHeader;

  MLINTERP_OPT opt_ML;
  } MISC_OPT;

extern TRAN_FILE_TYPE CheckFileType(const char * const str);

extern TRAN_FILE_TYPE CheckFileType(const std::string & str);

extern void SetAffineInvFlag(TRAN_OPT & opt, bool & set_current_affine_inv);

extern void DisplayOptQueue(const TRAN_OPT_QUEUE & opt_queue);

extern void DisplayOpt(const TRAN_OPT & opt);

// ##########################################################################

extern bool get_a_double_number(const char * const str, double & v);

// TODO: KENT:  These two functions have cross-platform-equivalent versions from kwSys and could be replaced.
extern void FilePartsWithgz(const std::string & filename, std::string & path, std::string & name, std::string & ext);

extern bool CheckFileExistence(const char * const str);

extern std::string GetPreferredTransformFileType(void);

extern void ConvertToLowerCase( std::string& str );

#endif // __antsUtilities_h__
