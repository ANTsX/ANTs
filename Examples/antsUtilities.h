/**
 * There were several functions that had been copied and
 * pasted over and over again in this library.  This
 * header files is contains a common definitoin for
 * those file.
 * \author Hans J. Johnson
 */
#ifndef __antsUtilities_h__
#define __antsUtilities_h__

#include "antscout.hxx"
#include <string>
#include <vector>
#include <algorithm>

#include <stdlib.h>
#include <errno.h>
#include <cmath>

#include "itkVector.h"
#include "itkBinaryThresholdImageFilter.h"

// We need to ensure that only one of these exists!
namespace ants
{
extern boost::iostreams::stream<ants_Sink> antscout;

// ##########################################################################
// ##########################################################################
// ##########################################################################
// ##########################################################################
// ##########################################################################
// ##########################################################################
// ##########################################################################
// Templates

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
                   const itk::Vector<TPixel, VDim> & v2 )
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

template <class ImageTypePointer, class AffineTransformPointer>
void GetAffineTransformFromImage(const ImageTypePointer& img, AffineTransformPointer & aff)
{
  typedef typename ImageTypePointer::ObjectType                        ImageType;
  typedef typename ImageType::DirectionType                            DirectionType;
  typedef typename ImageType::PointType                                PointType;
  typedef typename ImageType::SpacingType                              SpacingType;
  typedef typename AffineTransformPointer::ObjectType::TranslationType VectorType;

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

  antscout << "aff from image:" << aff << std::endl;
}

template <class WarperPointerType, class ImagePointerType, class SizeType, class PointType>
void GetLaregstSizeAfterWarp(WarperPointerType & warper, ImagePointerType & img, SizeType & largest_size,
                             PointType & origin_warped)
{
  typedef typename ImagePointerType::ObjectType ImageType;
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
  { ind[0] = 0; ind[1] = 0; ind[2] = 0; }
                                        break;
        case 1:
  { ind[0] = imgsz[0] - 1; ind[1] = 0; ind[2] = 0; }
                                                   break;
        case 2:
  { ind[0] = 0; ind[1] = imgsz[1] - 1; ind[2] = 0; }
                                                   break;
        case 3:
  { ind[0] = imgsz[0] - 1; ind[1] = imgsz[1] - 1; ind[2] = 0; }
                                                              break;
        case 4:
  { ind[0] = 0; ind[1] = 0; ind[2] = imgsz[2] - 1; }
                                                   break;
        case 5:
  { ind[0] = imgsz[0] - 1; ind[1] = 0; ind[2] = imgsz[2] - 1; }
                                                              break;
        case 6:
  { ind[0] = 0; ind[1] = imgsz[1] - 1; ind[2] = imgsz[2] - 1; }
                                                              break;
        case 7:
  { ind[0] = imgsz[0] - 1; ind[1] = imgsz[1] - 1; ind[2] = imgsz[2] - 1; }
                                                                         break;
        }
      PointType pt_orig, pt_warped;
      img->TransformIndexToPhysicalPoint(ind, pt_orig);
      if( warper->MultiInverseAffineOnlySinglePoint(pt_orig, pt_warped) == false )
        {
        antscout << "ERROR: outside of numeric boundary with affine transform." << std::endl;
        throw std::exception();
        }
      pts_warped.push_back(pt_warped);
      antscout << '[' << i << ']' << ind << ',' << pt_orig << "->" << pt_warped << std::endl;
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
  { ind[0] = 0; ind[1] = 0; }
                            break;
        case 1:
  { ind[0] = imgsz[0] - 1; ind[1] = 0; }
                                       break;
        case 2:
  { ind[0] = 0; ind[1] = imgsz[1] - 1; }
                                       break;
        case 3:
  { ind[0] = imgsz[0] - 1; ind[1] = imgsz[1] - 1; }
                                                  break;
        }
      PointType pt_orig, pt_warped;
      img->TransformIndexToPhysicalPoint(ind, pt_orig);
      if( warper->MultiInverseAffineOnlySinglePoint(pt_orig, pt_warped) == false )
        {
        antscout << "ERROR: outside of numeric boundary with affine transform." << std::endl;
        throw std::exception();
        }
      pts_warped.push_back(pt_warped);
      antscout << '[' << i << ']' << ind << ',' << pt_orig << "->" << pt_warped << std::endl;
      }
    }
  else
    {
    antscout << "could not determine the dimension after warping for non 2D/3D volumes" << std::endl;
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
  antscout << "origin_warped: " << origin_warped << std::endl;
  antscout << "pt_min: " << pt_min << " pt_max:" << pt_max << " largest_size:" << largest_size << std::endl;
}
}

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

extern TRAN_FILE_TYPE CheckFileType(const std::string str);

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
