#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>
#include "antsCommandLineOption.h"
#include "antsCommandLineParser.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"
#include <sstream>
#include <iostream>
#include <random>
#include <string>
#include <algorithm>
#include <vector>
#include <vnl/vnl_random.h>
#include <vnl/algo/vnl_qr.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_svd_economy.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_real_eigensystem.h>
#include <vnl/algo/vnl_generalized_eigensystem.h>
#include "antsSCCANObject.h"
#include "itkCSVNumericObjectFileWriter.h"
#include "itkCSVArray2DDataObject.h"
#include "itkCSVArray2DFileReader.h"
#include "itkExtractImageFilter.h"
#include "ReadWriteData.h"

namespace ants
{
// namespace antssccan {
template <typename TComp>
double
vnl_pearson_corr(vnl_vector<TComp> v1, vnl_vector<TComp> v2)
{
  double xysum = 0;

  for (unsigned int i = 0; i < v1.size(); i++)
  {
    xysum += v1(i) * v2(i);
  }
  double frac = 1.0 / (double)v1.size();
  double xsum = v1.sum(), ysum = v2.sum();
  double xsqr = v1.squared_magnitude();
  double ysqr = v2.squared_magnitude();
  double numer = xysum - frac * xsum * ysum;
  double denom = sqrt((xsqr - frac * xsum * xsum) * (ysqr - frac * ysum * ysum));
  if (denom <= 0)
  {
    return 0;
  }
  return numer / denom;
}

template <typename TImage, typename TComp>
void
WriteVectorToSpatialImage(std::string filename, std::string post, vnl_vector<TComp> w_p, typename TImage::Pointer mask)
{
  std::string::size_type pos = filename.rfind(".");
  std::string            filepre = std::string(filename, 0, pos);
  std::string            extension;
  if (pos != std::string::npos)
  {
    extension = std::string(filename, pos, filename.length() - 1);
    if (extension == std::string(".gz"))
    {
      pos = filepre.rfind(".");
      extension = std::string(filepre, pos, filepre.length() - 1) + extension;
      filepre = std::string(filepre, 0, pos);
    }
  }

  using SCCANType = itk::ants::antsSCCANObject<TImage, TComp>;
  typename SCCANType::Pointer sccanobj = SCCANType::New();
  typename TImage::Pointer    weights = sccanobj->ConvertVariateToSpatialImage(w_p, mask);
  std::string                 fn1 = filepre + post + extension;
  ANTs::WriteImage<TImage>(weights, fn1.c_str());
}

template <typename T>
inline std::string
sccan_to_string(const T & t)
{
  std::stringstream ss;

  ss << t;
  std::string stringout = ss.str();
  if (t < 100)
  {
    std::stringstream ss0;
    ss0 << 0;
    std::string extend = ss0.str();
    stringout = std::string(extend + stringout);
  }
  if (t < 10)
  {
    std::stringstream ss0;
    ss0 << 0;
    std::string extend = ss0.str();
    stringout = std::string(extend + stringout);
  }
  return stringout;
}

template <typename TImage, typename TComp>
void
WriteSortedVariatesToSpatialImage(std::string              filename,
                                  std::string              post,
                                  vnl_matrix<TComp>        varmat,
                                  typename TImage::Pointer mask,
                                  vnl_matrix<TComp>        data_mat,
                                  bool                     have_mask,
                                  vnl_vector<TComp>        l_array,
                                  vnl_matrix<TComp>        prior_mat)
{
  vnl_matrix<TComp> projections = data_mat * varmat;

  std::string::size_type pos = filename.rfind(".");
  std::string            filepre = std::string(filename, 0, pos);
  std::string            extension;
  if (pos != std::string::npos)
  {
    extension = std::string(filename, pos, filename.length() - 1);
    if (extension == std::string(".gz"))
    {
      pos = filepre.rfind(".");
      extension = std::string(filepre, pos, filepre.length() - 1) + extension;
      filepre = std::string(filepre, 0, pos);
    }
  }
  std::string   post2;
  std::ofstream myfile;
  std::string   fnmp1 = filepre + std::string("projections") + post + std::string(".csv");

  std::vector<std::string> ColumnHeaders1;
  for (unsigned int nv = 0; nv < projections.cols(); nv++)
  {
    std::string colname = std::string("Variate") + sccan_to_string<unsigned int>(nv);
    ColumnHeaders1.push_back(colname);
  }
  using WriterType = itk::CSVNumericObjectFileWriter<double>;
  WriterType::Pointer writer1 = WriterType::New();
  writer1->SetFileName(fnmp1.c_str());
  writer1->SetColumnHeaders(ColumnHeaders1);
  writer1->SetInput(&projections);
  try
  {
    writer1->Write();
  }
  catch (const itk::ExceptionObject & exp)
  {
    // std::cerr << "Exception caught!" << std::endl;
    // std::cerr << exp << std::endl;
    return;
  }

  // std::cout << data_mat.cols() << prior_mat.cols() << data_mat.rows() << prior_mat.rows() << std::endl;
  vnl_matrix<TComp> projectionsROI = data_mat * prior_mat.transpose();

  std::string::size_type posROI = filename.rfind(".");
  std::string            filepreROI = std::string(filename, 0, posROI);
  std::string            extensionROI;
  if (posROI != std::string::npos)
  {
    extensionROI = std::string(filename, posROI, filename.length() - 1);
    if (extension == std::string(".gz"))
    {
      posROI = filepreROI.rfind(".");
      extensionROI = std::string(filepreROI, posROI, filepreROI.length() - 1) + extensionROI;
      filepreROI = std::string(filepreROI, 0, posROI);
    }
  }
  // std::string post2;
  // std::ofstream myfile;
  std::string fnmp_prior = filepreROI + std::string("projectionsROI") + post + std::string(".csv");

  std::vector<std::string> ColumnHeadersROI;
  for (unsigned int nv = 0; nv < projectionsROI.cols(); nv++)
  {
    std::string colnameROI = std::string("Variate") + sccan_to_string<unsigned int>(nv);
    ColumnHeadersROI.push_back(colnameROI);
  }
  using WriterType = itk::CSVNumericObjectFileWriter<double>;
  WriterType::Pointer writerROI = WriterType::New();
  writerROI->SetFileName(fnmp_prior.c_str());
  writerROI->SetColumnHeaders(ColumnHeadersROI);
  writerROI->SetInput(&projectionsROI);
  try
  {
    writerROI->Write();
  }
  catch (const itk::ExceptionObject & exp)
  {
    // std::cerr << "Exception caught!" << std::endl;
    // std::cerr << exp << std::endl;
    return;
  }

  if (have_mask)
  {
    // std::cout << " have_mask WriteSortedVariatesToSpatialImage " << have_mask << std::endl;
    for (unsigned int vars = 0; vars < varmat.columns(); vars++)
    {
      post2 = post + sccan_to_string<unsigned int>(l_array.get(vars));
      vnl_vector<TComp> temp = varmat.get_column(l_array.get(vars));
      // // std::cout<<" File Name "<<filename<<" POST "<<l_array.get(vars)<<std::endl;
      WriteVectorToSpatialImage<TImage, TComp>(filename, post2, temp, mask);
    }
  }
  else
  {
    std::vector<std::string> ColumnHeaders;
    // write out the array2D object
    std::string fnmp = filepre + std::string("ViewVecs") + std::string(".csv");
    for (unsigned int nv = 0; nv < varmat.cols(); nv++)
    {
      std::string colname = std::string("Variate") + sccan_to_string<unsigned int>(nv);
      ColumnHeaders.push_back(colname);
    }
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(fnmp.c_str());
    writer->SetColumnHeaders(ColumnHeaders);
    writer->SetInput(&varmat);
    try
    {
      writer->Write();
    }
    catch (const itk::ExceptionObject & exp)
    {
      // std::cerr << "Exception caught!" << std::endl;
      // std::cerr << exp << std::endl;
      return;
    }
  }
}

template <typename TImage, typename TComp>
void
WriteVariatesToSpatialImage(std::string              filename,
                            std::string              post,
                            vnl_matrix<TComp>        varmat,
                            typename TImage::Pointer mask,
                            vnl_matrix<TComp>        data_mat,
                            bool                     have_mask,
                            vnl_matrix<TComp>        u_mat)
{
  vnl_matrix<TComp>      projections = data_mat * varmat;
  std::string::size_type pos = filename.rfind(".");
  std::string            filepre = std::string(filename, 0, pos);
  std::string            extension;
  if (pos != std::string::npos)
  {
    extension = std::string(filename, pos, filename.length() - 1);
    if (extension == std::string(".gz"))
    {
      pos = filepre.rfind(".");
      extension = std::string(filepre, pos, filepre.length() - 1) + extension;
      filepre = std::string(filepre, 0, pos);
    }
  }
  std::string              post2;
  std::ofstream            myfile;
  std::string              fnmp = filepre + std::string("projections") + post + std::string(".csv");
  std::vector<std::string> ColumnHeaders;
  for (unsigned int nv = 0; nv < projections.cols(); nv++)
  {
    std::string colname = std::string("Variate") + sccan_to_string<unsigned int>(nv);
    ColumnHeaders.push_back(colname);
  }
  using WriterType = itk::CSVNumericObjectFileWriter<double, 1, 1>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(fnmp.c_str());
  writer->SetColumnHeaders(ColumnHeaders);
  writer->SetInput(&projections);
  try
  {
    writer->Write();
  }
  catch (const itk::ExceptionObject & itkNotUsed(exp))
  {
    // std::cerr << "Exception caught!" << std::endl;
    // std::cerr << exp << std::endl;
    return;
  }
  if (have_mask)
  {
    // std::cerr << " have_mask " << have_mask << std::endl;
    for (unsigned int vars = 0; vars < varmat.columns(); vars++)
    {
      post2 = post + sccan_to_string<unsigned int>(vars);
      vnl_vector<TComp> temp = varmat.get_column(vars);
      WriteVectorToSpatialImage<TImage, TComp>(filename, post2, temp, mask);
    }
  }
  else
  {
    ColumnHeaders.clear();
    // write out the array2D object
    fnmp = filepre + std::string("_Variate_") + post + std::string(".csv");
    for (unsigned int nv = 0; nv < varmat.cols(); nv++)
    {
      std::string colname = std::string("Variate") + sccan_to_string<unsigned int>(nv);
      ColumnHeaders.push_back(colname);
    }
    writer = WriterType::New();
    writer->SetFileName(fnmp.c_str());
    writer->SetColumnHeaders(ColumnHeaders);
    writer->SetInput(&varmat);
    try
    {
      writer->Write();
    }
    catch (const itk::ExceptionObject & itkNotUsed(exp))
    {
      // std::cerr << "Exception caught!" << std::endl;
      // std::cerr << exp << std::endl;
      return;
    }
  }

  if (!u_mat.empty())
  {
    // write out the array2D object for U matrix
    ColumnHeaders.clear();
    fnmp = filepre + std::string("_Umatrix_") + post + std::string(".csv");
    for (unsigned int nv = 0; nv < varmat.cols(); nv++)
    {
      std::string colname = std::string("U") + sccan_to_string<unsigned int>(nv);
      ColumnHeaders.push_back(colname);
    }
    writer = WriterType::New();
    writer->SetFileName(fnmp.c_str());
    writer->SetColumnHeaders(ColumnHeaders);
    writer->SetInput(&u_mat);
    try
    {
      writer->Write();
    }
    catch (const itk::ExceptionObject & itkNotUsed(exp))
    {
      // std::cerr << "Exception caught!" << std::endl;
      // std::cerr << exp << std::endl;
      return;
    }
  }
}

template <typename TImage, typename TComp>
vnl_matrix<TComp>
CopyImageToVnlMatrix(typename TImage::Pointer p_img)
{
  using vMatrix = vnl_matrix<TComp>;

  typename TImage::SizeType pMatSize = p_img->GetLargestPossibleRegion().GetSize();
  vMatrix                   p(pMatSize[0], pMatSize[1]); // a (size)x(size+1)-matrix of int's
  for (unsigned long j = 0; j < p.columns(); ++j)        // loop over columns
  {
    for (unsigned long i = 0; i < p.rows(); ++i) // loop over rows
    {
      typename TImage::IndexType ind;
      ind[0] = i;
      ind[1] = j;
      TComp val = p_img->GetPixel(ind);
      p(i, j) = val; // to access matrix coefficients,
    }
  }
  return p;
}

template <typename TImage, typename TComp>
vnl_matrix<TComp>
DeleteRow(vnl_matrix<TComp> p_in, unsigned int row)
{
  using vMatrix = vnl_matrix<TComp>;
  unsigned int nrows = p_in.rows() - 1;
  if (row >= nrows)
  {
    nrows = p_in.rows();
  }
  vMatrix      p(nrows, p_in.columns());
  unsigned int rowct = 0;
  for (long i = 0; i < static_cast<long>(p.rows()); ++i) // loop over rows
  {
    if (i != static_cast<long>(row))
    {
      p.set_row(rowct, p_in.get_row(i));
      rowct++;
    }
  }
  return p;
}

int
sccanRandom(int n)
{
  return rand() % n;
}

template <typename TComp>
vnl_matrix<TComp>
PermuteMatrix(vnl_matrix<TComp> q, bool doperm = true)
{
  using vMatrix = vnl_matrix<TComp>;

  std::vector<unsigned long> permvec;
  for (unsigned long i = 0; i < q.rows(); i++)
  {
    permvec.push_back(i);
  }
// std::random_shuffle is deprecated since C++14,
// see: http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2014/n4190.htm
// 2023-08-30: MSVC 2019 on Github Actions runners having problems detecting C++ version,
// so use new std::shuffle
// #if __cplusplus >= 201402L
//  std::shuffle(permvec.begin(), permvec.end(), std::random_device());
// #else
//   std::random_shuffle(permvec.begin(), permvec.end(), sccanRandom);
// #endif
std::shuffle(permvec.begin(), permvec.end(), std::random_device());

  //    for (unsigned long i=0; i < q.rows(); i++)
  //  // std::cout << " permv " << i << " is " << permvec[i] << std::endl;
  // for (unsigned long i=0; i < q.rows(); i++)
  //  // std::cout << " permv " << i << " is " << permvec[i] << std::endl;
  // 1. permute q
  vMatrix q_perm(q.rows(), q.columns());
  for (unsigned long i = 0; i < q.rows(); i++)
  {
    unsigned long perm = permvec[i];
    if (doperm)
    {
      q_perm.set_row(i, q.get_row(perm));
    }
    else
    {
      q_perm.set_row(i, q.get_row(i));
    }
  }
  return q_perm;
}

template <unsigned int ImageDimension, typename PixelType>
int
matrixOperation(itk::ants::CommandLineParser::OptionType * option,
                itk::ants::CommandLineParser::OptionType * /* outputOption */ = nullptr)
{
  std::string funcName = std::string("matrixOperation");

  using ImageType = itk::Image<PixelType, ImageDimension>;
  typename ImageType::Pointer outputImage = nullptr;

  //   option->SetUsageOption( 2, "multires_matrix_invert[list.txt,maskhighres.nii.gz,masklowres.nii.gz,matrix.mhd]" );

  std::string value = option->GetFunction(0)->GetName();
  if (strcmp(value.c_str(), "multires_matrix_invert") == 0)
  {
    std::string listfn = option->GetFunction(0)->GetParameter(0);
    std::string maskhfn = option->GetFunction(0)->GetParameter(1);
    std::string masklfn = option->GetFunction(0)->GetParameter(2);
    //      vnl_matrix<matPixelType> matrixinv=MultiResMatrixInvert<ImageDimension,matPixelType>( listfn, maskhfn,
    // masklfn );
  }

  return EXIT_SUCCESS;
}

template <typename RealType>
int
CompareMatrixSizes(vnl_matrix<RealType> & p, vnl_matrix<RealType> & q)
{
  if (p.rows() != q.rows())
  {
    // std::cerr << " The number of rows must match !!" << std::endl;
    // std::cerr << " matrix-1 has " << p.rows() << " rows " << std::endl;
    // std::cerr << " matrix-2 has " << q.rows() << " rows " << std::endl;
    // std::cerr << " returning " << EXIT_FAILURE << std::endl;
    //    throw std::exception();
    return EXIT_FAILURE;
  }
  return 0;
}

template <typename PixelType>
void
ReadMatrixFromCSVorImageSet(std::string matname, vnl_matrix<PixelType> & p)
{
  using Scalar = PixelType;
  using MatrixImageType = itk::Image<PixelType, 2>;
  using matReaderType = itk::ImageFileReader<MatrixImageType>;
  std::string ext = itksys::SystemTools::GetFilenameExtension(matname);
  if (strcmp(ext.c_str(), ".csv") == 0)
  {
    using ReaderType = itk::CSVArray2DFileReader<double>;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(matname.c_str());
    reader->SetFieldDelimiterCharacter(',');
    reader->SetStringDelimiterCharacter('"');
    reader->HasColumnHeadersOn();
    reader->HasRowHeadersOff();
    reader->UseStringDelimiterCharacterOff();
    try
    {
      reader->Update();
    }
    catch (const itk::ExceptionObject & itkNotUsed(exp))
    {
      // std::cerr << "Exception caught!" << std::endl;
      // std::cerr << exp << std::endl;
    }
    using DataFrameObjectType = itk::CSVArray2DDataObject<double>;
    DataFrameObjectType::Pointer dfo = reader->GetOutput();
    p = dfo->GetMatrix();
    return;
  }
  else
  {
    typename matReaderType::Pointer matreader1 = matReaderType::New();
    matreader1->SetFileName(matname.c_str());
    matreader1->Update();
    p = CopyImageToVnlMatrix<MatrixImageType, Scalar>(matreader1->GetOutput());
  }
}

template <unsigned int ImageDimension, typename PixelType>
itk::Array2D<double>
ConvertImageListToMatrix(std::string imagelist, std::string maskfn, std::string outname)
{
  std::string ext = itksys::SystemTools::GetFilenameExtension(outname);

  using MatrixType = itk::Array2D<double>;
  std::vector<std::string> ColumnHeaders;
  MatrixType               zmat(1, 1);
  using ImageType = itk::Image<PixelType, ImageDimension>;
  using MatrixImageType = itk::Image<PixelType, 2>;
  using ReaderType = itk::ImageFileReader<ImageType>;
  typename ImageType::Pointer mask;
  ReadImage<ImageType>(mask, maskfn.c_str());
  unsigned long voxct = 0;
  using Iterator = itk::ImageRegionIteratorWithIndex<ImageType>;
  Iterator mIter(mask, mask->GetLargestPossibleRegion());
  for (mIter.GoToBegin(); !mIter.IsAtEnd(); ++mIter)
  {
    if (mIter.Get() >= 0.5)
    {
      voxct++;
    }
  }
  std::vector<std::string> image_fn_list;
  // first, count the number of files
  constexpr unsigned int maxChar = 512;
  char                   lineBuffer[maxChar];
  char                   filenm[maxChar];
  unsigned int           filecount = 0;
  {
    std::ifstream inputStreamA(imagelist.c_str(), std::ios::in);
    if (!inputStreamA.is_open())
    {
      // std::cout << "Can't open image list file: " << imagelist << std::endl;
      return zmat;
    }
    while (!inputStreamA.eof())
    {
      inputStreamA.getline(lineBuffer, maxChar, '\n');
      if (sscanf(lineBuffer, "%s ", filenm) != 1)
      {
        continue;
      }
      else
      {
        image_fn_list.emplace_back(filenm);
        filecount++;
      }
    }

    inputStreamA.close();
  }

  /** declare the output matrix image */
  unsigned long                      xsize = image_fn_list.size();
  unsigned long                      ysize = voxct;
  typename MatrixImageType::SizeType tilesize;
  tilesize[0] = xsize;
  tilesize[1] = ysize;

  //      // std::cout <<" have voxct " << voxct << " and nsub " << filecount << " or " << image_fn_list.size()<<
  // std::endl;

  MatrixType matrix(xsize, ysize);
  matrix.Fill(0);
  for (unsigned int j = 0; j < image_fn_list.size(); j++)
  {
    typename ReaderType::Pointer reader2 = ReaderType::New();
    reader2->SetFileName(image_fn_list[j]);
    reader2->Update();
    unsigned long xx = 0, yy = 0, tvoxct = 0;
    xx = j;
    for (mIter.GoToBegin(); !mIter.IsAtEnd(); ++mIter)
    {
      if (mIter.Get() >= 0.5)
      {
        yy = tvoxct;
        matrix[xx][yy] = reader2->GetOutput()->GetPixel(mIter.GetIndex());
        if (j == 0)
        {
          std::string colname = std::string("V") + sccan_to_string<unsigned long>(tvoxct);
          ColumnHeaders.push_back(colname);
        }
        tvoxct++;
      }
    }
  }

  if (strcmp(ext.c_str(), ".csv") == 0)
  {
    // write out the array2D object
    using WriterType = itk::CSVNumericObjectFileWriter<double, 1, 1>;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outname);
    writer->SetInput(&matrix);
    writer->SetColumnHeaders(ColumnHeaders);
    try
    {
      writer->Write();
    }
    catch (const itk::ExceptionObject & itkNotUsed(exp))
    {
      // std::cerr << "Exception caught!" << std::endl;
      // std::cerr << exp << std::endl;
      return matrix;
    }
    return matrix;
  }

  if (strcmp(ext.c_str(), ".mha") == 0 || strcmp(ext.c_str(), ".mhd") == 0)
  {
    typename MatrixImageType::RegionType region;
    region.SetSize(tilesize);
    typename MatrixImageType::Pointer matimage = AllocImage<MatrixImageType>(region);
    for (unsigned int j = 0; j < image_fn_list.size(); j++)
    {
      typename ReaderType::Pointer reader2 = ReaderType::New();
      reader2->SetFileName(image_fn_list[j]);
      reader2->Update();
      unsigned long xx = 0, yy = 0, tvoxct = 0;
      xx = j;
      typename MatrixImageType::IndexType mind;
      for (mIter.GoToBegin(); !mIter.IsAtEnd(); ++mIter)
      {
        if (mIter.Get() >= 0.5)
        {
          yy = tvoxct;
          mind[0] = xx;
          mind[1] = yy;
          matimage->SetPixel(mind, reader2->GetOutput()->GetPixel(mIter.GetIndex()));
          tvoxct++;
        }
      }
    }

    using WriterType = itk::ImageFileWriter<MatrixImageType>;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outname);
    writer->SetInput(matimage);
    writer->Update();
  }
  return matrix;
}

template <typename PixelType>
int
ConvertTimeSeriesImageToMatrix(std::string imagefn,
                               std::string maskfn,
                               std::string outname,
                               double      space_smoother,
                               double      time_smoother)
{
  constexpr unsigned int ImageDimension = 4;

  using ImageType = itk::Image<PixelType, ImageDimension>;
  using OutImageType = itk::Image<PixelType, ImageDimension - 1>;
  using OutIndexType = typename OutImageType::IndexType;
  using IndexType = typename ImageType::IndexType;

  using Scalar = double;
  std::string ext = itksys::SystemTools::GetFilenameExtension(outname);
  if (strcmp(ext.c_str(), ".csv") != 0)
  {
    // std::cerr << " must use .csv as output file extension " << std::endl;
    return EXIT_FAILURE;
  }
  typename ImageType::Pointer    image1 = nullptr;
  typename OutImageType::Pointer mask = nullptr;
  // std::cerr << " imagefn " << imagefn << std::endl;
  if (imagefn.length() > 3)
  {
    ReadImage<ImageType>(image1, imagefn.c_str());
  }
  else
  {
    // std::cerr << " cannot read image " << imagefn << std::endl;
    return 1;
  }

  if (space_smoother > 0)
  {
    typename ImageType::SpacingType spacing = image1->GetSpacing();
    typename ImageType::SpacingType spacing2 = image1->GetSpacing();
    // basically, don't do any dim-4 smoothing
    spacing2[3] = sqrt(spacing[0] * spacing[0] + spacing[1] * spacing[1] + spacing[2] * spacing[2]) * 1.e6;
    image1->SetSpacing(spacing2);
    using dgf = itk::DiscreteGaussianImageFilter<ImageType, ImageType>;
    typename dgf::Pointer filter = dgf::New();
    filter->SetVariance(space_smoother);
    filter->SetUseImageSpacing(true);
    filter->SetMaximumError(.01f);
    filter->SetInput(image1);
    filter->Update();
    image1 = filter->GetOutput();
    image1->SetSpacing(spacing);
  }

  if (time_smoother > 0)
  {
    typename ImageType::SpacingType spacing = image1->GetSpacing();
    typename ImageType::SpacingType spacing2 = image1->GetSpacing();
    // basically, don't do any dim-4 smoothing
    double bigspace = sqrt(spacing[0] * spacing[0] + spacing[1] * spacing[1] + spacing[2] * spacing[2]) * 1.e6;
    // basically no spatial smoothing
    spacing2.Fill(bigspace);
    spacing2[3] = 1;
    image1->SetSpacing(spacing2);
    using dgf = itk::DiscreteGaussianImageFilter<ImageType, ImageType>;
    typename dgf::Pointer filter = dgf::New();
    filter->SetVariance(time_smoother);
    filter->SetUseImageSpacing(true);
    filter->SetMaximumError(.01f);
    filter->SetInput(image1);
    filter->Update();
    image1 = filter->GetOutput();
    image1->SetSpacing(spacing);
  }

  if (maskfn.length() > 3)
  {
    ReadImage<OutImageType>(mask, maskfn.c_str());
  }
  else
  {
    // std::cerr << " cannot read mask " << maskfn << std::endl;
    return EXIT_FAILURE;
  }
  unsigned int  timedims = image1->GetLargestPossibleRegion().GetSize()[ImageDimension - 1];
  unsigned long voxct = 0;
  using ExtractFilterType = itk::ExtractImageFilter<ImageType, OutImageType>;
  using SliceIt = itk::ImageRegionIteratorWithIndex<OutImageType>;
  SliceIt mIter(mask, mask->GetLargestPossibleRegion());
  using LabelSetType = std::vector<unsigned int>;
  unsigned int maxlabel = 0;
  LabelSetType myLabelSet1;
  {
    /** count the labels in the mask image */
    for (mIter.GoToBegin(); !mIter.IsAtEnd(); ++mIter)
    {
      auto label = static_cast<unsigned int>(mIter.Get() + 0.5);
      if (label > 0)
      {
        if (find(myLabelSet1.begin(), myLabelSet1.end(), label) == myLabelSet1.end())
        {
          myLabelSet1.push_back(label);
          if (label > maxlabel)
          {
            maxlabel = label;
          }
        }
        voxct++;
      }
    }
  }
  if (maxlabel == 0)
  {
    // std::cerr << "FAILURE: Max label in input mask " << maskfn << " is 0 " << std::endl;
    return EXIT_FAILURE;
  }

  // std::cout << " timedims " << timedims << " n-Labels " << myLabelSet1.size() << std::endl;
  typename ImageType::RegionType extractRegion = image1->GetLargestPossibleRegion();
  extractRegion.SetSize(ImageDimension - 1, 0);
  unsigned int sub_vol = 0;
  extractRegion.SetIndex(ImageDimension - 1, sub_vol);
  typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
  extractFilter->SetInput(image1);
  //    extractFilter->SetDirectionCollapseToIdentity();
  extractFilter->SetDirectionCollapseToSubmatrix();
  extractFilter->SetExtractionRegion(extractRegion);
  extractFilter->Update();
  typename OutImageType::Pointer outimage = extractFilter->GetOutput();
  outimage->FillBuffer(0);

  using SliceIt = itk::ImageRegionIteratorWithIndex<OutImageType>;

  using timeVectorType = vnl_vector<Scalar>;
  timeVectorType mSample(timedims, 0);
  using MatrixType = itk::Array2D<double>;
  std::vector<std::string> ColumnHeaders;
  if (myLabelSet1.size() > 1)
  {
    /** sort the labels */
    std::sort(myLabelSet1.begin(), myLabelSet1.end());
    /** create a map between the roi and its matrix index */
    std::map<unsigned int, unsigned int> vectorindexMap;
    for (unsigned int i = 0; i < myLabelSet1.size(); i++)
    {
      std::string colname = std::string("Label") + sccan_to_string<unsigned int>(myLabelSet1[i]);
      ColumnHeaders.push_back(colname);
      vectorindexMap[myLabelSet1[i]] = i;
    }
    using countVectorType = vnl_vector<unsigned int>;
    countVectorType countVector(myLabelSet1.size(), 0);
    // std::cout << "Will map the image to its ROIs" << std::endl;
    MatrixType matrix(timedims, myLabelSet1.size());
    matrix.Fill(0);
    SliceIt vfIter2(outimage, outimage->GetLargestPossibleRegion());
    voxct = 0;
    for (vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2)
    {
      OutIndexType ind = vfIter2.GetIndex();
      auto         label = static_cast<unsigned int>(mask->GetPixel(ind) + 0.5);
      if (label > 0)
      {
        unsigned int vind = vectorindexMap[label];
        IndexType    tind;
        // first collect all samples for that location
        for (unsigned int i = 0; i < ImageDimension - 1; i++)
        {
          tind[i] = ind[i];
        }
        for (unsigned int t = 0; t < timedims; t++)
        {
          tind[ImageDimension - 1] = t;
          Scalar pix = image1->GetPixel(tind);
          mSample(t) = pix;
          matrix[t][vind] += pix;
          countVector[vind] += 1;
        }
      } // check mask
    }
    for (unsigned int i = 0; i < myLabelSet1.size(); i++)
    {
      matrix.set_column(i, matrix.get_column(i) / (double)countVector[i]);
    }
    using WriterType = itk::CSVNumericObjectFileWriter<double, 1, 1>;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outname);
    writer->SetInput(&matrix);
    writer->SetColumnHeaders(ColumnHeaders);
    try
    {
      writer->Write();
    }
    catch (const itk::ExceptionObject & itkNotUsed(exp))
    {
      // std::cerr << "Exception caught!" << std::endl;
      // std::cerr << exp << std::endl;
      return EXIT_FAILURE;
    }
    // std::cout << " done writing " << std::endl;
    return EXIT_SUCCESS;
  }
  MatrixType matrix(timedims, voxct);
  matrix.Fill(0);
  SliceIt vfIter2(outimage, outimage->GetLargestPossibleRegion());
  voxct = 0;
  for (vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2)
  {
    OutIndexType ind = vfIter2.GetIndex();
    if (mask->GetPixel(ind) >= 0.5)
    {
      IndexType tind;
      // first collect all samples for that location
      for (unsigned int i = 0; i < ImageDimension - 1; i++)
      {
        tind[i] = ind[i];
      }
      for (unsigned int t = 0; t < timedims; t++)
      {
        tind[ImageDimension - 1] = t;
        Scalar pix = image1->GetPixel(tind);
        mSample(t) = pix;
        matrix[t][voxct] = pix;
      }
      std::string colname = std::string("V") + sccan_to_string<unsigned int>(voxct);
      ColumnHeaders.push_back(colname);
      voxct++;
    } // check mask
  }

  using WriterType = itk::CSVNumericObjectFileWriter<double, 1, 1>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outname);
  writer->SetInput(&matrix);
  writer->SetColumnHeaders(ColumnHeaders);
  try
  {
    writer->Write();
  }
  catch (const itk::ExceptionObject & itkNotUsed(exp))
  {
    // std::cerr << "Exception caught!" << std::endl;
    // std::cerr << exp << std::endl;
    return EXIT_FAILURE;
  }
  // std::cout << " done writing " << std::endl;
  return EXIT_SUCCESS;
}

template <typename PixelType>
int
ConvertCSVVectorToImage(std::string csvfn, std::string maskfn, std::string outname, unsigned long rowOrCol)
{
  using Scalar = PixelType;
  using vMatrix = vnl_matrix<PixelType>;
  constexpr unsigned int ImageDimension = 3;
  using ImageType = itk::Image<PixelType, ImageDimension>;
  /** read the  images */
  typename ImageType::Pointer mask = nullptr;
  ReadImage<ImageType>(mask, maskfn.c_str());
  typename ImageType::Pointer outimage = nullptr;
  ReadImage<ImageType>(outimage, maskfn.c_str());
  outimage->FillBuffer(0);
  using Iterator = itk::ImageRegionIteratorWithIndex<ImageType>;
  unsigned long mct = 0;
  Iterator      mIter(mask, mask->GetLargestPossibleRegion());
  for (mIter.GoToBegin(); !mIter.IsAtEnd(); ++mIter)
  {
    if (mIter.Get() >= 0.5)
    {
      mct++;
    }
  }
  /** we refer to the two view matrices as P and Q */
  vMatrix p;
  p.fill(0);
  ReadMatrixFromCSVorImageSet<Scalar>(csvfn, p);
  if (mct != p.rows() && mct != p.cols())
  {
    // std::cout << " csv-vec rows " << p.rows() << " cols " << p.cols() << " mask non zero elements " << mct  <<
    // std::endl;
    throw std::exception();
  }

  if (mct == p.rows())
  {
    if (rowOrCol > p.cols() - 1)
    {
      // std::cout << " You are trying to select the " << rowOrCol << "th column but there are only " << p.cols()  << "
      // columns " << std::endl;
      throw std::exception();
    }
    mct = 0;
    for (mIter.GoToBegin(); !mIter.IsAtEnd(); ++mIter)
    {
      if (mIter.Get() >= 0.5)
      {
        PixelType val = p(mct, rowOrCol);
        outimage->SetPixel(mIter.GetIndex(), val);
        mct++;
      }
    }
  }
  else if (mct == p.cols()) // map the cols to the vector
  {
    if (rowOrCol > p.rows() - 1)
    {
      // std::cout << " You are trying to select the " << rowOrCol << "th row but there are only " << p.rows()  << "
      // rows " << std::endl;
      throw std::exception();
    }
    mct = 0;
    for (mIter.GoToBegin(); !mIter.IsAtEnd(); ++mIter)
    {
      if (mIter.Get() >= 0.5)
      {
        PixelType val = p(rowOrCol, mct);
        outimage->SetPixel(mIter.GetIndex(), val);
        mct++;
      }
    }
  }
  using WriterType = itk::ImageFileWriter<ImageType>;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outname);
  writer->SetInput(outimage);
  writer->Update();
  return EXIT_SUCCESS;
}

// p.d.
template <unsigned int ImageDimension, typename PixelType>
void
ConvertImageVecListToProjection(std::string veclist, std::string imagelist, std::string outname, bool average)
{
  // typedef itk::Image<PixelType,ImageDimension> ImageType;
  using ImageType = itk::Image<PixelType, ImageDimension>;
  using ReaderType = itk::ImageFileReader<ImageType>;

  std::vector<std::string> image_fn_list;
  std::vector<std::string> vec_fn_list;

  // first, count the number of files
  constexpr unsigned int maxChar = 512;
  char                   lineBuffer[maxChar], lineBufferVec[maxChar];
  char                   filenm[maxChar], filenmVec[maxChar];
  unsigned int           filecount = 0, filecountVec = 0;
  {
    std::ifstream inputStreamA(imagelist.c_str(), std::ios::in);
    if (!inputStreamA.is_open())
    {
      // std::cout << "Can't open image list file: " << imagelist << std::endl;
      return;
    }
    while (!inputStreamA.eof())
    {
      inputStreamA.getline(lineBuffer, maxChar, '\n');
      if (sscanf(lineBuffer, "%s ", filenm) != 1)
      {
        continue;
      }
      else
      {
        image_fn_list.emplace_back(filenm);
        filecount++;
      }
    }

    inputStreamA.close();
  }

  {
    std::ifstream inputStreamVec(veclist.c_str(), std::ios::in);
    if (!inputStreamVec.is_open())
    {
      // std::cout << "Can't open Vec list file: " << veclist << std::endl;
      return;
    }
    while (!inputStreamVec.eof())
    {
      inputStreamVec.getline(lineBufferVec, maxChar, '\n');
      if (sscanf(lineBufferVec, "%s ", filenmVec) != 1)
      {
        continue;
      }
      else
      {
        vec_fn_list.emplace_back(filenmVec);
        filecountVec++;
      }
    }

    inputStreamVec.close();
  }

  std::ofstream myfile;
  std::string   fnmp = outname + std::string(".csv");
  myfile.open(fnmp.c_str(), std::ios::out);
  using Iterator = itk::ImageRegionIteratorWithIndex<ImageType>;
  for (auto & j : image_fn_list)
  {
    for (unsigned int k = 0; k < vec_fn_list.size(); k++)
    {
      double                       proj = 0, dotSum = 0, dotCounter = 0, dotTotal = 0;
      typename ReaderType::Pointer reader1 = ReaderType::New();
      reader1->SetFileName(j);
      reader1->Update();
      typename ReaderType::Pointer reader2 = ReaderType::New();
      reader2->SetFileName(vec_fn_list[k]);
      reader2->Update();
      Iterator mIter(reader1->GetOutput(), reader1->GetOutput()->GetLargestPossibleRegion());
      Iterator mIter2(reader2->GetOutput(), reader2->GetOutput()->GetLargestPossibleRegion());
      for (mIter.GoToBegin(), mIter2.GoToBegin(); !mIter.IsAtEnd() && !mIter2.IsAtEnd(); ++mIter, ++mIter2)
      {
        proj = mIter.Get() * mIter2.Get();
        dotSum += proj;
        if (mIter2.Get() > 0)
        {
          dotCounter += mIter2.Get();
          dotTotal += mIter.Get() * mIter2.Get();
        }
      }
      if (average && dotCounter > 0)
      {
        dotSum = dotTotal / dotCounter;
      }
      if (k == vec_fn_list.size() - 1)
      {
        myfile << dotSum;
      }
      else
      {
        myfile << dotSum << " , ";
      }
    }
    myfile << std::endl;
  }
  myfile.close();
}

template <unsigned int ImageDimension, typename PixelType>
int
SVD_One_View(itk::ants::CommandLineParser * sccanparser,
             unsigned int                   permct,
             unsigned int                   n_evec,
             unsigned int                   robustify,
             unsigned int                   p_cluster_thresh,
             unsigned int                   iterct,
             unsigned int                   svd_option,
             PixelType                      usel1,
             PixelType                      row_sparseness,
             PixelType                      smoother,
             unsigned int                   covering,
             bool                           verbosity,
             bool                           getSmall = false)
{
  // std::cout << "SVD_One_View" << std::endl;

  using ImageType = itk::Image<PixelType, ImageDimension>;
  using Scalar = double;
  using SCCANType = itk::ants::antsSCCANObject<ImageType, Scalar>;
  using vMatrix = typename SCCANType::MatrixType;
  using vVector = typename SCCANType::VectorType;
  typename SCCANType::Pointer sccanobj = SCCANType::New();
  sccanobj->SetGetSmall(static_cast<bool>(getSmall));
  vMatrix priorROIMat;
  if (svd_option == 1)
  {
    // std::cout << " basic-svd " << std::endl;
  }
  else
  {
    // std::cout << " sparse-svd " << std::endl;   // note: 2 (in options) is for svd implementation
  }

  itk::ants::CommandLineParser::OptionType::Pointer initOpt = sccanparser->GetOption("initialization");
  itk::ants::CommandLineParser::OptionType::Pointer maskOpt = sccanparser->GetOption("mask");
  if (!initOpt || initOpt->GetNumberOfFunctions() == 0 || !maskOpt || maskOpt->GetNumberOfFunctions() == 0)
  {
    // std::cout << "Warning:  no initialization set, will use data-driven approach." << std::endl;
  }
  else
  {
    std::string maskfn = maskOpt->GetFunction(0)->GetName();
    std::string imagelistPrior = initOpt->GetFunction(0)->GetName();
    // std::cout << "you will initialize with " << imagelistPrior << std::endl;
    std::string outname = "none";
    priorROIMat = ConvertImageListToMatrix<ImageDimension, double>(imagelistPrior, maskfn, outname);
    // std::cout << priorROIMat.rows() << " " << priorROIMat.cols() << std::endl;
    sccanobj->SetMatrixPriorROI(priorROIMat);
  }

  itk::ants::CommandLineParser::OptionType::Pointer init2Opt = sccanparser->GetOption("initialization2");
  itk::ants::CommandLineParser::OptionType::Pointer mask2Opt = sccanparser->GetOption("mask2");
  if (!init2Opt || init2Opt->GetNumberOfFunctions() == 0 || !mask2Opt || mask2Opt->GetNumberOfFunctions() == 0)
  {
    // std::cout << "Warning:  no initialization set, will use data-driven approach." << std::endl;
  }
  else
  {
    std::string maskfn = mask2Opt->GetFunction(0)->GetName();
    std::string imagelistPrior = init2Opt->GetFunction(0)->GetName();
    // std::cout << "you will initialize Q with " << imagelistPrior << std::endl;
    std::string outname = "none";
    priorROIMat = ConvertImageListToMatrix<ImageDimension, double>(imagelistPrior, maskfn, outname);
    // std::cout << priorROIMat.rows() << " " << priorROIMat.cols() << std::endl;
    sccanobj->SetMatrixPriorROI2(priorROIMat);
  }

  itk::ants::CommandLineParser::OptionType::Pointer outputOption = sccanparser->GetOption("output");
  if (!outputOption || outputOption->GetNumberOfFunctions() == 0)
  {
    // std::cout << "Warning:  no output option set." << std::endl;
  }
  itk::ants::CommandLineParser::OptionType::Pointer option = sccanparser->GetOption("svd");
  PixelType                                         gradstep = itk::Math::abs(usel1);
  sccanobj->SetCovering(covering);
  sccanobj->SetSilent(!verbosity);
  if (usel1 > 0)
  {
    sccanobj->SetUseL1(true);
  }
  else
  {
    sccanobj->SetUseL1(false);
  }
  sccanobj->SetGradStep(gradstep);
  sccanobj->SetMaximumNumberOfIterations(iterct);
  sccanobj->SetRowSparseness(row_sparseness);
  sccanobj->SetSmoother(smoother);
  /** read the matrix images */
  /** we refer to the two view matrices as P and Q */

  std::string pmatname = std::string(option->GetFunction(0)->GetParameter(0));
  vMatrix     p;
  ReadMatrixFromCSVorImageSet<Scalar>(pmatname, p);
  typename ImageType::Pointer mask1 = nullptr;
  bool                        have_p_mask = false;
  have_p_mask = ReadImage<ImageType>(mask1, option->GetFunction(0)->GetParameter(1).c_str());
  auto    FracNonZero1 = sccanparser->Convert<double>(option->GetFunction(0)->GetParameter(2));
  vMatrix priorScaleMat;
  if (svd_option == 7)
  {
    FracNonZero1 = sccanparser->Convert<double>(option->GetFunction(0)->GetParameter(4));
    std::string imagelistPrior = option->GetFunction(0)->GetParameter(2);
    // std::string priorScaleFile = option->GetFunction( 0 )->GetParameter( 3 );
    std::string outname = "prior.mhd";
    ConvertImageListToMatrix<ImageDimension, double>(imagelistPrior, option->GetFunction(0)->GetParameter(1), outname);
    ReadMatrixFromCSVorImageSet<Scalar>(outname, priorROIMat);
    // ReadMatrixFromCSVorImageSet<Scalar>(priorScaleFile, priorScaleMat);
  }
  // std::cout << " frac nonzero " << FracNonZero1 << std::endl;

  if (robustify > 0)
  {
    // std::cout << " make robust " << std::endl;
    p = sccanobj->RankifyMatrixColumns(p);
  }

  /** the penalties define the fraction of non-zero values for each view */
  if (FracNonZero1 < 0)
  {
    FracNonZero1 = fabs(FracNonZero1);
    sccanobj->SetKeepPositiveP(false);
  }

  /** read the nuisance matrix image */
  vMatrix r;
  if (option->GetFunction(0)->GetNumberOfParameters() > 3)
  {
    std::string nuis_img = option->GetFunction(0)->GetParameter(3);
    if (nuis_img.length() > 3 && svd_option != 7)
    {
      // std::cout << " nuis_img " << nuis_img << std::endl;
      ReadMatrixFromCSVorImageSet<Scalar>(nuis_img, r);
      if (CompareMatrixSizes<Scalar>(p, r) == EXIT_FAILURE)
      {
        return EXIT_FAILURE;
      }
      itk::ants::CommandLineParser::OptionType::Pointer partialccaOpt = sccanparser->GetOption("partial-scca-option");
      std::string                                       partialccaoption = std::string("PQ");
      if (partialccaOpt)
      {
        //  enum SCCANFormulationType{ PQ , PminusRQ ,  PQminusR ,  PminusRQminusR , PQR  };
        if (partialccaOpt->GetNumberOfFunctions() > 0)
        {
          partialccaoption = sccanparser->Convert<std::string>(partialccaOpt->GetFunction()->GetName());
        }
        // std::cout << " Partial SCCA option " << partialccaoption << std::endl;
        if (!partialccaoption.compare(std::string("PQ")))
        {
          sccanobj->SetSCCANFormulation(SCCANType::PQ);
        }
        else if (!partialccaoption.compare(std::string("PminusRQ")))
        {
          sccanobj->SetSCCANFormulation(SCCANType::PminusRQ);
        }
        else if (!partialccaoption.compare(std::string("PQminusR")))
        {
          sccanobj->SetSCCANFormulation(SCCANType::PQminusR);
        }
        else if (!partialccaoption.compare(std::string("PminusRQminusR")))
        {
          sccanobj->SetSCCANFormulation(SCCANType::PminusRQminusR);
        }
      }
    }
  }
  else
  {
    // std::cout << " No nuisance parameters." << std::endl;
  }

  sccanobj->SetFractionNonZeroP(FracNonZero1);
  sccanobj->SetMinClusterSizeP(p_cluster_thresh);
  if (robustify > 0)
  {
    // std::cout << " make robust " << std::endl;
    p = sccanobj->RankifyMatrixColumns(p);
  }
  sccanobj->SetMatrixP(p);
  sccanobj->SetMatrixR(r);
  sccanobj->SetMaskImageP(mask1);

  double truecorr = 0;
  if (svd_option == 1)
  {
    truecorr = sccanobj->SparseReconHome(n_evec);
  }
  else if (svd_option == 3)
  {
    truecorr = sccanobj->SparseArnoldiSVD(n_evec); // cgsparse
  }
  else if (svd_option == 4)
  {
    truecorr = sccanobj->NetworkDecomposition(n_evec);
  }
  else if (svd_option == 5)
  {
    truecorr = sccanobj->LASSO(n_evec);
  }
  else if (svd_option == 2)
  {
    truecorr = sccanobj->CGSPCA(n_evec); // cgspca
  }
  else if (svd_option == 6)
  {
    truecorr = sccanobj->SparseRecon(n_evec); // sparse (default)
  }
  else if (svd_option == 7)
  {
    // sccanobj->SetPriorScaleMat( priorScaleMat);
    sccanobj->SetMatrixPriorROI(priorROIMat);
    sccanobj->SetFlagForSort();
    sccanobj->SetLambda(sccanparser->Convert<double>(option->GetFunction(0)->GetParameter(3)));
    truecorr = sccanobj->SparseReconPrior(n_evec, true); // Prior
  }
  else
  {
    truecorr = sccanobj->SparseArnoldiSVDGreedy(n_evec); // sparse (default)
  }
  vVector w_p = sccanobj->GetVariateP(0);
  // std::cout << " true-corr " << sccanobj->GetCanonicalCorrelations()  << std::endl;

  if (outputOption)
  {
    std::string filename = outputOption->GetFunction(0)->GetName();
    // std::cout << " write " << filename << std::endl;
    std::string::size_type pos = filename.rfind(".");
    std::string            filepre = std::string(filename, 0, pos);
    std::string            extension = std::string(filename, pos, filename.length() - 1);
    if (extension == std::string(".gz"))
    {
      pos = filepre.rfind(".");
      extension = std::string(filepre, pos, filepre.length() - 1) + extension;
      filepre = std::string(filepre, 0, pos);
    }
    std::string post = std::string("View1vec");
    WriteVariatesToSpatialImage<ImageType, Scalar>(
      filename, post, sccanobj->GetVariatesP(), mask1, sccanobj->GetMatrixP(), have_p_mask, sccanobj->GetMatrixU());

    // WriteSortedVariatesToSpatialImage<ImageType,Scalar>( filename, post, sccanobj->GetVariatesP() , mask1 ,
    // sccanobj->GetMatrixP() , have_p_mask,sccanobj->GetSortFinalLocArray(),priorROIMat );

    /** write the eigevalues to the csv file */
    std::string              fnmp = filepre + std::string("_eigenvalues.csv");
    std::vector<std::string> ColumnHeaders;
    std::string              colname = std::string("Eigenvalue");
    ColumnHeaders.push_back(colname);
    using CWriterType = itk::CSVNumericObjectFileWriter<double, 1, 1>;
    CWriterType::Pointer cwriter = CWriterType::New();
    cwriter->SetFileName(fnmp.c_str());
    cwriter->SetColumnHeaders(ColumnHeaders);
    vnl_matrix<double> evals;
    evals.set_size(sccanobj->GetCanonicalCorrelations().size(), 1);
    for (unsigned int i = 0; i < sccanobj->GetCanonicalCorrelations().size(); i++)
    {
      double evil = sccanobj->GetCanonicalCorrelations()(i);
      evals(i, 0) = evil;
    }
    cwriter->SetInput(&evals);
    cwriter->Write();
  }
  // permutation test
  if ((svd_option == 4 || svd_option == 5) && permct > 0)
  {
    // std::cout << "Begin" << permct << " permutations " << std::endl;
    unsigned long perm_exceed_ct = 0;
    for (unsigned long pct = 0; pct <= permct; pct++)
    {
      // 0. compute permutation for q ( switch around rows )
      vMatrix p_perm = PermuteMatrix<Scalar>(sccanobj->GetOriginalMatrixP());
      vMatrix r_perm = PermuteMatrix<Scalar>(sccanobj->GetOriginalMatrixR());
      sccanobj->SetMatrixP(p_perm);
      sccanobj->SetMatrixR(r_perm);
      double permcorr = 1.e9;
      // if ( pct > 76 && pct < 79 )
      if (svd_option == 4)
      {
        permcorr = sccanobj->NetworkDecomposition(n_evec); // cgsparse
      }
      if (svd_option == 5)
      {
        permcorr = sccanobj->LASSO(n_evec); // cgsparse
      }
      if (permcorr < truecorr)
      {
        perm_exceed_ct++;
      }
      // end solve cca permutation
      // std::cout << permcorr << " p-value " <<  (double)perm_exceed_ct
      //  / (pct + 1) << " ct " << pct << " true " << truecorr << " vs " << permcorr << std::endl;
    }
  }
  return EXIT_SUCCESS;
}

template <unsigned int ImageDimension, typename PixelType>
int
SCCA_vnl(itk::ants::CommandLineParser * sccanparser,
         unsigned int                   permct,
         unsigned int                   n_evec,
         unsigned int                   newimp,
         unsigned int                   robustify,
         unsigned int                   p_cluster_thresh,
         unsigned int                   q_cluster_thresh,
         unsigned int                   iterct,
         PixelType                      usel1,
         PixelType                      uselong,
         PixelType                      row_sparseness,
         PixelType                      smoother,
         unsigned int                   covering,
         PixelType                      priorWeight = 0,
         unsigned int                   verbosity = 0)
{
  itk::ants::CommandLineParser::OptionType::Pointer outputOption = sccanparser->GetOption("output");
  bool                                              writeoutput = true;

  if (!outputOption || outputOption->GetNumberOfFunctions() == 0)
  {
    // std::cout << "Warning:  no output option set." << std::endl;
    writeoutput = false;
  }
  if (newimp > 0)
  {
    // std::cout << "New imp irrelevant " << std::endl;
  }
  itk::ants::CommandLineParser::OptionType::Pointer option = sccanparser->GetOption("scca");
  using ImageType = itk::Image<PixelType, ImageDimension>;
  using Scalar = double;
  using SCCANType = itk::ants::antsSCCANObject<ImageType, Scalar>;
  using vMatrix = typename SCCANType::MatrixType;
  using vVector = typename SCCANType::VectorType;
  typename SCCANType::Pointer sccanobj = SCCANType::New();
  sccanobj->SetCovering(covering);
  sccanobj->SetSilent(!verbosity);
  sccanobj->SetPriorWeight(priorWeight);
  sccanobj->SetMaximumNumberOfIterations(iterct);
  if (uselong > 0)
  {
    sccanobj->SetUseLongitudinalFormulation(uselong);
  }
  PixelType gradstep = itk::Math::abs(usel1);
  if (usel1 > 0)
  {
    sccanobj->SetUseL1(true);
  }
  else
  {
    sccanobj->SetUseL1(false);
  }
  vMatrix                                           priorROIMat;
  vMatrix                                           priorROIMat2;
  itk::ants::CommandLineParser::OptionType::Pointer initOpt = sccanparser->GetOption("initialization");
  itk::ants::CommandLineParser::OptionType::Pointer maskOpt = sccanparser->GetOption("mask");
  if (!initOpt || initOpt->GetNumberOfFunctions() == 0 || !maskOpt || maskOpt->GetNumberOfFunctions() == 0)
  {
  }
  else
  {
    std::string maskfn = maskOpt->GetFunction(0)->GetName();
    std::string imagelistPrior = initOpt->GetFunction(0)->GetName();
    // std::cout << "you will initialize P with " << imagelistPrior << " and " << maskfn << std::endl;
    std::string outname = "none";
    priorROIMat = ConvertImageListToMatrix<ImageDimension, double>(imagelistPrior, maskfn, outname);
    // std::cout << priorROIMat.rows() << " " << priorROIMat.cols() << std::endl;
    sccanobj->SetMatrixPriorROI(priorROIMat);
  }

  itk::ants::CommandLineParser::OptionType::Pointer init2Opt = sccanparser->GetOption("initialization2");
  itk::ants::CommandLineParser::OptionType::Pointer mask2Opt = sccanparser->GetOption("mask2");
  if (!init2Opt || init2Opt->GetNumberOfFunctions() == 0 || !mask2Opt || mask2Opt->GetNumberOfFunctions() == 0)
  {
    // itkDebugStatement( std::cerr << "Warning:  no Q initialization set, will use data-driven approach." << std::endl
    // );
  }
  else
  {
    std::string maskfn = mask2Opt->GetFunction(0)->GetName();
    std::string imagelistPrior = init2Opt->GetFunction(0)->GetName();
    // std::cout << "you will initialize Q with " << imagelistPrior << " and " << maskfn << std::endl;
    std::string outname = "none";
    priorROIMat2 = ConvertImageListToMatrix<ImageDimension, double>(imagelistPrior, maskfn, outname);
    // std::cout << priorROIMat2.rows() << " " << priorROIMat2.cols() << std::endl;
    sccanobj->SetMatrixPriorROI2(priorROIMat2);
  }
  sccanobj->SetGradStep(gradstep);
  sccanobj->SetSmoother(smoother);
  sccanobj->SetRowSparseness(row_sparseness);
  Scalar pinvtoler = 1.e-6;
  /** read the matrix images */
  /** we refer to the two view matrices as P and Q */
  std::string pmatname = std::string(option->GetFunction(0)->GetParameter(0));
  vMatrix     p;
  // // std::cout <<" read-p "<< std::endl;
  ReadMatrixFromCSVorImageSet<Scalar>(pmatname, p);
  std::string qmatname = std::string(option->GetFunction(0)->GetParameter(1));
  vMatrix     q;
  // // std::cout <<" read-q "<< std::endl;
  ReadMatrixFromCSVorImageSet<Scalar>(qmatname, q);
  //  // std::cout << q.get_row(0) << std::endl;
  //  // std::cout << q.mean() << std::endl;
  if (CompareMatrixSizes<Scalar>(p, q) == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  }

  typename ImageType::Pointer mask1 = nullptr;
  std::string                 mask1fn = option->GetFunction(0)->GetParameter(2);
  bool                        have_p_mask = false;
  if (mask1fn.length() > 5)
  {
    have_p_mask = ReadImage<ImageType>(mask1, mask1fn.c_str());
  }

  typename ImageType::Pointer mask2 = nullptr;
  std::string                 mask2fn = option->GetFunction(0)->GetParameter(3);
  bool                        have_q_mask = false;
  if (mask2fn.length() > 5)
  {
    have_q_mask = ReadImage<ImageType>(mask2, mask2fn.c_str());
  }

  /** the penalties define the fraction of non-zero values for each view */
  auto FracNonZero1 = sccanparser->Convert<double>(option->GetFunction(0)->GetParameter(4));
  if (FracNonZero1 < 0)
  {
    FracNonZero1 = fabs(FracNonZero1);
    sccanobj->SetKeepPositiveP(false); // true if P sparsity > 0
  }
  auto FracNonZero2 = sccanparser->Convert<double>(option->GetFunction(0)->GetParameter(5));
  if (FracNonZero2 < 0)
  {
    FracNonZero2 = fabs(FracNonZero2);
    sccanobj->SetKeepPositiveQ(false); // true if Q sparsity > 0
  }

  sccanobj->SetFractionNonZeroP(FracNonZero1);
  sccanobj->SetFractionNonZeroQ(FracNonZero2);
  sccanobj->SetMinClusterSizeP(p_cluster_thresh);
  sccanobj->SetMinClusterSizeQ(q_cluster_thresh);
  if (robustify > 0)
  {
    // std::cout << " make robust " << std::endl;
    p = sccanobj->RankifyMatrixColumns(p);
    q = sccanobj->RankifyMatrixColumns(q);
  }
  sccanobj->SetMatrixP(p);
  sccanobj->SetMatrixQ(q);
  sccanobj->SetMaskImageP(mask1);
  sccanobj->SetMaskImageQ(mask2);
  sccanobj->SparsePartialArnoldiCCA(n_evec);
  vVector w_p = sccanobj->GetVariateP(0);
  vVector w_q = sccanobj->GetVariateQ(0);
  vVector sccancorrs = sccanobj->GetCanonicalCorrelations();
  // std::cout << " true-corr " << sccancorrs << std::endl;

  std::string            filename = outputOption->GetFunction(0)->GetName();
  std::string::size_type pos = filename.rfind(".");
  std::string            filepre = std::string(filename, 0, pos);
  std::string            extension = std::string(filename, pos, filename.length() - 1);
  if (extension == std::string(".gz"))
  {
    pos = filepre.rfind(".");
    extension = std::string(filepre, pos, filepre.length() - 1) + extension;
    filepre = std::string(filepre, 0, pos);
  }
  std::string post = std::string("View1vec");
  if (writeoutput)
  {
    WriteVariatesToSpatialImage<ImageType, Scalar>(
      filename, post, sccanobj->GetVariatesP(), mask1, sccanobj->GetMatrixP(), have_p_mask, sccanobj->GetMatrixU());
    post = std::string("View2vec");
    WriteVariatesToSpatialImage<ImageType, Scalar>(
      filename, post, sccanobj->GetVariatesQ(), mask2, sccanobj->GetMatrixQ(), have_q_mask, sccanobj->GetMatrixU());
  }

  /** begin permutation 1. q_pvMatrix CqqInv=vnl_svd_inverse<Scalar>(Cqq);
   q=q*CqqInv;
  sermuted ;  2. scca ;  3. test corrs and weights significance */
  if (permct > 0)
  {
    vnl_vector<unsigned long> perm_exceed_ct(sccancorrs.size(), 0);
    vVector                   w_p_signif_ct(w_p.size(), 0);
    vVector                   w_q_signif_ct(w_q.size(), 0);
    for (unsigned long pct = 0; pct <= permct; pct++)
    {
      // 0. compute permutation for q ( switch around rows )
      sccanobj->SetFractionNonZeroP(FracNonZero1);
      sccanobj->SetFractionNonZeroQ(FracNonZero2);
      sccanobj->SetGradStep(gradstep);
      sccanobj->SetMatrixP(p);
      ReadMatrixFromCSVorImageSet<Scalar>(qmatname, q);
      vMatrix q_perm = PermuteMatrix<Scalar>(q);
      sccanobj->SetMatrixQ(q_perm);
      sccanobj->SparsePartialArnoldiCCA(n_evec);
      vVector permcorrs = sccanobj->GetCanonicalCorrelations();
      // std::cout << " perm-corr " << permcorrs << " ct " << pct << " p-values ";
      for (unsigned int kk = 0; kk < permcorrs.size(); kk++)
      {
        if (permcorrs[kk] > sccancorrs[kk])
        {
          perm_exceed_ct[kk]++;
        }
        // std::cout <<  ( double ) perm_exceed_ct[kk] / (pct + 1) << " ";
      }
      // std::cout << std::endl;
      vVector w_p_perm = sccanobj->GetVariateP(0);
      vVector w_q_perm = sccanobj->GetVariateQ(0);
      for (unsigned long j = 0; j < w_p.size(); j++)
      {
        if (w_p_perm(j) > w_p(j))
        {
          w_p_signif_ct(j) = w_p_signif_ct(j)++;
        }
      }
      for (unsigned long j = 0; j < w_q.size(); j++)
      {
        if (w_q_perm(j) > w_q(j))
        {
          w_q_signif_ct(j) = w_q_signif_ct(j)++;
        }
      }
      // end solve cca permutation
      if (pct == permct)
      {
        // std::cout << "final_p_values" << ",";
        for (unsigned int kk = 0; kk < permcorrs.size(); kk++)
        {
          // std::cout << ( double ) perm_exceed_ct[kk] / (pct + 1) << ",";
        }
        // std::cout << "x" << std::endl;

        std::ofstream myfile;
        std::string   fnmp = filepre + std::string("_summary.csv");
        myfile.open(fnmp.c_str(), std::ios::out);
        myfile << "TypeOfMeasure"
               << ",";
        for (unsigned int kk = 0; kk < permcorrs.size(); kk++)
        {
          std::string colname = std::string("Variate") + sccan_to_string<unsigned int>(kk);
          myfile << colname << ",";
        }
        myfile << "x" << std::endl;
        myfile << "final_p_values"
               << ",";
        for (unsigned int kk = 0; kk < permcorrs.size(); kk++)
        {
          myfile << (double)perm_exceed_ct[kk] / (pct + 1) << ",";
        }
        myfile << "x" << std::endl;
        myfile << "corrs"
               << ",";
        for (unsigned int kk = 0; kk < permcorrs.size(); kk++)
        {
          myfile << sccancorrs[kk] << ",";
        }
        myfile << "x" << std::endl;
        myfile.close();
      }
    }
    unsigned long psigct = 0, qsigct = 0;
    for (unsigned long j = 0; j < w_p.size(); j++)
    {
      if (w_p(j) > pinvtoler)
      {
        w_p_signif_ct(j) = 1.0 - (double)w_p_signif_ct(j) / (double)(permct);
        if (w_p_signif_ct(j) > 0.949)
        {
          psigct++;
        }
      }
      else
      {
        w_p_signif_ct(j) = 0;
      }
    }
    for (unsigned long j = 0; j < w_q.size(); j++)
    {
      if (w_q(j) > pinvtoler)
      {
        w_q_signif_ct(j) = 1.0 - (double)w_q_signif_ct(j) / (double)(permct);
        if (w_q_signif_ct(j) > 0.949)
        {
          qsigct++;
        }
      }
      else
      {
        w_q_signif_ct(j) = 0;
      }
    }
    if (writeoutput)
    {
      post = std::string("View1pval");
      if (have_p_mask)
      {
        WriteVectorToSpatialImage<ImageType, Scalar>(filename, post, w_p_signif_ct, mask1);
      }
      post = std::string("View2pval");
      if (have_q_mask)
      {
        WriteVectorToSpatialImage<ImageType, Scalar>(filename, post, w_q_signif_ct, mask2);
      }
    }
  }
  return EXIT_SUCCESS;
}

template <unsigned int ImageDimension, typename PixelType>
int
mSCCA_vnl(itk::ants::CommandLineParser * sccanparser,
          unsigned int                   permct,
          bool                           run_partial_scca = false,
          unsigned int                   n_e_vecs = 3,
          unsigned int                   newimp = 0,
          unsigned int                   robustify = 0,
          unsigned int                   p_cluster_thresh = 100,
          unsigned int                   q_cluster_thresh = 1,
          unsigned int                   iterct = 20)
{
  // std::cout << " Entering MSCCA --- computing " << n_e_vecs << " canonical variates by default. " << std::endl;
  itk::ants::CommandLineParser::OptionType::Pointer outputOption = sccanparser->GetOption("output");

  if (!outputOption || outputOption->GetNumberOfFunctions() == 0)
  {
    // std::cout << "Warning:  no output option set." << std::endl;
  }
  // std::cout << " newimp " << newimp << std::endl;
  itk::ants::CommandLineParser::OptionType::Pointer option = sccanparser->GetOption("scca");
  using ImageType = itk::Image<PixelType, ImageDimension>;
  using Scalar = double;
  using SCCANType = itk::ants::antsSCCANObject<ImageType, Scalar>;
  using MatrixImageType = itk::Image<Scalar, 2>;
  typename SCCANType::Pointer sccanobj = SCCANType::New();
  sccanobj->SetMaximumNumberOfIterations(iterct);
  using vMatrix = typename SCCANType::MatrixType;
  using vVector = typename SCCANType::VectorType;

  /** we refer to the two view matrices as P and Q */
  using ImageType = itk::Image<PixelType, ImageDimension>;
  using Scalar = double;
  using MatrixImageType = itk::Image<Scalar, 2>;

  /** read the matrix images */
  std::string pmatname = std::string(option->GetFunction(0)->GetParameter(0));
  vMatrix     pin;
  ReadMatrixFromCSVorImageSet<Scalar>(pmatname, pin);
  std::string qmatname = std::string(option->GetFunction(0)->GetParameter(1));
  vMatrix     qin;
  ReadMatrixFromCSVorImageSet<Scalar>(qmatname, qin);
  std::string rmatname = std::string(option->GetFunction(0)->GetParameter(2));
  vMatrix     rin;
  ReadMatrixFromCSVorImageSet<Scalar>(rmatname, rin);
  if (CompareMatrixSizes<Scalar>(pin, qin) == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  }
  if (CompareMatrixSizes<Scalar>(qin, rin) == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  }
  if (CompareMatrixSizes<Scalar>(pin, rin) == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  }

  typename ImageType::Pointer mask1 = nullptr;
  bool have_p_mask = ReadImage<ImageType>(mask1, option->GetFunction(0)->GetParameter(3).c_str());
  typename ImageType::Pointer mask2 = nullptr;
  bool have_q_mask = ReadImage<ImageType>(mask2, option->GetFunction(0)->GetParameter(4).c_str());
  typename ImageType::Pointer mask3 = nullptr;

  /** the penalties define the fraction of non-zero values for each view */
  auto FracNonZero1 = sccanparser->Convert<double>(option->GetFunction(0)->GetParameter(6));
  if (FracNonZero1 < 0)
  {
    FracNonZero1 = fabs(FracNonZero1);
    sccanobj->SetKeepPositiveP(false);
  }
  auto FracNonZero2 = sccanparser->Convert<double>(option->GetFunction(0)->GetParameter(7));
  if (FracNonZero2 < 0)
  {
    FracNonZero2 = fabs(FracNonZero2);
    sccanobj->SetKeepPositiveQ(false);
  }
  auto FracNonZero3 = sccanparser->Convert<double>(option->GetFunction(0)->GetParameter(8));
  if (FracNonZero3 < 0)
  {
    FracNonZero3 = fabs(FracNonZero3);
    sccanobj->SetKeepPositiveR(false);
  }

  sccanobj->SetFractionNonZeroP(FracNonZero1);
  sccanobj->SetFractionNonZeroQ(FracNonZero2);
  sccanobj->SetFractionNonZeroR(FracNonZero3);
  for (unsigned int leave_out = pin.rows(); leave_out <= pin.rows(); leave_out++)
  {
    // std::cout << " Leaving Out " << leave_out << std::endl;
    vVector p_leave_out;
    vVector q_leave_out;
    if (leave_out < pin.rows())
    {
      p_leave_out = pin.get_row(leave_out);
      q_leave_out = qin.get_row(leave_out);
    }
    vMatrix p = DeleteRow<MatrixImageType, Scalar>(pin, leave_out);
    vMatrix q = DeleteRow<MatrixImageType, Scalar>(qin, leave_out);
    vMatrix r = DeleteRow<MatrixImageType, Scalar>(rin, leave_out);
    sccanobj->SetMinClusterSizeP(p_cluster_thresh);
    sccanobj->SetMinClusterSizeQ(q_cluster_thresh);
    if (robustify > 0)
    {
      // std::cout << " make robust " << std::endl;
      p = sccanobj->RankifyMatrixColumns(p);
      q = sccanobj->RankifyMatrixColumns(q);
      r = sccanobj->RankifyMatrixColumns(r);
    }
    double truecorr = 0;
    if (run_partial_scca)
    {
      // std::cout << " begin partial PQ " << std::endl;
      typename SCCANType::Pointer sccanobjCovar = SCCANType::New();
      sccanobjCovar->SetMaximumNumberOfIterations(iterct);
      sccanobjCovar->SetMatrixP(p);
      sccanobjCovar->SetMatrixQ(q);
      sccanobjCovar->SetMatrixR(r);
      sccanobjCovar->SetMinClusterSizeP(p_cluster_thresh);
      sccanobjCovar->SetMinClusterSizeQ(q_cluster_thresh);

      itk::ants::CommandLineParser::OptionType::Pointer partialccaOpt = sccanparser->GetOption("partial-scca-option");
      std::string                                       partialccaoption = std::string("PQ");
      if (partialccaOpt)
      {
        //  enum SCCANFormulationType{ PQ , PminusRQ ,  PQminusR ,  PminusRQminusR , PQR  };
        if (partialccaOpt->GetNumberOfFunctions() > 0)
        {
          partialccaoption = sccanparser->Convert<std::string>(partialccaOpt->GetFunction()->GetName());
        }
        // std::cout << " Partial SCCA option " << partialccaoption << std::endl;
        if (!partialccaoption.compare(std::string("PQ")))
        {
          sccanobjCovar->SetSCCANFormulation(SCCANType::PQ);
        }
        else if (!partialccaoption.compare(std::string("PminusRQ")))
        {
          sccanobjCovar->SetSCCANFormulation(SCCANType::PminusRQ);
        }
        else if (!partialccaoption.compare(std::string("PQminusR")))
        {
          sccanobjCovar->SetSCCANFormulation(SCCANType::PQminusR);
        }
        else if (!partialccaoption.compare(std::string("PminusRQminusR")))
        {
          sccanobjCovar->SetSCCANFormulation(SCCANType::PminusRQminusR);
        }
      }
      sccanobjCovar->SetFractionNonZeroP(FracNonZero1);
      sccanobjCovar->SetKeepPositiveP(sccanobj->GetKeepPositiveP());
      sccanobjCovar->SetFractionNonZeroQ(FracNonZero2);
      sccanobjCovar->SetKeepPositiveQ(sccanobj->GetKeepPositiveQ());
      sccanobjCovar->SetMaskImageP(mask1);
      sccanobjCovar->SetMaskImageQ(mask2);
      if (newimp == 1)
      {
        truecorr = sccanobjCovar->SparsePartialCCA(n_e_vecs);
      }
      else if (newimp == 0)
      {
        truecorr = sccanobjCovar->SparsePartialArnoldiCCA(n_e_vecs);
      }
      //  truecorr=sccanobjCovar->RunSCCAN2multiple(n_e_vecs );
      // std::cout << " partialed out corr ";
      for (unsigned int ff = 0; ff < sccanobjCovar->GetCanonicalCorrelations().size(); ff++)
      {
        // std::cout << " " << sccanobjCovar->GetCanonicalCorrelations()[ff];
      }
      // std::cout << std::endl;

      if (outputOption)
      {
        std::string            filename = outputOption->GetFunction(0)->GetName();
        std::string::size_type pos = filename.rfind(".");
        std::string            filepre = std::string(filename, 0, pos);
        std::string            extension = std::string(filename, pos, filename.length() - 1);
        if (extension == std::string(".gz"))
        {
          pos = filepre.rfind(".");
          extension = std::string(filepre, pos, filepre.length() - 1) + extension;
          filepre = std::string(filepre, 0, pos);
        }
        std::string post = std::string("View1vec");
        WriteVariatesToSpatialImage<ImageType, Scalar>(filename,
                                                       post,
                                                       sccanobjCovar->GetVariatesP(),
                                                       mask1,
                                                       sccanobjCovar->GetMatrixP(),
                                                       have_p_mask,
                                                       sccanobj->GetMatrixU());
        post = std::string("View2vec");
        WriteVariatesToSpatialImage<ImageType, Scalar>(filename,
                                                       post,
                                                       sccanobjCovar->GetVariatesQ(),
                                                       mask2,
                                                       sccanobjCovar->GetMatrixQ(),
                                                       have_q_mask,
                                                       sccanobj->GetMatrixU());
      }

      /** begin permutation 1. q_pvMatrix CqqInv=vnl_svd_inverse<Scalar>(Cqq);
       q=q*CqqInv;
      sermuted ;  2. scca ;  3. test corrs and weights significance */
      if (permct > 0)
      {
        unsigned long perm_exceed_ct = 0;
        vVector       w_p_signif_ct(p.cols(), 0);
        vVector       w_q_signif_ct(q.cols(), 0);
        for (unsigned long pct = 0; pct <= permct; pct++)
        {
          /** both the new object and the copy object should produce the same results - verified in 1 example!*/
          //      typename SCCANType::Pointer sccanobjPerm=sccanobjCovar;
          typename SCCANType::Pointer sccanobjPerm = SCCANType::New();
          sccanobjPerm->SetMaximumNumberOfIterations(iterct);
          sccanobjPerm->SetMinClusterSizeP(p_cluster_thresh);
          sccanobjPerm->SetMinClusterSizeQ(q_cluster_thresh);
          sccanobjPerm->SetFractionNonZeroP(FracNonZero1);
          sccanobjPerm->SetKeepPositiveP(sccanobj->GetKeepPositiveP());
          sccanobjPerm->SetKeepPositiveQ(sccanobj->GetKeepPositiveQ());
          sccanobjPerm->SetFractionNonZeroQ(FracNonZero2);
          sccanobjPerm->SetMaskImageP(mask1);
          sccanobjPerm->SetMaskImageQ(mask2);
          // 0. compute permutation for q ( switch around rows )
          vMatrix p_perm = PermuteMatrix<Scalar>(sccanobjCovar->GetMatrixP());
          vMatrix q_perm = PermuteMatrix<Scalar>(sccanobjCovar->GetMatrixQ());
          vMatrix r_perm = PermuteMatrix<Scalar>(r);
          sccanobjPerm->SetMatrixP(p_perm);
          sccanobjPerm->SetMatrixQ(q_perm);
          sccanobjPerm->SetMatrixR(r_perm);
          //      double permcorr=sccanobjPerm->RunSCCAN2();
          sccanobjPerm->SetSCCANFormulation(sccanobjCovar->GetSCCANFormulation());
          sccanobjPerm->SetAlreadyWhitened(false);
          double permcorr = 0;
          if (newimp == 0)
          {
            permcorr = sccanobjPerm->SparsePartialArnoldiCCA(n_e_vecs);
          }
          else if (newimp == 1)
          {
            permcorr = sccanobjPerm->SparsePartialCCA(n_e_vecs);
          }
          // permcorr=sccanobjPerm->RunSCCAN2multiple(n_e_vecs);//
          // else permcorr=
          // std::cout << " partialed out corr ";
          for (unsigned int ff = 0; ff < sccanobjPerm->GetCanonicalCorrelations().size(); ff++)
          {
            // std::cout << " " << sccanobjPerm->GetCanonicalCorrelations()[ff];
          }
          // std::cout << std::endl;
          if (permcorr > truecorr)
          {
            perm_exceed_ct++;
          }
          /*
               vVector w_p_perm=sccanobjPerm->GetVariateP(0);
               vVector w_q_perm=sccanobjPerm->GetVariateQ(0);
               for (unsigned long j=0; j<p.cols(); j++)
             if ( w_p_perm(j) > sccanobjCovar->GetVariateP(0)(j))
               {
                 w_p_signif_ct(j)=w_p_signif_ct(j)++;
               }
               for (unsigned long j=0; j<q.cols(); j++)
             if ( w_q_perm(j) >  sccanobjCovar->GetVariateQ(0)(j) )
               {
                 w_q_signif_ct(j)=w_q_signif_ct(j)++;
               }
               // end solve cca permutation
               */
          // std::cout << permcorr << " p-value " <<  (double)perm_exceed_ct
          //  / (pct + 1) << " ct " << pct << " true " << truecorr << std::endl;
        }
        unsigned long psigct = 0, qsigct = 0;
        Scalar        pinvtoler = 1.e-6;
        for (unsigned long j = 0; j < sccanobjCovar->GetVariateP(0).size(); j++)
        {
          if (sccanobjCovar->GetVariateP(0)(j) > pinvtoler)
          {
            w_p_signif_ct(j) = 1.0 - (double)w_p_signif_ct(j) / (double)(permct);
            if (w_p_signif_ct(j) > 0.949)
            {
              psigct++;
            }
          }
          else
          {
            w_p_signif_ct(j) = 0;
          }
        }
        for (unsigned long j = 0; j < sccanobjCovar->GetVariateQ(0).size(); j++)
        {
          if (sccanobjCovar->GetVariateQ(0)(j) > pinvtoler)
          {
            w_q_signif_ct(j) = 1.0 - (double)w_q_signif_ct(j) / (double)(permct);
            if (w_q_signif_ct(j) > 0.949)
            {
              qsigct++;
            }
          }
          else
          {
            w_q_signif_ct(j) = 0;
          }
        }
        // std::cout <<  " p-value " <<  (double)perm_exceed_ct / (permct) << " ct " << permct << std::endl;
        // std::cout << " p-vox " <<  (double)psigct / sccanobjCovar->GetVariateP(0).size() << " ct " << permct
        //          << std::endl;
        // std::cout << " q-vox " <<  (double)qsigct / sccanobjCovar->GetVariateP(0).size() << " ct " << permct
        //          << std::endl;
        //	// std::cout << "Here in sccan after printing "<<pct<<std::endl;
      }

      return EXIT_SUCCESS;
    } // run_partial_scca
    // std::cout << " VNL mSCCA " << std::endl;
    sccanobj->SetMatrixP(p);
    sccanobj->SetMatrixQ(q);
    sccanobj->SetMatrixR(r);
    sccanobj->SetMaskImageP(mask1);
    sccanobj->SetMaskImageQ(mask2);
    sccanobj->SetMaskImageR(mask3);
    truecorr = sccanobj->RunSCCAN3();
    vVector w_p = sccanobj->GetPWeights();
    vVector w_q = sccanobj->GetQWeights();
    vVector w_r = sccanobj->GetRWeights();
    // std::cout << " final correlation  " << truecorr  << std::endl;
    //  // std::cout << " Projection-P " << p*w_p << std::endl;
    // // std::cout << " Projection-Q " << q*w_q << std::endl;
    if (leave_out < pin.rows())
    {
      // std::cout << " Projection-leave-P " << dot_product(p_leave_out, w_p) << std::endl;
      // std::cout << " Projection-leave-Q " << dot_product(q_leave_out, w_q) << std::endl;
    }
    //  // std::cout <<  " r weights " << w_r << std::endl;
    for (unsigned long j = 0; j < w_r.size(); j++)
    {
      if (w_r(j) > 0)
      {
        // std::cout << " r-weight " << j << "," << w_r(j) << std::endl;
      }
    }
    if (outputOption)
    {
      std::string filename = outputOption->GetFunction(0)->GetName();
      // std::cout << " write " << filename << std::endl;
      std::string::size_type pos = filename.rfind(".");
      std::string            filepre = std::string(filename, 0, pos);
      std::string            extension = std::string(filename, pos, filename.length() - 1);
      if (extension == std::string(".gz"))
      {
        pos = filepre.rfind(".");
        extension = std::string(filepre, pos, filepre.length() - 1) + extension;
        filepre = std::string(filepre, 0, pos);
      }
      std::string post = std::string("View1vec");
      WriteVariatesToSpatialImage<ImageType, Scalar>(
        filename, post, sccanobj->GetVariatesP(), mask1, sccanobj->GetMatrixP(), have_p_mask, sccanobj->GetMatrixU());
      post = std::string("View2vec");
      WriteVariatesToSpatialImage<ImageType, Scalar>(
        filename, post, sccanobj->GetVariatesQ(), mask2, sccanobj->GetMatrixQ(), have_q_mask, sccanobj->GetMatrixU());
    }

    /** begin permutation 1. q_pvMatrix CqqInv=vnl_svd_inverse<Scalar>(Cqq);
     q=q*CqqInv;
    sermuted ;  2. scca ;  3. test corrs and weights significance */
    unsigned long perm_exceed_ct = 0;
    if (permct > 0)
    {
      vVector w_p_signif_ct(w_p.size(), 0);
      vVector w_q_signif_ct(w_q.size(), 0);
      vVector w_r_signif_ct(w_r.size(), 0);
      for (unsigned long pct = 0; pct <= permct; pct++)
      {
        // 0. compute permutation for q ( switch around rows )
        // // std::cout << " dont permute q " << std::endl;
        vMatrix q_perm = PermuteMatrix<Scalar>(sccanobj->GetMatrixQ());
        vMatrix r_perm = PermuteMatrix<Scalar>(sccanobj->GetMatrixR());
        sccanobj->SetMatrixQ(q_perm);
        sccanobj->SetMatrixR(r_perm);
        double permcorr = sccanobj->RunSCCAN3();
        if (permcorr > truecorr)
        {
          perm_exceed_ct++;
        }
        vVector w_p_perm = sccanobj->GetPWeights();
        vVector w_q_perm = sccanobj->GetQWeights();
        vVector w_r_perm = sccanobj->GetRWeights();
        for (unsigned long j = 0; j < w_r.size(); j++)
        {
          if (w_r_perm(j) > w_r(j))
          {
            w_r_signif_ct(j) = w_r_signif_ct(j)++;
          }
        }
        //      // std::cout << " only testing correlation with biserial predictions " << std::endl;
        // end solve cca permutation
        // std::cout << permcorr << " p-value " <<  (double)perm_exceed_ct
        //  / (pct + 1) << " ct " << pct << " true " << truecorr << std::endl;
        for (unsigned long j = 0; j < w_r.size(); j++)
        {
          if (w_r(j) > 0)
          {
            // std::cout << " r entry " << j << " signif " <<  (double)w_r_signif_ct(j) / (double)(pct + 1) <<
            // std::endl;
          }
        }
      }
    }
    //  // std::cout <<  " p-value " <<  (double)perm_exceed_ct/(permct+1) << " ct " << permct << std::endl;
  }
  return EXIT_SUCCESS;
}

int
sccan(itk::ants::CommandLineParser * sccanparser)
{
  // Define dimensionality
  constexpr unsigned int ImageDimension = 3;
  using matPixelType = double;

  itk::ants::CommandLineParser::OptionType::Pointer outputOption = sccanparser->GetOption("output");
  if (!outputOption || outputOption->GetNumberOfFunctions() == 0)
  {
    // std::cout << "Warning:  no output option set." << std::endl;
  }
  unsigned int                                      permct = 0;
  itk::ants::CommandLineParser::OptionType::Pointer permoption = sccanparser->GetOption("n_permutations");
  if (!permoption || permoption->GetNumberOfFunctions() == 0)
  {
  }
  else
  {
    permct = sccanparser->Convert<unsigned int>(permoption->GetFunction()->GetName());
  }

  unsigned int iterct = 20;
  permoption = sccanparser->GetOption("iterations");
  if (permoption && permoption->GetNumberOfFunctions() > 0)
  {
    iterct = sccanparser->Convert<unsigned int>(permoption->GetFunction()->GetName());
  }

  unsigned int                                      evec_ct = 1;
  itk::ants::CommandLineParser::OptionType::Pointer evec_option = sccanparser->GetOption("n_eigenvectors");
  if (!evec_option || evec_option->GetNumberOfFunctions() == 0)
  {
  }
  else
  {
    evec_ct = sccanparser->Convert<unsigned int>(evec_option->GetFunction()->GetName());
  }

  matPixelType                                      uselong = 0;
  itk::ants::CommandLineParser::OptionType::Pointer long_option = sccanparser->GetOption("uselong");
  if (!long_option || long_option->GetNumberOfFunctions() == 0)
  {
  }
  else
  {
    uselong = sccanparser->Convert<matPixelType>(long_option->GetFunction()->GetName());
  }

  unsigned int                                      covering = 1;
  itk::ants::CommandLineParser::OptionType::Pointer covering_option = sccanparser->GetOption("covering");
  if (!covering_option || covering_option->GetNumberOfFunctions() == 0)
  {
  }
  else
  {
    covering = sccanparser->Convert<unsigned int>(covering_option->GetFunction()->GetName());
  }

  matPixelType                                      usel1 = 0.1;
  itk::ants::CommandLineParser::OptionType::Pointer l1_option = sccanparser->GetOption("l1");
  if (!l1_option || l1_option->GetNumberOfFunctions() == 0)
  {
  }
  else
  {
    usel1 = sccanparser->Convert<matPixelType>(l1_option->GetFunction()->GetName());
  }

  unsigned int                                      robustify = 0;
  itk::ants::CommandLineParser::OptionType::Pointer robust_option = sccanparser->GetOption("robustify");
  if (!robust_option || robust_option->GetNumberOfFunctions() == 0)
  {
  }
  else
  {
    robustify = sccanparser->Convert<unsigned int>(robust_option->GetFunction()->GetName());
  }

  unsigned int                                      verbosity = 0;
  itk::ants::CommandLineParser::OptionType::Pointer verbopt = sccanparser->GetOption("verbose");
  if (!verbopt || verbopt->GetNumberOfFunctions() == 0)
  {
  }
  else
  {
    verbosity = sccanparser->Convert<unsigned int>(verbopt->GetFunction()->GetName());
  }


  itk::ants::CommandLineParser::OptionType::Pointer evecg_option = sccanparser->GetOption("EvecGradPenalty");
  if (!evecg_option || evecg_option->GetNumberOfFunctions() == 0)
  {
  }
  else
  {
    sccanparser->Convert<matPixelType>(evecg_option->GetFunction()->GetName());
  }

  matPixelType                                      smoother = 0;
  itk::ants::CommandLineParser::OptionType::Pointer smooth_option = sccanparser->GetOption("smoother");
  if (!smooth_option || smooth_option->GetNumberOfFunctions() == 0)
  {
  }
  else
  {
    smoother = sccanparser->Convert<matPixelType>(smooth_option->GetFunction()->GetName());
  }

  matPixelType                                      priorWeight = 0;
  itk::ants::CommandLineParser::OptionType::Pointer pwoption = sccanparser->GetOption("prior-weight");
  if (!pwoption || pwoption->GetNumberOfFunctions() == 0)
  {
  }
  else
  {
    priorWeight = sccanparser->Convert<matPixelType>(pwoption->GetFunction()->GetName());
  }

  bool                                              getSmall = false;
  itk::ants::CommandLineParser::OptionType::Pointer getsmallopt = sccanparser->GetOption("get-small");
  if (getsmallopt && getsmallopt->GetNumberOfFunctions() > 0)
  {
    getSmall = sccanparser->Convert<bool>(getsmallopt->GetFunction()->GetName());
  }

  matPixelType                                      row_sparseness = 0;
  itk::ants::CommandLineParser::OptionType::Pointer row_option = sccanparser->GetOption("row_sparseness");
  if (!row_option || row_option->GetNumberOfFunctions() == 0)
  {
  }
  else
  {
    row_sparseness = sccanparser->Convert<matPixelType>(row_option->GetFunction()->GetName());
  }

  unsigned int                                      p_cluster_thresh = 1;
  itk::ants::CommandLineParser::OptionType::Pointer clust_option = sccanparser->GetOption("PClusterThresh");
  if (!clust_option || clust_option->GetNumberOfFunctions() == 0)
  {
  }
  else
  {
    p_cluster_thresh = sccanparser->Convert<unsigned int>(clust_option->GetFunction()->GetName());
  }

  unsigned int q_cluster_thresh = 1;
  clust_option = sccanparser->GetOption("QClusterThresh");
  if (!clust_option || clust_option->GetNumberOfFunctions() == 0)
  {
  }
  else
  {
    q_cluster_thresh = sccanparser->Convert<unsigned int>(clust_option->GetFunction()->GetName());
  }

  itk::ants::CommandLineParser::OptionType::Pointer positivity_option = sccanparser->GetOption("PositivityConstraint");
  if (!positivity_option || positivity_option->GetNumberOfFunctions() == 0)
  {
  }
  else
  {
    sccanparser->Convert<unsigned int>(positivity_option->GetFunction()->GetName());
  }

  bool                                              eigen_imp = false;
  itk::ants::CommandLineParser::OptionType::Pointer eigen_option = sccanparser->GetOption("ridge_cca");
  if (!eigen_option || eigen_option->GetNumberOfFunctions() == 0)
  {
  }
  else
  {
    eigen_imp = sccanparser->Convert<bool>(eigen_option->GetFunction()->GetName());
  }

  //  operations on individual matrices
  itk::ants::CommandLineParser::OptionType::Pointer matrixOption = sccanparser->GetOption("imageset-to-matrix");
  if (matrixOption && matrixOption->GetNumberOfFunctions() > 0)
  {
    std::string outname = outputOption->GetFunction(0)->GetName();
    std::string imagelist = matrixOption->GetFunction(0)->GetParameter(0);
    std::string maskfn = matrixOption->GetFunction(0)->GetParameter(1);
    ConvertImageListToMatrix<ImageDimension, double>(imagelist, maskfn, outname);
    return EXIT_SUCCESS;
  }

  //  operations on individual matrices
  itk::ants::CommandLineParser::OptionType::Pointer matrixOptionTimeSeries =
    sccanparser->GetOption("timeseriesimage-to-matrix");
  if (matrixOptionTimeSeries && matrixOptionTimeSeries->GetNumberOfFunctions() > 0)
  {
    std::string outname = outputOption->GetFunction(0)->GetName();
    std::string imagefn = matrixOptionTimeSeries->GetFunction(0)->GetParameter(0);
    std::string maskfn = matrixOptionTimeSeries->GetFunction(0)->GetParameter(1);
    double      smoother_space = 0;
    if (matrixOptionTimeSeries->GetFunction(0)->GetNumberOfParameters() > 2)
    {
      smoother_space = sccanparser->Convert<double>(matrixOptionTimeSeries->GetFunction(0)->GetParameter(2));
    }
    double smoother_time = 0;
    if (matrixOptionTimeSeries->GetFunction(0)->GetNumberOfParameters() > 3)
    {
      smoother_time = sccanparser->Convert<double>(matrixOptionTimeSeries->GetFunction(0)->GetParameter(3));
    }
    ConvertTimeSeriesImageToMatrix<double>(imagefn, maskfn, outname, smoother_space, smoother_time);
    // std::cout << " outname done " << outname << std::endl;
    return EXIT_SUCCESS;
  }

  //  operations on individual matrices
  itk::ants::CommandLineParser::OptionType::Pointer matrixOptionV2I = sccanparser->GetOption("vector-to-image");
  if (matrixOptionV2I && matrixOptionV2I->GetNumberOfFunctions() > 0)
  {
    std::string outname = outputOption->GetFunction(0)->GetName();
    std::string csvfn = matrixOptionV2I->GetFunction(0)->GetParameter(0);
    std::string maskfn = matrixOptionV2I->GetFunction(0)->GetParameter(1);
    auto        rowOrCol = sccanparser->Convert<unsigned long>(matrixOptionV2I->GetFunction(0)->GetParameter(2));
    ConvertCSVVectorToImage<double>(csvfn, maskfn, outname, rowOrCol);
    // std::cout << " V2I done " << outname << std::endl;
    return EXIT_SUCCESS;
  }

  // p.d.
  itk::ants::CommandLineParser::OptionType::Pointer matrixProjectionOption =
    sccanparser->GetOption("imageset-to-projections");
  if (matrixProjectionOption && matrixProjectionOption->GetNumberOfFunctions() > 0)
  {
    std::string outFilename = outputOption->GetFunction(0)->GetName();
    std::string vecList = matrixProjectionOption->GetFunction(0)->GetParameter(0);
    std::string imageList = matrixProjectionOption->GetFunction(0)->GetParameter(1);
    bool        average = sccanparser->Convert<bool>(matrixProjectionOption->GetFunction(0)->GetParameter(2));
    // // std::cout <<"here" << outFilename << " " << vecList << " " <<imageList << std::endl;
    if (average)
    {
      // std::cout << " doing average instead of dot product " << std::endl;
    }
    ConvertImageVecListToProjection<ImageDimension, double>(vecList, imageList, outFilename, average);
    return EXIT_SUCCESS;
  }

  itk::ants::CommandLineParser::OptionType::Pointer svdOption = sccanparser->GetOption("svd");
  if (svdOption && svdOption->GetNumberOfFunctions() > 0)
  {
    std::string initializationStrategy = svdOption->GetFunction()->GetName();
    if (!initializationStrategy.compare(std::string("sparse")))
    {
      SVD_One_View<ImageDimension, double>(sccanparser,
                                           permct,
                                           evec_ct,
                                           robustify,
                                           p_cluster_thresh,
                                           iterct,
                                           1,
                                           usel1,
                                           row_sparseness,
                                           smoother,
                                           covering,
                                           verbosity);
      return EXIT_SUCCESS;
    }
    if (!initializationStrategy.compare(std::string("cgspca")))
    {
      SVD_One_View<ImageDimension, double>(sccanparser,
                                           permct,
                                           evec_ct,
                                           robustify,
                                           p_cluster_thresh,
                                           iterct,
                                           2,
                                           usel1,
                                           row_sparseness,
                                           smoother,
                                           covering,
                                           verbosity);
      return EXIT_SUCCESS;
    }
    if (!initializationStrategy.compare(std::string("network")))
    {
      SVD_One_View<ImageDimension, double>(sccanparser,
                                           permct,
                                           evec_ct,
                                           robustify,
                                           p_cluster_thresh,
                                           iterct,
                                           4,
                                           usel1,
                                           row_sparseness,
                                           smoother,
                                           covering,
                                           verbosity);
      return EXIT_SUCCESS;
    }
    if (!initializationStrategy.compare(std::string("lasso")))
    {
      SVD_One_View<ImageDimension, double>(sccanparser,
                                           permct,
                                           evec_ct,
                                           robustify,
                                           p_cluster_thresh,
                                           iterct,
                                           5,
                                           usel1,
                                           row_sparseness,
                                           smoother,
                                           covering,
                                           verbosity);
      return EXIT_SUCCESS;
    }
    if (!initializationStrategy.compare(std::string("recon")))
    {
      SVD_One_View<ImageDimension, double>(sccanparser,
                                           permct,
                                           evec_ct,
                                           robustify,
                                           p_cluster_thresh,
                                           iterct,
                                           6,
                                           usel1,
                                           row_sparseness,
                                           smoother,
                                           covering,
                                           verbosity,
                                           getSmall);
      return EXIT_SUCCESS;
    }
    if (!initializationStrategy.compare(std::string("recon4d")))
    {
      SVD_One_View<ImageDimension + 1, double>(sccanparser,
                                               permct,
                                               evec_ct,
                                               robustify,
                                               p_cluster_thresh,
                                               iterct,
                                               6,
                                               usel1,
                                               row_sparseness,
                                               smoother,
                                               covering,
                                               verbosity);
      return EXIT_SUCCESS;
    }
    if (!initializationStrategy.compare(std::string("prior")))
    {
      SVD_One_View<ImageDimension, double>(sccanparser,
                                           permct,
                                           evec_ct,
                                           robustify,
                                           p_cluster_thresh,
                                           iterct,
                                           7,
                                           usel1,
                                           row_sparseness,
                                           smoother,
                                           covering,
                                           verbosity);
      return EXIT_SUCCESS;
    }
    SVD_One_View<ImageDimension, double>(sccanparser,
                                         permct,
                                         evec_ct,
                                         robustify,
                                         p_cluster_thresh,
                                         iterct,
                                         1,
                                         usel1,
                                         row_sparseness,
                                         smoother,
                                         covering,
                                         verbosity);
    return EXIT_SUCCESS;
  }

  // std::cout << " scca-max-iterations " << iterct << " you will assess significance with " << permct << "
  // permutations." << std::endl;
  //  operations on pairs of matrices
  itk::ants::CommandLineParser::OptionType::Pointer matrixPairOption = sccanparser->GetOption("scca");
  if (matrixPairOption && matrixPairOption->GetNumberOfFunctions() > 0)
  {
    if (matrixPairOption && matrixPairOption->GetFunction(0)->GetNumberOfParameters() < 2)
    {
      std::cerr << "  Incorrect number of parameters." << std::endl;
      return EXIT_FAILURE;
    }
    std::string initializationStrategy = matrixPairOption->GetFunction()->GetName();
    // call RCCA_eigen or RCCA_vnl
    unsigned int exitvalue = EXIT_FAILURE;
    if (!initializationStrategy.compare(std::string("two-view")))
    {
      exitvalue = SCCA_vnl<ImageDimension, double>(sccanparser,
                                                   permct,
                                                   evec_ct,
                                                   eigen_imp,
                                                   robustify,
                                                   p_cluster_thresh,
                                                   q_cluster_thresh,
                                                   iterct,
                                                   usel1,
                                                   uselong,
                                                   row_sparseness,
                                                   smoother,
                                                   covering,
                                                   priorWeight,
                                                   verbosity);
    }
    else if (!initializationStrategy.compare(std::string("three-view")))
    {
      // std::cout << " mscca 3-view " << std::endl;
      exitvalue = mSCCA_vnl<ImageDimension, double>(
        sccanparser, permct, false, evec_ct, eigen_imp, robustify, p_cluster_thresh, q_cluster_thresh, iterct);
    }
    else if (!initializationStrategy.compare(std::string("partial")))
    {
      // std::cout << " pscca " << std::endl;
      exitvalue = mSCCA_vnl<ImageDimension, double>(
        sccanparser, permct, true, evec_ct, eigen_imp, robustify, p_cluster_thresh, q_cluster_thresh, iterct);
    }
    else if (!initializationStrategy.compare(std::string("dynsccan")))
    {
      // std::cout << " tscca " << std::endl;
      exitvalue = SCCA_vnl<ImageDimension + 1, double>(sccanparser,
                                                       permct,
                                                       evec_ct,
                                                       eigen_imp,
                                                       robustify,
                                                       p_cluster_thresh,
                                                       q_cluster_thresh,
                                                       iterct,
                                                       usel1,
                                                       uselong,
                                                       row_sparseness,
                                                       smoother,
                                                       covering,
                                                       verbosity);
    }
    else
    {
      std::cerr << " unrecognized option in matrixPairOperation " << std::endl;
      return exitvalue;
    }
    // std::cout << " exit value " << exitvalue << std::endl;
    return exitvalue;
  } // matrixPairOption
  else
  {
    std::cerr << " no option specified " << std::endl;
  }
  return EXIT_FAILURE;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
sccan(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the sccanparser should handle
  args.insert(args.begin(), "sccan");

  std::ignore = std::remove(args.begin(), args.end(), std::string(""));
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

  itk::ants::CommandLineParser::Pointer sccanparser = itk::ants::CommandLineParser::New();

  sccanparser->SetCommand(argv[0]);

  std::string commandDescription =
    std::string("A tool for sparse statistical analysis on images : ") +
    std::string(
      " scca, pscca (with options), mscca.  Can also convert an imagelist/mask pair to a binary matrix image.  ");

  sccanparser->SetCommandDescription(commandDescription);

  /** in this function, list all the operations you will perform */

  using OptionType = itk::ants::CommandLineParser::OptionType;
  {
    std::string         description = std::string("Print the help menu (short version).");
    OptionType::Pointer option = OptionType::New();
    option->SetShortName('h');
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string         description = std::string("Print the help menu (long version).");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("help");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string         description = std::string("Output dependent on which option is called.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("output");
    option->SetShortName('o');
    option->SetUsageOption(0, "outputImage");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string         description = std::string("Number of permutations to use in scca.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("n_permutations");
    option->SetShortName('p');
    option->SetUsageOption(0, "500");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string         description = std::string("Smoothing function for variates");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("smoother");
    option->SetShortName('s');
    option->SetUsageOption(0, "0");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string description =
      std::string("Row sparseness - if (+) then keep values (+) otherwise allow +/- values --- always L1");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("row_sparseness");
    option->SetShortName('z');
    option->SetUsageOption(0, "0");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string         description = std::string("Max iterations for scca optimization (min 20).");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("iterations");
    option->SetShortName('i');
    option->SetUsageOption(0, "20");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string         description = std::string("Number of eigenvectors to compute in scca/spca.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("n_eigenvectors");
    option->SetShortName('n');
    option->SetUsageOption(0, "2");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string         description = std::string("rank-based scca");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("robustify");
    option->SetShortName('r');
    option->SetUsageOption(0, "0");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string         description = std::string("try to make the decomposition cover the whole domain, if possible ");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("covering");
    option->SetShortName('c');
    option->SetUsageOption(0, "0");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string         description = std::string("use longitudinal formulation ( > 0 ) or not ( <= 0 ) ");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("uselong");
    option->SetShortName('g');
    option->SetUsageOption(0, "0");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string description =
      std::string("use l1 ( > 0 ) or l0 ( < 0 ) penalty, also sets gradient step size e.g. -l 0.5 ( L1 ) , -l -0.5 "
                  "(L0)  will set 0.5 grad descent step for either penalty");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("l1");
    option->SetShortName('l');
    option->SetUsageOption(0, "0");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string         description = std::string("cluster threshold on view P");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("PClusterThresh");
    option->SetUsageOption(0, "1");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }
  {
    std::string         description = std::string("cluster threshold on view Q");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("QClusterThresh");
    option->SetUsageOption(0, "1");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string         description = std::string("Ridge cca.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("ridge_cca");
    option->SetShortName('e');
    option->SetUsageOption(0, "0");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string description = std::string("Initialization file list for Eigenanatomy - must also pass mask option");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("initialization");
    option->SetUsageOption(0, "NA");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string description =
      std::string("Initialization file list for SCCAN-Eigenanatomy - must also pass mask option");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("initialization2");
    option->SetUsageOption(0, "NA");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string         description = std::string("Mask file for Eigenanatomy initialization");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("mask");
    option->SetUsageOption(0, "NA");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string         description = std::string("Mask file for Eigenanatomy initialization 2");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("mask2");
    option->SetUsageOption(0, "NA");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string         description = std::string("Choices for pscca: PQ, PminusRQ, PQminusR, PminusRQminusR ");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("partial-scca-option");
    option->SetUsageOption(0, "PminusRQ");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string description = std::string("Scalar value weight on prior between 0 (prior is weak) and 1 (prior is "
                                          "strong).  Only engaged if initialization is used.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("prior-weight");
    option->SetUsageOption(0, "0.0");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string         description = std::string("Find smallest eigenvectors");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("get-small");
    option->SetUsageOption(0, "0.0");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string         description = std::string("set whether output is verbose");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("verbose");
    option->SetShortName('v');
    option->SetUsageOption(0, "0");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string description = std::string("takes a list of image files names (one per line) ") +
                              std::string("and converts it to a 2D matrix / image in binary or csv format depending on "
                                          "the filetype used to define the output.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("imageset-to-matrix");
    option->SetUsageOption(0, "[list.txt,mask.nii.gz]");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string description = std::string("takes a timeseries (4D) image ") +
                              std::string("and converts it to a 2D matrix csv format as output.") +
                              std::string("If the mask has multiple labels ( more the one ) then the average time "
                                          "series in each label will be computed and put in the csv.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("timeseriesimage-to-matrix");
    option->SetUsageOption(
      0,
      "[four_d_image.nii.gz,three_d_mask.nii.gz, optional-spatial-smoothing-param-in-spacing-units-default-zero, "
      "optional-temporal-smoothing-param-in-time-series-units-default-zero  ]");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string description = std::string(
      "converts the 1st column vector in a csv file back to an image --- currently needs the csv file to have > 1 "
      "columns.  if the number of entries in the column does not equal the number of entries in the mask but the "
      "number of rows does equal the number of entries in the mask, then it will convert the row vector to an image. ");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("vector-to-image");
    option->SetUsageOption(0, "[vector.csv,three_d_mask.nii.gz, which-row-or-col ]");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  // p.d.
  {
    std::string description = std::string("takes a list of image and projection files names (one per line) ") +
                              std::string("and writes them to a  csv file --- basically computing X*Y (matrices).");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("imageset-to-projections");
    option->SetUsageOption(0, "[list_projections.txt,list_images.txt, bool do-average-not-real-projection ]");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string description =
      std::string("Matrix-based scca operations for 2 and 3 views.") +
      std::string("For all these options, the FracNonZero terms set the fraction of variables to use in the estimate. "
                  "E.g. if one sets 0.5 then half of the variables will have non-zero values.  If the FracNonZero is "
                  "(+) then the weight vectors must be positive.  If they are negative, weights can be (+) or (-).  "
                  "partial does partial scca for 2 views while partialing out the 3rd view. ");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("scca");
    option->SetUsageOption(0, "two-view[matrix-view1.mhd,matrix-view2.mhd,mask1,mask2,FracNonZero1,FracNonZero2] ");
    option->SetUsageOption(1,
                           "three-view[matrix-view1.mhd,matrix-view2.mhd,matrix-view3.mhd,mask1,mask2,mask3,"
                           "FracNonZero1,FracNonZero2,FracNonZero3]");
    option->SetUsageOption(2,
                           "partial[matrix-view1.mhd,matrix-view2.mhd,matrix-view3.mhd,mask1,mask2,mask3,FracNonZero1,"
                           "FracNonZero2,FracNonZero3]");
    option->SetUsageOption(3, "dynsccan[matrix-view1.mhd,matrix-view2.mhd,mask1,mask2,FracNonZero1,FracNonZero2] ");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  {
    std::string description = std::string("a sparse svd implementation --- will report correlation of eigenvector with "
                                          "original data columns averaged over columns with non-zero weights.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("svd");
    option->SetUsageOption(0,
                           "sparse[matrix-view1.mhd,mask1,FracNonZero1,nuisance-matrix] --- will only use view1 ... "
                           "unless nuisance matrix is specified.");
    option->SetUsageOption(1,
                           "classic[matrix-view1.mhd,mask1,FracNonZero1,nuisance-matrix] --- will only use view1 ... "
                           "unless nuisance matrix is specified.");
    option->SetUsageOption(
      2,
      "cgspca[matrix-view1.mhd,mask1,FracNonZero1,nuisance-matrix] --- will only use view1 ... unless nuisance matrix "
      "is specified, -i controls the number of sparse approximations per eigenvector, -n controls the number of "
      "eigenvectors.  total output will then be  i*n sparse eigenvectors.");
    option->SetUsageOption(3,
                           "prior[ matrix.mha , mask.nii.gz , PriorList.txt , PriorScale.csv , PriorWeightIn0to1 , "
                           "sparseness ] ... if sparseness is set to zero, we take sparseness from the priors.");
    option->SetUsageOption(4, "network[matrix-view1.mhd,mask1,FracNonZero1,guidance-matrix]");
    option->SetUsageOption(5, "lasso[matrix-view1.mhd,mask1,Lambda,guidance-matrix]");
    option->SetUsageOption(6, "recon[matrix-view1.mhd,mask1,FracNonZero1,nuisance-matrix]");
    option->SetUsageOption(7, "recon4d[matrix-view1.mhd,mask1,FracNonZero1,nuisance-matrix]");
    option->SetDescription(description);
    sccanparser->AddOption(option);
  }

  if (sccanparser->Parse(argc, argv) == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  }

  // Print the entire help menu
  itk::ants::CommandLineParser::OptionType::Pointer shortHelpOption = sccanparser->GetOption('h');
  itk::ants::CommandLineParser::OptionType::Pointer longHelpOption = sccanparser->GetOption("help");
  if (argc < 2 || (shortHelpOption->GetFunction() &&
                   sccanparser->Convert<unsigned int>(shortHelpOption->GetFunction()->GetName()) == 1))
  {
    sccanparser->PrintMenu(std::cout, 5, true);
    if (argc < 2)
    {
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }
  if (longHelpOption->GetFunction() &&
      sccanparser->Convert<unsigned int>(longHelpOption->GetFunction()->GetName()) == 1)
  {
    sccanparser->PrintMenu(std::cout, 5, false);
    return EXIT_SUCCESS;
  }

  // Print the long help menu for specific items
  if (longHelpOption && longHelpOption->GetNumberOfFunctions() > 0 &&
      sccanparser->Convert<unsigned int>(longHelpOption->GetFunction()->GetName()) != 0)
  {
    itk::ants::CommandLineParser::OptionListType options = sccanparser->GetOptions();
    for (unsigned int n = 0; n < longHelpOption->GetNumberOfFunctions(); n++)
    {
      std::string                                                  value = longHelpOption->GetFunction(n)->GetName();
      itk::ants::CommandLineParser::OptionListType::const_iterator it;
      for (it = options.begin(); it != options.end(); ++it)
      {
        if ((*it)->GetLongName().find(value) == 0)
        {
          sccanparser->PrintMenu(std::cout, 5, false);
        }
      }
    }
    return EXIT_FAILURE;
  }

  // Call main routine
  return sccan(sccanparser);
}

// now compute covariance matrices
// covariance matrix --- Cov(X, Y) = E[XY] - E[X].E[Y]
/* input matrix
ta<-matrix(c(-1,1,2,2,-2,3,1,1,4,0,3,4),nrow=3,ncol=4)
ta<-matrix(c(-1,1,2,-2,3,1,4,0,3),nrow=3,ncol=3)

> ta
   [,1] [,2] [,3]
[1,]   -1    1    2
[2,]   -2    3    1
[3,]    4    0    3

// cov(ta,ta)
        [,1]      [,2] [,3]
[1,] 10.333333 -4.166667  3.0
[2,] -4.166667  2.333333 -1.5
[3,]  3.000000 -1.500000  1.0

> cov(a[1,]-mean(a[1,]),a[1,]-mean(a[1,]))
[1] 10.33333
v<-a[1,]-mean(a[1,])
> v*v
[1]  1.777778  5.444444 13.444444
> sum(v*v)
[1] 20.66667
> sum(v*v)/(2)  <-  sum( v*v ) / ( n - 1 ) = covariance of two vectors ...
[1] 10.33333
*/

/** try  q.colwise.sum() to compute mean ... */
/** try  q.rowwise.sum() to compute mean ...

eMatrix testM(3,3);
testM(0,0)=-1;testM(1,0)=1; testM(2,0)=2;
testM(0,1)=-2;testM(1,1)=3; testM(2,1)=1;
testM(0,2)= 4;testM(1,2)=0; testM(2,2)=3;
p=testM;
pMatSize[0]=3;
pMatSize[1]=3;
*/

/*
            //1.0/(double)q.columns(); //randgen.drand32();
for (unsigned int it=0; it<4; it++)
{
  //    // std::cout << " 2norm(v0) " << v_0.two_norm() << std::endl;
  vVector v_1=(q)*v_0;
  double vnorm=v_1.two_norm();
  // std::cout << " von " << vnorm << std::endl;
  v_0=v_1/(vnorm);
  // std::cout << " vo " << v_0 << std::endl;
// check if power method works ....
vVector Xv=q*v_0;
Scalar vdotXv = dot_product(v_0,Xv);
// std::cout << " vdotXv " << vdotXv << std::endl;
vVector Xv2=Xv-v_0*vdotXv;
// this value should be small -- i.e. v_0 is an eigenvector of X
// std::cout << " init eigenvector result " << Xv2.squared_magnitude() << std::endl;}
*/
//  // std::cout << v_0 << std::endl;
/*

function [Up,Sp,Vp] = rank_one_svd_update( U, S, V, a, b, force_orth )
% function [Up,Sp,Vp] = rank_one_svd_update( U, S, V, a, b, force_orth )
%
% Given the SVD of
%
%   X = U*S*V'
%
% update it to be the SVD of
%
%   X + ab' = Up*Sp*Vp'
%
% that is, implement a rank-one update to the SVD of X.
%
% Depending on a,b there may be considerable structure that could
% be exploited, but which this code does not.
%
% The subspace rotations involved may not preserve orthogonality due
% to numerical round-off errors.  To compensate, you can set the
% "force_orth" flag, which will force orthogonality via a QR plus
% another SVD.  In a long loop, you may want to force orthogonality
% every so often.
%
% See Matthew Brand, "Fast low-rank modifications of the thin
% singular value decomposition".
%
% D. Wingate 8/17/2007
%

current_rank = size( U, 2 );

% P is an orthogonal basis of the column-space
% of (I-UU')a, which is the component of "a" that is
% orthogonal to U.
m = U' * a;
p = a - U*m;
Ra = sqrt(p'*p);
P = (1/Ra)*p;

% XXX this has problems if a is already in the column space of U!
% I don't know what to do in that case.
if ( Ra < 1e-13 )
 fprintf('------> Whoa! No orthogonal component of m!\n');
end;

% Q is an orthogonal basis of the column-space
% of (I-VV')b.
n = V' * b;
q = b - V*n;
Rb = sqrt(q'*q);
Q = (1/Rb)*q;

if ( Rb < 1e-13 )
 fprintf('------> Whoa! No orthogonal component of n!\n');
end;

%
% Diagonalize K, maintaining rank
%

% XXX note that this diagonal-plus-rank-one, so we should be able
% to take advantage of the structure!
z = zeros( size(m) );

K = [ S z ; z' 0 ] + [ m; Ra ]*[ n; Rb ]';

[tUp,tSp,tVp] = svds( K, current_rank );

%
% Now update our matrices!
%

Sp = tSp;

Up = [ U P ] * tUp;
Vp = [ V Q ] * tVp;

% The above rotations may not preserve orthogonality, so we explicitly
% deal with that via a QR plus another SVD.  In a long loop, you may
% want to force orthogonality every so often.

if ( force_orth )
 [UQ,UR] = qr( Up, 0 );
 [VQ,VR] = qr( Vp, 0 );
 [tUp,tSp,tVp] = svds( UR * Sp * VR', current_rank );
 Up = UQ * tUp;
 Vp = VQ * tVp;
 Sp = tSp;
end;

return;
*/
// } // namespace antssccan
} // namespace ants
