#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "vnl/algo/vnl_determinant.h"

#include "ReadWriteData.h"

#include "vnl/algo/vnl_determinant.h"
#include "itkDiscreteGaussianImageFilter.h"

#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkVectorCurvatureAnisotropicDiffusionImageFilter.h"

#include "itkMatrixOffsetTransformBase.h"
#include "itkWarpImageMultiTransformFilter.h"

namespace ants
{
template <typename TImage>
typename TImage::Pointer
VectorAniDiff(typename TImage::Pointer img, unsigned int iters)
{
  double timeStep = 0.065;

  using VectorImageType = TImage;
  using FilterType = itk::VectorCurvatureAnisotropicDiffusionImageFilter<VectorImageType, VectorImageType>;

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(img);
  filter->SetNumberOfIterations(iters);
  filter->SetTimeStep(timeStep);
  filter->SetConductanceParameter(1.0);
  filter->Update();
  // Software Guide : EndCodeSnippet

  return filter->GetOutput();
}

template <typename TImage>
typename TImage::Pointer
GenerateGridImage(TImage * img, unsigned int gridsize)
{
  using ImageType = TImage;
  enum
  {
    ImageDimension = TImage::ImageDimension
  };

  itk::ImageRegionIteratorWithIndex<ImageType> wimIter(img, img->GetLargestPossibleRegion());
  wimIter.GoToBegin();
  for (; !wimIter.IsAtEnd(); ++wimIter)
  {
    wimIter.Set(2);
  }

  wimIter.GoToBegin();
  for (; !wimIter.IsAtEnd(); ++wimIter)
  {
    typename ImageType::IndexType ind = wimIter.GetIndex();
    for (int i = 0; i < 2; i++)
    {
      if (ind[i] % gridsize == 0)
      {
        wimIter.Set(0);
      }
      //          if (ind[i] % (gridsize+1) == 0) wimIter.Set(0);// tartan
    }
  }

  wimIter.GoToBegin();
  for (; !wimIter.IsAtEnd(); ++wimIter)
  {
    typename ImageType::IndexType ind = wimIter.GetIndex();
    typename ImageType::IndexType ind2 = wimIter.GetIndex();
    for (int i = 0; i < 2; i++)
    {
      ind2[i] = ind[i] - 1;
      if (ind2[i] < 0)
      {
        ind2[i] = 0;
      }
    }
    for (int i = 0; i < 2; i++)
    {
      // this creates a 3-d effect
      // if (ind[i] % gridsize == 0) img->SetPixel(ind2,2);
      // this gives double thickness
      if (ind[i] % gridsize == 0)
      {
        img->SetPixel(ind2, 0);
      }
    }
  }

  return img;
}

template <typename ImageType>
typename ImageType::Pointer
ReadAnImage(char * fn)
{
  // Read the image files begin
  using FileSourceType = itk::ImageFileReader<ImageType>;
  typename FileSourceType::Pointer reffilter = FileSourceType::New();
  reffilter->SetFileName(fn);
  try
  {
    reffilter->Update();
  }
  catch (...)
  {
    return nullptr;
  }
  return reffilter->GetOutput();
}

template <typename TImage, typename TDisplacementField>
typename TDisplacementField::PixelType
TransformVector(TDisplacementField * field, typename TImage::IndexType index)
{
  enum
  {
    ImageDimension = TImage::ImageDimension
  };
  typename TDisplacementField::PixelType vec = field->GetPixel(index);
  /* buggy code from before */
  typename TDisplacementField::PixelType newvec;
  newvec.Fill(0);
  for (unsigned int row = 0; row < ImageDimension; row++)
  {
    for (unsigned int col = 0; col < ImageDimension; col++)
    {
      newvec[row] += vec[col] * static_cast<float>(field->GetDirection()[row][col]);
    }
  }
  return newvec;
}

template <typename TImage, typename TDisplacementField>
typename TDisplacementField::PixelType
ProjectVector(typename TDisplacementField::PixelType invec, typename TDisplacementField::PixelType projvec)
{
  enum
  {
    ImageDimension = TImage::ImageDimension
  };
  typename TDisplacementField::PixelType newvec;
  double                                 ip = 0;
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    ip += static_cast<double>(invec[i] * projvec[i]);
  }
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    newvec[i] = static_cast<float>(ip) * projvec[i];
  }
  return newvec;
}

void
antsjacobiansplit(const std::string & s, char c, std::vector<std::string> & v)
{
  std::string::size_type i = 0;
  std::string::size_type j = s.find(c);

  while (j != std::string::npos)
  {
    v.push_back(s.substr(i, j - i));
    i = ++j;
    j = s.find(c, j);
    if (j == std::string::npos)
    {
      v.push_back(s.substr(i, s.length()));
    }
  }
}

template <typename TImage, typename TDisplacementField>
void
ComputeJacobian(TDisplacementField * field,
                char *               fnm,
                char *               maskfn,
                bool                 uselog = false,
                bool                 norm = false,
                std::string          projvec = "")
{
  std::vector<std::string> v;

  if (projvec.length() > 2)
  {
    antsjacobiansplit(projvec, 'x', v);
    for (const auto & i : v)
    {
      std::cout << i << '\n';
    }
  }

  using ImageType = TImage;
  using FieldType = TDisplacementField;
  enum
  {
    ImageDimension = TImage::ImageDimension
  };
  using FloatImageType = itk::Image<float, ImageDimension>;
  typename FloatImageType::Pointer mask = nullptr;
  typename FieldType::PixelType    pvec;
  if (!v.empty())
  {
    for (unsigned int i = 0; i < ImageDimension; i++)
    {
      pvec[i] = atof(v[i].c_str());
    }
    pvec = pvec / pvec.GetNorm();
    std::cout << " using projection vector " << pvec << std::endl;
  }
  mask = ReadAnImage<FloatImageType>(maskfn);

  if (!field)
  {
    return;
  }
  typename TImage::SizeType    s = field->GetLargestPossibleRegion().GetSize();
  typename TImage::SpacingType sp = field->GetSpacing();

  typename FloatImageType::Pointer m_FloatImage = AllocImage<FloatImageType>(field, 0);

  if (false)
  {
    using TransformType = itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension>;
    using WarperType = itk::WarpImageMultiTransformFilter<ImageType, ImageType, FieldType, TransformType>;
    typename WarperType::Pointer warper = WarperType::New();
    warper->SetInput(nullptr);
    warper->SetEdgePaddingValue(0);
    warper->SetSmoothScale(1);
    warper->PushBackDisplacementFieldTransform(field);
    warper->SetOutputParametersFromImage(field);
    warper->Update();
    //    grid=warper->GetOutput();
    using writertype = itk::ImageFileWriter<ImageType>;
    typename writertype::Pointer writer = writertype::New();
    std::string                  fng = std::string(fnm) + "grid.nii.gz";
    writer->SetFileName(fng.c_str());
    writer->SetInput(nullptr);
    writer->Write();
    std::cout << " Grid done ";
  }
  typename FloatImageType::SizeType m_FieldSize = field->GetLargestPossibleRegion().GetSize();

  using Iterator = itk::ImageRegionIteratorWithIndex<FloatImageType>;
  Iterator wimIter(m_FloatImage, m_FloatImage->GetLargestPossibleRegion());
  wimIter.GoToBegin();
  for (; !wimIter.IsAtEnd(); ++wimIter)
  {
    wimIter.Set(1.0);
  }

  using MatrixType = vnl_matrix<double>;
  MatrixType jMatrix, idMatrix, avgMatrix;
  jMatrix.set_size(ImageDimension, ImageDimension);
  avgMatrix.set_size(ImageDimension, ImageDimension);
  avgMatrix.fill(0);
  itk::ImageRegionIteratorWithIndex<TDisplacementField> m_FieldIter(field, field->GetLargestPossibleRegion());
  typename TImage::IndexType                            rindex;
  typename TImage::IndexType                            ddrindex;
  typename TImage::IndexType                            ddlindex;

  typename TImage::IndexType difIndex[ImageDimension][2];

  double       det = 0.0;
  unsigned int posoff = 1;
  float        space = 1.0;


  typename FieldType::PixelType dPix;
  typename FieldType::PixelType lpix;
  typename FieldType::PixelType rpix;
  typename FieldType::PixelType cpix;

  //   double totaljac=0.0;

  // /the finite difference equations
  // float wC, wLL, wL, wR, wRR;
  // 3rd deriv - 4th order
  // wC = 0.0;
  // wLL = 1.; wL = -2.0; wR =  2.0; wRR = -1.0;
  // 4th deriv - 4th order
  // wC = -6.0;
  // wLL = 1.; wL = -4.0; wR = -4.0; wRR = 1.0;
  // 2nd deriv - 4th order
  // wC = 30.0;
  // wLL = -1.0; wL = 16.0; wR = 16.0; wRR = -1.0;

  unsigned long ct = 0;
  for (m_FieldIter.GoToBegin(); !m_FieldIter.IsAtEnd(); ++m_FieldIter)
  {
    rindex = m_FieldIter.GetIndex();
    float mindist = 1.0;
    bool  oktosample = true;
    float dist = 100.0;
    for (unsigned int row = 0; row < ImageDimension; row++)
    {
      dist = fabs((float)rindex[row]);
      if (dist < mindist)
      {
        oktosample = false;
      }
      dist = fabs((float)s[row] - (float)rindex[row]);
      if (dist < mindist)
      {
        oktosample = false;
      }
    }
    if (oktosample)
    {
      ct++;
      cpix = TransformVector<ImageType, FieldType>(field, rindex);
      if (!v.empty())
      {
        cpix = ProjectVector<ImageType, FieldType>(cpix, pvec);
      }
      for (unsigned int row = 0; row < ImageDimension; row++)
      {
        difIndex[row][0] = rindex;
        difIndex[row][1] = rindex;
        ddrindex = rindex;
        ddlindex = rindex;
        if ((unsigned int)rindex[row] < (unsigned int)m_FieldSize[row] - 2)
        {
          difIndex[row][0][row] = rindex[row] + posoff;
          ddrindex[row] = rindex[row] + posoff * 2;
        }
        if (rindex[row] > 1)
        {
          difIndex[row][1][row] = rindex[row] - 1;
          ddlindex[row] = rindex[row] - 2;
        }

        float h = 1;
        space = 1.0; // should use image spacing here?

        rpix = TransformVector<ImageType, FieldType>(field, difIndex[row][1]);
        if (!v.empty())
        {
          rpix = ProjectVector<ImageType, FieldType>(rpix, pvec);
        }
        rpix = rpix * h + cpix * (1.f - h);
        lpix = TransformVector<ImageType, FieldType>(field, difIndex[row][0]);
        if (!v.empty())
        {
          lpix = ProjectVector<ImageType, FieldType>(lpix, pvec);
        }
        lpix = lpix * h + cpix * (1.f - h);
        //    dPix = ( rpix - lpix)*(1.0)/(2.0);

        //    rrpix = TransformVector<ImageType,FieldType>(field,ddrindex);
        // rrpix = rrpix*h+rpix*(1.-h);
        // llpix = TransformVector<ImageType,FieldType>(field,ddlindex);
        // llpix = llpix*h+lpix*(1.-h);
        //      dPix=( rrpix*(-1.0) + rpix*8.0 - lpix*8.0 + lpix )*(-1.0)*space/(12.0*h); //4th order centered
        // difference
        dPix = (lpix - rpix) * space / (2.0f * h); // 4th order centered difference
        for (unsigned int col = 0; col < ImageDimension; col++)
        {
          float val;
          if (row == col)
          {
            val = static_cast<float>(dPix[col]) / static_cast<float>(sp[col]) + 1.0f;
          }
          else
          {
            val = static_cast<float>(dPix[col]) / static_cast<float>(sp[col]);
          }
          //        std::cout << " row " << row << " col " << col << " val " << val << std::endl;
          jMatrix.put(col, row, val);
          avgMatrix.put(col, row, avgMatrix.get(col, row) + static_cast<double>(val));
        }
      }

      // the determinant of the jacobian matrix
      // std::cout << " get det " << std::endl;
      det = vnl_determinant(jMatrix);
      //    float prodval = m_FloatImage->GetPixel(rindex);
      if (det < 0.0)
      {
        det = 0;
      }

      m_FloatImage->SetPixel(rindex, det);

      // totaljac+=det;
    } // oktosample if
  }

  if (norm && mask)
  {
    std::cout << " using mask && normalizing " << std::endl;
    /*
    typedef itk::DiscreteGaussianImageFilter<TImage, TImage> dgf;
    float sig=2.0;
    typename FloatImageType::Pointer temp;
    {
    typename dgf::Pointer filter = dgf::New();
    filter->SetVariance(sig);
    filter->SetUseImageSpacing(false);
    filter->SetMaximumError(.01f);
    filter->SetInput(m_FloatImage);
    filter->Update();
    //  m_FloatImage=filter->GetOutput();
    temp=filter->GetOutput();
    }
    {
    typename dgf::Pointer filter = dgf::New();
    filter->SetVariance(sig);
    filter->SetUseImageSpacing(false);
    filter->SetMaximumError(.01f);
    filter->SetInput(temp);
    filter->Update();
    //  m_FloatImage=filter->GetOutput();
    temp=filter->GetOutput();
    } */

    double        total = 0.0;
    unsigned long _ct = 0;
    for (m_FieldIter.GoToBegin(); !m_FieldIter.IsAtEnd(); ++m_FieldIter)
    {
      rindex = m_FieldIter.GetIndex();
      if (mask->GetPixel(rindex) > 0)
      {
        total += static_cast<double>(m_FloatImage->GetPixel(rindex));
        _ct++;
      }
      else
      {
        m_FloatImage->SetPixel(rindex, 0);
      }
    }
    total /= (double)_ct;
    for (m_FieldIter.GoToBegin(); !m_FieldIter.IsAtEnd(); ++m_FieldIter)
    {
      rindex = m_FieldIter.GetIndex();
      double val = static_cast<double>(m_FloatImage->GetPixel(rindex)) / total;
      if (mask->GetPixel(rindex) > 0)
      {
        m_FloatImage->SetPixel(rindex, val);
      }
      else
      {
        m_FloatImage->SetPixel(rindex, 0);
      }
    }
  }
  for (m_FieldIter.GoToBegin(); !m_FieldIter.IsAtEnd(); ++m_FieldIter)
  {
    rindex = m_FieldIter.GetIndex();
    double val = m_FloatImage->GetPixel(rindex);
    if (uselog && val > 0)
    {
      val = log(val);
    }
    else if (uselog && val < 0)
    {
      val = log(0.01);
    }
    if (uselog)
    {
      m_FloatImage->SetPixel(rindex, val);
    }
  }

  using writertype = itk::ImageFileWriter<TImage>;
  typename writertype::Pointer writer = writertype::New();
  std::string                  fn = std::string(fnm) + "jacobian.nii.gz";
  if (uselog)
  {
    fn = std::string(fnm) + "logjacobian.nii.gz";
  }
  writer->SetFileName(fn.c_str());
  writer->SetInput(m_FloatImage);
  writer->Write();
}

template <unsigned int ImageDimension>
int
Jacobian(int argc, char * argv[])
{
  //  std::cout << " enter " << ImageDimension << std::endl;
  if (argc < 3)
  {
    std::cout << "Usage:   Jacobian gWarp outfile uselog maskfn normbytotalbool VectorToProjectWarpAgainst "
              << std::endl;
    std::cout << " VectorToProjectWarpAgainst should be in the form 1.0x0.0x0.0 where x separates vector components "
              << std::endl;
    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }
  using PixelType = float;
  using VectorType = itk::Vector<float, ImageDimension>;
  using FieldType = itk::Image<VectorType, ImageDimension>;
  using ImageType = itk::Image<PixelType, ImageDimension>;

  using ReaderType = itk::ImageFileReader<FieldType>;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[1]);
  reader->Update();
  typename FieldType::Pointer gWarp = reader->GetOutput();
  bool                        uselog = false;
  if (argc > 3)
  {
    uselog = (bool)atoi(argv[3]);
  }
  bool norm = false;
  if (argc > 5)
  {
    norm = (bool)atoi(argv[5]);
  }
  std::string projvec;
  if (argc > 6)
  {
    projvec = std::string(argv[6]);
  }
  std::string maskfn = std::string("ThugLifeIsDeadToMe---2pac");
  if (argc > 4)
  {
    maskfn = std::string(argv[4]);
  }

  //  std::cout << " name "<< argv[2] <<  " mask " << argv[4] << " norm " << norm << " Log " << uselog << std::endl;
  ComputeJacobian<ImageType, FieldType>(gWarp, argv[2], const_cast<char *>(maskfn.c_str()), uselog, norm, projvec);
  //  DiffeomorphicJacobian<ImageType,ImageType,FieldType>(gWarp,1,argv[2]);
  //  if (argc > 3) DiffeomorphicMetric<ImageType,ImageType,FieldType>(gWarp,argv[2]);

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
ANTSJacobian(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "ANTSJacobian");
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
  std::cout << "WARNING! " << argv[0]
            << " may not be working correctly, see CreateJacobianDeterminantImage for an alternative method  "
            << std::endl;
  // std::cout << "Please use CreateJacobianDeterminantImage " << std::endl;
  // return EXIT_SUCCESS;

  if (argc < 3)
  {
    std::cout << "Usage: " << argv[0] << " ImageDim gWarp outfile uselog maskfn normbytotalbool projectionvector "
              << std::endl;
    std::cout << " for example " << std::endl
              << " ANTSJacobian 3  myWarp.nii   Output  1   templatebrainmask.nii   1 1x0 " << std::endl;
    std::cout << " the last 1 normalizes the jacobian by the total in the mask.  use this to adjust for head size. 1x0 "
                 "will project the warp along direction 1,0 --- don't add this option if you dont want to do this "
              << std::endl;
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
      Jacobian<2>(argc - 1, argv + 1);
    }
    break;
    case 3:
    {
      Jacobian<3>(argc - 1, argv + 1);
    }
    break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
} // namespace ants
