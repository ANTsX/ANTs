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
#include "antsAllocImage.h"
#include <algorithm>

#include "antsCommandLineOption.h"
#include "antsCommandLineParser.h"
#include "antsSCCANObject.h"
#include "ReadWriteData.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include "itkVariableLengthVector.h"
#include "itkVectorContainer.h"

namespace ants
{

template <typename TComp>
double
vnl_pearson_corr(vnl_vector<TComp> v1, vnl_vector<TComp> v2)
{
  double xysum = 0;

  for (unsigned int i = 0; i < v1.size(); i++)
  {
    xysum += static_cast<double>(v1(i) * v2(i));
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

template <typename NetworkType>
bool
RegionSCCA(typename NetworkType::Pointer network,
           typename NetworkType::Pointer time,
           typename NetworkType::Pointer labels,
           unsigned int                  nLabels,
           unsigned int                  minRegionSize,
           unsigned int                  n_evec,
           unsigned int                  iterct,
           float                         sparsity,
           bool                          robust,
           bool                          useL1,
           float                         gradstep,
           bool                          keepPositive,
           unsigned int                  minClusterSize)
{
  using SCCANType = itk::ants::antsSCCANObject<NetworkType, double>;
  using MatrixType = typename SCCANType::MatrixType;
  using VectorType = typename SCCANType::VectorType;

  // Determine the number of regions to examine
  std::set<unsigned int>                         labelset;
  itk::ImageRegionIteratorWithIndex<NetworkType> it(labels, labels->GetLargestPossibleRegion());

  if (nLabels == 0)
  {

    while (!it.IsAtEnd())
    {
      if (it.Value() > 0)
      {
        labelset.insert(it.Value());
      }
      ++it;
    }
  }
  else
  {
    for (unsigned int i = 0; i < nLabels; i++)
    {
      labelset.insert(i + 1);
    }
  }

  unsigned int N = labelset.size();
  std::cout << "Network Size = " << N << " x " << N << std::endl;

  typename NetworkType::RegionType           region;
  typename NetworkType::RegionType::SizeType size;
  size[0] = N;
  size[1] = N;
  region.SetSize(size);
  network->SetRegions(region);
  network->AllocateInitialized();

  unsigned int nVoxels = labels->GetLargestPossibleRegion().GetSize()[0];
  unsigned int nTimes = time->GetLargestPossibleRegion().GetSize()[0];

  if (nVoxels != time->GetLargestPossibleRegion().GetSize()[1])
  {
    std::cout << "number of labels does not match number of voxels" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Examining " << N << " regions, covering " << nVoxels << " voxels with " << nTimes << " time points each"
            << std::endl;

  // unsigned int labelCounts[N];
  auto * labelCounts = new unsigned int[N];

  for (unsigned int i = 0; i < N; i++)
  {
    typename NetworkType::IndexType idx;
    idx[1] = 0;

    labelCounts[i] = 0;
    for (unsigned int v = 0; v < nVoxels; v++)
    {
      idx[0] = v;
      if (itk::Math::FloatAlmostEqual(static_cast<float>(labels->GetPixel(idx)), static_cast<float>(i + 1)))
      {
        ++labelCounts[i];
      }
    }
  }

  // used to rankify matrices if using robust
  typename SCCANType::Pointer cca_rankify = SCCANType::New();

  for (unsigned int i = 0; i < N; i++)
  {
    typename NetworkType::IndexType idx;
    idx[1] = 0;

    MatrixType P(nTimes, labelCounts[i], 0.0);

    unsigned int iCount = 0;
    for (unsigned int v = 0; v < nVoxels; v++)
    {
      idx[0] = v;
      typename NetworkType::IndexType timeIdx;
      timeIdx[1] = v;

      if (itk::Math::FloatAlmostEqual(static_cast<float>(labels->GetPixel(idx)), static_cast<float>(i + 1)))
      {
        for (unsigned int t = 0; t < nTimes; t++)
        {
          timeIdx[0] = t;
          P(t, iCount) = time->GetPixel(timeIdx);
        }
        ++iCount;
      }
    }

    if (robust && (labelCounts[i] >= minRegionSize))
    {
      P = cca_rankify->RankifyMatrixColumns(P);
    }


    if (labelCounts[i] >= minRegionSize)
    {
      for (unsigned int j = i + 1; j < N; j++)
      {
        MatrixType                      Q(nTimes, labelCounts[j], 0.0);
        typename NetworkType::IndexType idx2;
        idx2[1] = 0;

        unsigned int jCount = 0;
        for (unsigned int v2 = 0; v2 < nVoxels; v2++)
        {
          idx2[0] = v2;
          typename NetworkType::IndexType timeIdx2;
          timeIdx2[1] = v2;

          if (itk::Math::FloatAlmostEqual(static_cast<float>(labels->GetPixel(idx2)), static_cast<float>(j + 1)))
          {
            for (unsigned int t2 = 0; t2 < nTimes; t2++)
            {
              timeIdx2[0] = t2;
              Q(t2, jCount) = time->GetPixel(timeIdx2);
            }
            ++jCount;
          }
        }

        if (robust)
        {
          Q = cca_rankify->RankifyMatrixColumns(Q);
        }

        if (labelCounts[j] >= minRegionSize)
        {

          // Correlation magic goes here
          typename SCCANType::Pointer cca = SCCANType::New();
          cca->SetSilent(true);
          cca->SetMaximumNumberOfIterations(iterct);
          cca->SetUseL1(useL1);
          cca->SetGradStep(gradstep);
          cca->SetKeepPositiveP(keepPositive);
          cca->SetKeepPositiveQ(keepPositive);
          cca->SetFractionNonZeroP(sparsity);
          cca->SetFractionNonZeroQ(sparsity);
          cca->SetMinClusterSizeP(minClusterSize);
          cca->SetMinClusterSizeQ(minClusterSize);
          cca->SetMatrixP(P);
          cca->SetMatrixQ(Q);

          // is truecorr just sccancorrs[0]?
          cca->SparsePartialArnoldiCCA(n_evec);
          VectorType sccancorrs = cca->GetCanonicalCorrelations();

          VectorType pVec = cca->GetVariateP();
          for (unsigned int ip = 0; ip < pVec.size(); ip++)
          {
            pVec[ip] = itk::Math::abs(pVec[ip]);
          }
          // pVec = pVec.normalize();
          pVec = P * pVec;

          VectorType qVec = cca->GetVariateQ();
          for (unsigned int iq = 0; iq < qVec.size(); iq++)
          {
            qVec[iq] = itk::Math::abs(qVec[iq]);
          }
          // qVec = qVec.normalize();
          qVec = Q * qVec;

          double final_corr = vnl_pearson_corr(pVec, qVec);
          if (!std::isfinite(final_corr))
          {
            final_corr = 0.0;
          }

          typename NetworkType::IndexType connIdx;
          connIdx[0] = i;
          connIdx[1] = j;
          network->SetPixel(connIdx, final_corr);
          connIdx[0] = j;
          connIdx[1] = i;
          network->SetPixel(connIdx, final_corr);
        }
      }
    }
  }

  delete[] labelCounts;

  return network;
}


template <typename NetworkType>
bool
RegionAveraging(typename NetworkType::Pointer network,
                typename NetworkType::Pointer time,
                typename NetworkType::Pointer labels,
                unsigned int                  nLabels,
                unsigned int                  minSize)
{

  using VectorType = vnl_vector<float>;
  using MatrixType = vnl_matrix<float>;

  // Determine the number of regions to examine
  std::set<unsigned int>                         labelset;
  itk::ImageRegionIteratorWithIndex<NetworkType> it(labels, labels->GetLargestPossibleRegion());

  if (nLabels == 0)
  {

    while (!it.IsAtEnd())
    {
      if (it.Value() > 0)
      {
        labelset.insert(it.Value());
      }
      ++it;
    }
  }
  else
  {
    for (unsigned int i = 0; i < nLabels; i++)
    {
      labelset.insert(i + 1);
    }
  }

  unsigned int N = labelset.size();
  std::cout << "Network Size = " << N << " x " << N << std::endl;

  typename NetworkType::RegionType           region;
  typename NetworkType::RegionType::SizeType size;
  size[0] = N;
  size[1] = N;
  region.SetSize(size);
  network->SetRegions(region);
  network->AllocateInitialized();

  unsigned int nVoxels = labels->GetLargestPossibleRegion().GetSize()[0];
  unsigned int nTimes = time->GetLargestPossibleRegion().GetSize()[0];

  if (nVoxels != time->GetLargestPossibleRegion().GetSize()[1])
  {
    std::cout << "number of labels does not match number of voxels" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Examining " << N << " regions, covering " << nVoxels << " voxels with " << nTimes << " time points each"
            << std::endl;

  VectorType labelCounts(N, 0);
  MatrixType timeSig(N, nTimes, 0.0);
  for (unsigned int i = 0; i < N; i++)
  {
    typename NetworkType::IndexType idx;
    idx[1] = 0;

    for (unsigned int v = 0; v < nVoxels; v++)
    {
      idx[0] = v;
      if (itk::Math::FloatAlmostEqual(static_cast<float>(labels->GetPixel(idx)), static_cast<float>(i + 1)))
      {
        labelCounts[i]++;

        typename NetworkType::IndexType timeIdx;
        timeIdx[1] = v;
        for (unsigned int t = 0; t < nTimes; t++)
        {
          timeIdx[0] = t;
          timeSig(i, t) += time->GetPixel(timeIdx);
        }
      }
    }
  }

  for (unsigned int i = 0; i < N; i++)
  {
    for (unsigned int j = 0; j < nTimes; j++)
    {
      if (labelCounts[i] > 0)
      {
        timeSig(i, j) /= labelCounts[i];
      }
    }
  }

  for (unsigned int i = 0; i < N; i++)
  {
    for (unsigned int j = (i + 1); j < N; j++)
    {

      if ((labelCounts[i] > minSize) && (labelCounts[j] > minSize))
      {
        VectorType p = timeSig.get_row(i);
        VectorType q = timeSig.get_row(j);

        double corr = vnl_pearson_corr(p, q);

        if (!std::isfinite(corr))
        {
          corr = 0.0;
        }

        typename NetworkType::IndexType connIdx;
        connIdx[0] = i;
        connIdx[1] = j;
        network->SetPixel(connIdx, corr);
        connIdx[0] = j;
        connIdx[1] = i;
        network->SetPixel(connIdx, corr);
      }
    }
  }

  return network;
}

int
timesccan(itk::ants::CommandLineParser * parser)
{

  using NetworkType = itk::Image<float, 2>;


  std::string                                       outname = "output.nii.gz";
  itk::ants::CommandLineParser::OptionType::Pointer outputOption = parser->GetOption("output");
  if (!outputOption || outputOption->GetNumberOfFunctions() == 0)
  {
    std::cout << "Warning:  no output option set." << std::endl;
  }
  else
  {
    outname = outputOption->GetFunction()->GetName();
    std::cout << "Writing output to: " << outname << std::endl;
  }

  unsigned int                                      nLabels = 0;
  itk::ants::CommandLineParser::OptionType::Pointer labels_option = parser->GetOption("number-consecutive-labels");
  if (!labels_option || labels_option->GetNumberOfFunctions() == 0)
  {
    //    std::cout << "Warning:  no permutation option set." << std::endl;
  }
  else
  {
    nLabels = parser->Convert<unsigned int>(labels_option->GetFunction()->GetName());
  }

  unsigned int                                      roiSize = 1;
  itk::ants::CommandLineParser::OptionType::Pointer size_option = parser->GetOption("minimum-region-size");
  if (!size_option || size_option->GetNumberOfFunctions() == 0)
  {
    //    std::cout << "Warning:  no permutation option set." << std::endl;
  }
  else
  {
    roiSize = parser->Convert<unsigned int>(size_option->GetFunction()->GetName());
  }

  unsigned int                                      clusterSize = 1;
  itk::ants::CommandLineParser::OptionType::Pointer clust_option = parser->GetOption("minimum-cluster-size");
  if (!clust_option || clust_option->GetNumberOfFunctions() == 0)
  {
    //    std::cout << "Warning:  no permutation option set." << std::endl;
  }
  else
  {
    clusterSize = parser->Convert<unsigned int>(clust_option->GetFunction()->GetName());
  }

  unsigned int                                      evec_ct = 5;
  itk::ants::CommandLineParser::OptionType::Pointer evec_option = parser->GetOption("n_eigenvectors");
  if (!evec_option || evec_option->GetNumberOfFunctions() == 0)
  {
    //    std::cout << "Warning:  no permutation option set." << std::endl;
  }
  else
  {
    evec_ct = parser->Convert<unsigned int>(evec_option->GetFunction()->GetName());
  }

  unsigned int                                      iterations = 20;
  itk::ants::CommandLineParser::OptionType::Pointer iter_option = parser->GetOption("iterations");
  if (!iter_option || iter_option->GetNumberOfFunctions() == 0)
  {
    //    std::cout << "Warning:  no permutation option set." << std::endl;
  }
  else
  {
    iterations = parser->Convert<unsigned int>(iter_option->GetFunction()->GetName());
  }

  float                                             sparsity = 0.1;
  itk::ants::CommandLineParser::OptionType::Pointer sparse_option = parser->GetOption("sparsity");
  if (!sparse_option || sparse_option->GetNumberOfFunctions() == 0)
  {
    //    std::cout << "Warning:  no permutation option set." << std::endl;
  }
  else
  {
    sparsity = parser->Convert<unsigned int>(sparse_option->GetFunction()->GetName());
  }

  unsigned int                                      keepPositive = 1;
  itk::ants::CommandLineParser::OptionType::Pointer pos_option = parser->GetOption("keep-positive");
  if (!pos_option || pos_option->GetNumberOfFunctions() == 0)
  {
    //    std::cout << "Warning:  no permutation option set." << std::endl;
  }
  else
  {
    keepPositive = parser->Convert<unsigned int>(pos_option->GetFunction()->GetName());
  }

  unsigned int                                      usel1 = 1;
  itk::ants::CommandLineParser::OptionType::Pointer l1_option = parser->GetOption("l1");
  if (!l1_option || l1_option->GetNumberOfFunctions() == 0)
  {
    //    std::cout << "Warning:  no permutation option set." << std::endl;
  }
  else
  {
    usel1 = parser->Convert<unsigned int>(l1_option->GetFunction()->GetName());
  }

  unsigned int                                      robustify = 0;
  itk::ants::CommandLineParser::OptionType::Pointer robust_option = parser->GetOption("robustify");
  if (!robust_option || robust_option->GetNumberOfFunctions() == 0)
  {
    //    std::cout << "Warning:  no permutation option set." << std::endl;
  }
  else
  {
    robustify = parser->Convert<unsigned int>(robust_option->GetFunction()->GetName());
  }

  itk::ants::CommandLineParser::OptionType::Pointer evecg_option = parser->GetOption("EvecGradPenalty");
  if (evecg_option && evecg_option->GetNumberOfFunctions() != 0)
  {
    parser->Convert<unsigned int>(evecg_option->GetFunction()->GetName());
  }

  itk::ants::CommandLineParser::OptionType::Pointer eigen_option = parser->GetOption("ridge_cca");
  if (eigen_option && eigen_option->GetNumberOfFunctions() != 0)
  {
    parser->Convert<bool>(eigen_option->GetFunction()->GetName());
  }

  NetworkType::Pointer network = NetworkType::New();

  itk::ants::CommandLineParser::OptionType::Pointer netOption = parser->GetOption("network");
  if (netOption && netOption->GetNumberOfFunctions() > 0)
  {
    std::cout << "Build network" << std::endl;
    if (netOption && netOption->GetFunction(0)->GetNumberOfParameters() < 2)
    {
      std::cout << "  Incorrect number of parameters." << std::endl;
      return EXIT_FAILURE;
    }
    std::string connectivityStrategy = netOption->GetFunction()->GetName();
    std::string timeMatrixName = std::string(netOption->GetFunction(0)->GetParameter(0));
    std::string labelMatrixName = std::string(netOption->GetFunction(0)->GetParameter(1));

    std::cout << "Method: " << connectivityStrategy << std::endl;

    if (connectivityStrategy == "scca")
    {
      std::cout << "Time Series Data: " << timeMatrixName << std::endl;
      std::cout << "Time Series Labels: " << labelMatrixName << std::endl;

      NetworkType::Pointer timeMat = nullptr;
      ReadImage<NetworkType>(timeMat, timeMatrixName.c_str());

      NetworkType::Pointer labelMat = nullptr;
      ReadImage<NetworkType>(labelMat, labelMatrixName.c_str());

      float gradstep = -0.5 + itk::Math::abs(usel1);

      RegionSCCA<NetworkType>(network,
                              timeMat,
                              labelMat,
                              nLabels,
                              roiSize,
                              evec_ct,
                              iterations,
                              sparsity,
                              robustify,
                              usel1,
                              gradstep,
                              keepPositive,
                              clusterSize);
    }
    else if (connectivityStrategy == "region-averaging")
    {
      std::cout << "Time Series Data: " << timeMatrixName << std::endl;
      std::cout << "Time Series Labels: " << labelMatrixName << std::endl;

      NetworkType::Pointer timeMat = nullptr;
      ReadImage<NetworkType>(timeMat, timeMatrixName.c_str());

      NetworkType::Pointer labelMat = nullptr;
      ReadImage<NetworkType>(labelMat, labelMatrixName.c_str());

      RegionAveraging<NetworkType>(network, timeMat, labelMat, nLabels, roiSize);
    }
    else
    {
      std::cout << "Unknown method:" << connectivityStrategy << std::endl;
      return EXIT_FAILURE;
    }
  }

  ANTs::WriteImage<NetworkType>(network, outname.c_str());

  return 0;
}

void
InitializeCommandLineOptions(itk::ants::CommandLineParser * parser)
{
  /** in this function, list all the operations you will perform */

  using OptionType = itk::ants::CommandLineParser::OptionType;

  {
    std::string         description = std::string("Print the help menu (short version).");
    OptionType::Pointer option = OptionType::New();
    option->SetShortName('h');
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("Print the help menu (long version).");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("help");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("Output is a 2D correlation matrix.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("output");
    option->SetShortName('o');
    option->SetUsageOption(0, "outputImage");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("Number of consecutive labels in data");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("number-consecutive-labels");
    option->SetShortName('l');
    option->SetUsageOption(0, "0");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("Minimum size of a region: regions below this size are given a 0.0 connectivity value");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("minimum-region-size");
    option->SetShortName('R');
    option->SetUsageOption(0, "1");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("Number of iterations");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("iterations");
    option->SetShortName('i');
    option->SetUsageOption(0, "20");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Sparsity - a float from (0,1] indicating what fraction of the data to use");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("sparsity");
    option->SetShortName('s');
    option->SetUsageOption(0, "0.10");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("Number of permutations to use in scca.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("n_eigenvectors");
    option->SetShortName('n');
    option->SetUsageOption(0, "2");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("rank-based scca");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("robustify");
    option->SetShortName('r');
    option->SetUsageOption(0, "0");
    option->SetDescription(description);
    parser->AddOption(option);
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
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("cluster threshold on view P");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("ClusterThresh");
    option->SetUsageOption(0, "1");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("Number of permutations to use in scca.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("ridge_cca");
    option->SetShortName('e');
    option->SetUsageOption(0, "0");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("Choices for pscca: PQ, PminusRQ, PQminusR, PminusRQminusR ");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("partial-scca-option");
    option->SetUsageOption(0, "PminusRQ");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("takes a timeseries (4D) image ") +
                              std::string("and converts it to a 2D matrix csv format as output.") +
                              std::string("If the mask has multiple labels ( more the one ) then the average time "
                                          "series in each label will be computed and put in the csv.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("timeseriesimage-to-matrix");
    option->SetUsageOption(0, "[four_d_image.nii.gz,three_d_mask.nii.gz]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("takes a labeled (3D) image ") + std::string("and converts it to a 2D matrix csv format as output.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("labelsimage-to-matrix");
    option->SetUsageOption(0, "[three_d_mask.nii.gz]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("Build the network connectivity matrix");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("network");
    option->SetUsageOption(0, "scca[time-matrix.mhd,label-matrix.mhd]");
    option->SetUsageOption(1, "region-averaging[time-matrix.mhd,label-matrix.mhd]");
    option->SetDescription(description);
    parser->AddOption(option);
  }
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
TimeSCCAN(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "TimeSCCAN");

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

  itk::ants::CommandLineParser::Pointer parser = itk::ants::CommandLineParser::New();

  parser->SetCommand(argv[0]);

  std::string commandDescription =
    std::string("A tool for sparse statistical analysis on connectivity within a subject : ");

  parser->SetCommandDescription(commandDescription);
  InitializeCommandLineOptions(parser);

  if (parser->Parse(argc, argv) == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  }

  // Print the entire help menu
  itk::ants::CommandLineParser::OptionType::Pointer shortHelpOption = parser->GetOption('h');
  itk::ants::CommandLineParser::OptionType::Pointer longHelpOption = parser->GetOption("help");
  if (argc < 2 ||
      (shortHelpOption->GetFunction() && parser->Convert<unsigned int>(shortHelpOption->GetFunction()->GetName()) == 1))
  {
    parser->PrintMenu(std::cout, 5, true);
    if (argc < 2)
    {
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }
  if (longHelpOption->GetFunction() && parser->Convert<unsigned int>(longHelpOption->GetFunction()->GetName()) == 1)
  {
    parser->PrintMenu(std::cout, 5, false);
    return EXIT_SUCCESS;
  }

  // Print the long help menu for specific items
  if (longHelpOption && longHelpOption->GetNumberOfFunctions() > 0 &&
      parser->Convert<unsigned int>(longHelpOption->GetFunction()->GetName()) != 0)
  {
    itk::ants::CommandLineParser::OptionListType options = parser->GetOptions();
    for (unsigned int n = 0; n < longHelpOption->GetNumberOfFunctions(); n++)
    {
      const std::string & value = longHelpOption->GetFunction(n)->GetName();
      for (auto it = options.cbegin(); it != options.cend(); ++it)
      {
        const std::string & longname = (*it)->GetLongName();
        if (longname.rfind(value, 0) == 0) // determining if `longname` starts with `value`
        {
          parser->PrintMenu(std::cout, 5, false);
        }
      }
    }
    return EXIT_FAILURE;
  }

  // Call main routine
  return timesccan(parser);
}

} // namespace ants
