/*===================================================================

  Program:   ASHS (Automatic Segmentation of Hippocampal Subfields)
  Module:    $Id: WeightedVotingLabelFusionImageFilter.txx,v 1.1 2012/09/05 14:39:51 zhuzhu771 Exp $
  Language:  C++ program
  Copyright (c) 2012 Paul A. Yushkevich, University of Pennsylvania

  This file is part of ASHS

  ASHS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  ===================================================================

  CITATION:

    This program implements the method described in the paper

    H. Wang, J. W. Suh, S. Das, J. Pluta, C. Craige, P. Yushkevich,
    "Multi-atlas segmentation with joint label fusion," IEEE Trans.
    on Pattern Analysis and Machine Intelligence, 35(3), 611-623, 2013

  =================================================================== */

#include <itkNeighborhoodIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_svd.h>

#include <set>
#include <map>
#include <vector>

template <class TInputImage, class TOutputImage>
void
WeightedVotingLabelFusionImageFilter<TInputImage, TOutputImage>
::UpdateInputs()
{
  // Set all the inputs
  this->SetNumberOfIndexedInputs(1 + 2 * m_Atlases.size() + m_Exclusions.size());

  size_t kInput = 0;
  this->SetNthInput(kInput++, m_Target);
  for(size_t i = 0; i < m_Atlases.size(); i++)
    {
    this->SetNthInput(kInput++, m_Atlases[i]);
    this->SetNthInput(kInput++, m_AtlasSegs[i]);
    }

  for(typename ExclusionMap::iterator it = m_Exclusions.begin(); it != m_Exclusions.end(); ++it)
    {
    this->SetNthInput(kInput++, it->second);
    }
}




template <class TInputImage, class TOutputImage>
void
WeightedVotingLabelFusionImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();

  // Get the output requested region
  RegionType outRegion = this->GetOutput()->GetRequestedRegion();

  // Pad this region by the search window and patch size
  outRegion.PadByRadius(m_SearchRadius);
  outRegion.PadByRadius(m_PatchRadius);

  // Iterate over all the inputs to this filter
  for(size_t i = 0; i < this->GetNumberOfInputs(); i++)
    {
    // Get i-th input
    InputImageType *input = const_cast<InputImageType *>(this->GetInput(i));
    RegionType region = outRegion;
    region.Crop(input->GetLargestPossibleRegion());
    input->SetRequestedRegion(region);
    }
}

template<class TInputImage, class TOutputImage>
void
WeightedVotingLabelFusionImageFilter<TInputImage, TOutputImage>
::ComputeOffsetTable(
  const InputImageType *image,
  const SizeType &radius,
  int **offset,
  size_t &nPatch,
  int **manhattan)
{
  // Use iterators to construct offset tables
  RegionType r = image->GetBufferedRegion();
  NIter itTempPatch(radius, image, r);

  // Position the iterator in the middle to avoid problems with boundary conditions
  IndexType iCenter;
  for(size_t i = 0; i < InputImageDimension; i++)
    iCenter[i] = r.GetIndex(i) + r.GetSize(i)/2;
  itTempPatch.SetLocation(iCenter);

  // Compute the offsets
  nPatch = itTempPatch.Size();
  (*offset) = new int[nPatch];
  if(manhattan)
    (*manhattan) = new int[nPatch];
  for(size_t i = 0; i < nPatch; i++)
  {
    (*offset)[i] = itTempPatch[i] - itTempPatch.GetCenterPointer();
    if(manhattan)
      {
      typename NIter::OffsetType off = itTempPatch.GetOffset(i);
      (*manhattan)[i] = 0;
      for(int d = 0; d < InputImageDimension; d++)
        (*manhattan)[i] += abs(off[d]);
      }
  }
}



template <class TInputImage, class TOutputImage>
void
WeightedVotingLabelFusionImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  // Allocate the output
  this->GetOutput()->SetBufferedRegion( this->GetOutput()->GetRequestedRegion() );
  this->GetOutput()->Allocate();

  // Get the target image
  InputImageType *target = m_Target;

  // Create a neighborhood iterator for the target image
  NIter itTarget(m_PatchRadius, target, this->GetOutput()->GetRequestedRegion());

  // Get the number of atlases
  int n = m_Atlases.size();

  // Construct offset tables for all the images (these can be different because they
  // depend on the buffered region)
  size_t nPatch, nSearch;
  int *offPatchTarget,
      **offPatchAtlas = new int *[n],
      **offPatchSeg = new int *[n],
      **offSearchAtlas = new int *[n],
      *manhattan;

  // Collect search statistics
  std::vector<int> searchHisto(100, 0);

  // Compute the offset table for the target image
  ComputeOffsetTable(target, m_PatchRadius, &offPatchTarget, nPatch);

  // Find all unique labels in the requested region
  std::set<InputImagePixelType> labelset;
  for(int i = 0; i < n; i++)
    {
    // Find all the labels
    const InputImageType *seg = m_AtlasSegs[i];
    itk::ImageRegionConstIteratorWithIndex<InputImageType> it(seg, seg->GetRequestedRegion());
    for(; !it.IsAtEnd(); ++it)
      {
      InputImagePixelType label = it.Get();
      labelset.insert(label);
      }

    // Compute the offset table for that atlas
    ComputeOffsetTable(m_Atlases[i], m_PatchRadius, offPatchAtlas+i, nPatch);
    ComputeOffsetTable(m_AtlasSegs[i], m_PatchRadius, offPatchSeg+i, nPatch);
    ComputeOffsetTable(m_Atlases[i], m_SearchRadius, offSearchAtlas+i, nSearch, &manhattan);
    }

  // Initialize the posterior maps
  m_PosteriorMap.clear();

  // Allocate posterior images for the different labels
  for(typename std::set<InputImagePixelType>::iterator sit = labelset.begin();
    sit != labelset.end(); ++sit)
    {
    m_PosteriorMap[*sit] = PosteriorImage::New();
    m_PosteriorMap[*sit]->SetLargestPossibleRegion(target->GetLargestPossibleRegion());
    m_PosteriorMap[*sit]->SetRequestedRegion(this->GetOutput()->GetRequestedRegion());
    m_PosteriorMap[*sit]->SetBufferedRegion(this->GetOutput()->GetRequestedRegion());
    m_PosteriorMap[*sit]->Allocate();
    m_PosteriorMap[*sit]->FillBuffer(0.0f);
    }

  int iter = 0;


  // We need an array of absolute patch differences between target image and atlases
  // (apd - atlas patch difference)
  InputImagePixelType **apd = new InputImagePixelType*[n];
  for(int i = 0; i < n; i++)
    apd[i] = new InputImagePixelType[nPatch];

  // Also an array of pointers to the segmentations of different atlases
  const InputImagePixelType **patchSeg = new const InputImagePixelType*[n];

  // Create an array for storing the normalized target patch to save more time
  InputImagePixelType *xNormTargetPatch = new InputImagePixelType[nPatch];

  // Iterate over voxels in the output region
  typedef itk::ImageRegionIteratorWithIndex<TOutputImage> OutIter;
  for(OutIter it(this->GetOutput(), this->GetOutput()->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    // Point the target iterator to the output location
    itTarget.SetLocation(it.GetIndex());
    InputImagePixelType *pTargetCurrent = target->GetBufferPointer() + target->ComputeOffset(it.GetIndex());

    // Compute stats for the target patch
    InputImagePixelType mu, sigma;
    PatchStats(pTargetCurrent, nPatch, offPatchTarget, mu, sigma);
    for(unsigned int i = 0; i < nPatch; i++)
      xNormTargetPatch[i] = (*(pTargetCurrent + offPatchTarget[i]) - mu) / sigma;

    // In each atlas, search for a patch that matches our patch
    for(int i = 0; i < n; i++)
      {
      const InputImageType *atlas = m_Atlases[i];
      const InputImageType *seg = m_AtlasSegs[i];
      int *offPatch = offPatchAtlas[i], *offSearch = offSearchAtlas[i];

      // Get the requested region for the atlas
      RegionType rr = atlas->GetRequestedRegion();

      // Define the search region
      RegionType rSearch;
      for(int j = 0; j < InputImageDimension; j++)
        {
        rr.SetIndex(j, rr.GetIndex(j) + m_SearchRadius[j]);
        rr.SetSize(j, rr.GetSize(j) - 2 * m_SearchRadius[j]);
        rSearch.SetIndex(j, it.GetIndex()[j] - m_SearchRadius[j]);
        rSearch.SetSize(j, 2 * m_SearchRadius[j]+1);
        }
      rSearch.Crop(rr);

      // Search over neighborhood
      const InputImagePixelType *pAtlasCurrent = atlas->GetBufferPointer() + atlas->ComputeOffset(it.GetIndex());
      double bestMatch = 1e100;
      const InputImagePixelType *bestMatchPtr = NULL;
      InputImagePixelType bestMatchSum = 0, bestMatchSSQ = 0;
      int bestK = 0;
      for(unsigned int k = 0; k < nSearch; k++)
        {
        // Pointer to the voxel at the center of the search
        const InputImagePixelType *pSearchCenter = pAtlasCurrent + offSearch[k];
        InputImagePixelType matchSum = 0, matchSSQ = 0;
        double match = this->PatchSimilarity(pSearchCenter, xNormTargetPatch, nPatch, offPatch,
                                             matchSum, matchSSQ);
        /*
        if(it.GetIndex()[0] == 326 && it.GetIndex()[1] == 244 && it.GetIndex()[2] == 12 && i == 9)
          {
          NIter itt(m_SearchRadius, atlas, atlas->GetBufferedRegion());
          itt.SetLocation(it.GetIndex());
          if(itt.GetIndex(k)[0] == 323 && itt.GetIndex(k)[1] == 247 && itt.GetIndex(k)[2]==12)
            {
            std::cout << "BAD PATCH" << std::endl;
            for(int z = 0; z < nPatch; z++)
              std::cout << *(pSearchCenter + offPatch[z]) <<
                "\t" << xNormTargetPatch[z] << std::endl;
            printf("match = %f, sum = %f, ssq = %f\n", match, matchSum, matchSSQ);
            }
          }
          */
        if(match < bestMatch)
          {
          bestMatch = match;
          bestMatchPtr = pSearchCenter;
          bestMatchSum = matchSum;
          bestMatchSSQ = matchSSQ;
          bestK = k;
          }
        }

      // Update the manhattan distance histogram
      searchHisto[manhattan[bestK]]++;

      /*
      if(it.GetIndex()[0] == 326 && it.GetIndex()[1] == 244 && it.GetIndex()[2] == 12)
        {
        NIter itt(m_SearchRadius, atlas, atlas->GetBufferedRegion());
        itt.SetLocation(it.GetIndex());
        std::cout << "Best match " << i << ":" <<
          itt.GetOffset(bestK)[0] + 326 << ", " <<
          itt.GetOffset(bestK)[1] + 244 << ", " <<
          itt.GetOffset(bestK)[2] + 12 << "; " <<
          bestMatch << std::endl;

        double bM = 1e100; int bxs, bys, bzs;
        for(int xs = -3; xs <= 3; xs++)
          {
          for(int ys = -3; ys <= 3; ys++)
            {
            for(int zs = -1; zs <= 1; zs++)
              {
              vnl_vector<double> Xt(147), Xa(147);
              int g = 0;
              for(int xp = -3; xp <= 3; xp++)
                {
                for(int yp = -3; yp <= 3; yp++)
                  {
                  for(int zp = -1; zp <= 1; zp++)
                    {
                    itk::Index<InputImageDimension> ia, it;
                    it[0] = 326 + xp; it[1] = 244 + yp; it[2] = 12 + zp;
                    ia[0] = 326 + xp + xs; ia[1] = 244 + yp + ys; ia[2] = 12 + zp + zs;
                    Xt[g] = target->GetPixel(it);
                    Xa[g] = atlas->GetPixel(ia);
                    g++;
                    }
                  }
                }
              if(xs == -3 && ys == 3 && zs == 0 && i==9)
                {
                std::cout << "BAD PATCH" << std::endl;
                for(int z = 0; z < nPatch; z++)
                  std::cout << Xa[z] << std::endl;
                }
              Xt = (Xt - Xt.mean()).normalize();
              Xa = (Xa - Xa.mean()).normalize();
              double match = (Xt - Xa).magnitude();
              // if(i == 9)
              //   printf("{%d, %d, %d} => %f\n", 326+xs, 244+ys, 12+zs, match);
              if(match < bM)
                {
                bM = match;
                bxs=xs; bys=ys; bzs=zs;
                }
              }
            }
          }
        std::cout << " --- " << 326+bxs << " " << 244+bys << 12 + bzs << std::endl;
        }
        */

      // Once the patch has been found, compute the absolute difference with target image
      InputImagePixelType bestMatchMean = bestMatchSum / nPatch;
      InputImagePixelType bestMatchVar =
        (bestMatchSSQ - nPatch * bestMatchMean * bestMatchMean) / (nPatch - 1);
      if(bestMatchVar < 1.0e-12)
        bestMatchVar = 1.0e-12;
      InputImagePixelType bestMatchSD = sqrt(bestMatchVar);

      for(unsigned int m = 0; m < nPatch; m++)
        {
        InputImagePixelType x = *(bestMatchPtr + offPatch[m]);
        apd[i][m] = fabs(xNormTargetPatch[m] - (x - bestMatchMean) / bestMatchSD);
        }

      // Store the best found neighborhood
      patchSeg[i] = (bestMatchPtr - atlas->GetBufferPointer()) + seg->GetBufferPointer();
      }

    // Allocate Mx
    typedef vnl_matrix<double> MatrixType;
    MatrixType Mx(n, n);

    // Now we can compute Mx
    for(int i = 0; i < n; i++)
      {
      for(int k = 0; k <= i; k++)
        {
        // Multiply through the apd arrays
        //InputImagePixelType mxval = 0.0;
        double mxval = 0.0;
        for(unsigned int m = 0; m < nPatch; m++)
          mxval += apd[i][m] * apd[k][m];

        mxval /= (nPatch - 1);

        if(m_Beta == 2)
          mxval *= mxval;
        else
          mxval = pow(mxval, m_Beta);

        Mx(i,k) = Mx(k,i) = (InputImagePixelType)mxval;
        }
      }

    /*
    if(it.GetIndex()[0] == 326 && it.GetIndex()[1] == 244 && it.GetIndex()[2] == 12)
      {
      std::cout << Mx << std::endl;

      // Check the computation
      const InputImageType *atlas1 = this->GetInput(1 + 2 * 9);
      const InputImageType *atlas2 = this->GetInput(1 + 2 * 9);
      vnl_vector<double> xt(nPatch), x1(nPatch), x2(nPatch);
      int g=0;
      for(int xp = -3; xp <= 3; xp++)
        {
        for(int yp = -3; yp <= 3; yp++)
          {
          for(int zp = -1; zp <= 1; zp++)
            {
            itk::Index<InputImageDimension> i1, i2, it;
            i1[0] = 324 + xp; i1[1] = 245 + yp; i1[2] = 11 + zp;
            i2[0] = 325 + xp; i2[1] = 244 + yp; i2[2] = 12 + zp;
            it[0] = 326 + xp; it[1] = 244 + yp; it[2] = 12 + zp;
            xt[g] = target->GetPixel(it);
            x1[g] = atlas1->GetPixel(i1);
            x2[g] = atlas2->GetPixel(i2);
            g++;
            }
          }
        }

      // Compute the entry
      xt = (xt - xt.mean()).normalize();
      x1 = (x1 - x1.mean()).normalize();
      x2 = (x2 - x2.mean()).normalize();

      double mxval = dot_product((x1 - xt).apply(fabs), (x2 - xt).apply(fabs));
      mxval = pow(mxval, 2);
      std::cout << "Mx = " << mxval << " vs " << Mx(9,9) << std::endl;
      }
      */

    // Now we can compute the weights by solving for the inverse of Mx
    MatrixType Mx_bar(n, n, 0.0);
    Mx_bar.fill_diagonal(m_Alpha);
    Mx_bar += Mx;

    // Define a vector of all ones
    vnl_vector<double> ones(n, 1.0);

    // Solve for the weights
    vnl_vector<double> W = vnl_svd<double>(Mx_bar).solve(ones);

    // Normalize the weights
    W *= 1.0 / dot_product(W, ones);

    /*
    if(it.GetIndex()[0] == 326 && it.GetIndex()[1] == 244 && it.GetIndex()[2] == 12)
      {
      std::cout << W << std::endl;
      }
      */

    // Perform voting using Hongzhi's averaging scheme. Iterate over all segmentation patches
    for(unsigned int ni = 0; ni < nPatch; ni++)
      {
      IndexType idx = itTarget.GetIndex(ni);
      if(this->GetOutput()->GetRequestedRegion().IsInside(idx))
        {
        for(int i = 0; i < n; i++)
          {
          // The segmentation at the corresponding patch location in atlas i
          InputImagePixelType label = *(patchSeg[i] + offPatchSeg[i][ni]);

          // Add that weight the posterior map for voxel at idx
          m_PosteriorMap[label]->SetPixel(idx, m_PosteriorMap[label]->GetPixel(idx) + W[i]);
          }
        }
      }

      if(++iter % 1000 == 0)
        std::cout << "." << std::flush;
    }

  std::cout << std::endl << "Search Manhattan Distance Histogram " << std::endl;
  for(size_t i = 0; i < searchHisto.size() && searchHisto[i] > 0; i++)
    std::cout << "    " << i << "\t" << searchHisto[i] << std::endl;

  std::cout << std::endl << "VOTING " << std::endl;

  // Perform voting at each voxel
  for(OutIter it(this->GetOutput(), this->GetOutput()->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    double wmax = 0;
    InputImagePixelType winner = 0;

    for(typename std::set<InputImagePixelType>::iterator sit = labelset.begin();
      sit != labelset.end(); ++sit)
      {
      double posterior = m_PosteriorMap[*sit]->GetPixel(it.GetIndex());

      // check if the label is excluded
      typename ExclusionMap::iterator xit = m_Exclusions.find(*sit);
      bool excluded = (xit != m_Exclusions.end() && xit->second->GetPixel(it.GetIndex()) != 0);

      // Vote!
      if (wmax < posterior && !excluded)
        {
        wmax = posterior;
        winner = *sit;
        }
      }

    it.Set(winner);
    }

  // Clear posterior maps
  if(!m_RetainPosteriorMaps)
    m_PosteriorMap.clear();
}


template <class TInputImage, class TOutputImage>
void
WeightedVotingLabelFusionImageFilter<TInputImage, TOutputImage>
::PatchStats(const InputImagePixelType *p, size_t n, int *offsets, InputImagePixelType &mean, InputImagePixelType &std)
{
  InputImagePixelType sum = 0, ssq = 0;
  for(unsigned int i = 0; i < n; i++)
    {
    InputImagePixelType v = *(p + offsets[i]);
    sum += v;
    ssq += v * v;
    }

  mean = sum / n;
  std = (ssq - n * mean * mean) / (n - 1);
  if(std < 1e-6)
    std = 1e-6;
  std = sqrt(std);
}

/**
 * This function computes similarity between a normalized patch (normtrg) and a patch
 * that has not been normalized (psearch). It can be shown that the sum of squared
 * differences between a normalized patch u and a unnormalized patch v is equal to
 *
 * 2 [ (n-1) - (\Sum u_i v_i ) / \sigma_v ]
 *
 * Since we are only interested in finding the patch with the smallest SSD, we can simplify
 * this expression further to minimizing - [ (\Sum u_i v_i ) / \sigma_v ] ^ 2. To further
 * simpify computation, we return
 *
 *        - (\Sum u_i v_i)^2 / z,   where z = sigma_v^2 * (n-1)
 */
template <class TInputImage, class TOutputImage>
double
WeightedVotingLabelFusionImageFilter<TInputImage, TOutputImage>
::PatchSimilarity(
  const InputImagePixelType *psearch,
  const InputImagePixelType *normtrg,
  size_t n,
  int *offsets,
  InputImagePixelType &sum_psearch,
  InputImagePixelType &ssq_psearch)
{
  // Here the patch normtrg should already be normalized.
  // We simultaneously compute the patch stats and solve the problem
  InputImagePixelType sum_uv = 0;
  for(unsigned int i = 0; i < n; i++)
    {
    InputImagePixelType u = *(psearch + offsets[i]);
    InputImagePixelType v = normtrg[i];
    sum_psearch += u;
    ssq_psearch += u * u;
    sum_uv += u * v;
    }

  InputImagePixelType var_u_unnorm = ssq_psearch - sum_psearch * sum_psearch / n;
  if(var_u_unnorm < 1.0e-6)
    var_u_unnorm = 1.0e-6;

  if(sum_uv > 0)
    return - (sum_uv * sum_uv) / var_u_unnorm;
  else
    return (sum_uv * sum_uv) / var_u_unnorm;

  // InputImagePixelType sd_u = std::max(1e-6, sqrt((ssq_u - sum_u * sum_u / n) / (n - 1)));
  // return 2 * ((n - 1) - sum_uv / sd_u);
}

template <class TInputImage, class TOutputImage>
double
WeightedVotingLabelFusionImageFilter<TInputImage, TOutputImage>
::JointErrorEstimate(const InputImagePixelType *t, const InputImagePixelType *a1, const InputImagePixelType *a2, size_t n, int *offsets)
{
  InputImagePixelType mu_t, sigma_t, mu1, sigma1, mu2, sigma2;
  PatchStats(t, n, offsets, mu_t, sigma_t);
  PatchStats(a1, n, offsets, mu1, sigma1);
  PatchStats(a2, n, offsets, mu2, sigma2);

  // What should we return when patches have zero variance?
  if(sigma1 == 0 || sigma2 == 0 || sigma_t == 0)
    return 0;

  double Mxval = 0.0;
  for(unsigned int i = 0; i < n; i++)
    {
    int off = offsets[i];
    InputImagePixelType ft = (*(t + off) - mu_t) / sigma_t;
    InputImagePixelType f1 = (*(a1 + off) - mu1) / sigma1;
    InputImagePixelType f2 = (*(a2 + off) - mu2) / sigma2;

    Mxval += fabs(ft-f1) * fabs(ft-f2);
    }

  return pow(Mxval, -m_Beta);
}

