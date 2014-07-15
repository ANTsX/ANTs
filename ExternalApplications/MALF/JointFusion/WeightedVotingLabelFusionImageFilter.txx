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
using namespace std;

template <class TInputImage, class TOutputImage>
void
WeightedVotingLabelFusionImageFilter<TInputImage, TOutputImage>
::UpdateInputs()
{
  // Set all the inputs
  this->SetNumberOfIndexedInputs(m_Target.size() + 1 + (m_Target.size()+2) * m_Atlases.size() + m_Exclusions.size());
//  this->SetNumberOfInputs(1 + 2 * m_Atlases.size() + m_Exclusions.size());

  size_t kInput = 0;
  for(size_t i = 0; i < m_Target.size(); i++)
    {
    this->SetNthInput(kInput++, m_Target[i]);
    }
//  this->SetNthInput(kInput++, m_Target);

  for(size_t i = 0; i < m_Atlases.size(); i++)
    {
    this->SetNthInput(kInput++, m_Atlases[i]);
    }

  for(size_t i = 0; i < m_AtlasSegs.size(); i++)
    {
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
  InputImageList target = m_Target;

  // Get the number of atlases
  int n = m_AtlasSegs.size();

  typedef vector<NIter> NItList;
  NItList itTargetList;
  for (int j=0;j<m_Modality;j++)
  {
    NIter aitPatch(m_PatchRadius, target[j], target[j]->GetRequestedRegion());
    itTargetList.push_back(aitPatch);
  }


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
  ComputeOffsetTable(target[0], m_PatchRadius, &offPatchTarget, nPatch);

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
    ComputeOffsetTable(m_Atlases[i*m_Modality], m_PatchRadius, offPatchAtlas+i, nPatch);
    ComputeOffsetTable(m_AtlasSegs[i], m_PatchRadius, offPatchSeg+i, nPatch);
    ComputeOffsetTable(m_Atlases[i*m_Modality], m_SearchRadius, offSearchAtlas+i, nSearch, &manhattan);
    }

  // Initialize the posterior maps
  m_PosteriorMap.clear();

  // Allocate posterior images for the different labels
  for(typename std::set<InputImagePixelType>::iterator sit = labelset.begin();
    sit != labelset.end(); ++sit)
    {
    m_PosteriorMap[*sit] = PosteriorImage::New();
    m_PosteriorMap[*sit]->SetLargestPossibleRegion(target[0]->GetLargestPossibleRegion());
    m_PosteriorMap[*sit]->SetRequestedRegion(this->GetOutput()->GetRequestedRegion());
    m_PosteriorMap[*sit]->SetBufferedRegion(this->GetOutput()->GetRequestedRegion());
    m_PosteriorMap[*sit]->Allocate();
    m_PosteriorMap[*sit]->FillBuffer(0.0f);
    }

  std::cout<<m_RetainVotingWeight<<endl;
  if (m_RetainVotingWeight)
  {
    m_VotingWeight.clear();
    for (int i=0;i<n;i++)
    {
      m_VotingWeight[i] = PosteriorImage::New();
      m_VotingWeight[i]->SetLargestPossibleRegion(target[0]->GetLargestPossibleRegion());
      m_VotingWeight[i]->SetRequestedRegion(this->GetOutput()->GetRequestedRegion());
      m_VotingWeight[i]->SetBufferedRegion(this->GetOutput()->GetRequestedRegion());
      m_VotingWeight[i]->Allocate();
      m_VotingWeight[i]->FillBuffer(0.0f);
    }
  }
  PosteriorImagePtr countermap = PosteriorImage::New();
  countermap->SetLargestPossibleRegion(target[0]->GetLargestPossibleRegion());
  countermap->SetRequestedRegion(this->GetOutput()->GetRequestedRegion());
  countermap->SetBufferedRegion(this->GetOutput()->GetRequestedRegion());
  countermap->Allocate();
  countermap->FillBuffer(0.0f);

  typedef vnl_matrix<double> MatrixType;
  MatrixType GroupMx(n, m_GroupWeight.size(), 0.0);
  for (size_t i=0;i<m_GroupWeight.size();i++)
  {
    for (int j=0;j<n;j++)
      if (m_GroupID[j]==(int)i)
        GroupMx[j][i]=1;
  }

  MatrixType GroupWeightMx(m_GroupWeight.size(), 1, 0.0);
  for (size_t i=0;i<m_GroupWeight.size();i++)
    GroupWeightMx[i][0]=m_GroupWeight[i];

  int iter = 0;

  // We need an array of absolute patch differences between target image and atlases
  // (apd - atlas patch difference)
  InputImagePixelType **apd = new InputImagePixelType*[n];
  for(int i = 0; i < n; i++)
    apd[i] = new InputImagePixelType[nPatch * m_Modality];

  // Also an array of pointers to the segmentations of different atlases
  const InputImagePixelType **patchSeg = new const InputImagePixelType*[n];

  // Create an array for storing the normalized target patch to save more time
  InputImagePixelType *xNormTargetPatch = new InputImagePixelType[nPatch * m_Modality];

  // Iterate over voxels in the output region
  typedef itk::ImageRegionIteratorWithIndex<TOutputImage> OutIter;
  for(OutIter it(this->GetOutput(), this->GetOutput()->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    // Point the target iterator to the output location
    for (int i1=0;i1<m_Modality;i1++)
      itTargetList[i1].SetLocation(it.GetIndex());

    int flag=0;
    for(int i = 0; i < n; i++)
    {
      const InputImageType *seg = m_AtlasSegs[i];
      if (seg->GetPixel(it.GetIndex()))
        flag++;
    }
    if (flag==0)
    {
      continue;
    }

    // Compute stats for the target patch
    InputImagePixelType mu, sigma;
    for (int i1=0;i1<m_Modality;i1++)
    {
      int toff=i1 * nPatch;
      InputImagePixelType *pTargetCurrent = target[i1]->GetBufferPointer() + target[i1]->ComputeOffset(it.GetIndex());
      PatchStats(pTargetCurrent, nPatch, offPatchTarget, mu, sigma);
      for(unsigned int i = 0; i < nPatch; i++)
        xNormTargetPatch[toff + i] = (*(pTargetCurrent + offPatchTarget[i]) - mu) / sigma;
    }

    // In each atlas, search for a patch that matches our patch
    for(int i = 0; i < n; i++)
      {
      const InputImageType *seg = m_AtlasSegs[i];
      int toff=i*m_Modality;
      int *offPatch = offPatchAtlas[i], *offSearch = offSearchAtlas[i];

      const InputImageType *tatlas = m_Atlases[toff];

      // Get the requested region for the tatlas
      RegionType rr = tatlas->GetRequestedRegion();

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
      double bestMatch = 1e100;
      const InputImagePixelType **bestMatchPtr = new const InputImagePixelType *[m_Modality];
      const InputImagePixelType **tMatchPtr = new const InputImagePixelType *[m_Modality];

      InputImagePixelType *bestMatchSum =  new InputImagePixelType[m_Modality];
      InputImagePixelType *bestMatchSSQ =  new InputImagePixelType[m_Modality];
      InputImagePixelType *bestMatchMean = new InputImagePixelType[m_Modality];
      InputImagePixelType *bestMatchVar =  new InputImagePixelType[m_Modality];
      InputImagePixelType *bestMatchSD =   new InputImagePixelType[m_Modality];
      InputImagePixelType *MatchSum =      new InputImagePixelType[m_Modality];
      InputImagePixelType *MatchSSQ =      new InputImagePixelType[m_Modality];

      int bestK = 0;
      for(unsigned int k = 0; k < nSearch; k++)
        {
        double tmatch=0;
        for (int i1=0; i1<m_Modality; i1++)
        {
          int toff1= i1 * nPatch;
          const InputImageType *atlas = m_Atlases[toff+i1];
          const InputImagePixelType *pAtlasCurrent = atlas->GetBufferPointer() + atlas->ComputeOffset(it.GetIndex());

          // Pointer to the voxel at the center of the search
          const InputImagePixelType *pSearchCenter = pAtlasCurrent + offSearch[k];
          InputImagePixelType matchSum = 0, matchSSQ = 0;
          double match = this->PatchSimilarity(pSearchCenter, xNormTargetPatch, toff1, nPatch, offPatch,
                                             matchSum, matchSSQ);
          tMatchPtr[i1] = pSearchCenter;
          MatchSum[i1] = matchSum;
          MatchSSQ[i1] = matchSSQ;
          tmatch += match;
        }
        if(tmatch < bestMatch)
          {
          bestMatch = tmatch;
          for (int j1=0;j1<m_Modality;j1++)
          {
            bestMatchSum[j1] = MatchSum[j1];
            bestMatchSSQ[j1] = MatchSSQ[j1];
            bestMatchPtr[j1] = tMatchPtr[j1];
          }
          bestK = k;
          }
        }

      // Update the manhattan distance histogram
      searchHisto[manhattan[bestK]]++;

      // Once the patch has been found, compute the absolute difference with target image
      for (int j1=0;j1<m_Modality;j1++)
      {
        bestMatchMean[j1] = bestMatchSum[j1] / nPatch;
        bestMatchVar[j1]  = (bestMatchSSQ[j1] - nPatch * bestMatchMean[j1] * bestMatchMean[j1]) / (nPatch - 1);
        if(bestMatchVar[j1] < 1.0e-12)
          bestMatchVar[j1] = 1.0e-12;
        bestMatchSD[j1] = sqrt(bestMatchVar[j1]);
      }

      for (int i1=0;i1<m_Modality;i1++)
      {
        int loff = i1*nPatch;
        for(unsigned int m = 0; m < nPatch; m++)
        {
        InputImagePixelType x = *(bestMatchPtr[i1] + offPatch[m]);
        apd[i][loff+m] = fabs(xNormTargetPatch[loff+m] - (x - bestMatchMean[i1]) / bestMatchSD[i1]);
        }
      }

      // Store the best found neighborhood
      patchSeg[i] = (bestMatchPtr[0] - m_Atlases[toff]->GetBufferPointer()) + seg->GetBufferPointer();
      delete [] bestMatchSum;
      delete [] bestMatchSSQ;
      delete [] bestMatchMean;
      delete [] bestMatchVar;
      delete [] bestMatchSD;
      delete [] MatchSum;
      delete [] MatchSSQ;
      }

    // Allocate Mx
    MatrixType Mx(n, n);

    // Now we can compute Mx
    for(int i = 0; i < n; i++)
      {
      for(int k = 0; k <= i; k++)
        {
        // Multiply through the apd arrays
        InputImagePixelType mxval = 0.0;
        for(unsigned int m = 0; m < nPatch*m_Modality; m++)
          mxval += apd[i][m] * apd[k][m];

        mxval /= (nPatch - 1);

        if(m_Beta == 2)
          mxval *= mxval;
        else
          mxval = pow(mxval, float(m_Beta));

        Mx(i,k) = Mx(k,i) = mxval;
        }
      }

    // Now we can compute the weights by solving for the inverse of Mx
    MatrixType Mx_bar(n, n, 0.0);
    Mx_bar.fill_diagonal(m_Alpha);
    Mx_bar += Mx;

    // Define a vector of all ones
    vnl_vector<double> ones(n, 1.0);

    // Solve for the weights
    MatrixType InvMx = vnl_svd<double>(Mx_bar).inverse();
    MatrixType GroupMxT=vnl_matrix<double>(GroupMx).transpose();
    MatrixType Mx1 = GroupMxT * InvMx * GroupMx;
    MatrixType InvMx1 = vnl_svd<double>(Mx1).inverse();
    MatrixType tv= InvMx1 * GroupWeightMx;
    MatrixType tw(n, 1, 0.0);
    for (int i=0;i<n;i++)
      tw[i][0] = tv[m_GroupID[i]][0];

    MatrixType W1 = InvMx * tw;

    vnl_vector<double> W(n);
    for (int i=0;i<n;i++)
      W[i]=W1[i][0];

//    vnl_vector<double> W0 = vnl_svd<double>(Mx_bar).solve(ones);

//    // Normalize the weights
//    W0 *= 1.0 / dot_product(W0, ones);
//    // Solve for the weights
//    W = vnl_svd<double>(Mx_bar).solve(ones);
//
//    // Normalize the weights
//    W *= 1.0 / dot_product(W, ones);

    // Perform voting using Hongzhi's averaging scheme. Iterate over all segmentation patches
    for(unsigned int ni = 0; ni < nPatch; ni++)
      {
      IndexType idx = itTargetList[0].GetIndex(ni);
      if(this->GetOutput()->GetRequestedRegion().IsInside(idx))
        {
        for(int i = 0; i < n; i++)
          {
          // The segmentation at the corresponding patch location in atlas i
          InputImagePixelType label = *(patchSeg[i] + offPatchSeg[i][ni]);

          // Add that weight the posterior map for voxel at idx
          m_PosteriorMap[label]->SetPixel(idx, m_PosteriorMap[label]->GetPixel(idx) + W[i]);

          if (m_RetainVotingWeight)
            m_VotingWeight[i]->SetPixel(idx, m_VotingWeight[i]->GetPixel(idx) + W[i]);

          countermap->SetPixel(idx, countermap->GetPixel(idx) + W[i]);
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
  std::cout << std::endl << "VOTING finished" << std::endl;
  // Clear posterior maps
  if(!m_RetainPosteriorMaps)
    m_PosteriorMap.clear();

  for(OutIter it(this->GetOutput(), this->GetOutput()->GetBufferedRegion()); !it.IsAtEnd(); ++it)
  {
    IndexType idx = it.GetIndex();
    if (countermap->GetPixel(idx)<0.1)
      continue;

    if (m_RetainPosteriorMaps)
    {
      for(typename std::set<InputImagePixelType>::iterator sit = labelset.begin();
        sit != labelset.end(); ++sit)
      {
        m_PosteriorMap[*sit]->SetPixel(idx,m_PosteriorMap[*sit]->GetPixel(idx)/countermap->GetPixel(idx));
      }
    }
    if(m_RetainVotingWeight)
    {
      for (int i=0; i < n; i++)
        m_VotingWeight[i]->SetPixel(idx,m_VotingWeight[i]->GetPixel(idx)/countermap->GetPixel(idx));
    }
  }
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
  int offset,
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
    InputImagePixelType v = normtrg[i+offset];
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

