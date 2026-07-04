/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/ANTsX/ANTs/blob/main/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkN3SharpenFFT_h
#define itkN3SharpenFFT_h

#include "vnl/vnl_vector.h"
#include <complex>

// itk::PocketFFTCommon (ITKFFT) supersedes the deprecated vnl_fft_1d, but only
// exists in ITK versions that ship itkPocketFFTCommon.h. Fall back to vnl_fft_1d
// on older ITK so ANTs builds against either.
#if defined(__has_include) && __has_include("itkPocketFFTCommon.h")
#  include "itkPocketFFTCommon.h"
#  define ITK_ANTS_HAS_POCKETFFT_COMMON 1
#else
#  include "vnl/algo/vnl_fft_1d.h"
#endif

namespace itk
{
namespace ants
{
// In-place unnormalized forward 1-D transform (vnl fwd_transform convention:
// PocketFFT forward=false, scale=1).
template <typename TReal>
inline void
N3SharpenFFTForward(vnl_vector<std::complex<TReal>> & v)
{
#if defined(ITK_ANTS_HAS_POCKETFFT_COMMON)
  PocketFFTCommon::Transform1D(v.data_block(), v.size(), false, TReal{ 1 });
#else
  vnl_fft_1d<TReal> fft(v.size());
  fft.fwd_transform(v);
#endif
}

// In-place unnormalized inverse 1-D transform (vnl bwd_transform convention:
// PocketFFT forward=true, scale=1).
template <typename TReal>
inline void
N3SharpenFFTInverse(vnl_vector<std::complex<TReal>> & v)
{
#if defined(ITK_ANTS_HAS_POCKETFFT_COMMON)
  PocketFFTCommon::Transform1D(v.data_block(), v.size(), true, TReal{ 1 });
#else
  vnl_fft_1d<TReal> fft(v.size());
  fft.bwd_transform(v);
#endif
}
} // namespace ants
} // namespace itk

#endif
