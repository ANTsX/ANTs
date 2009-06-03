/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPreservationOfPrincipalDirectionTensorReorientationImageFilter.cxx,v $
  Language:  C++
  Date:      $Date: 2009/03/17 19:01:36 $
  Version:   $Revision: 1.2 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkPreservationOfPrincipalDirectionTensorReorientationImageFilter_txx
#define _itkPreservationOfPrincipalDirectionTensorReorientationImageFilter_txx

#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkOffset.h"
#include "itkProgressReporter.h"
#include "itkObjectFactory.h"
#include "itkVector.h"
#include "itkPreservationOfPrincipalDirectionTensorReorientationImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkNumericTraitsFixedArrayPixel.h"

#include "itkVariableSizeMatrix.h"
#include "itkDecomposeTensorFunction.h"

#include <vnl/vnl_cross.h>
#include <vnl/vnl_inverse.h>
#include "vnl/algo/vnl_qr.h"
#include "vnl/algo/vnl_svd.h"
// #include <vnl/vnl_inverse_transpose.h>

namespace itk
{
template <typename TTensorImage, typename TVectorImage>
PreservationOfPrincipalDirectionTensorReorientationImageFilter<TTensorImage, TVectorImage>
::PreservationOfPrincipalDirectionTensorReorientationImageFilter()
{
  m_DeformationField = NULL;
  m_ReferenceImage = NULL;
}

template <typename TTensorImage, typename TVectorImage>
vnl_matrix<float>
PreservationOfPrincipalDirectionTensorReorientationImageFilter<TTensorImage, TVectorImage>
::FindRotation(  vnl_vector<float>  from,   vnl_vector<float> to)
{
#define CROSS(dest, v1, v2) { \
    dest[0] = v1[1] * v2[2] - v1[2] * v2[1]; \
    dest[1] = v1[2] * v2[0] - v1[0] * v2[2]; \
    dest[2] = v1[0] * v2[1] - v1[1] * v2[0]; \
    }

#define DOT(v1, v2) (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])

#define SUB(dest, v1, v2) { \
    dest[0] = v1[0] - v2[0]; \
    dest[1] = v1[1] - v2[1]; \
    dest[2] = v1[2] - v2[2]; }
/*
* A function for creating a rotation matrix that rotates a vector called
* "f rom" into another vector called "to".
* Input : from[3], to[3] which both must be *normalized* non-zero vectors
* Output: mtx[3][3] -- a 3x3 matrix in colum-major form
* Author: Tomas Moller, 1999
*/

  VnlMatrixType M(3, 3);
  M.fill(0);
  float v[3];
  float e, h;
  CROSS(v, from, to);
  e = DOT(from, to);
  float EPSILON = 1.e-5;
  if( e > 1.0 - EPSILON ) /* "from" almost or equal to "to"-vector? */
    {
    /* return identity */
    M(0, 0) = 1.0; M(0, 1) = 0.0; M(0, 2) = 0.0;
    M(1, 0) = 0.0; M(1, 1) = 1.0; M(1, 2) = 0.0;
    M(2, 0) = 0.0; M(2, 1) = 0.0; M(2, 2) = 1.0;
    }
  else if( e < -1.0 + EPSILON ) /* "from" almost or equal to negated "to"? */
    {
    float up[3], left[3];
    float invlen;
    float fxx, fyy, fzz, fxy, fxz, fyz;
    float uxx, uyy, uzz, uxy, uxz, uyz;
    float lxx, lyy, lzz, lxy, lxz, lyz;

    /* left=CROSS(from, (1,0,0)) */
    left[0] = 0.0; left[1] = from[2]; left[2] = -from[1];
    if( DOT(left, left) < EPSILON ) /* was left=CROSS(from,(1,0,0)) a good choice? */
      {
      /* here we now that left = CROSS(from, (1,0,0)) will be a good choice */
      left[0] = -from[2]; left[1] = 0.0; left[2] = from[0];
      }

    /* normalize "left" */
    invlen = 1.0 / sqrt(DOT(left, left) );
    left[0] *= invlen;
    left[1] *= invlen;
    left[2] *= invlen;
    CROSS(up, left, from);

    /* now we have a coordinate system, i.e., a basis; */
    /* M=(from, up, left), and we want to rotate to: */
    /* N=(-from, up, -left). This is done with the matrix:*/
    /* N*M^T where M^T is the transpose of M */
    fxx = -from[0] * from[0]; fyy = -from[1] * from[1]; fzz = -from[2] * from[2];
    fxy = -from[0] * from[1]; fxz = -from[0] * from[2]; fyz = -from[1] * from[2];

    uxx = up[0] * up[0]; uyy = up[1] * up[1]; uzz = up[2] * up[2];
    uxy = up[0] * up[1]; uxz = up[0] * up[2]; uyz = up[1] * up[2];

    lxx = -left[0] * left[0]; lyy = -left[1] * left[1]; lzz = -left[2] * left[2];
    lxy = -left[0] * left[1]; lxz = -left[0] * left[2]; lyz = -left[1] * left[2];

    /* symmetric matrix */
    M(0, 0) = fxx + uxx + lxx; M(0, 1) = fxy + uxy + lxy; M(0, 2) = fxz + uxz + lxz;
    M(1, 0) = M(0, 1); M(1, 1) = fyy + uyy + lyy; M(1, 2) = fyz + uyz + lyz;
    M(2, 0) = M(0, 2); M(2, 1) = M(1, 2); M(2, 2) = fzz + uzz + lzz;
    }
  else /* the most common case, unless "from"="to", or "from"=-"to" */
    {
#if 0
    /* unoptimized version - a good compiler will optimize this. */
    h = (1.0 - e) / DOT(v, v);
    M(0, 0) = e + h * v[0] * v[0]; M(0, 1) = h * v[0] * v[1] - v[2]; M(0,
                                                                       2) = h * v[0] * v[2] + v[1];
    M(1, 0) = h * v[0] * v[1] + v[2]; M(1, 1) = e + h * v[1] * v[1]; M(1,
                                                                       2) h * v[1] * v[2] - v[0];
    M(2, 0) = h * v[0] * v[2] - v[1]; M(2, 1) = h * v[1] * v[2] + v[0]; M(2,
                                                                          2) = e + h * v[2] * v[2];
#else
    /* ...otherwise use this hand optimized version (9 mults less) */
    float hvx, hvz, hvxy, hvxz, hvyz;
    h = (1.0 - e) / DOT(v, v);
    hvx = h * v[0];
    hvz = h * v[2];
    hvxy = hvx * v[1];
    hvxz = hvx * v[2];
    hvyz = hvz * v[1];
    M(0, 0) = e + hvx * v[0]; M(0, 1) = hvxy - v[2]; M(0, 2) = hvxz + v[1];
    M(1, 0) = hvxy + v[2]; M(1, 1) = e + h * v[1] * v[1]; M(1, 2) = hvyz - v[0];
    M(2, 0) = hvxz - v[1]; M(2, 1) = hvyz + v[0]; M(2, 2) = e + hvz * v[2];
#endif

    return M;
    }

  return M;
}

template <typename TTensorImage, typename TVectorImage>
Vector<float, 6>
PreservationOfPrincipalDirectionTensorReorientationImageFilter<TTensorImage, TVectorImage>
::ReorientTensors( Vector<float, 6>  fixedTens, Vector<float, 6> movingTens,  typename TTensorImage::IndexType index)
{
// get from and to vectors

  VnlMatrixType fDT(3, 3);

  fDT.fill(0);
  fDT(0, 0) = fixedTens[0];
  fDT(1, 1) = fixedTens[3];
  fDT(2, 2) = fixedTens[5];
  fDT(1, 0) = fDT(0, 1) = fixedTens[1];
  fDT(2, 0) = fDT(0, 2) = fixedTens[2];
  fDT(2, 1) = fDT(1, 2) = fixedTens[4];

  VnlMatrixType mDT(3, 3);
  mDT.fill(0);
  mDT(0, 0) = movingTens[0];
  mDT(1, 1) = movingTens[3];
  mDT(2, 2) = movingTens[5];
  mDT(1, 0) = mDT(0, 1) = movingTens[1];
  mDT(2, 0) = mDT(0, 2) = movingTens[2];
  mDT(2, 1) = mDT(1, 2) = movingTens[4];

// find rotation that maps 1st eigvecs - apply
  vnl_symmetric_eigensystem<float> eig(fDT);
  vnl_symmetric_eigensystem<float> meig(mDT);
// get angle between vecs

  vvec  from = meig.get_eigenvector(2);
  vvec  to = eig.get_eigenvector(2);
  float dp = dot_product(from, to);
  float angletheta = 0;
  if( dp >= 0.9999 )
    {
    angletheta = 0;
    }
  else
    {
    angletheta = acos(dp ) * 180 / 3.1614;
    }
  if(  angletheta  > 90 )
    {
    from = from * (-1.0);  angletheta = 180 - 90;
    }
  //  this->m_ThetaE+=theta;
//  if ( theta >= 45 ) { from = from *  (90-theta)/90 +  to * theta/90;  from = from / sqrt(dot_product(from,from)); }
// from = from *  0.5 +  to * 0.5;  from = from / sqrt(dot_product(from,from));
  VnlMatrixType M = this->FindRotation(from, to );

//  VnlMatrixType rDT = M*mDT*M.transpose();
  vvec rvec = M * meig.get_eigenvector(1);

// find rotation that maps 2nd eigvecs - apply
  from = rvec;
  to = eig.get_eigenvector(1);
  dp = dot_product(from, to);
  float anglepsi = 0;
  if( dp >= 0.9999 )
    {
    anglepsi = 0;
    }
  else
    {
    anglepsi = acos(dp ) * 180 / 3.1614;
    }
  if(  anglepsi  > 90 )
    {
    from = from * (-1.0);  anglepsi = 180 - 90;
    }

  if( this->m_AngleEnergyImage )
    {
    float mag = angletheta; // +anglepsi*anglepsi);
    this->m_AngleEnergyImage->SetPixel(index, mag);
    }

  VnlMatrixType M2 = this->FindRotation(from, to );
  M2 = M2 * M;
  VnlMatrixType rDT = M2 * mDT * M2.transpose();

  TensorType rmovingTens;
  rmovingTens[0] = rDT[0][0];
  rmovingTens[1] = rDT[1][0];
  rmovingTens[2] = rDT[2][0];
  rmovingTens[3] = rDT[1][1];
  rmovingTens[4] = rDT[1][2];
  rmovingTens[5] = rDT[2][2];

  return rmovingTens;
}

template <typename TTensorImage, typename TVectorImage>
typename TTensorImage::PixelType
PreservationOfPrincipalDirectionTensorReorientationImageFilter<TTensorImage, TVectorImage>
::PPD(typename TTensorImage::IndexType index, vnl_matrix<float> jMatrix  )
{
//  typedef vnl_matrix<float>        VMatrixType;
  InputImagePointer input = this->GetInput();
  InputPixelType    dtv = input->GetPixel(index);

  itk::Vector<float, 6> dtv2;

  VnlMatrixType DT(3, 3);
  DT.fill(0);
  DT(0, 0) = dtv[0];
  DT(1, 1) = dtv[3];
  DT(2, 2) = dtv[5];
  DT(1, 0) = DT(0, 1) = dtv[1];
  DT(2, 0) = DT(0, 2) = dtv[2];
  DT(2, 1) = DT(1, 2) = dtv[4];
  vnl_symmetric_eigensystem<float> eig(DT);
  typedef vnl_vector<float> vvec;

  vvec e1 = eig.get_eigenvector(2);
  vvec e1p =  jMatrix * e1;
  e1p = e1p / sqrt(dot_product(e1p, e1p) );

  vvec e2 = jMatrix * eig.get_eigenvector(1);
  vvec e2p = e2 - dot_product( e2, e1p) * e1p;
  e2p = e2p / sqrt(dot_product(e2p, e2p) );

  vvec e3p = vnl_cross_3d<float>(e1p, e2p);
  e3p = e3p / sqrt(dot_product(e3p, e3p) );

  VnlMatrixType DTrec = eig.get_eigenvalue(2) * outer_product(e1p, e1p) + eig.get_eigenvalue(1) * outer_product(e2p,
                                                                                                                e2p)
    + eig.get_eigenvalue(0) * outer_product(e3p, e3p);

//  if (  (eig.get_eigenvalue(2) -   eig.get_eigenvalue(1)) > 0.5 )
    {
    // std::cout << " DT IN " << DT << std::endl;
    // std::cout << " DT OU " << DTrec << std::endl;
    // std::cout << " eig1 " <<  eig.get_eigenvalue(0) << " 2 " << eig.get_eigenvalue(1) << std::endl;
    }

  dtv2[0] = DTrec(0, 0);
  dtv2[3] = DTrec(1, 1);
  dtv2[5] = DTrec(2, 2);
  dtv2[1] = DTrec(0, 1);
  dtv2[2] = DTrec(0, 2);
  dtv2[4] = DTrec(1, 2);

  return dtv2;
}

template <typename TTensorImage, typename TVectorImage>
void
PreservationOfPrincipalDirectionTensorReorientationImageFilter<TTensorImage, TVectorImage>
::GenerateData()
{
  typedef VariableSizeMatrix<float> MatrixType2;

  // get input and output images
  InputImagePointer  input = this->GetInput();
  OutputImagePointer output = this->GetOutput();
  InputSizeType      inputSize = input->GetLargestPossibleRegion().GetSize();

  output->SetRegions( m_DeformationField->GetLargestPossibleRegion() );
  output->SetSpacing( m_DeformationField->GetSpacing() );
  output->SetOrigin( m_DeformationField->GetOrigin() );
  output->SetDirection( m_DeformationField->GetDirection() );
  output->Allocate();

  this->m_AngleEnergyImage = FloatImageType::New();
  this->m_AngleEnergyImage->SetRegions( m_DeformationField->GetLargestPossibleRegion() );
  this->m_AngleEnergyImage->SetSpacing(m_DeformationField->GetSpacing() );
  this->m_AngleEnergyImage->SetOrigin(m_DeformationField->GetOrigin() );
  this->m_AngleEnergyImage->Allocate();
  this->m_AngleEnergyImage->FillBuffer(0);

  ImageRegionConstIterator<InputImageType>      inputIt( input, input->GetLargestPossibleRegion() );
  ImageRegionIteratorWithIndex<OutputImageType> outputIt( output, output->GetLargestPossibleRegion() );

  typedef itk::VectorLinearInterpolateImageFunction<InputImageType> VectorInterpType;
  typename VectorInterpType::Pointer inputInterp = VectorInterpType::New();
  inputInterp->SetInputImage( input );

  // DeformationFieldPointer tempField=DeformationFieldType::New();

  // copy field to TempField
  // tempField->SetSpacing(m_DeformationField->GetSpacing() );
  // tempField->SetOrigin(m_DeformationField->GetOrigin() );
  // tempField->SetLargestPossibleRegion(
  // m_DeformationField->GetLargestPossibleRegion() );
  // tempField->SetRequestedRegion(
  // m_DeformationField->GetRequestedRegion() );
  // tempField->SetBufferedRegion(m_DeformationField->GetBufferedRegion() );
  // tempField->Allocate();
  //  tempField->FillBuffer(  );

  // for ( inputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt )
  // {
  // tempField->SetPixel( inputIt.GetIndex(), m_DeformationField->GetPixel( inputIt.GetIndex()));
  // }

  ImageRegionIteratorWithIndex<DeformationFieldType> dispIt( m_DeformationField,
                                                             m_DeformationField->GetLargestPossibleRegion() );

  // double det=0.0;
  unsigned int posoff = 1;
  // float difspace=1.0;
  float         space = 1.0;
  unsigned long ct = 0;
  // for all voxels
  for( outputIt.GoToBegin(), dispIt.GoToBegin();
       !outputIt.IsAtEnd(), !dispIt.IsAtEnd();
       ++outputIt, ++dispIt )
    {
    MatrixType2 jMatrix;
    jMatrix.SetSize(ImageDimension, ImageDimension);
    jMatrix.Fill(0.0);

    // Get Jacobian of displacement field
    typename DeformationFieldType::IndexType rindex = dispIt.GetIndex();
    float mindist = 1.0;
    bool  oktosample = true;
    float dist = 100.0;
    // ignore outer edges of warpfield
    for( unsigned int row = 0; row < ImageDimension; row++ )
      {
      dist = fabs( (float)rindex[row] );
      if( dist < mindist )
        {
        oktosample = false;
        }

      dist = fabs( (float)m_DeformationField->GetLargestPossibleRegion().GetSize()[row] - (float)rindex[row] );
      if( dist < mindist )
        {
        oktosample = false;
        }
      }

    typename OutputImageType::PointType outputPt;
    output->TransformIndexToPhysicalPoint(outputIt.GetIndex(), outputPt);
    if( !inputInterp->IsInsideBuffer(outputPt) )
      {
      oktosample = false;
      }

    if( oktosample )
      {
      typename DeformationFieldType::IndexType temp = rindex;
      typename DeformationFieldType::PixelType cpix = dispIt.Value();
      typename DeformationFieldType::IndexType difIndex[ImageDimension][2];
      for( unsigned int row = 0; row < ImageDimension; row++ )
        {
        difIndex[row][0] = rindex;
        difIndex[row][1] = rindex;
        typename DeformationFieldType::IndexType ddrindex = rindex;
        typename DeformationFieldType::IndexType ddlindex = rindex;

        if( rindex[row] < (int)m_DeformationField->GetLargestPossibleRegion().GetSize()[row] - 2 )
          {
          difIndex[row][0][row] = rindex[row] + posoff;
          ddrindex[row] = rindex[row] + posoff * 2;
          }
        if( rindex[row] > 1 )
          {
          difIndex[row][1][row] = rindex[row] - 1;
          ddlindex[row] = rindex[row] - 2;
          }

        float h = 0.25;
        space = 1.0; // should use image spacing here?

        typename DeformationFieldType::PixelType rpix = m_DeformationField->GetPixel(difIndex[row][1]);
        typename DeformationFieldType::PixelType lpix = m_DeformationField->GetPixel(difIndex[row][0]);

        typename DeformationFieldType::PixelType rrpix = m_DeformationField->GetPixel(ddrindex);
        typename DeformationFieldType::PixelType llpix = m_DeformationField->GetPixel(ddlindex);

        rpix = rpix * h + cpix * (1. - h);
        lpix = lpix * h + cpix * (1. - h);
        rrpix = rrpix * h + rpix * (1. - h);
        llpix = llpix * h + lpix * (1. - h);

        // 4th order centered difference
        typename DeformationFieldType::PixelType dPix;
        dPix = ( lpix * 8.0 - rpix * 8.0 + rrpix - llpix ) * space / (12.0); // 4th order centered differenceX
        for( unsigned int col = 0; col < ImageDimension; col++ )
          {
          float val;
          if( row == col )
            {
            val = dPix[col] / m_DeformationField->GetSpacing()[col] + 1.0;
            }
          else
            {
            val = dPix[col] / m_DeformationField->GetSpacing()[col];
            }
          jMatrix(col, row) = val;
          }
        }

      bool           finitestrain = false;
      bool           ppd = true;
      InputPixelType outTensor;

      if( finitestrain )
        {
        // Copy tensor to matrix form
        // InputPixelType inTensor = inputIt.Value();
        InputPixelType inTensor = inputInterp->Evaluate(outputPt);
        MatrixType2    tensor;
        tensor.SetSize(3, 3);
        tensor[0][0] = inTensor[0];
        tensor[0][1] = tensor[1][0] = inTensor[1];
        tensor[0][2] = tensor[2][0] = inTensor[2];
        tensor[1][1] = inTensor[3];
        tensor[1][2] = tensor[2][1] = inTensor[4];
        tensor[2][2] = inTensor[5];

        // replace with Polar Decomposition
        MatrixType2 rotationMatrix;
        MatrixType2 rotationMatrixT;

        VnlMatrixType      M = jMatrix.GetVnlMatrix();
        VnlMatrixType      PQ = M;
        VnlMatrixType      NQ = M;
        VnlMatrixType      PQNQDiff;
        const unsigned int maximumIterations = 100;
        for( unsigned int ni = 0; ni < maximumIterations; ni++ )
          {
          // Average current Qi with its inverse transpose
          NQ = ( PQ + vnl_inverse_transpose( PQ ) ) / 2.0;
          PQNQDiff = NQ - PQ;
          if( PQNQDiff.frobenius_norm() < 1e-4 )
            {
            // std::cout << "Polar decomposition used "      << ni << " iterations " << std::endl;
            break;
            }
          else
            {
            PQ = NQ;
            }
          }

        rotationMatrix = NQ;
        rotationMatrixT = rotationMatrix.GetTranspose();
        // Rotate tensor
        MatrixType2 rTensor = rotationMatrixT * tensor * rotationMatrix;

        // Copy tensor back to vector form and update output image
        outTensor[0] = rTensor[0][0];
        outTensor[1] = rTensor[1][0];
        outTensor[2] = rTensor[2][0];
        outTensor[3] = rTensor[1][1];
        outTensor[4] = rTensor[1][2];
        outTensor[5] = rTensor[2][2];
        }

      else if( ppd )
        {
        // InputPixelType inTensor = inputIt.Value();
        InputPixelType inTensor = inputInterp->Evaluate(outputPt);
        InputPixelType dtv = inTensor;

        // valid values?
        bool docompute = true;
        for( unsigned int jj = 0; jj < 6; jj++ )
          {
          float ff = dtv[jj];
          if( vnl_math_isnan(ff) )
            {
            docompute = false;
            }
          }

        if( docompute )
          {
          itk::Vector<float, 6> dtv2;
          bool                  verbose = false;
          VnlMatrixType         DT(3, 3);
          DT.fill(0);
          DT(0, 0) = dtv[0];
          DT(1, 1) = dtv[3];
          DT(2, 2) = dtv[5];
          DT(1, 0) = DT(0, 1) = dtv[1];
          DT(2, 0) = DT(0, 2) = dtv[2];
          DT(2, 1) = DT(1, 2) = dtv[4];

          if( verbose )
            {
            std::cout << " DT " << DT << std::endl;
            }
          if( verbose )
            {
            std::cout << " jM " << jMatrix  << std::endl;
            }

          vnl_symmetric_eigensystem<float> eig(DT);
          typedef vnl_vector<float> vvec;

          vvec e1 = eig.get_eigenvector(2);
          vvec e1p =  jMatrix.GetTranspose() * e1;
          e1p = e1p / sqrt(dot_product(e1p, e1p) );

          vvec e2 = jMatrix.GetTranspose() * eig.get_eigenvector(1);
          vvec e2p = e2 - dot_product( e2, e1p) * e1p;
          e2p = e2p / sqrt(dot_product(e2p, e2p) );

          vvec e3p = vnl_cross_3d<float>(e1p, e2p);
          e3p = e3p / sqrt(dot_product(e3p, e3p) );

          VnlMatrixType DTrec = eig.get_eigenvalue(2) * outer_product(e1p, e1p) + eig.get_eigenvalue(1) * outer_product(
              e2p, e2p) + eig.get_eigenvalue(0) * outer_product(e3p, e3p);

          if( verbose )
            {
            std::cout << " DTrec " << DTrec  << " ind "  <<  rindex << std::endl;
            }

          dtv2[0] = DTrec(0, 0);
          dtv2[3] = DTrec(1, 1);
          dtv2[5] = DTrec(2, 2);
          dtv2[1] = DTrec(0, 1);
          dtv2[2] = DTrec(0, 2);
          dtv2[4] = DTrec(1, 2);

          outTensor = dtv2;
          }
        else
          {
          outTensor = inTensor;
          }
        }

      if( this->m_ReferenceImage )
        {
        typename InputImageType::PixelType refTensor = this->m_ReferenceImage->GetPixel(outputIt.GetIndex() );
        outTensor = this->ReorientTensors(refTensor, outTensor, outputIt.GetIndex() );
        }
      outputIt.Set( outTensor );

      // std::cout << inputIt.Value() << " - " << outputIt.Value() << " ";
      } // oktosample if
    else
      {
      outputIt.Set( inputIt.Value() );
      }

    ct++;
    } // pixelIterations
}

/**
 * Standard "PrintSelf" method
 */
template <typename TTensorImage, typename TVectorImage>
void
PreservationOfPrincipalDirectionTensorReorientationImageFilter<TTensorImage, TVectorImage>
::PrintSelf(
  std::ostream& os,
  Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}
} // end namespace itk

#endif
