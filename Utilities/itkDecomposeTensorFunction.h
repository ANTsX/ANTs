/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: itkDecomposeTensorFunction.h,v $
  Language:  C++
  Date:      $Date: 2008/11/15 23:46:06 $
  Version:   $Revision: 1.16 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDecomposeTensorFunction_h
#define __itkDecomposeTensorFunction_h

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "itkProcessObject.h"

#include "itkVariableSizeMatrix.h"

namespace itk
{
/** \class DecomposeTensorFunction
 *
 */
template <typename TInput,
          typename TRealType = float,
          typename TOutput = itk::VariableSizeMatrix<TRealType>
          >
class DecomposeTensorFunction : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef DecomposeTensorFunction  Self;
  typedef ProcessObject            Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( DecomposeTensorFunction, ProcessObject );

  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
  typedef TInput  InputMatrixType;
  typedef TOutput OutputMatrixType;

  /** Define the data type and the vector of data type used in calculations. */
  typedef TRealType RealType;

  // Wrappers for vnl routines
  void EvaluateEigenDecomposition( InputMatrixType &, OutputMatrixType &, OutputMatrixType & );

  void EvaluateSymmetricEigenDecomposition( InputMatrixType &, OutputMatrixType &, OutputMatrixType & );

  void EvaluateQRDecomposition( InputMatrixType &, OutputMatrixType &, OutputMatrixType & );

  void EvaluateSVDDecomposition( InputMatrixType &, OutputMatrixType &, OutputMatrixType &, OutputMatrixType & );

  void EvaluateSVDEconomyDecomposition( InputMatrixType &, OutputMatrixType &, OutputMatrixType & );

  void EvaluateLeftPolarDecomposition( InputMatrixType &, OutputMatrixType &, OutputMatrixType & );

  void EvaluateRightPolarDecomposition( InputMatrixType &, OutputMatrixType &, OutputMatrixType & );

  void EvaluateCholeskyDecomposition( InputMatrixType &, OutputMatrixType & );

  RealType EvaluateDeterminant( InputMatrixType & );

  DecomposeTensorFunction();
  virtual ~DecomposeTensorFunction()
  {
  }

protected:

  void PrintSelf( std::ostream& os, Indent indent ) const;

private:

  DecomposeTensorFunction(const Self &); // purposely not implemented
  void operator=(const Self &);          // purposely not implemented
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDecomposeTensorFunction.hxx"
#endif

#endif
