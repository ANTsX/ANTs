/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkFEMElement3DMembrane1DOF_hxx
#define __itkFEMElement3DMembrane1DOF_hxx


namespace itk
{
namespace fem
{
template <typename TBaseClass>
Element3DMembrane1DOF<TBaseClass>::Element3DMembrane1DOF()
  : Superclass()
  , m_mat(0)
{}

// ////////////////////////////////////////////////////////////////////////
/*
 * Methods related to the physics of the problem.
 */

template <typename TBaseClass>
void
Element3DMembrane1DOF<TBaseClass>::GetStrainDisplacementMatrix(MatrixType & B, const MatrixType & shapeDgl) const
{}

template <typename TBaseClass>
void
Element3DMembrane1DOF<TBaseClass>::GetMassMatrix(MatrixType & Me) const
{
  // Call the parent's get matrix function
  Superclass::GetMassMatrix(Me);

  // Since parent class doesn't have the material properties,
  // we need to adjust Me matrix here for the density of the element.
  Me = Me * m_mat->RhoC;
}

template <typename TBaseClass>
void
Element3DMembrane1DOF<TBaseClass>::GetMaterialMatrix(MatrixType & D) const
{
  unsigned int d = 3;

  D.set_size(d, d);

  D.fill(0.0);

  // This is the main difference from the linear elasticity problem.
  /* Material properties matrix.  Simpler than linear elasticity. */
  Float disot = m_mat->E;
  for (unsigned int i = 0; i < d; i++)
  {
    D[i][i] = disot;
  }
}

template <typename TBaseClass>
void
Element3DMembrane1DOF<TBaseClass>::GetStiffnessMatrix(MatrixType & Ke) const
{
  Superclass::GetStiffnessMatrix(Ke);
}

template <typename TBaseClass>
void
Element3DMembrane1DOF<TBaseClass>::Read(std::istream & f, void * info)
{
  int n;
  /*
   * Convert the info pointer to a usable objects
   */
  ReadInfoType::MaterialArrayPointer mats = static_cast<ReadInfoType *>(info)->m_mat;

  /* first call the parent's read function */
  Superclass::Read(f, info);

  try
  {
    /*
     * Read and set the material pointer
     */
    FEMLightObject::SkipWhiteSpace(f);
    f >> n;
    if (!f)
    {
      goto out;
    }
    m_mat = dynamic_cast<const MaterialLinearElasticity *>(&*mats->Find(n));
  }
  catch (const FEMExceptionObjectNotFound & e)
  {
    throw FEMExceptionObjectNotFound(__FILE__, __LINE__, "Element3DMembrane1DOF::Read()", e.m_baseClassName, e.m_GN);
  }

  // Check if the material object was of correct class
  if (!m_mat)
  {
    throw FEMExceptionWrongClass(__FILE__, __LINE__, "Element3DMembrane1DOF::Read()");
  }

out:

  if (!f)
  {
    throw FEMExceptionIO(__FILE__, __LINE__, "Element3DMembrane1DOF::Read()", "Error reading FEM element!");
  }
}

/*
 * Write the element to the output stream.
 */
template <typename TBaseClass>
void
Element3DMembrane1DOF<TBaseClass>::Write(std::ostream & f) const
{
  // First call the parent's write function
  Superclass::Write(f);

  /*
   * then write the actual data (material number)
   * We also add some comments in the output file
   */
  f << "\t" << m_mat->GN << "\t% MaterialLinearElasticity ID\n";

  // check for errors
  if (!f)
  {
    throw FEMExceptionIO(__FILE__, __LINE__, "Element3DMembrane1DOF::Write()", "Error writing FEM element!");
  }
}

#ifdef _MSC_VER
// Declare a static dummy function to prevent a MSVC 6.0 SP5 from crashing.
// I have no idea why things don't work when this is not declared, but it
// looks like this declaration makes compiler forget about some of the
// troubles it has with templates.
static void
Dummy();

#endif // #ifdef _MSC_VER
} // namespace fem
} // namespace itk

#endif // #ifndef __itkFEMElement3DMembrane1DOF_hxx
