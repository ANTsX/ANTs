#ifndef itkantsReadWriteTransform_h
#define itkantsReadWriteTransform_h

#include "itkDisplacementFieldTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"

#include "itkCompositeTransform.h"

namespace itk
{
namespace ants
{
template <class T, unsigned VImageDimension>
typename itk::Transform<T, VImageDimension, VImageDimension>::Pointer
ReadTransform(const std::string & filename,
              const bool useStaticCastForR = false) // This parameter changes to true by the programs that use R, so this code
                                                     // returns a different output for them.
{
  // We must explicitly check for file existance because failed reading is an acceptable
  // state for non-displacment feilds.
  if( !itksys::SystemTools::FileExists( filename.c_str() ) )
    {
    std::cerr << "Transform file does not exist: " << filename << std::endl;
    return ITK_NULLPTR;
    }

  bool hasTransformBeenRead = false;

  typedef typename itk::DisplacementFieldTransform<T, VImageDimension> DisplacementFieldTransformType;
  typedef typename DisplacementFieldTransformType::DisplacementFieldType    DisplacementFieldType;
  typedef itk::ImageFileReader<DisplacementFieldType>                       DisplacementFieldReaderType;
  typename DisplacementFieldReaderType::Pointer fieldReader = DisplacementFieldReaderType::New();

  // There are known tranform type extentions that should not be considered as imaging files
  // That would be used as deformatino feilds
  // If file is an hdf5 file, assume it is a tranform instead of an image.
  if( filename.find(".h5") == std::string::npos
      && filename.find(".hdf5") == std::string::npos
      && filename.find(".hdf4") == std::string::npos
      && filename.find(".mat") == std::string::npos
      && filename.find(".txt") == std::string::npos
      )
    {
    try
      {
      fieldReader->SetFileName( filename.c_str() );
      fieldReader->Update();
      hasTransformBeenRead = true;
      }
    catch( ... )
      {
      hasTransformBeenRead = false;
      }
    }

  typedef typename itk::Transform<T, VImageDimension, VImageDimension> TransformType;
  typename TransformType::Pointer transform;
  if( hasTransformBeenRead )
    {
    typename DisplacementFieldTransformType::Pointer displacementFieldTransform =
      DisplacementFieldTransformType::New();
    displacementFieldTransform->SetDisplacementField( fieldReader->GetOutput() );
    transform = dynamic_cast<TransformType *>( displacementFieldTransform.GetPointer() );
    }
  else
    {
    typename itk::TransformFileReaderTemplate<T>::Pointer transformReader
      = itk::TransformFileReaderTemplate<T>::New();

    transformReader->SetFileName( filename.c_str() );
    try
      {
      transformReader->Update();
      }
    catch( const itk::ExceptionObject & e )
      {
      std::cerr << "Transform reader for "
                << filename << " caught an ITK exception:\n";
      e.Print( std::cerr );
      return transform;
      }
    catch( const std::exception & e )
      {
      std::cerr << "Transform reader for "
                << filename << " caught an exception:\n";
      std::cerr << e.what() << std::endl;
      return transform;
      }
    catch( ... )
      {
      std::cerr << "Transform reader for "
                << filename << " caught an unknown exception!!!\n";
      return transform;
      }

    const typename itk::TransformFileReaderTemplate<T>::TransformListType * const listOfTransforms =
      transformReader->GetTransformList();
    transform = dynamic_cast<TransformType *>( listOfTransforms->front().GetPointer() );

    /** below is a bad thing but it's the only temporary fix i could find for ANTsR on unix --- B.A. */
    if ( transform.IsNull() && ( useStaticCastForR == true ) )
       {
       transform = static_cast<TransformType *>( listOfTransforms->front().GetPointer() );
       }
    }
  return transform;
}

template <class T, unsigned int VImageDimension>
int
WriteTransform(typename itk::Transform<T, VImageDimension, VImageDimension>::Pointer & xfrm,
               const std::string & filename)
{
  typedef typename itk::DisplacementFieldTransform<T, VImageDimension>      DisplacementFieldTransformType;
  typedef typename DisplacementFieldTransformType::DisplacementFieldType    DisplacementFieldType;
  typedef typename itk::ImageFileWriter<DisplacementFieldType>              DisplacementFieldWriter;
  typedef itk::TransformFileWriterTemplate<T>                               TransformWriterType;

  DisplacementFieldTransformType *dispXfrm =
    dynamic_cast<DisplacementFieldTransformType *>(xfrm.GetPointer() );

  // if it's a displacement field transform
  try
    {
    if( dispXfrm != ITK_NULLPTR )
      {
      typename DisplacementFieldType::Pointer dispField = dispXfrm->GetModifiableDisplacementField();
      typename DisplacementFieldWriter::Pointer writer = DisplacementFieldWriter::New();
      writer->SetInput(dispField);
      writer->SetFileName(filename.c_str() );
      writer->Update();
      }
    else
    // regular transform
      {
      typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
      transformWriter->SetInput(xfrm);
      transformWriter->SetFileName(filename.c_str() );
      transformWriter->Update();
      }
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Can't write transform file " << filename << std::endl;
    std::cerr << "Exception Object caught: " << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
} // namespace itk
#endif // itkantsReadWriteTransform_h
