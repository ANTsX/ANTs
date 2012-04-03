#ifndef itkantsReadWriteTransform_h
#define itkantsReadWriteTransform_h

#include "itkDisplacementFieldTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
namespace itk
{
namespace ants
{
template <unsigned VImageDimension>
typename itk::Transform<double, VImageDimension, VImageDimension>::Pointer
ReadTransform(const std::string & filename)
{
  bool hasTransformBeenRead = false;

  typedef typename itk::DisplacementFieldTransform<double, VImageDimension> DisplacementFieldTransformType;
  typedef typename DisplacementFieldTransformType::DisplacementFieldType    DisplacementFieldType;
  typedef itk::ImageFileReader<DisplacementFieldType>                       DisplacementFieldReaderType;
  typename DisplacementFieldReaderType::Pointer fieldReader = DisplacementFieldReaderType::New();
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

  typedef typename itk::Transform<double, VImageDimension, VImageDimension> TransformType;
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
    typename itk::TransformFileReader::Pointer transformReader
      = itk::TransformFileReader::New();

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

    transform =
      dynamic_cast<TransformType *>( transformReader->GetTransformList()->front().GetPointer() );
    }
  return transform;
}

template <unsigned int VImageDimension>
int
WriteTransform(typename itk::Transform<double, VImageDimension, VImageDimension>::Pointer & xfrm,
               const std::string & filename)
{
  typedef typename itk::DisplacementFieldTransform<double, VImageDimension> DisplacementFieldTransformType;
  typedef typename DisplacementFieldTransformType::DisplacementFieldType    DisplacementFieldType;
  typedef typename itk::ImageFileWriter<DisplacementFieldType>              DisplacementFieldWriter;
  typedef itk::TransformFileWriter                                          TransformWriterType;

  DisplacementFieldTransformType *dispXfrm =
    dynamic_cast<DisplacementFieldTransformType *>(xfrm.GetPointer() );

  // if it's a displacement field transform
  try
    {
    if( dispXfrm != 0 )
      {
      typename DisplacementFieldType::Pointer dispField =
        dispXfrm->GetDisplacementField();
      typename DisplacementFieldWriter::Pointer writer =
        DisplacementFieldWriter::New();
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
