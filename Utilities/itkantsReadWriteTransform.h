#ifndef itkantsReadWriteTransform_h
#define itkantsReadWriteTransform_h

#include "itkDisplacementFieldTransform.h"
#include "itkGaussianExponentialDiffeomorphicTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itksys/SystemTools.hxx"
#include "itkCompositeTransform.h"

namespace itk
{
namespace ants
{
template <typename T, unsigned VImageDimension>
typename itk::Transform<T, VImageDimension, VImageDimension>::Pointer
ReadTransform(const std::string & filename,
              const bool useStaticCastForR = false) // This parameter changes to true by the programs that use R, so
                                                    // this code returns a different output for them.
{
  // We must explicitly check for file existance because failed reading is an acceptable
  // state for non-displacment feilds.
  if (!itksys::SystemTools::FileExists(filename.c_str()))
  {
    std::cerr << "Transform file does not exist: " << filename << std::endl;
    return nullptr;
  }

  bool hasTransformBeenRead = false;

  typedef typename itk::DisplacementFieldTransform<T, VImageDimension>   DisplacementFieldTransformType;
  typedef typename DisplacementFieldTransformType::DisplacementFieldType DisplacementFieldType;
  typedef itk::ImageFileReader<DisplacementFieldType>                    DisplacementFieldReaderType;

  typename DisplacementFieldReaderType::Pointer                fieldReader = DisplacementFieldReaderType::New();
  typedef typename itk::CompositeTransform<T, VImageDimension> CompositeTransformType;

  // There are known transform type extentions that should not be considered as imaging files
  // That would be used as deformatino feilds
  // If file is an hdf5 file, assume it is a transform instead of an image.
  bool recognizedExtension = false;
  recognizedExtension |= (itksys::SystemTools::GetFilenameLastExtension(filename) == ".h5");
  recognizedExtension |= (itksys::SystemTools::GetFilenameLastExtension(filename) == ".hdf5");
  recognizedExtension |= (itksys::SystemTools::GetFilenameLastExtension(filename) == ".hdf4");
  recognizedExtension |= (itksys::SystemTools::GetFilenameLastExtension(filename) == ".mat");
  recognizedExtension |= (itksys::SystemTools::GetFilenameLastExtension(filename) == ".txt");
  recognizedExtension |= (itksys::SystemTools::GetFilenameLastExtension(filename) == ".xfm");

  if (!recognizedExtension)
  {
    try
    {
      fieldReader->SetFileName(filename.c_str());
      fieldReader->Update();
      hasTransformBeenRead = true;
    }
    catch (...)
    {
      hasTransformBeenRead = false;
    }
  }

  typedef typename itk::Transform<T, VImageDimension, VImageDimension> TransformType;
  typename TransformType::Pointer                                      transform;
  if (hasTransformBeenRead)
  {
    typename DisplacementFieldTransformType::Pointer displacementFieldTransform = DisplacementFieldTransformType::New();
    displacementFieldTransform->SetDisplacementField(fieldReader->GetOutput());
    transform = dynamic_cast<TransformType *>(displacementFieldTransform.GetPointer());
  }
  else
  {
    typename itk::TransformFileReaderTemplate<T>::Pointer transformReader = itk::TransformFileReaderTemplate<T>::New();

    transformReader->SetFileName(filename.c_str());
    try
    {
      transformReader->Update();
    }
    catch (const itk::ExceptionObject & e)
    {
      std::cerr << "Transform reader for " << filename << " caught an ITK exception:\n";
      e.Print(std::cerr);
      return transform;
    }
    catch (const std::exception & e)
    {
      std::cerr << "Transform reader for " << filename << " caught an exception:\n";
      std::cerr << e.what() << std::endl;
      return transform;
    }
    catch (...)
    {
      std::cerr << "Transform reader for " << filename << " caught an unknown exception!!!\n";
      return transform;
    }

    const typename itk::TransformFileReaderTemplate<T>::TransformListType * const listOfTransforms =
      transformReader->GetTransformList();

    if (listOfTransforms->size() > 1)
    {
      typename CompositeTransformType::Pointer comp_transform = CompositeTransformType::New();
      for (typename itk::TransformFileReaderTemplate<T>::TransformListType::const_iterator i =
             listOfTransforms->begin();
           i != listOfTransforms->end();
           ++i)
      {
        comp_transform->AddTransform(dynamic_cast<TransformType *>(i->GetPointer()));
      }
      transform = dynamic_cast<TransformType *>(comp_transform.GetPointer());
    }
    else
    {

      transform = dynamic_cast<TransformType *>(listOfTransforms->front().GetPointer());
    }

    /** below is a bad thing but it's the only temporary fix i could find for ANTsR on unix --- B.A. */
    if (transform.IsNull() && (useStaticCastForR == true))
    {
      transform = static_cast<TransformType *>(listOfTransforms->front().GetPointer());
    }
  }
  return transform;
}

template <typename T, unsigned int VImageDimension>
int
WriteTransform(typename itk::Transform<T, VImageDimension, VImageDimension>::Pointer & xfrm,
               const std::string &                                                     filename)
{
  typedef typename itk::Transform<T, VImageDimension, VImageDimension>   GenericTransformType;
  typedef typename itk::DisplacementFieldTransform<T, VImageDimension>   DisplacementFieldTransformType;
  typedef typename DisplacementFieldTransformType::DisplacementFieldType DisplacementFieldType;
  typedef typename itk::ImageFileWriter<DisplacementFieldType>           DisplacementFieldWriter;
  typedef itk::TransformFileWriterTemplate<T>                            TransformWriterType;

  DisplacementFieldTransformType * dispXfrm = dynamic_cast<DisplacementFieldTransformType *>(xfrm.GetPointer());

  // if it's a displacement field transform or output file indicates it should be a transform
  try
  {
    if (dispXfrm != nullptr && filename.find(".mat") == std::string::npos && filename.find(".txt") == std::string::npos)
    {
      typename DisplacementFieldType::Pointer dispField = dispXfrm->GetModifiableDisplacementField();
      if (filename.find(".xfm") == std::string::npos && filename.find(".h5") == std::string::npos &&
          filename.find(".hdf5") == std::string::npos && filename.find(".hdf4") == std::string::npos)
      {
        typename DisplacementFieldWriter::Pointer writer = DisplacementFieldWriter::New();
        writer->SetInput(dispField);
        writer->SetFileName(filename.c_str());
        writer->Update();
      }
      else // creating a DisplacementFieldTransformType object to make everybody happy
      {
        typename DisplacementFieldTransformType::Pointer tmp_xfrm = DisplacementFieldTransformType::New();
        tmp_xfrm->SetDisplacementField(dispField);
        typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
        transformWriter->SetInput(tmp_xfrm);
        transformWriter->SetFileName(filename.c_str());
#if ITK_VERSION_MAJOR >= 5
        transformWriter->SetUseCompression(true);
#endif
        transformWriter->Update();
      }
    }
    else
    // regular transform, hope that everything works as expected!
    {
      typedef itk::CompositeTransform<T, VImageDimension> CompositeTransformType;
      typedef typename CompositeTransformType::Pointer    CompositeTransformPointer;
      typename TransformWriterType::Pointer               transformWriter = TransformWriterType::New();

      CompositeTransformType * comp_xfm = dynamic_cast<CompositeTransformType *>(xfrm.GetPointer());
      if (comp_xfm != nullptr)
      { // this is a composite transform, make sure it doesn't contain wiered stuff
        CompositeTransformPointer tmp_comp_xfm = CompositeTransformType::New();

        size_t numTransforms = comp_xfm->GetNumberOfTransforms();
        for (size_t i = 0; i < numTransforms; i++)
        {
          GenericTransformType *           _xfm = comp_xfm->GetNthTransform(i);
          DisplacementFieldTransformType * _dispXfrm = dynamic_cast<DisplacementFieldTransformType *>(_xfm);

          if (_dispXfrm != nullptr)
          { // assume that we have to make it DisplacementFieldTransform
            typename DisplacementFieldTransformType::Pointer _xfm_disp = DisplacementFieldTransformType::New();
            _xfm_disp->SetDisplacementField(_dispXfrm->GetModifiableDisplacementField());
            tmp_comp_xfm->AddTransform(_xfm_disp);
          }
          else
          { // asume we just pass it on
            tmp_comp_xfm->AddTransform(_xfm);
          }
        }
        transformWriter->SetInput(tmp_comp_xfm);
      }
      else
      { // this is  a simple transform
        transformWriter->SetInput(xfrm);
      }
      transformWriter->SetFileName(filename.c_str());
#if ITK_VERSION_MAJOR >= 5
      transformWriter->SetUseCompression(true);
#endif
      transformWriter->Update();
    }
  }
  catch (const itk::ExceptionObject & err)
  {
    std::cerr << "Can't write transform file " << filename << std::endl;
    std::cerr << "Exception Object caught: " << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

template <typename T, unsigned int VImageDimension>
int
WriteInverseTransform(typename itk::DisplacementFieldTransform<T, VImageDimension>::Pointer & xfrm,
                      const std::string &                                                     filename)
{
  typedef typename itk::DisplacementFieldTransform<T, VImageDimension>   DisplacementFieldTransformType;
  typedef typename DisplacementFieldTransformType::DisplacementFieldType DisplacementFieldType;
  typedef typename itk::ImageFileWriter<DisplacementFieldType>           DisplacementFieldWriter;
  typedef itk::TransformFileWriterTemplate<T>                            TransformWriterType;

  typename DisplacementFieldType::Pointer inverseDispField = xfrm->GetModifiableInverseDisplacementField();
  try
  {
    if (filename.find(".xfm") == std::string::npos && filename.find(".h5") == std::string::npos &&
        filename.find(".hdf5") == std::string::npos && filename.find(".hdf4") == std::string::npos)
    {
      typename DisplacementFieldWriter::Pointer writer = DisplacementFieldWriter::New();
      writer->SetInput(inverseDispField);
      writer->SetFileName(filename.c_str());
      writer->Update();
    }
    else
    // regular transform, but need to create inverse of the right type
    {
      typename DisplacementFieldTransformType::Pointer inv_xfrm = DisplacementFieldTransformType::New();
      inv_xfrm->SetDisplacementField(inverseDispField);
      typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
      transformWriter->SetInput(inv_xfrm);
      transformWriter->SetFileName(filename.c_str());
#if ITK_VERSION_MAJOR >= 5
      transformWriter->SetUseCompression(true);
#endif
      transformWriter->Update();
    }
  }
  catch (const itk::ExceptionObject & err)
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
