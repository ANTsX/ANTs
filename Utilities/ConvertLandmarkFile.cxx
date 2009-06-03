#include "itkLabeledPointSetFileReader.h"
#include "itkLabeledPointSetFileWriter.h"

#include "global.h"
#include "string.h"
#include <fstream.h>

int main( unsigned int argc, char *argv[] )
{
  if( argc < 2 )
    {
    std::cout << "Usage: " << argv[0] << " inputFile outputFile [percentage]" << std::endl;
    exit( 1 );
    }

  typedef double RealType;

  typedef itk::PointSet<long, ImageDimension> PointSetType;

  typedef itk::LabeledPointSetFileReader<PointSetType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  if( argc > 3 )
    {
    reader->SetRandomPercentage( atof( argv[3] ) );
    }
  reader->Update();

  std::cout << "Number of labels: " << reader->GetNumberOfLabels() << std::endl;
  std::cout << "Labels: ";
  for( unsigned int i = 0; i < reader->GetNumberOfLabels(); i++ )
    {
    std::cout << reader->GetLabelSet()->operator[](i) << " ";
    }
  std::cout << std::endl;

  typedef itk::LabeledPointSetFileWriter<PointSetType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( reader->GetOutput() );
  writer->Update();

  return 0;
}
