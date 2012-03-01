#ifndef _TensorFunctions_cxx_
#define _TensorFunctions_cxx_

// #include "itkSmallStrainDiffusionTensorReorientationImageFilter.h"
// #include "itkVectorIndexSelectionCastImageFilter.h"

#include "vnl/algo/vnl_matrix_inverse.h"
/*
template <class TensorType>
float  GetTensorFA( TensorType dtv)
{

// inner product
  float xx = dtv[0];
  float xy = dtv[1];
  float xz = dtv[2];
  float yy = dtv[3];
  float yz = dtv[4];
  float zz = dtv[5];
  float isp = ( xx*xx + yy*yy + zz*zz + 2.0*(xy*xy + xz*xz + yz*yz) );

  // Computed as
  // FA = vcl_sqrt(1.5*sum(sum(N.*N))/sum((sum(D.*D))))
  // where N = D - ((1/3)*trace(D)*eye(3,3))
  // equation (28) in http://lmi.bwh.harvard.edu/papers/pdfs/2002/westinMEDIA02.pdf

  if( isp > 0.0 )
    {
      float trace = dtv[0];
      trace += dtv[3];
      trace += dtv[5];

      float anisotropy = 3.0 * isp - trace * trace;
      float fractionalAnisotropy =( vcl_sqrt(anisotropy / ( 2.0 * isp ) ) );
      return fractionalAnisotropy;
    }

   return 0.0 ;



}

template <class TVectorType, class TTensorType>
float  GetMetricTensorCost(  TVectorType dpath,  TTensorType dtv , unsigned int matrixpower)
{

  typedef vnl_matrix<double>        MatrixType;
  MatrixType DT(3,3);
  DT.fill(0);
  DT(0,0)=dtv[0];
  DT(1,1)=dtv[3];
  DT(2,2)=dtv[5];
  DT(1,0)=DT(0,1)=dtv[1];
  DT(2,0)=DT(0,2)=dtv[2];
  DT(2,1)=DT(1,2)=dtv[4];
  vnl_symmetric_eigensystem< double > eig(DT);
  double e1 = (eig.D(0,0));
  double e2 = (eig.D(1,1));
  double e3 = (eig.D(2,2));
  double etot=e1+e2+e3;
  if (etot==0) etot=1;

  MatrixType  vec(3,1);
  vec(0,0)=dpath[0];
  vec(1,0)=dpath[1];
  vec(2,0)=dpath[2];
  MatrixType inv =  vnl_matrix_inverse<double>(DT);

  for (unsigned int lo=1; lo<matrixpower; lo++)   inv = inv*inv;

MatrixType sol= vec.transpose()*inv*vec;
float cost = sol(0,0);///etot;

  return sqrt(cost);

}

template <class TVectorType, class TTensorType>
TVectorType ChangeTensorByVector(  TVectorType dpath,  TTensorType dtv, float epsilon)
{

  typedef TVectorType VectorType;

  typedef vnl_matrix<double>        MatrixType;
  MatrixType DT(3,3);
  DT.fill(0);
  DT(0,0)=dtv[0];
  DT(1,1)=dtv[3];
  DT(2,2)=dtv[5];
  DT(1,0)=DT(0,1)=dtv[1];
  DT(2,0)=DT(0,2)=dtv[2];
  DT(2,1)=DT(1,2)=dtv[4];
  vnl_symmetric_eigensystem< double > eig(DT);
  double e3 = (eig.D(0,0));
  double e2 = (eig.D(1,1));
  double e1 = (eig.D(2,2));
  double etot=e1+e2+e3;
  if (etot==0) etot=1;

  MatrixType  vec(3,1);
  vec(0,0)=dpath[0];
  vec(1,0)=dpath[1];
  vec(2,0)=dpath[2];

  MatrixType  evec1(3,1);//biggest
  evec1(0,0)=eig.V(0,2);
  evec1(1,0)=eig.V(1,2);
  evec1(2,0)=eig.V(2,2);
  MatrixType  evec2(3,1);// middle
  evec2(0,0)=eig.V(0,1);
  evec2(1,0)=eig.V(1,1);
  evec2(2,0)=eig.V(2,1);
  MatrixType  evec3(3,1);//smallest
  evec3(0,0)=eig.V(0,0);
  evec3(1,0)=eig.V(1,0);
  evec3(2,0)=eig.V(2,0);

  float temp;
  temp=(vec.transpose()*evec1)(0,0);
  temp=sqrt(temp*temp);
  e1 *=( 1.0 - epsilon*temp);
  if (e1 < 1.e-11) e1=1.e-11;

  temp=(vec.transpose()*evec2)(0,0);
  temp=sqrt(temp*temp);
  e2 *=( 1.0 - epsilon*temp);
  if (e2 < 1.e-11) e2=1.e-11;

  temp=(vec.transpose()*evec3)(0,0);
  temp=sqrt(temp*temp);
  e3 *=( 1.0 - epsilon*temp);
   if (e3 < 1.e-11) e3=1.e-11;

  DT =  (evec3*evec3.transpose())*e3 +   (evec2*evec2.transpose())*e2  +  (evec1*evec1.transpose())*e1;

  itk::Vector<float,6> newtens;
  newtens[0]=DT(0,0);
  newtens[3]=DT(1,1);
  newtens[5]=DT(2,2);
  newtens[1]=DT(0,1);
  newtens[2]=DT(0,2);
  newtens[4]=DT(2,1);

  return newtens;

}

template<class TTensorType>
float  GetTensorADC(  TensorType dtv,  unsigned int opt = 0)
{

  float eps=1.e-9,mag=0;
  for (unsigned int jj=0; jj<6; jj++)
  {
    float ff=dtv[jj];
    mag+=ff*ff;
    if ( vnl_math_isnan( ff ) || vnl_math_isinf(ff)   )
    {
      return 0;
    }
  }
  mag=sqrt(mag);

  if (  dtv[1]==0 && dtv[2] == 0 && dtv[4]==0) return 0;
  if (mag < eps) { return 0; }


  //  typedef itk::Vector<float, 6> TensorTypeIn;
  typedef itk::Vector<float, 6> TensorType;
  typedef itk::Vector<float, 6> TensorTypeOut;
  typedef itk::Vector<float, 6> TensorTypeIn;
  const unsigned int ImageDimension = 3;
  typedef itk::Image<TensorTypeIn, ImageDimension> InTensorImageType;
  typedef itk::Image<TensorTypeOut, ImageDimension> OutTensorImageType;
  typedef itk::Image<TensorTypeOut, ImageDimension> TensorImageType;
  typedef     InTensorImageType::Pointer TensorImagePointer;
  typedef float  PixelType;
  typedef itk::Vector<float,ImageDimension>         VectorType1;
  typedef itk::Image<VectorType1,ImageDimension>     FieldType;
  typedef itk::Image<PixelType,ImageDimension> ImageType;
  typedef itk::ImageFileReader<InTensorImageType> readertype;
  typedef itk::ImageFileWriter<OutTensorImageType> writertype;
  typedef ImageType::IndexType IndexType;
  typedef ImageType::SizeType SizeType;
  typedef ImageType::SpacingType SpacingType;
//   typedef itk::AffineTransform<double,ImageDimension>   AffineTransformType;
//   typedef itk::LinearInterpolateImageFunction<ImageType,double>  InterpolatorType1;
//   typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double>  InterpolatorType2;

  itk::Vector<float, 6> dtv2;
  typedef vnl_matrix<double>        MatrixType;
  MatrixType DT(3,3);
  DT.fill(0);
  DT(0,0)=dtv[0];
  DT(1,1)=dtv[3];
  DT(2,2)=dtv[5];
  DT(1,0)=DT(0,1)=dtv[1];
  DT(2,0)=DT(0,2)=dtv[2];
  DT(2,1)=DT(1,2)=dtv[4];
  //  if (takelog ) std::cout << " TAKING LOG " << std::endl;  else std::cout << "TAKING EXP " << std::endl;
  //  std::cout << " dtv " << dtv << std::endl;
  vnl_symmetric_eigensystem< double > eig(DT);
  double e1 = (eig.D(0,0));
  double e2 = (eig.D(1,1));
  double e3 = (eig.D(2,2));

  if (opt <= 1 )  return (e1+e1+e3)/3.0;
  //  else if (opt == 4 ) return e2;
  else if (opt == 3 ) return (e2+e1)/2.0;
  else if (opt == 2 ) return e3;

}

template<class TTensorType>
itk::RGBPixel< unsigned char >   GetTensorRGB( TTensorType dtv )
{

  itk::RGBPixel< unsigned char > zero;
  zero.Fill(0);
  float eps=1.e-9,mag=0;
  for (unsigned int jj=0; jj<6; jj++)
  {
    float ff=dtv[jj];
    mag+=ff*ff;
    if ( vnl_math_isnan( ff ) || vnl_math_isinf(ff)   )
    {
      return zero;
    }
  }
  mag=sqrt(mag);

  if (  dtv[1]==0 && dtv[2] == 0 && dtv[4]==0) return zero;
  if (mag < eps) { return zero; }

  //  typedef itk::Vector<float, 6> TensorTypeIn;
  typedef itk::Vector<float, 6> TensorType;
  typedef itk::Vector<float, 6> TensorTypeOut;
  typedef itk::Vector<float, 6> TensorTypeIn;
  const unsigned int ImageDimension = 3;
  typedef itk::Image<TensorTypeIn, ImageDimension> InTensorImageType;
  typedef itk::Image<TensorTypeOut, ImageDimension> OutTensorImageType;
  typedef itk::Image<TensorTypeOut, ImageDimension> TensorImageType;
  typedef     InTensorImageType::Pointer TensorImagePointer;
  typedef float  PixelType;
  typedef itk::Vector<float,ImageDimension>         VectorType1;
  typedef itk::Image<VectorType1,ImageDimension>     FieldType;
  typedef itk::Image<PixelType,ImageDimension> ImageType;
  typedef itk::ImageFileReader<InTensorImageType> readertype;
  typedef itk::ImageFileWriter<OutTensorImageType> writertype;
  typedef ImageType::IndexType IndexType;
  typedef ImageType::SizeType SizeType;
  typedef ImageType::SpacingType SpacingType;
//   typedef itk::AffineTransform<double,ImageDimension>   AffineTransformType;
//   typedef itk::LinearInterpolateImageFunction<ImageType,double>  InterpolatorType1;
//   typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double>  InterpolatorType2;

  itk::Vector<float, 6> dtv2;
  typedef vnl_matrix<double>        MatrixType;
  MatrixType DT(3,3);
  DT.fill(0);
  DT(0,0)=dtv[0];
  DT(1,1)=dtv[3];
  DT(2,2)=dtv[5];
  DT(1,0)=DT(0,1)=dtv[1];
  DT(2,0)=DT(0,2)=dtv[2];
  DT(2,1)=DT(1,2)=dtv[4];
  //  if (takelog ) std::cout << " TAKING LOG " << std::endl;  else std::cout << "TAKING EXP " << std::endl;
  //  std::cout << " dtv " << dtv << std::endl;
  vnl_symmetric_eigensystem< double > eig(DT);


  double e1 = (eig.D(0,0));
  double e2 = (eig.D(1,1));
  double e3 = (eig.D(2,2));

  itk::RGBPixel< float > rgb;

  float xx = dtv[0];
  float xy = dtv[1];
  float xz = dtv[2];
  float yy = dtv[3];
  float yz = dtv[4];
  float zz = dtv[5];
  float isp = ( xx*xx + yy*yy + zz*zz + 2.0*(xy*xy + xz*xz + yz*yz) );

  float fa=0.0;
  if( isp > 0.0 )
    {
      float trace = dtv[0];
      trace += dtv[3];
      trace += dtv[5];
            float anisotropy = 3.0 * isp - trace * trace;
      fa=( vcl_sqrt(anisotropy / ( 2.0 * isp ) ) );
    }

  //rgb[0]=eig.V(2,0)*fa*255;//+eig.V(1,0)*e2;
  //rgb[1]=eig.V(2,1)*fa*255;//+eig.V(1,1)*e2;
  //  rgb[2]=eig.V(2,2)*fa*255;//+eig.V(1,2)*e2;

  rgb[0]=eig.V(0,2)*fa*255;//+eig.V(1,0)*e2;
  rgb[1]=eig.V(1,2)*fa*255;//+eig.V(1,1)*e2;
  rgb[2]=eig.V(2,2)*fa*255;//+eig.V(1,2)*e2;

  return rgb;

}

template <class TTensorType>
itk::RGBPixel< float >   GetTensorPrincipalEigenvector( TTensorType dtv )
{

  itk::RGBPixel< float > zero;
  zero.Fill(0);
  float eps=1.e-9,mag=0;
  for (unsigned int jj=0; jj<6; jj++)
  {
    float ff=dtv[jj];
    mag+=ff*ff;
    if ( vnl_math_isnan( ff ) || vnl_math_isinf(ff)   )
    {
      return zero;
    }
  }
  mag=sqrt(mag);

  if (  dtv[1]==0 && dtv[2] == 0 && dtv[4]==0) return zero;
  if (mag < eps) { return zero; }

  //  typedef itk::Vector<float, 6> TensorTypeIn;
  typedef itk::Vector<float, 6> TensorType;
  typedef itk::Vector<float, 6> TensorTypeOut;
  typedef itk::Vector<float, 6> TensorTypeIn;
  const unsigned int ImageDimension = 3;
  typedef itk::Image<TensorTypeIn, ImageDimension> InTensorImageType;
  typedef itk::Image<TensorTypeOut, ImageDimension> OutTensorImageType;
  typedef itk::Image<TensorTypeOut, ImageDimension> TensorImageType;
  typedef     InTensorImageType::Pointer TensorImagePointer;
  typedef float  PixelType;
  typedef itk::Vector<float,ImageDimension>         VectorType1;
  typedef itk::Image<VectorType1,ImageDimension>     FieldType;
  typedef itk::Image<PixelType,ImageDimension> ImageType;
  typedef itk::ImageFileReader<InTensorImageType> readertype;
  typedef itk::ImageFileWriter<OutTensorImageType> writertype;
  typedef ImageType::IndexType IndexType;
  typedef ImageType::SizeType SizeType;
  typedef ImageType::SpacingType SpacingType;
//   typedef itk::AffineTransform<double,ImageDimension>   AffineTransformType;
//   typedef itk::LinearInterpolateImageFunction<ImageType,double>  InterpolatorType1;
//   typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double>  InterpolatorType2;

  itk::Vector<float, 6> dtv2;
  typedef vnl_matrix<double>        MatrixType;
  MatrixType DT(3,3);
  DT.fill(0);
  DT(0,0)=dtv[0];
  DT(1,1)=dtv[3];
  DT(2,2)=dtv[5];
  DT(1,0)=DT(0,1)=dtv[1];
  DT(2,0)=DT(0,2)=dtv[2];
  DT(2,1)=DT(1,2)=dtv[4];
  //  if (takelog ) std::cout << " TAKING LOG " << std::endl;  else std::cout << "TAKING EXP " << std::endl;
  //  std::cout << " dtv " << dtv << std::endl;
  vnl_symmetric_eigensystem< double > eig(DT);


  double e1 = (eig.D(0,0));
  double e2 = (eig.D(1,1));
  double e3 = (eig.D(2,2));

  itk::RGBPixel< float > rgb;

  float xx = dtv[0];
  float xy = dtv[1];
  float xz = dtv[2];
  float yy = dtv[3];
  float yz = dtv[4];
  float zz = dtv[5];
  float isp = ( xx*xx + yy*yy + zz*zz + 2.0*(xy*xy + xz*xz + yz*yz) );

  float fa=0.0;
  if( isp > 0.0 )
    {
      float trace = dtv[0];
      trace += dtv[3];
      trace += dtv[5];
            float anisotropy = 3.0 * isp - trace * trace;
      fa=( vcl_sqrt(anisotropy / ( 2.0 * isp ) ) );
    }

  //rgb[0]=eig.V(2,0)*fa*255;//+eig.V(1,0)*e2;
  //rgb[1]=eig.V(2,1)*fa*255;//+eig.V(1,1)*e2;
  //  rgb[2]=eig.V(2,2)*fa*255;//+eig.V(1,2)*e2;

  // biggest evec
  rgb[0]=eig.V(0,2);//+eig.V(1,0)*e2;
  rgb[1]=eig.V(1,2);//+eig.V(1,1)*e2;
  rgb[2]=eig.V(2,2);//+eig.V(1,2)*e2;


  return rgb;
  mag=rgb[0]*rgb[0]+rgb[1]*rgb[1]+rgb[2]*rgb[2];

  mag=sqrt(mag);
  rgb[0]=rgb[0]/mag;
  rgb[1]=rgb[1]/mag;
  rgb[2]=rgb[2]/mag;

  return rgb;

}

template <class TTensorType>
TTensorType TensorLogAndExp( TTensorType dtv, bool takelog , bool &success)
{
  success = true;
  float eps=1.e-9,mag=0;
  for (unsigned int jj=0; jj<6; jj++)
  {
    float ff=dtv[jj];
    mag+=ff*ff;
    if ( vnl_math_isnan( ff ) || vnl_math_isinf(ff)   )
    {
      dtv.Fill(0); //dtv[0]=eps;   dtv[3]=eps;  dtv[5]=eps;
      return dtv;
    }
  }
  mag=sqrt(mag);

  if (  dtv[1]==0 && dtv[2] == 0 && dtv[4]==0) return dtv;
  if (mag < eps)
    {
    success = false; return dtv;
    }


  //  typedef itk::Vector<float, 6> TensorTypeIn;
  typedef TTensorType TensorType;
  typedef TTensorType TensorTypeOut;
  typedef TTensorType TensorTypeIn;
  const unsigned int ImageDimension = 3;
  typedef itk::Image<TensorTypeIn, ImageDimension> InTensorImageType;
  typedef itk::Image<TensorTypeOut, ImageDimension> OutTensorImageType;
  typedef itk::Image<TensorTypeOut, ImageDimension> TensorImageType;
  typedef     InTensorImageType::Pointer TensorImagePointer;
  typedef float  PixelType;
  typedef itk::Vector<float,ImageDimension>         VectorType1;
  typedef itk::Image<VectorType1,ImageDimension>     FieldType;
  typedef itk::Image<PixelType,ImageDimension> ImageType;
  typedef itk::ImageFileReader<InTensorImageType> readertype;
  typedef itk::ImageFileWriter<OutTensorImageType> writertype;
  typedef ImageType::IndexType IndexType;
  typedef ImageType::SizeType SizeType;
  typedef ImageType::SpacingType SpacingType;
//   typedef itk::AffineTransform<double,ImageDimension>   AffineTransformType;
//   typedef itk::LinearInterpolateImageFunction<ImageType,double>  InterpolatorType1;
//   typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double>  InterpolatorType2;

  itk::Vector<float, 6> dtv2;
  typedef vnl_matrix<double>        MatrixType;
  MatrixType DT(3,3);
  DT.fill(0);
  DT(0,0)=dtv[0];
  DT(1,1)=dtv[2];
  DT(2,2)=dtv[5];
  DT(1,0)=DT(0,1)=dtv[1];
  DT(2,0)=DT(0,2)=dtv[3];
  DT(2,1)=DT(1,2)=dtv[4];
  //  if (takelog ) std::cout << " TAKING LOG " << std::endl;  else std::cout << "TAKING EXP " << std::endl;
  //  std::cout << " dtv " << dtv << std::endl;
  vnl_symmetric_eigensystem< double > eig(DT);

  //  std::cout << " eig.D " << eig.D << std::endl;
  //  std::cout << " eig V " << eig.V << std::endl;

  double e1 = (eig.D(0,0));
  double e2 = (eig.D(1,1));
  double e3 = (eig.D(2,2));

  float peigeps=1.e-2;
  float eigeps=1.e-6;
  float eigmx=10000;
  if ( e3 < peigeps )
    {
    success=false;
    return dtv;
    }

  MatrixType eigmat(3,3);
  eigmat.fill(0);
      if (takelog)
    {
      eigmat(0,0)=log(e1);
      eigmat(1,1)=log(e2);
      eigmat(2,2)=log(e3);
    }
      else //take exp
    {
      eigmat(0,0)=exp(e1);
      eigmat(1,1)=exp(e2);
      eigmat(2,2)=exp(e3);
    }
      //      std::cout << " e1 " << e1 <<  " e2 " << e2 << " e3 " << std::endl;

      MatrixType DTrec=eig.V*eigmat*eig.V.transpose();

      dtv2[0]=DTrec(0,0);
      dtv2[2]=DTrec(1,1);
      dtv2[5]=DTrec(2,2);
      dtv2[1]=DTrec(0,1);
      dtv2[3]=DTrec(0,2);
      dtv2[4]=DTrec(1,2);

      return dtv2;

}


static float  GetMetricTensorCost(  itk::Vector<float, 3> dpath,  itk::Vector<float, 6> dtv )
{

  typedef vnl_matrix<double>        MatrixType;
  MatrixType DT(3,3);
  DT.fill(0);
  DT(0,0)=dtv[0];
  DT(1,1)=dtv[3];
  DT(2,2)=dtv[5];
  DT(1,0)=DT(0,1)=dtv[1];
  DT(2,0)=DT(0,2)=dtv[2];
  DT(2,1)=DT(1,2)=dtv[4];
  vnl_symmetric_eigensystem< double > eig(DT);
  double e1 = (eig.D(0,0));
  double e2 = (eig.D(1,1));
  double e3 = (eig.D(2,2));
  double etot=e1+e2+e3;
  if (etot==0) etot=1;

  MatrixType  vec(3,1);
  vec(0,0)=dpath[0];
  vec(1,0)=dpath[1];
  vec(2,0)=dpath[2];
  MatrixType inv =  vnl_matrix_inverse<double>(DT);

  MatrixType sol= vec.transpose()*inv*vec;
  float cost = sol(0,0)/etot;

  return cost;

}
*/

#endif
