/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "iMathFunctions.h"
#include "ReadWriteData.h"
#include "antsUtilities.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCannyEdgeDetectionImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkGrayscaleMorphologicalClosingImageFilter.h"
#include "itkGrayscaleMorphologicalOpeningImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkLaplacianSharpeningImageFilter.h"
#include "itkMultiScaleLaplacianBlobDetectorImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"


namespace ants
{

template <class ImageType>
typename ImageType::Pointer
BlobDetector( typename ImageType::Pointer image, unsigned int nBlobs,
              typename ImageType::Pointer itkNotUsed(image2),
              double itkNotUsed(corrThresh), double itkNotUsed(radius), double itkNotUsed(distanceThresh) )
{
  typedef float RealType;

  // FIXME - wtf goes on here

  // sensitive parameters are set here - begin
  //RealType     gradsig = 1.0;      // sigma for gradient filter
  unsigned int stepsperoctave = 10; // number of steps between doubling of scale
  RealType     minscale = vcl_pow( 1.0, 1.0 );
  RealType     maxscale = vcl_pow( 2.0, 10.0 );
  //RealType     uniqfeat_thresh = 0.01;
  //RealType     smallval = 1.e-2; // assumes images are normalizes in [ 0, 1 ]
  //bool         dosinkhorn = false;
  //RealType     maxradiusdiffallowed = 0.25; // IMPORTANT feature size difference
  //RealType     kneighborhoodval = 3;        // IMPORTANT - defines how many nhood nodes to use in k-hood definition
  //unsigned int radval = 20;                 // IMPORTANT radius for correlation
  //RealType     dthresh = 0.02;              // IMPORTANT distance preservation threshold
  // sensitive parameters are set here - end


  typedef itk::MultiScaleLaplacianBlobDetectorImageFilter<ImageType> BlobFilterType;
  typename BlobFilterType::Pointer blobFilter = BlobFilterType::New();
  typedef typename BlobFilterType::BlobPointer BlobPointer;
  blobFilter->SetStartT( minscale );
  blobFilter->SetEndT( maxscale );
  blobFilter->SetStepsPerOctave( stepsperoctave );
  blobFilter->SetNumberOfBlobs( nBlobs );
  blobFilter->SetInput( image );
  blobFilter->Update();

  //typedef typename BlobFilterType::BlobRadiusImageType BlobRadiusImageType;
  //typename BlobRadiusImageType::Pointer labimg = blobFilter->GetBlobRadiusImage();
  //typename BlobRadiusImageType::Pointer labimg2;

  return blobFilter->GetBlobRadiusImage();

  typedef typename BlobFilterType::BlobsListType BlobsListType;
  BlobsListType blobs1 =  blobFilter->GetBlobs();

  vnl_matrix<RealType> correspondencematrix1;
  correspondencematrix1.set_size( blobs1.size(), blobs1.size() );
  correspondencematrix1.fill( 1 );
  vnl_matrix<RealType> correspondencematrix2;
  vnl_matrix<RealType> correspondencematrix;

  BlobsListType blobs2;
/*
  if( fn2.length() > 3 )
    {
    ReadImage<ImageType>( image2, fn2.c_str() );
    typename BlobFilterType::Pointer blobFilter2 = BlobFilterType::New();
    blobFilter2->SetStartT( minscale );
    blobFilter2->SetEndT( maxscale );
    blobFilter2->SetStepsPerOctave( stepsperoctave );
    blobFilter2->SetNumberOfBlobs( nblobs );
    blobFilter2->SetInput( image2 );
    blobFilter2->Update();
    labimg2 = blobFilter2->GetBlobRadiusImage();
    WriteImage<BlobRadiusImageType>( labimg2, outname2.c_str() );
    labimg->FillBuffer( 0 );
    labimg2->FillBuffer( 0 );
    blobs2 =  blobFilter2->GetBlobs();
    correspondencematrix2.set_size( blobs2.size(), blobs2.size() );
    correspondencematrix2.fill( 1 );
    }
  else
    {
    return EXIT_SUCCESS;
    }
    */

/*

  // now compute some feature characteristics in each blob
  typedef typename ImageType::IndexType IndexType;
  IndexType   zeroind;  zeroind.Fill( radval );
  BlobPointer bestblob = ITK_NULLPTR;
  if( ( !blobs2.empty() ) && ( !blobs1.empty() ) )
    {
    getBlobCorrespondenceMatrix<ImageDimension, ImageType, BlobsListType>
      ( radval, image, image2, correspondencematrix, blobs1, blobs2, gradsig, dosinkhorn );
    // vnl_matrix<RealType> diagcorr = outer_product( correspondencematrix1.get_diagonal(),
    //  correspondencematrix2.get_diagonal() );
    // Sinkhorn<RealType>( diagcorr );
    // for ( unsigned int row = 0; row < correspondencematrix.rows(); row++ )
    //  for ( unsigned int col = 0; col < correspondencematrix.cols(); col++ )
    //	correspondencematrix( row, col ) *= diagcorr( row, col );
    //    Sinkhorn<RealType>( correspondencematrix );
    unsigned int matchpt = 1;
    std::cout << " now compute pairwise matching " << correspondencematrix.max_value() << " reducing to "
              << corrthresh << std::endl;
    unsigned int count1 = 0;
    typedef std::pair<BlobPointer, BlobPointer> BlobPairType;
    std::vector<BlobPairType> blobpairs;
    while( ( matchpt < ( corrthresh + 1 ) ) && ( count1 <  blobs1.size() ) )
      {
      unsigned int maxpair = correspondencematrix.arg_max();
      unsigned int maxrow = ( unsigned int )  maxpair / correspondencematrix.cols();
      unsigned int maxcol = maxpair - maxrow * correspondencematrix.cols();
      BlobPointer  blob1 = blobs1[maxrow];
      bestblob = blobs2[maxcol];
      if( bestblob &&  bestblob->GetObjectRadius() > 1 )
        {
        if( fabs( bestblob->GetObjectRadius() - blob1->GetObjectRadius() ) < maxradiusdiffallowed )
          {
          if( bestblob && ( image->GetPixel( blob1->GetCenter() ) > smallval )  &&
              ( image2->GetPixel( bestblob->GetCenter() )  > smallval ) &&
              ( correspondencematrix1(maxrow, maxrow) > uniqfeat_thresh ) &&
              ( correspondencematrix2(maxcol, maxcol) > uniqfeat_thresh  )
              )
            {
            BlobPairType blobpairing = std::make_pair( blob1, bestblob );
            blobpairs.push_back( blobpairing );
            std::cout << " best correlation " << correspondencematrix.absolute_value_max() << " rad1 "
                      << blob1->GetObjectRadius() << " rad2 " << bestblob->GetObjectRadius() << " : " << matchpt
                      << std::endl;
            labimg->SetPixel(     blob1->GetCenter(), matchpt ); // ( int ) ( 0.5 +   ( *i )->GetObjectRadius() ) );
            labimg2->SetPixel( bestblob->GetCenter(), matchpt ); // ( int ) ( 0.5 + bestblob->GetObjectRadius() ) );
            matchpt++;
            }
          }
        }
      correspondencematrix.set_row( maxrow, correspondencematrix.get_row( 0 ).fill( 0 ) );
      correspondencematrix.set_column( maxcol, correspondencematrix.get_column( 0 ).fill( 0 ) );
      count1++;
      }

    // For every blob, compute the distance to its neighbors before and after matching
    vnl_matrix<RealType> distmatpre( blobpairs.size(), blobpairs.size() );
    distmatpre.fill( 0 );
    vnl_matrix<RealType> distmatpost( blobpairs.size(), blobpairs.size() );
    distmatpost.fill( 0 );
    vnl_matrix<RealType> distratiomat( blobpairs.size(), blobpairs.size() );
    distratiomat.fill( 0 );
    if( true )
      {
      for( unsigned int bp = 0; bp < blobpairs.size(); bp++ )
        {
        IndexType             blobind = blobpairs[bp].first->GetCenter();
        IndexType             blobpairind  = blobpairs[bp].second->GetCenter();
        std::vector<RealType> distspre;
        std::vector<RealType> distspost;
        std::vector<size_t>   distspreind;
        std::vector<size_t>   distspostind;
        for( unsigned int bp2 = 0; bp2 < blobpairs.size(); bp2++ )
          {
          IndexType blobneighborind = blobpairs[bp2].first->GetCenter();
          IndexType blobpairneighborind  = blobpairs[bp2].second->GetCenter();
          RealType  dist1 = 0;
          RealType  dist2 = 0;
          for( unsigned int dim = 0; dim < ImageDimension; dim++ )
            {
            RealType delta1 = blobind[dim] - blobneighborind[dim];
            RealType delta2 = blobpairind[dim]  - blobpairneighborind[dim];
            dist1 += delta1 * delta1;
            dist2 += delta2 * delta2;
            }
          RealType drat = 0;
          if( dist1 > 0 )
            {
            drat = dist2 / dist1;
            }
          distspre.push_back( dist1 );
          distspost.push_back( dist2 );
          distspreind.push_back(  bp2 );
          distspostind.push_back( bp2 );
          //        // std::cout << " blob " << bp << " vs " << bp2 << sqrt( dist1 ) << " v " << sqrt( dist2 ) <<
          // std::endl;
          distmatpre(   bp, bp2 ) = distmatpre(   bp2, bp ) = dist1;
          distmatpost(  bp, bp2 ) = distmatpost(  bp2, bp ) = dist2;
          distratiomat( bp, bp2 ) = distratiomat( bp2, bp ) = drat;
          }
        }
      // now we have the distance ratio matrix - let's find a cluster of nodes with values near 1
      // count the k neighborhood for each blobpair possibility
      for( unsigned int bp = 0; bp < blobpairs.size(); bp++ )
        {
        IndexType    blobind = blobpairs[bp].first->GetCenter();
        IndexType    blobpairind  = blobpairs[bp].second->GetCenter();
        unsigned int kct = 0;
        for( unsigned int bp2 = 0; bp2 < blobpairs.size(); bp2++ )
          {
          if( ( bp2 != bp ) && ( vnl_math_abs( distratiomat( bp2, bp ) - 1 ) <  dthresh ) )
            {
            kct++;
            }
          }
        if( kct < kneighborhoodval )
          {
          labimg->SetPixel(  blobind, 0 );     // ( int ) ( 0.5 +   ( *i )->GetObjectRadius() ) );
          labimg2->SetPixel( blobpairind, 0 ); // ( int ) ( 0.5 + bestblob->GetObjectRadius() ) );
          }
        else
          {
          // std::cout << " blob " << bp << " keep " <<  distratiomat.get_row( bp ) << std::endl;
          }
        }
      } // if false
    if(  ( correspondencematrix1.rows() > 0 ) && ( false ) )
      {
      typedef itk::CSVNumericObjectFileWriter<RealType, 1, 1> WriterType;
      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( "corrmat1.csv" );
      writer->SetInput( &correspondencematrix1 );
      writer->SetInput( &distmatpre );
      writer->Write();
      }
    if(  ( correspondencematrix2.rows() > 0 ) && ( false ) )
      {
      typedef itk::CSVNumericObjectFileWriter<RealType, 1, 1> WriterType;
      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( "corrmat2.csv" );
      writer->SetInput( &correspondencematrix2 );
      writer->SetInput( &distratiomat );
      writer->Write();
      }
    WriteImage<BlobRadiusImageType>( labimg, outname.c_str() );
    WriteImage<BlobRadiusImageType>( labimg2, outname2.c_str() );
    // std::cout << " Matched " << matchpt << " blobs " << std::endl;
    }
  return EXIT_SUCCESS;
  */
}

template <class ImageType>
typename ImageType::Pointer
iMathCanny( typename ImageType::Pointer image,
            double sigma,
            double lowerThreshold,
            double upperThreshold )
{

  typedef typename ImageType::PixelType            PixelType;
  typedef itk::CannyEdgeDetectionImageFilter< ImageType, ImageType >  FilterType;

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image );
  filter->SetVariance( sigma );
  filter->SetUpperThreshold( (PixelType) upperThreshold );
  filter->SetLowerThreshold( (PixelType) lowerThreshold );
  filter->Update();

  return filter->GetOutput();

}

template <class ImageType>
typename ImageType::Pointer
iMathDistanceMap( typename ImageType::Pointer image, bool useSpacing )
{
  typedef itk::DanielssonDistanceMapImageFilter<ImageType, ImageType> FilterType;

  typename  FilterType::Pointer filter = FilterType::New();
  filter->InputIsBinaryOff();
  filter->SetUseImageSpacing(useSpacing);
  filter->SetInput(image);
  filter->Update();

  return filter->GetOutput();
}


// algorithm :
// 1. get distance map of object
// 2. threshold
// 3. label connected components
// 4. label surface
// 5. if everywhere on surface is next to object then it's a hole
// 6. make sure it's not the background
template <class ImageType>
typename ImageType::Pointer
iMathFillHoles( typename ImageType::Pointer image, double holeParam )
{

  if ( (holeParam < 0) || (holeParam > 2) )
    {
    //itk::itkExceptionMacro("FillHoles: holeParam value must lie in [0,2]");
    }

  typedef typename ImageType::Pointer                ImagePointerType;
  typedef itk::Image<int, ImageType::ImageDimension> MaskType;
  typedef typename ImageType::PixelType              PixelType;
  typedef typename MaskType::PixelType               LabelType;

  const PixelType imageMax = itk::NumericTraits<PixelType>::max();
  const LabelType labelMax = itk::NumericTraits<LabelType>::max();
  const PixelType objectMin = 0.5;
  const PixelType distanceMin = 0.001;

  typedef itk::CastImageFilter<MaskType,ImageType>            MaskToImage;
  typedef itk::BinaryThresholdImageFilter<ImageType,MaskType> ThresholdFilterType;
  typedef itk::BinaryThresholdImageFilter<MaskType,MaskType>  ThresholdMaskFilterType;

  typename ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();
  threshold->SetInput( image );
  threshold->SetInsideValue(1);
  threshold->SetOutsideValue(0);
  threshold->SetLowerThreshold(objectMin);
  threshold->SetUpperThreshold(imageMax);

  typedef itk::DanielssonDistanceMapImageFilter<MaskType, ImageType> FilterType;
  typename  FilterType::Pointer distance = FilterType::New();
  distance->InputIsBinaryOff();
  distance->SetUseImageSpacing(false);
  distance->SetInput(threshold->GetOutput());

  typename ThresholdFilterType::Pointer dThreshold = ThresholdFilterType::New();
  dThreshold->SetInput( distance->GetOutput() );
  dThreshold->SetInsideValue(1);
  dThreshold->SetOutsideValue(0);
  dThreshold->SetLowerThreshold(distanceMin);
  dThreshold->SetUpperThreshold(imageMax);
  dThreshold->Update();

  typedef itk::ConnectedComponentImageFilter<MaskType,MaskType> ConnectedFilterType;
  typename ConnectedFilterType::Pointer connected = ConnectedFilterType::New();
  connected->SetInput( dThreshold->GetOutput() );
  connected->SetFullyConnected( false );

  typedef itk::RelabelComponentImageFilter<MaskType, MaskType>  RelabelFilterType;
  typename RelabelFilterType::Pointer relabel = RelabelFilterType::New();
  relabel->SetInput( connected->GetOutput() );
  relabel->SetMinimumObjectSize( 0 );
  relabel->Update();

  if( holeParam == 2 )
    {
    typename ThresholdMaskFilterType::Pointer oThreshold = ThresholdMaskFilterType::New();
    oThreshold->SetInput( relabel->GetOutput() );
    oThreshold->SetInsideValue(1);
    oThreshold->SetOutsideValue(0);
    oThreshold->SetLowerThreshold(2);
    oThreshold->SetUpperThreshold(labelMax);

    typedef itk::AddImageFilter<MaskType,MaskType> AddFilterType;
    typename AddFilterType::Pointer add = AddFilterType::New();
    add->SetInput1( threshold->GetOutput() );
    add->SetInput2( oThreshold->GetOutput() );

    typename MaskToImage::Pointer maskToImage = MaskToImage::New();
    maskToImage->SetInput( add->GetOutput() );
    maskToImage->Update();

    return maskToImage->GetOutput();
    }

  // FIXME - add filter for below -- avoid iterators in these functions
  typename MaskToImage::Pointer caster = MaskToImage::New();
  caster->SetInput( threshold->GetOutput() );
  caster->Update();
  ImagePointerType imageout = caster->GetOutput();

  typedef itk::NeighborhoodIterator<MaskType> iteratorType;
  typename iteratorType::RadiusType rad;
  for( unsigned int j = 0; j < ImageType::ImageDimension; j++ )
    {
    rad[j] = 1;
    }
  iteratorType GHood(rad, relabel->GetOutput(), relabel->GetOutput()->GetLargestPossibleRegion() );

  float maximum = relabel->GetNumberOfObjects();
  // now we have the exact number of objects labeled independently
  for( int lab = 2; lab <= maximum; lab++ )
    {
    float erat = 2;
    if( holeParam <= 1 )
      {
      GHood.GoToBegin();

      unsigned long objectedge = 0;
      unsigned long backgroundedge = 0;
      unsigned long totaledge = 0;
      unsigned long volume = 0;

      while( !GHood.IsAtEnd() )
        {
        typename ImageType::PixelType p = GHood.GetCenterPixel();
        typename ImageType::IndexType ind2;
        if( p == lab )
          {
          volume++;
          for( unsigned int i = 0; i < GHood.Size(); i++ )
            {
            ind2 = GHood.GetIndex(i);
            float val2 = threshold->GetOutput()->GetPixel(ind2);
            if( (val2 == 1) && GHood.GetPixel(i) != lab )
              {
              objectedge++;
              totaledge++;
              }
            else if( (val2 == 1) && GHood.GetPixel(i) != lab )
              {
              backgroundedge++;
              totaledge++;
              }
            }
          }
        ++GHood;
        }

      erat = (float)objectedge / (float)totaledge;
      }

    if( erat > holeParam ) // fill the hole
      {
      // std::cout << " Filling " << lab << " of " << maximum <<  std::endl;
      typedef itk::ImageRegionIteratorWithIndex<MaskType> RelabelIterator;
      RelabelIterator vfIter( relabel->GetOutput(),
                              relabel->GetOutput()->GetLargestPossibleRegion() );
      for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
        {
        if( vfIter.Get() == lab )
          {
          imageout->SetPixel(vfIter.GetIndex(), 1);
          }
        }
      }
    }

  return imageout;
}


template <class ImageType>
typename ImageType::Pointer
iMathGC(typename ImageType::Pointer image, unsigned long radius)
{

  const unsigned int ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::PixelType            PixelType;

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension>
    StructuringElementType;

  typedef itk::GrayscaleMorphologicalClosingImageFilter< ImageType, ImageType, StructuringElementType >  FilterType;

  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image );
  filter->SetKernel( structuringElement );
  filter->Update();

  return filter->GetOutput();

}

template <class ImageType>
typename ImageType::Pointer
iMathGD(typename ImageType::Pointer image, unsigned long radius)
{

  const unsigned int ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::PixelType                         PixelType;

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension>
    StructuringElementType;

  typedef itk::GrayscaleDilateImageFilter< ImageType, ImageType, StructuringElementType >  FilterType;

  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image );
  filter->SetKernel( structuringElement );
  filter->Update();

  return filter->GetOutput();

}

template <class ImageType>
typename ImageType::Pointer
iMathGE( typename ImageType::Pointer image, unsigned long radius)
{
  const unsigned int ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::PixelType                         PixelType;

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension>
    StructuringElementType;

  typedef itk::GrayscaleErodeImageFilter< ImageType, ImageType, StructuringElementType >   FilterType;

  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image );
  filter->SetKernel( structuringElement );
  filter->Update();

  return filter->GetOutput();
}

template <class ImageType>
typename ImageType::Pointer
iMathGO( typename ImageType::Pointer image, unsigned long radius)
{
  const unsigned int ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::PixelType                         PixelType;

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension>
    StructuringElementType;

  typedef itk::GrayscaleMorphologicalOpeningImageFilter< ImageType, ImageType, StructuringElementType >  FilterType;

  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image );
  filter->SetKernel( structuringElement );
  filter->Update();

  return filter->GetOutput();
}

template <class ImageType>
typename ImageType::Pointer
iMathGetLargestComponent( typename ImageType::Pointer image,
                     unsigned long smallest )
{
  const unsigned int ImageDimension = ImageType::ImageDimension;

  if ( image->GetNumberOfComponentsPerPixel() != 1 )
    {
    // NOPE
    }

  typedef typename ImageType::PixelType                 PixelType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>  Iterator;

  // compute the voxel volume
  typename ImageType::SpacingType spacing = image->GetSpacing();
  float volumeelement = 1.0;
  for( unsigned int i = 0;  i < spacing.Size(); i++ )
    {
    volumeelement *= spacing[i];
    }

  typedef itk::Image<unsigned long, ImageDimension>                          LabelImageType;
  typedef itk::BinaryThresholdImageFilter<ImageType, LabelImageType>         ThresholdFilterType;
  typedef itk::ConnectedComponentImageFilter<LabelImageType, LabelImageType> FilterType;
  typedef itk::RelabelComponentImageFilter<LabelImageType, ImageType>        RelabelType;

  typename ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();
  typename FilterType::Pointer filter = FilterType::New();
  typename RelabelType::Pointer relabel = RelabelType::New();

  threshold->SetInput(image);
  threshold->SetInsideValue(1);
  threshold->SetOutsideValue(0);
  threshold->SetLowerThreshold(0.25);  //FIXME - why these values?
  threshold->SetUpperThreshold(1.e9);
  threshold->Update();

  filter->SetInput(threshold->GetOutput() );
  filter->SetFullyConnected( 0 );
  filter->Update();
  relabel->SetInput( filter->GetOutput() );
  relabel->SetMinimumObjectSize( smallest );
  //    relabel->SetUseHistograms(true);

  try
    {
    relabel->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    // std::cout << "Relabel: exception caught !" << std::endl;
    // std::cout << excep << std::endl;
    }

  //  WriteImage<ImageType>(relabel->GetOutput(),outname.c_str());
  //  return 0;
  typename ImageType::Pointer Clusters = MakeNewImage<ImageType>(relabel->GetOutput(), 0);
  // typename ImageType::Pointer Clusters=relabel->GetOutput();
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  Iterator vfIter( relabel->GetOutput(),  relabel->GetOutput()->GetLargestPossibleRegion() );

  float maximum = relabel->GetNumberOfObjects();
  float                     maxtstat = 0;
  std::vector<unsigned int> histogram( (int)maximum + 1);
  std::vector<float>        clustersum( (int)maximum + 1);
  for( int i = 0; i <= maximum; i++ )
    {
    histogram[i] = 0;
    clustersum[i] = 0;
    }
  for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    if( vfIter.Get() > 0 )
      {
      float vox = image->GetPixel(vfIter.GetIndex() );
      histogram[(unsigned int)vfIter.Get()] = histogram[(unsigned int)vfIter.Get()] + 1;
      clustersum[(unsigned int)vfIter.Get()] += vox;
      if( vox > maxtstat )
        {
        maxtstat = vox;
        }
      }
    }
  for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    if( vfIter.Get() > 0 )
      {
      Clusters->SetPixel( vfIter.GetIndex(), histogram[(unsigned int)vfIter.Get()]  );
      //  if ( Clusters->GetPixel( vfIter.GetIndex() ) > maximgval )
      //    maximgval=Clusters->GetPixel( vfIter.GetIndex());
      }
    else
      {
      Clusters->SetPixel(vfIter.GetIndex(), 0);
      }
    }

  float maximgval = 0;
  for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    if( Clusters->GetPixel( vfIter.GetIndex() ) > maximgval )
      {
      maximgval = Clusters->GetPixel( vfIter.GetIndex() );
      }
    }

  for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    if( Clusters->GetPixel( vfIter.GetIndex() ) >= maximgval )
      {
      image->SetPixel( vfIter.GetIndex(), 1);
      }
    else
      {
      image->SetPixel( vfIter.GetIndex(), 0);
      }
    }

  return image;
}

template <class ImageType>
typename ImageType::Pointer
iMathGrad(typename ImageType::Pointer image, double sigma, bool normalize )
{

  typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType,ImageType> FilterType;
  typename FilterType::Pointer grad = FilterType::New();
  grad->SetInput( image );
  grad->SetSigma( sigma );
  grad->Update();

  typename ImageType::Pointer output = grad->GetOutput();
  if ( normalize )
    {
    typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescaleFilterType;
    typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
    rescaler->SetOutputMinimum( 0 );
    rescaler->SetOutputMaximum( 1 );
    rescaler->SetInput( grad->GetOutput() );
    rescaler->Update();
    output = rescaler->GetOutput();
    }

  return output;
}

template <class ImageType>
typename ImageType::Pointer
iMathLaplacian(typename ImageType::Pointer image, double sigma, bool normalize )
{

  typedef itk::LaplacianRecursiveGaussianImageFilter<ImageType,ImageType> FilterType;
  typename FilterType::Pointer laplacian = FilterType::New();
  laplacian->SetInput( image );
  laplacian->SetSigma( sigma );
  laplacian->Update();

  typename ImageType::Pointer output = laplacian->GetOutput();
  if ( normalize )
    {
    typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescaleFilterType;
    typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
    rescaler->SetOutputMinimum( 0 );
    rescaler->SetOutputMaximum( 1 );
    rescaler->SetInput( laplacian->GetOutput() );
    rescaler->Update();
    output = rescaler->GetOutput();
    }

  return output;
}

template <class ImageType>
typename ImageType::Pointer
iMathMaurerDistance(typename ImageType::Pointer image,
                    typename ImageType::PixelType foreground )
{

  typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( image);
  thresholder->SetLowerThreshold( foreground );
  thresholder->SetUpperThreshold( foreground );
  thresholder->SetInsideValue( 1 );
  thresholder->SetOutsideValue( 0 );

  typedef itk::SignedMaurerDistanceMapImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( thresholder->GetOutput() );
  filter->SetSquaredDistance( false );
  filter->SetUseImageSpacing( true );
  filter->SetInsideIsPositive( false );
  filter->Update();

  return filter->GetOutput();
}

template <class ImageType>
typename ImageType::Pointer
iMathMC(typename ImageType::Pointer image, unsigned long radius,
        typename ImageType::PixelType closeValue)
{

  const unsigned int ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::PixelType                         PixelType;

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension>
    StructuringElementType;

  typedef itk::BinaryMorphologicalClosingImageFilter< ImageType, ImageType, StructuringElementType >  FilterType;

  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image );
  filter->SetKernel( structuringElement );
  filter->SetForegroundValue( closeValue );
  //filter->SetBackgroundValue(0);
  filter->Update();

  return filter->GetOutput();

}

template <class ImageType>
typename ImageType::Pointer
iMathMD(typename ImageType::Pointer image, unsigned long radius,
        typename ImageType::PixelType dilateValue)
{

  const unsigned int ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::PixelType                         PixelType;

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension>
    StructuringElementType;

  typedef itk::BinaryDilateImageFilter< ImageType, ImageType, StructuringElementType >  FilterType;

  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image );
  filter->SetKernel( structuringElement );
  filter->SetDilateValue( dilateValue );
  filter->SetBackgroundValue(0);
  filter->Update();

  return filter->GetOutput();

}

template <class ImageType>
typename ImageType::Pointer
iMathME( typename ImageType::Pointer image, unsigned long radius,
         typename ImageType::PixelType erodeValue )
{
  const unsigned int ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::PixelType                         PixelType;

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension>
    StructuringElementType;

  typedef itk::BinaryErodeImageFilter< ImageType, ImageType, StructuringElementType >   FilterType;

  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image );
  filter->SetKernel( structuringElement );
  filter->SetErodeValue( erodeValue );
  filter->SetBackgroundValue(0);
  filter->Update();

  return filter->GetOutput();
}

template <class ImageType>
typename ImageType::Pointer
iMathMO( typename ImageType::Pointer image, unsigned long radius,
         typename ImageType::PixelType openValue )
{
  const unsigned int ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::PixelType                         PixelType;

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension>
    StructuringElementType;

  typedef itk::BinaryMorphologicalOpeningImageFilter< ImageType, ImageType, StructuringElementType >  FilterType;

  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image );
  filter->SetKernel( structuringElement );
  filter->SetForegroundValue( openValue );
  filter->SetBackgroundValue( 0 );
  filter->Update();

  return filter->GetOutput();
}

template <class ImageType>
typename ImageType::Pointer
iMathNormalize( typename ImageType::Pointer image )
{
  if ( image->GetNumberOfComponentsPerPixel() != 1 )
    {
    // NOPE
    }

  typedef typename ImageType::PixelType                 PixelType;
  typedef typename ImageType::Pointer                   ImagePointerType;

  typedef itk::RescaleIntensityImageFilter<ImageType,ImageType> NormFilterType;
  typename NormFilterType::Pointer normFilter = NormFilterType::New();
  normFilter->SetInput( image );
  normFilter->SetOutputMinimum( itk::NumericTraits<PixelType>::ZeroValue() );
  normFilter->SetOutputMaximum( itk::NumericTraits<PixelType>::OneValue() );
  normFilter->Update();

  return normFilter->GetOutput();
}

template <class ImageType>
typename ImageType::Pointer
iMathPad( typename ImageType::Pointer image, int padding )
{
  typedef typename ImageType::PixelType                 PixelType;
  typedef typename ImageType::Pointer                   ImagePointerType;

  typename ImageType::PointType origin = image->GetOrigin();

  typename ImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
  for (unsigned int i=0; i<ImageType::ImageDimension; i++)
    {
    size[i] += 2*padding;
    origin[i] -= (padding * image->GetSpacing()[i]);
    }

  typedef itk::IdentityTransform<double,ImageType::ImageDimension> TransformType;
  typename TransformType::Pointer id = TransformType::New();
  id->SetIdentity();

  typedef itk::NearestNeighborInterpolateImageFunction<ImageType> InterpType;
  typename InterpType::Pointer interp = InterpType::New();

  typedef itk::ResampleImageFilter<ImageType,ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image );
  filter->SetOutputSpacing( image->GetSpacing() );
  filter->SetOutputOrigin( origin );
  filter->SetOutputDirection( image->GetDirection() );
  filter->SetDefaultPixelValue( 0 );
  filter->SetSize( size );
  filter->SetTransform( id );
  filter->SetInterpolator( interp );
  filter->Update();

  return filter->GetOutput();
}


template <class ImageType>
typename ImageType::Pointer
iMathPeronaMalik( typename ImageType::Pointer image, unsigned long nIterations,
  double conductance )
{
  if ( image->GetNumberOfComponentsPerPixel() != 1 )
    {
    // NOPE
    }

  typedef typename ImageType::PixelType                 PixelType;
  typedef typename ImageType::Pointer                   ImagePointerType;

  typedef itk::GradientAnisotropicDiffusionImageFilter< ImageType, ImageType >
    FilterType;
  typedef typename FilterType::TimeStepType             TimeStepType;

  // Select time step size.
  TimeStepType  spacingsize = 0;
  for( unsigned int d = 0; d < ImageType::ImageDimension; d++ )
    {
    TimeStepType sp = image->GetSpacing()[d];
    spacingsize += sp * sp;
    }
  spacingsize = sqrt( spacingsize );

  // FIXME - cite reason for this step
  double dimPlusOne = ImageType::ImageDimension + 1;
  TimeStepType mytimestep = spacingsize / vcl_pow( 2.0 , dimPlusOne );
  TimeStepType reftimestep = 0.4 / vcl_pow( 2.0 , dimPlusOne );
  if ( mytimestep > reftimestep )
    {
    mytimestep = reftimestep;
    }

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image );
  filter->SetConductanceParameter( conductance ); // might need to change this
  filter->SetNumberOfIterations( nIterations );
  filter->SetTimeStep( mytimestep );

  filter->Update();
  return filter->GetOutput();
}

template <class ImageType>
typename ImageType::Pointer
iMathSharpen( typename ImageType::Pointer image )
{
  if ( image->GetNumberOfComponentsPerPixel() != 1 )
    {
    // NOPE
    }

  typedef typename ImageType::PixelType                 PixelType;
  typedef typename ImageType::Pointer                   ImagePointerType;

  typedef itk::LaplacianSharpeningImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer sharpenFilter = FilterType::New();
  sharpenFilter->SetInput( image );
  sharpenFilter->Update();

  return sharpenFilter->GetOutput();
}

template <class ImageType>
typename ImageType::Pointer
iMathTruncateIntensity( typename ImageType::Pointer image, double lowerQ, double upperQ, int nBins,
                        typename itk::Image<unsigned int, ImageType::ImageDimension>::Pointer mask )
{

  typedef typename ImageType::PixelType                     PixelType;
  typedef typename ImageType::Pointer                       ImagePointerType;
  typedef unsigned int                                      LabelType;
  typedef itk::Image<LabelType, ImageType::ImageDimension>  MaskType;

  if( mask.IsNull() )
    {
    typedef itk::BinaryThresholdImageFilter<ImageType, MaskType> ThresholdFilterType;
    typename ThresholdFilterType::Pointer thresh = ThresholdFilterType::New();
    thresh->SetInput( image );
    thresh->SetLowerThreshold( 1e-6 );
    thresh->SetUpperThreshold( itk::NumericTraits<PixelType>::max() );
    thresh->SetInsideValue(1);
    thresh->SetOutsideValue(0);
    thresh->Update();
    mask = thresh->GetOutput();
    }
  typedef itk::LabelStatisticsImageFilter<ImageType, MaskType> HistogramFilterType;
  typename HistogramFilterType::Pointer stats = HistogramFilterType::New();

  stats->SetInput( image );
  stats->SetLabelInput( mask );
  stats->Update();
  PixelType minValue = stats->GetMinimum( 1 );
  PixelType maxValue = stats->GetMaximum( 1 );

  // Hack increment by delta
  if (minValue == 0)
    {
    minValue = (PixelType) (minValue + 1e-6);
    }
  if (minValue == 0)
    {
    minValue++;
    }

  stats->SetUseHistograms( true );
  stats->SetHistogramParameters( nBins, minValue, maxValue );
  stats->Update();

  typedef typename HistogramFilterType::HistogramPointer HistogramPointer;
  HistogramPointer histogram = stats->GetHistogram( 1 );

  PixelType lowerQuantile = histogram->Quantile( 0, lowerQ );
  PixelType upperQuantile = histogram->Quantile( 0, upperQ );

  typedef itk::IntensityWindowingImageFilter<ImageType,ImageType> WindowFilterType;
  typename WindowFilterType::Pointer windowFilter = WindowFilterType::New();
  windowFilter->SetInput( image );
  windowFilter->SetWindowMinimum( lowerQuantile );
  windowFilter->SetOutputMinimum( lowerQuantile );
  windowFilter->SetWindowMaximum( upperQuantile );
  windowFilter->SetOutputMaximum( upperQuantile );
  windowFilter->Update();

  return windowFilter->GetOutput();

}


} // namespace ants
