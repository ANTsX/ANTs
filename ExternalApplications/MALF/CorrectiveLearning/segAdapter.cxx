/*
segAdaptor: this algorithm corrects segmentation errors produced by a host segmentation method.
The learning is performed by the program, bl (BiasLearn.cxx)

Hongzhi Wang, 07/18/2010
*/

#include "util.h"
#include "AdaBoost.h"
using namespace std;

int main( int argc, char * * argv )
{
  if( argc < 5 )
    {
    std::cerr << "Missing parameters. " << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << " inputImage inputSegmentation AdaBoostOutPutPrefix outputSegmentation [-x label image.nii]" << endl
              << endl
              << " Apply corrections to segmentations produced by the host segmentation method." << endl << endl
              << " Meanings of the parameters:" << endl
              << " inputImage:           the image file" << endl
              << " inputSegmentation:    the segmentation produced by the host segmentation method for the image"
              << endl
              << " AdaBoostOutPutPrefix: the path and prefix to the learned AdaBoost files used in ./bl" << endl
              << " outputSegmentation:   the path and file name of the corrected segmentation." << endl
              << " -x label image.nii    Specify an exclusion region for the given label. " << endl
              << "                       If a voxel has non-zero value in an exclusion image," << endl
              << "                       the corresponding label is not allowed at that voxel." << endl
              << std::endl;
    return -1;
    }
  WriterType::Pointer writer = WriterType::New();

  int                 r, c, d;
  ReaderType::Pointer autoseg = ReaderType::New();
  ReaderType::Pointer im = ReaderType::New();
  ReaderType::Pointer exclusionMask = ReaderType::New();
  im->SetFileName( argv[1] );
  autoseg->SetFileName( argv[2] );
  try
    {
    autoseg->Update();
    im->Update();
    r = im->GetOutput()->GetRequestedRegion().GetSize()[0];
    c = im->GetOutput()->GetRequestedRegion().GetSize()[1];
    d = im->GetOutput()->GetRequestedRegion().GetSize()[2];
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return -1;
    }
  cout << r << ":" << c << ":" << d << endl;
  ImageType::Pointer nseg = ImageType::New();
  ImageType::Pointer cseg = ImageType::New();
  ImageType::Pointer emask = ImageType::New();

  cseg->SetRegions(autoseg->GetOutput()->GetRequestedRegion() );
  cseg->SetSpacing( autoseg->GetOutput()->GetSpacing() );
  cseg->SetOrigin( autoseg->GetOutput()->GetOrigin() );
  cseg->SetDirection(autoseg->GetOutput()->GetDirection() );
  cseg->Allocate();

  nseg->SetRegions(autoseg->GetOutput()->GetRequestedRegion() );
  nseg->SetSpacing( autoseg->GetOutput()->GetSpacing() );
  nseg->SetOrigin( autoseg->GetOutput()->GetOrigin() );
  nseg->SetDirection(autoseg->GetOutput()->GetDirection() );
  nseg->Allocate();

  emask->SetRegions(autoseg->GetOutput()->GetRequestedRegion() );
  emask->SetSpacing( autoseg->GetOutput()->GetSpacing() );
  emask->SetOrigin( autoseg->GetOutput()->GetOrigin() );
  emask->SetDirection(autoseg->GetOutput()->GetDirection() );
  emask->Allocate();

  ImageType::Pointer mask = ImageType::New();
  ImageType::Pointer segmask = ImageType::New();
  ImageType::Pointer nim = ImageType::New();

  nim->SetRegions(im->GetOutput()->GetRequestedRegion() );
  nim->SetSpacing(im->GetOutput()->GetSpacing() );
  nim->SetOrigin(im->GetOutput()->GetOrigin() );
  nim->SetDirection(im->GetOutput()->GetDirection() );
  nim->Allocate();

  segmask->SetRegions(autoseg->GetOutput()->GetRequestedRegion() );
  segmask->SetSpacing( autoseg->GetOutput()->GetSpacing() );
  segmask->SetOrigin( autoseg->GetOutput()->GetOrigin() );
  segmask->SetDirection(autoseg->GetOutput()->GetDirection() );
  segmask->Allocate();

  mask->SetRegions(autoseg->GetOutput()->GetRequestedRegion() );
  mask->SetSpacing( autoseg->GetOutput()->GetSpacing() );
  mask->SetOrigin( autoseg->GetOutput()->GetOrigin() );
  mask->SetDirection(autoseg->GetOutput()->GetDirection() );
  mask->Allocate();

  IteratorType imit0(im->GetOutput(), im->GetOutput()->GetRequestedRegion() );
  IteratorType nimit0(nim, nim->GetRequestedRegion() );
  IteratorType maskit0(mask, mask->GetRequestedRegion() );
  IteratorType emaskit0(emask, emask->GetRequestedRegion() );
  IteratorType autoit0(autoseg->GetOutput(), autoseg->GetOutput()->GetRequestedRegion() );

  map<int, string> fnExclusion;
  int              j, k;
  for( j = 5; j < argc - 2; j++ )
    {
    string arg = argv[j];
    if( arg == "-x" )
      {
      int    label = atoi(argv[++j]);
      string image = argv[++j];
      fnExclusion[label] = image;
      cout << "excluded: " << label << " " << image << endl;
      }
    }
  IndexIteratorType asegit(autoseg->GetOutput(), autoseg->GetOutput()->GetRequestedRegion() );
  IndexIteratorType csegit(cseg, cseg->GetRequestedRegion() );
  IndexIteratorType nsegit(nseg, nseg->GetRequestedRegion() );
  IndexIteratorType segmaskit(segmask, segmask->GetRequestedRegion() );

  char tfn[1024];
  sprintf(tfn, "%s-AdaBoostResults-param-Tlabel0", argv[3]);
  double   initt;
    {
    ifstream initt_ifs( tfn, ifstream::in );
    if( !initt_ifs.good() )
      {
      initt = 0.5;
      }
    else
      {
      initt = 0;
      initt_ifs.close();
      }
    }
  for( csegit.GoToBegin(), nsegit.GoToBegin(), asegit.GoToBegin(); !csegit.IsAtEnd(); ++csegit, ++nsegit, ++asegit )
    {
    csegit.Set(initt);
    nsegit.Set(0);
    }

  int ML = 0;;
  for( autoit0.GoToBegin(); !autoit0.IsAtEnd(); ++autoit0 )
    {
    if( autoit0.Value() > ML )
      {
      ML = int(autoit0.Value() );
      }
    }
  cout << "label #: " << ML << endl;

  int LC;

  double t;
  int    DX, DY, DZ, MD, dilateR;
  time_t second1 = time(NULL), second2;
  for( int Tlabel = 0; Tlabel <= ML; Tlabel++ )
    {
    if( fnExclusion[Tlabel] != "" )
      {
      exclusionMask->SetFileName( fnExclusion[Tlabel] );
      try
        {
        exclusionMask->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        return -1;
        }
      cout << "using exclusion mask: " << fnExclusion[Tlabel] << endl;
      IteratorType it0(exclusionMask->GetOutput(), exclusionMask->GetOutput()->GetRequestedRegion() );
      for( emaskit0.GoToBegin(), it0.GoToBegin(); !emaskit0.IsAtEnd(); ++emaskit0, ++it0 )
        {
        emaskit0.Set(it0.Value() );
        }
      }
    else
      {
      for( emaskit0.GoToBegin(); !emaskit0.IsAtEnd(); ++emaskit0 )
        {
        emaskit0.Set(0);
        }
      }
    IteratorType exclusionMaskit(emask, emask->GetRequestedRegion() );

    vector<int>    sign;
    vector<double> alpha;
    vector<double> threshold;
    vector<int>    featureID;

    // load AdaBoost learning parameters
    cout << "Tlabel  " << Tlabel << endl;
    sprintf(tfn, "%s-AdaBoostResults-param-Tlabel%d", argv[3], Tlabel);
    ifstream ifs( tfn, ifstream::in );
    if( !ifs.good() )
      {
      continue;
      }
    cout << tfn << endl;
    ifs >> dilateR >> DX >> DY >> DZ;
    ifs.close();
    sprintf(tfn, "%s-AdaBoostResults-Tlabel%d", argv[3], Tlabel);
    ifs.open( tfn, ifstream::in );
    if( !ifs.good() )
      {
      continue;
      }
    cout << tfn << endl;

    while( ifs.good() )
      {
      ifs >> t;
      if( !ifs.good() )
        {
        break;
        }
      ifs >> t;
      alpha.push_back(t);
      ifs >> t;
      ifs >> t;
      featureID.push_back(int(t) );
      ifs >> t;
      sign.push_back(int(t) );
      ifs >> t;
      threshold.push_back(t);
      ifs >> t;
      }
    ifs.close();
    LC = alpha.size();
    cout << "#weak learners: " << LC << endl;

    // compute features
    int NFeature = (2 * (DX * 2 + 1) * (DY * 2 + 1) * (DZ * 2 + 1) + 3) * 4 - 3;
    int tc = 0;;

    cout << "preparing data..." << endl;
    myDilate(autoseg->GetOutput(), mask, -1, dilateR);
    myDilate(autoseg->GetOutput(), segmask, Tlabel, dilateR);

    double mim = 0;
    for( imit0.GoToBegin(), maskit0.GoToBegin(); !imit0.IsAtEnd(); ++imit0, ++maskit0 )
      {
      if( maskit0.Value() > 0 )
        {
        mim += imit0.Value();
        tc++;
        }
      }
    mim /= tc;
    for( imit0.GoToBegin(), nimit0.GoToBegin(); !imit0.IsAtEnd(); ++imit0, ++nimit0 )
      {
      nimit0.Set(imit0.Value() - mim);
      }
    cout << "mean intensity : sample #" << mim << " : " << tc << endl;

    MD = DX;
    if( DY > MD )
      {
      MD = DY;
      }
    if( DZ > MD )
      {
      MD = DZ;
      }
    NeighborhoodIteratorType::RadiusType radius;
    radius.Fill(MD);

    tc = 0;
    double mx = 0, my = 0, mz = 0;
    for( segmaskit.GoToBegin(); !segmaskit.IsAtEnd(); ++segmaskit )
      {
      ImageType::IndexType idx = segmaskit.GetIndex();
      if( segmask->GetPixel(idx) > 0 && (mask->GetPixel(idx) > 0) )
        {
        mx += idx[0];
        my += idx[1];
        mz += idx[2];
        tc++;
        }
      }
    if( tc < 1 )
      {
      continue;
      }
    mx = mx / tc;
    my = my / tc;
    mz = mz / tc;
    cout << "center location:(" << mx << "," << my << "," << mz << ") feature #:" << NFeature << endl;

    NeighborhoodIteratorType imnit( radius, nim, nim->GetRequestedRegion() );
    NeighborhoodIteratorType masknit( radius, mask, mask->GetRequestedRegion() );

    cout << "making corrections..." << endl;
    int                      l, C;
    double                   H, cH;
    NeighborhoodIteratorType autosegnit(radius, autoseg->GetOutput(), autoseg->GetOutput()->GetRequestedRegion() );
    double *                 X = new double[NFeature];
    for( exclusionMaskit.GoToBegin(), segmaskit.GoToBegin(), csegit.GoToBegin(), nsegit.GoToBegin();
         !segmaskit.IsAtEnd(); ++segmaskit, ++csegit, ++nsegit, ++exclusionMaskit )
      {
      ImageType::IndexType idx = segmaskit.GetIndex();
      // excluding voxels which are too close to the edge of the image from consideration
      if( idx[0] < DY || idx[0] >= r - DY || idx[1] < DX || idx[1] >= c - DX || idx[2] < DZ || idx[2] >= d - DZ )
        {
        continue;
        }

      if( exclusionMaskit.Value() == 0 && mask->GetPixel(idx) > 0 && segmask->GetPixel(idx) > 0 && cseg->GetPixel(idx) <
          1 )
        {
        imnit.SetLocation(idx);
        autosegnit.SetLocation(idx);

        X[0] = idx[1] - my;
        X[1] = idx[0] - mx;
        X[2] = idx[2] - mz;
        C = 3;
        for( j = -DX; j < DX + 1; j++ )
          {
          for( k = -DY; k < DY + 1; k++ )
            {
            for( l = -DZ; l < DZ + 1; l++ )
              {
              X[C++] = imnit.GetPixel( (NeighborhoodIteratorType::OffsetType) {{k, j, l}});
              X[C++] = autosegnit.GetPixel( (NeighborhoodIteratorType::OffsetType) {{k, j, l}});
//		  cout<<X[C]<<" ";
              }
            }
          }
        tc = C;
        for( j = 0; j < 3; j++ )
          {
          for( k = j; k < C; k++ )
            {
            X[tc++] = X[j] * X[k];
            }
          }

        // Apply AdaBoost classifier on the features to obtain classifier response in probability, stored in cH.
        AdaBoostClassify(X, 1, NFeature, LC, &featureID[0], &alpha[0], &sign[0], &threshold[0], &H, &cH);
        if( cseg->GetPixel(idx) < cH )
          {
          csegit.Set(cH);
          nsegit.Set(Tlabel);
          }
        }
      }
    delete [] X;
    second2 = time(NULL);
    cout << "Label " << Tlabel << " is processed in " << (second2 - second1) << " seconds." << endl;
    second1 = second2;
    }

  cout << "saving corrected segmentation to " << argv[4] << endl;
  // write the corrected segmentation into the output file
  writer->SetFileName( argv[4] );
  writer->SetInput(nseg);
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return -1;
    }

  return 0;
}
