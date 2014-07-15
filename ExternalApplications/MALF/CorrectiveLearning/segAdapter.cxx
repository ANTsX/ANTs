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
  if( argc < 4 )
    {
    std::cerr << "segAdapter: Apply corrections to segmentations produced by the host segmentation method. "
              << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << " input_Segmentation AdaBoost_OutPutPrefix output_Segmentation [options]" << endl << endl
              << " Meanings of the parameters:" << endl
              << " input_Segmentation:    Segmentation produced by the host segmentation method for the image" << endl
              << " AdaBoost_OutPutPrefix: Path and prefix to the AdaBoost files specified in ./bl" << endl
              << " output_Segmentation:   Path and file name of the output corrected segmentation." << endl
              << " options: " << endl
              << "  -f feature1.nii feature2.nii ... " << endl
              << "                        Feature images for the target subject. " << endl
              << "                        The name pattern could be in C printf format, e.g. feature%04d.nii" << endl
              <<
      "                        Then feature0000.nii will be used for label 0 and feature0001.nii for label 1, etc."
              << endl
              <<
      "  -m mask:              Specify ROI for the training images. Should be in C printf format, e.g. mask%04d.nii"
              << endl
              << "                        ROI will be derived by performing dilation on this mask." << endl
              << "  -x label image.nii    Specify an exclusion region for the given label. " << endl
              << "                        If a voxel has non-zero value in an exclusion image," << endl
              << "                        the corresponding label is not allowed at that voxel." << endl
              << "  -p filenamePattern    Save the posterior maps (probability that each " << endl
              << "                        voxel belongs to each label) as images. The number of " << endl
              << "                        images saved equals the number of labels. The filename" << endl
              << "                        pattern must be in C printf format, e.g. posterior%04d.nii.gz" << endl
              << "  -mrf <method> [parameters]  Apply Markov Random Field prior to derive the segmentation." << endl
              << "                              Options: ICM (Iterated Conditional Modes) " << endl
              <<
      "                              May be followed by optional parameters in brackets, e.g., -mrf ICM[beta,iter]."
              << endl
              <<
      "                              beta: weight for the MRF prior, must be a non-negative number. Default: 0.1"
              << endl
              << "                              iter: max iteration for ICM optimization. Default: 10" << endl
              << std::endl;
    return -1;
    }
  WriterType::Pointer writer = WriterType::New();

  float               ICM_beta = 0.1;
  int                 ICM_iter = 10;
  int                 MRFFlag = 0;
  int                 r = 0, c = 0, d = 0, maskFlag = 0, j, k;
  ReaderType::Pointer autoseg = ReaderType::New();
  ReaderType::Pointer exclusionMask = ReaderType::New();
  ReaderType::Pointer ROIMask = ReaderType::New();

  typedef itk::Image<float, 3>               PosteriorImage;
  typedef PosteriorImage::Pointer            PosteriorImagePtr;
  typedef std::map<float, PosteriorImagePtr> PosteriorMap;

  PosteriorMap m_PosteriorMap;
  // Initialize the posterior maps
  m_PosteriorMap.clear();

  string           maskFn;
  string           posteriorFn = "";
  vector<string>   featurefn;
  vector<int>      featurefnflag;
  map<int, string> fnExclusion;
  char             buffer[4096];
  cout << argv[1] << endl;
  autoseg->SetFileName( argv[1] );
  try
    {
    autoseg->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return -1;
    }
  r = autoseg->GetOutput()->GetRequestedRegion().GetSize()[0];
  c = autoseg->GetOutput()->GetRequestedRegion().GetSize()[1];
  d = autoseg->GetOutput()->GetRequestedRegion().GetSize()[2];

  cout << r << ":" << c << ":" << d << endl;
  for( j = 4; j < argc; j++ )
    {
    string arg = argv[j];

    if( arg == "-f" )
      {
      featurefn.push_back(argv[++j]);
      while( j + 1 < argc && argv[j + 1][0] != '-' )
        {
        featurefn.push_back(argv[++j]);
        }
      }
    else if( arg == "-x" )
      {
      int    label = atoi(argv[++j]);
      string image = argv[++j];
      fnExclusion[label] = image;
      cout << "excluded: " << label << " " << image << endl;
      }
    else if( arg == "-m" )
      {
      maskFn = argv[++j];
      maskFlag = 1;
      }
    else if( arg == "-p" )
      {
      posteriorFn = argv[++j];

      }
    else if( arg == "-mrf" && j < argc - 1 )
      {
      char *parm = argv[++j];
      if( !strncmp(parm, "ICM", 3) )
        {
        sscanf(parm, "ICM[%f,%d]", &ICM_beta, &ICM_iter);
        cout << "ICM parameters: " << ICM_beta << " , " << ICM_iter << endl;
        if( ICM_beta < 0 || ICM_iter < 0 )
          {
          cout << "the ICM parameters cannot be negative!..";
          return -1;
          }
        MRFFlag = 1;
        }
      }
    }
  cout << "feature #:" << featurefn.size() << endl;
  cout << "maskFlag:" << maskFlag << "  " << maskFn << endl;

  ImageType::Pointer nseg = ImageType::New();
  ImageType::Pointer cseg = ImageType::New();
  ImageType::Pointer emask = ImageType::New();
  ImageType::Pointer posterior = ImageType::New();

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

  posterior->SetRegions(autoseg->GetOutput()->GetRequestedRegion() );
  posterior->SetSpacing( autoseg->GetOutput()->GetSpacing() );
  posterior->SetOrigin( autoseg->GetOutput()->GetOrigin() );
  posterior->SetDirection(autoseg->GetOutput()->GetDirection() );
  posterior->Allocate();

  ImageType::Pointer mask = ImageType::New();
  ImageType::Pointer segmask = ImageType::New();

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

//    IteratorType imit0(im->GetOutput(), im->GetOutput()->GetRequestedRegion());
  IteratorType maskit0(mask, mask->GetRequestedRegion() );
  IteratorType emaskit0(emask, emask->GetRequestedRegion() );
  IteratorType autoit0(autoseg->GetOutput(), autoseg->GetOutput()->GetRequestedRegion() );

  IndexIteratorType asegit(autoseg->GetOutput(), autoseg->GetOutput()->GetRequestedRegion() );
  IndexIteratorType csegit(cseg, cseg->GetRequestedRegion() );
  IndexIteratorType nsegit(nseg, nseg->GetRequestedRegion() );
  IndexIteratorType segmaskit(segmask, segmask->GetRequestedRegion() );
//    IndexIteratorType posteriorit(posterior, posterior->GetRequestedRegion());

  char tfn[1024];
  sprintf(tfn, "%s-AdaBoostResults-param-Tlabel0", argv[2]);
  ifstream ifs( tfn, ifstream::in );
  double   initt;
  if( !ifs.good() )
    {
    initt = 0.5;
    }
  else
    {
    initt = 0;
    ifs.close();
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

  double      t;
  int         DX, DY, DZ, MD, dilateR, featureChannel;
  vector<int> labelset;
  time_t      second1 = time(NULL), second2;
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
    sprintf(tfn, "%s-AdaBoostResults-param-Tlabel%d", argv[2], Tlabel);
    ifs.open( tfn, ifstream::in );
    if( !ifs.good() )
      {
      continue;
      }
    cout << tfn << endl;
    ifs >> dilateR >> DX >> DY >> DZ >> featureChannel;
    ifs.close();
    if( (size_t)featureChannel != featurefn.size() )
      {
      cout << "feature images donot match training data!" << endl
           << "feature images used in training:" << featureChannel << endl
           << "feature images used now        :" << featurefn.size() << endl;
      }
    sprintf(tfn, "%s-AdaBoostResults-Tlabel%d", argv[2], Tlabel);
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
    labelset.push_back(Tlabel);
    std::vector<ImageType::Pointer> ims;
    for( size_t i0 = 0; i0 < featurefn.size(); i0++ )
      {
      std::vector<int>::iterator it;
      sprintf(buffer, featurefn[i0].data(), Tlabel);
      cout << "loading " << buffer << endl;
      // Read the following options as images
      ReaderType::Pointer im = ReaderType::New();
      im->SetFileName( buffer );
      try
        {
        im->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        return -1;
        }
      ims.push_back(im->GetOutput() );
      }
    cout << ims.size() << endl;
    if( maskFlag )
      {
      sprintf(buffer, maskFn.data(), Tlabel);
      cout << buffer << endl;
      ROIMask->SetFileName( buffer );
      try
        {
        ROIMask->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        return -1;
        }
      }

    m_PosteriorMap[Tlabel] = PosteriorImage::New();
    m_PosteriorMap[Tlabel]->SetRegions(autoseg->GetOutput()->GetRequestedRegion() );
    m_PosteriorMap[Tlabel]->SetRequestedRegion(autoseg->GetOutput()->GetRequestedRegion() );
    m_PosteriorMap[Tlabel]->SetBufferedRegion(autoseg->GetOutput()->GetRequestedRegion() );
    m_PosteriorMap[Tlabel]->SetSpacing( autoseg->GetOutput()->GetSpacing() );
    m_PosteriorMap[Tlabel]->SetOrigin( autoseg->GetOutput()->GetOrigin() );
    m_PosteriorMap[Tlabel]->SetDirection(autoseg->GetOutput()->GetDirection() );
    m_PosteriorMap[Tlabel]->Allocate();
    m_PosteriorMap[Tlabel]->FillBuffer(0.0f);

    IndexIteratorType posteriorit(m_PosteriorMap[Tlabel], m_PosteriorMap[Tlabel]->GetRequestedRegion() );

    // compute features
    int NFeature;
    NFeature = ( (1 + ims.size() ) * (DX * 2 + 1) * (DY * 2 + 1) * (DZ * 2 + 1) + 3) * 4 - 3;
    cout << "NFeature: " << NFeature << endl;

    cout << "preparing data..." << endl;
    IndexIteratorType imaskit(mask, mask->GetRequestedRegion() );
    IndexIteratorType autosegit(autoseg->GetOutput(), autoseg->GetOutput()->GetRequestedRegion() );

    if( maskFlag )
      {
      myDilate(ROIMask->GetOutput(), segmask, -1, dilateR);
      myDilate(autoseg->GetOutput(), mask, -1, dilateR);
      }
    else
      {
      cout << "dilateR: " << dilateR << endl;
      // The dilated ROI for the specific label. Note that this label can be the background label, which is 0.
      myDilate(autoseg->GetOutput(), segmask, Tlabel, dilateR);
      myDilate(autoseg->GetOutput(), mask, -1, dilateR);
      }

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
    std::vector<IteratorType>             imits;
    std::vector<NeighborhoodIteratorType> imnits;
    for( size_t it = 0; it < ims.size(); it++ )
      {
      double                   mim = 0;
      int                      tc = 0;
      IteratorType             imit0(ims[it], ims[it]->GetRequestedRegion() );
      NeighborhoodIteratorType imnit(radius, ims[it], ims[it]->GetRequestedRegion() );
      imits.push_back(imit0);
      imnits.push_back(imnit);
      // normalize image intensity by the mean intensity of the whole figure region
      for( imit0.GoToBegin(), maskit0.GoToBegin(); !imit0.IsAtEnd(); ++imit0, ++maskit0 )
        {
        if( maskit0.Value() > 0 )
          {
          mim += imit0.Value();
          tc++;
          }
        }
      mim /= tc;
      cout << "mean intensity : sample #" << mim << ":" << tc << endl;
      for( imit0.GoToBegin(); !imit0.IsAtEnd(); ++imit0 )
        {
        imit0.Set(imit0.Value() - mim);
        }
      }

    int    tc = 0;
    double mx = 0, my = 0, mz = 0, totalval = 0;
    for( segmaskit.GoToBegin(), autosegit.GoToBegin(); !segmaskit.IsAtEnd(); ++segmaskit, ++autosegit )
      {
      ImageType::IndexType idx = segmaskit.GetIndex();
      // excluding voxels which are too close to the edge of the image from consideration
      if( idx[0] < DX || idx[0] >= r - DX || idx[1] < DY || idx[1] >= c - DY || idx[2] < DZ || idx[2] >= d - DZ )
        {
        continue;
        }
      if( segmask->GetPixel(idx) > 0 && (mask->GetPixel(idx) > 0) )
        {
        mx += idx[0];
        my += idx[1];
        mz += idx[2];
        totalval += 1;
        tc++;
        }
      }
    if( tc < 1 )
      {
      continue;
      }
    mx = mx / totalval;
    my = my / totalval;
    mz = mz / totalval;
    cout << "center location: (" << mx << "," << my << "," << mz << ") feature #:" << NFeature << endl;

    NeighborhoodIteratorType masknit( radius, mask, mask->GetRequestedRegion() );
    NeighborhoodIteratorType autosegnit(radius, autoseg->GetOutput(), autoseg->GetOutput()->GetRequestedRegion() );

//        if (posteriorFn.size())
      {
      for( posteriorit.GoToBegin(); !posteriorit.IsAtEnd(); ++posteriorit )
        {
        posteriorit.Set(0);
        }
      }

    cout << "making corrections..." << endl;
    int      l, C;
    double   H, cH;
    double * X = new double[NFeature];
    for( exclusionMaskit.GoToBegin(),
         segmaskit.GoToBegin(), csegit.GoToBegin(), nsegit.GoToBegin(), posteriorit.GoToBegin();
         !segmaskit.IsAtEnd(); ++segmaskit, ++csegit, ++nsegit, ++exclusionMaskit, ++posteriorit )
      {
      ImageType::IndexType idx = segmaskit.GetIndex();
      // excluding voxels which are too close to the edge of the image from consideration
      if( idx[0] < DX || idx[0] >= r - DX || idx[1] < DY || idx[1] >= c - DY || idx[2] < DZ || idx[2] >= d - DZ )
        {
        continue;
        }
      if( exclusionMaskit.Value() == 0 && mask->GetPixel(idx) > 0 && segmask->GetPixel(idx) > 0 && cseg->GetPixel(idx) <
          1 )
        {
        autosegnit.SetLocation(idx);
        for( size_t it = 0; it < ims.size(); it++ )
          {
          imnits[it].SetLocation(idx);
          }

        X[0] = idx[0] - mx;
        X[1] = idx[1] - my;
        X[2] = idx[2] - mz;
        C = 3;
        for( j = -DX; j < DX + 1; j++ )
          {
          for( k = -DY; k < DY + 1; k++ )
            {
            for( l = -DZ; l < DZ + 1; l++ )
              {
              NeighborhoodIteratorType::OffsetType offset;
              offset[0] = j;
              offset[1] = k;
              offset[2] = l;
              X[C++] = autosegnit.GetPixel(offset);
              for( size_t it = 0; it < ims.size(); it++ )
                {
                X[C++] = imnits[it].GetPixel(offset);
                }
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
        posteriorit.Set(cH);
        }
      }
    delete X;
    if( posteriorFn.size() )
      {
      cout << "saving posterior map for label " << Tlabel << " .... " << endl;
      sprintf(buffer, posteriorFn.c_str(), Tlabel);
      writer->SetFileName( buffer );
      writer->SetInput(m_PosteriorMap[Tlabel]);
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
      }
    second2 = time(NULL);
    cout << "Label " << Tlabel << " is processed in " << (second2 - second1) << " seconds." << endl;
    second1 = second2;
    }

  if( MRFFlag )
    {
    cout << "start ICM MRF optimization..." << endl;
    MRF_ICM(nseg, m_PosteriorMap, ML, labelset, ICM_beta, ICM_iter);
    }

  cout << "saving corrected segmentation to " << argv[3] << endl;
  // write the corrected segmentation into the output file
  writer->SetFileName( argv[3] );
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
  cout << "finished.." << endl;
  return 0;
}
