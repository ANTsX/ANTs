/*
Learns segmentation bias produced by one host segmentation method using AdaBoost and context features.
Hongzhi Wang, 07/18/2010
*/
#include "util.h"
#include "AdaBoost.h"
using namespace std;

int main( int argc, char * * argv )
{
  if( argc != 10 )
    {
    std::cerr << "Missing parameters. " << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              <<
      " inputImageFile manualSegmentationFile autoSegmentationFile Tlabel dilationR featureRadius sampleRatio"
              << " iteration AdaBoostOutputPrefix" << endl << endl
              << " Learns to make corrections for a host segmentation method." << endl << endl
              << " Meanings of the parameters: " << endl
              << " inputImageFile:         A text file saves the locations of training images." << endl
              << "                         Each row of the file stores the location of one training image." << endl
              << "                         For example, a training set with two training images may have an" << endl
              << "                         inputImageFile with the following content:" << endl
              << "                         /home/subject01.nii.gz" << endl
              << "                         /home/subject02.nii.gz" << endl << endl
              << " manualSegmentationFile: A text file saves the locations of the corresponding manual segmentationds"
              << endl
              << "                         for the training images." << endl
              << "                         Each row of the file stores the location of one training image" << endl
              << "                         For example, the manualSegmentationFile for  the above example may have"
              << endl
              << "                         the following content:" << endl
              << "                         /home/subject01_manualSeg.nii.gz" << endl
              << "                         /home/subject02_manulSeg.nii.gz" << endl << endl
              <<
      " autoSegmentationFile:   A text file saves the locations of the corresponding automatic segmentationds" << endl
              << "                         for the training images produced by the host segmentation method." << endl
              << "                         Each row of the file stores  the location of one training image." << endl
              << "                         For example, the autoSegmentationFile for the above example may have the "
              << endl
              << "                         following content:" << endl
              << "                         /home/subject01_autoSeg.nii.gz" << endl
              << "                         /home/subject02_autoSeg.nii.gz" << endl << endl
              << " Tlabel:                 Our method learns to correct one label at a time. This parameter specifies"
              << endl
              << "                         the target label that is going to be learned by this learning task."
              << endl << endl
              << " dilationR:              This parameter specifies the region of interest for this learning task."
              << endl
              << "                         The ROI is defined as the region assigned with the Tlabel by the host"
              << endl
              <<
      "                         segmentation plus the some dilation. dilationR specifies the dilation radius." << endl
              << "                         The dilated ROI should cover most or all voxels that are manually labeled"
              << endl
              << "                         with Tlabel." << endl
              << " featureRadius:          This parameter specifies the neighborhood region for constructing features"
              << endl
              << "                         for learning. An example of this parameter can be 2x3x4, which defines the"
              << endl
              << "                         radius of the region to be 2, 3, 4 for the y-axis,x-axis, and z-axis" << endl
              << "                         respectively. So the neighbohood region has dimension 5x7x9." << endl << endl
              <<
      " sampleRate:             0<= sampleRate <=1. When there are too many training data, loading every voxel in ROI as a"
              << endl
              << "                         training sample may be impossible due to the memory limit. To address this "
              << endl
              <<
      "                         problem, sampleRate specifies the percentage voxels from ROI that will be used" << endl
              << "                         as training samples. If sampleRate=0.01, 1 percent voxels will be used."
              << endl << endl
              << " iteration:              We use AdaBoost learning. This parameter specifies the number of iterations"
              << endl
              << "                         for this learning task. Typically, hundreds of iterations are needed."
              << endl << endl
              << " AdaBoostOutputPrefix:   Our program stores the learned results in two text files. This parameter"
              << endl
              << "                         specifies the path and the prefix of these two file names. For instance, if"
              << endl
              << "                         this parameter is given as /home/mytest, then the output files are" << endl
              << "                         /home/mytest-AdaBoostResults-param-Tlabel?" << endl
              << "                         /home/mytest-AdaBoostResults-Tlabel?" << endl
              << "                         where ? represents the target label. For one host segmentation method, the"
              << endl
              <<
      "                         learning task for each segmentation label should have the same AdaBoostOutputPrefix."
              << endl
              << std::endl;
    return -1;
    }
  WriterType::Pointer writer = WriterType::New();
  ReaderType::Pointer autoseg = ReaderType::New();
  ReaderType::Pointer im = ReaderType::New();

  ifstream imfile;
  string   imstr;
  string   manualstr;
  string   autostr;
  imfile.open(argv[1]);
  int    subjectN = 0;
  string t;
  while( !imfile.eof() )
    {
    getline(imfile, t);
    if( t.length() > 0 )
      {
      subjectN++;
      }
    }

  imfile.close();
  imfile.open(argv[1]);
  cout << "subject #: " << subjectN << endl;
  ifstream manualfile;
  manualfile.open(argv[2]);
  ifstream autofile;
  autofile.open(argv[3]);
  int    Tlabel = atoi(argv[4]);
  int    dilateR = atoi(argv[5]);
  string rstr(argv[6]);
  double sampleRatio = atof(argv[7]);
  cout << "sampleRatio: " << sampleRatio << endl;
  if( sampleRatio > 1 )
    {
    std::cerr << "sampleRatio can not be greater than 1!" << endl;
    return -1;
    }
  int sampleRate = int(1.0 / sampleRatio);
  int iterN = atoi(argv[8]);
  cout << rstr << " : " << rstr.length() << endl;

  int DX, DY, DZ;
  const int p = rstr.find("x");
  const int p1 = rstr.find("x", p + 1);
  t.assign(rstr, 0, p);
  DX = atoi(t.c_str() );
  t.assign(rstr, p + 1, p1 - p);
  DY = atoi(t.c_str() );
  t.assign(rstr, p1 + 1, rstr.length() - 1 - p1);
  DZ = atoi(t.c_str() );
  cout << p << " " << p1 << " " << DX << " " << DY << " " << DZ << endl;
  int NFeature = (2 * (DX * 2 + 1) * (DY * 2 + 1) * (DZ * 2 + 1) + 3) * 4 - 3;
  int MD = DX;
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
  int fileInd = 0;

  vector<double> X;
  vector<int>    Y;
  int            j, k, l, r, c, d,  iNFeature = 0, iNFeature1 = 0;
  int            totalSample = 0;

  while( !imfile.eof() )
    {
    getline(imfile, imstr);
    cout << imstr << endl;
    getline(manualfile, manualstr);
    cout << manualstr << endl;
    getline(autofile, autostr);
    cout << autostr << endl;

    if( imstr.length() && manualstr.length() && autostr.length() )
      {
      cout << ++fileInd << endl;
      ReaderType::Pointer im1 = ReaderType::New();
      ReaderType::Pointer manualSeg = ReaderType::New();
      ReaderType::Pointer autoSeg = ReaderType::New();
      cout << "loading training subject: " << fileInd << " ... " << endl;
      im1->SetFileName(imstr);
      try
        {
        im1->Update();
        cout << im1->GetOutput()->GetRequestedRegion().GetSize()[0] << ":"
             << im1->GetOutput()->GetRequestedRegion().GetSize()[1] << ":"
             << im1->GetOutput()->GetRequestedRegion().GetSize()[2] << endl;
        }
      catch( itk::ExceptionObject & err )
        {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        return -1;
        }

      manualSeg->SetFileName(manualstr);
      try
        {
        manualSeg->Update();
        r = manualSeg->GetOutput()->GetRequestedRegion().GetSize()[0];
        c = manualSeg->GetOutput()->GetRequestedRegion().GetSize()[1];
        d = manualSeg->GetOutput()->GetRequestedRegion().GetSize()[2];

        cout << manualSeg->GetOutput()->GetRequestedRegion().GetSize()[0] << ":"
             << manualSeg->GetOutput()->GetRequestedRegion().GetSize()[1] << ":"
             << manualSeg->GetOutput()->GetRequestedRegion().GetSize()[2] << endl;
        }
      catch( itk::ExceptionObject & err )
        {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        return -1;
        }

      autoSeg->SetFileName(autostr);
      try
        {
        autoSeg->Update();
        cout << autoSeg->GetOutput()->GetRequestedRegion().GetSize()[0] << ":"
             << autoSeg->GetOutput()->GetRequestedRegion().GetSize()[1] << ":"
             << autoSeg->GetOutput()->GetRequestedRegion().GetSize()[2] << endl;
        }
      catch( itk::ExceptionObject & err )
        {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        return -1;
        }

      ImageType::Pointer mask = ImageType::New();
      ImageType::Pointer segmask = ImageType::New();

      segmask->SetRegions(autoSeg->GetOutput()->GetRequestedRegion() );
      segmask->SetSpacing( autoSeg->GetOutput()->GetSpacing() );
      segmask->SetOrigin( autoSeg->GetOutput()->GetOrigin() );
      segmask->SetDirection(autoSeg->GetOutput()->GetDirection() );
      segmask->Allocate();

      mask->SetRegions(autoSeg->GetOutput()->GetRequestedRegion() );
      mask->SetSpacing( autoSeg->GetOutput()->GetSpacing() );
      mask->SetOrigin( autoSeg->GetOutput()->GetOrigin() );
      mask->SetDirection(autoSeg->GetOutput()->GetDirection() );
      mask->Allocate();

      cout << "preparing the data..." << endl;
      // compute ROI,  mask, using the dilation radius
      // The dilated ROI for the whole non-background segmentation, i.e. regions with non-zero labels.
      myDilate(autoSeg->GetOutput(), mask, -1, dilateR);
      // The dilated ROI for the specific label. Note that this label can be the background label, which is 0.
      myDilate(autoSeg->GetOutput(), segmask, Tlabel, dilateR);
      IteratorType      imit0(im1->GetOutput(), im1->GetOutput()->GetRequestedRegion() );
      IteratorType      maskit0(mask, mask->GetRequestedRegion() );
      IteratorType      segmaskit0(segmask, segmask->GetRequestedRegion() );
      IndexIteratorType imaskit(mask, mask->GetRequestedRegion() );
      IndexIteratorType isegmaskit(segmask, segmask->GetRequestedRegion() );

      double mim = 0;
      int    tc = 0;
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

      // compute the spatial center of the specific label
      IndexIteratorType segmaskit(segmask, segmask->GetRequestedRegion() );
      double            mx = 0, my = 0, mz = 0;
      tc = 0;
      int counter = 0;
      for( segmaskit.GoToBegin(); !segmaskit.IsAtEnd(); ++segmaskit )
        {
        ImageType::IndexType idx = segmaskit.GetIndex();
        // excluding voxels which are too close to the edge of the image from consideration
        if( idx[0] < DY || idx[0] >= r - DY || idx[1] < DX || idx[1] >= c - DX || idx[2] < DZ || idx[2] >= d - DZ )
          {
          continue;
          }
        if( segmask->GetPixel(idx) > 0 && (mask->GetPixel(idx) > 0) )
          {
          if( counter++ % sampleRate != 0 )
            {
            continue;
            }

          mx += idx[0];
          my += idx[1];
          mz += idx[2];
          tc++;
          // class label  1 voxels are actually labeled with the same label by human
          //            -1 voxels are not labeled with the same label by human
          if( manualSeg->GetOutput()->GetPixel(idx) == Tlabel )
            {
            Y.push_back(1);
            }
          else
            {
            Y.push_back(-1);
            }
          }
        }
      if( tc < 1 )
        {
        continue;
        }
      mx = mx / tc;
      my = my / tc;
      mz = mz / tc;
      cout << "center location: (" << mx << "," << my << "," << mz << ") feature #:" << NFeature << endl;

      NeighborhoodIteratorType imnit(radius, im1->GetOutput(), im1->GetOutput()->GetRequestedRegion() );
      NeighborhoodIteratorType autosegnit(radius, autoSeg->GetOutput(),
                                          autoSeg->GetOutput()->GetRequestedRegion() );
      NeighborhoodIteratorType::OffsetType offset;
      counter = 0;
      cout << "computing features..." << endl;
      for( segmaskit.GoToBegin(); !segmaskit.IsAtEnd(); ++segmaskit )
        {
        ImageType::IndexType idx = segmaskit.GetIndex();
        // excluding voxels which are too close to the edge of the image from consideration
        if( idx[0] < DY || idx[0] >= r - DY || idx[1] < DX || idx[1] >= c - DX || idx[2] < DZ || idx[2] >= d - DZ )
          {
          continue;
          }

        if( mask->GetPixel(idx) > 0 && segmask->GetPixel(idx) > 0 )
          {
          if( counter++ % sampleRate != 0 )
            {
            continue;
            }
          imnit.SetLocation(idx);
          autosegnit.SetLocation(idx);

          iNFeature = X.size();
          // spatial features: xyz coordinates normalized by the center
          X.push_back(idx[1] - my);
          X.push_back(idx[0] - mx);
          X.push_back(idx[2] - mz);
          // apearance and segmentation features
          for( j = -DX; j < DX + 1; j++ )
            {
            for( k = -DY; k < DY + 1; k++ )
              {
              for( l = -DZ; l < DZ + 1; l++ )
                {
                offset[0] = k;
                offset[1] = j;
                offset[2] = l;
                X.push_back(imnit.GetPixel(offset) );
                X.push_back(autosegnit.GetPixel(offset) );
                }
              }
            }
          // enhance spatial correlations by multiple spatial features with appearnace and segmentation features
          iNFeature1 = X.size();
          for( j = iNFeature; j < iNFeature + 3; j++ )
            {
            for( k = j; k < iNFeature1; k++ )
              {
              X.push_back(X[j] * X[k]);
              }
            }
          }
        }
      totalSample += tc;
      cout << "Total samples:" << tc << "  feature #:" << NFeature << "  #of training data:" << totalSample << endl;
      }
    }

  imfile.close();
  manualfile.close();
  autofile.close();

  if( totalSample < 50 )
    {
    cout << "Too few training samples...";
    return -1;
    }

  char fileName[1024];
  sprintf(fileName, "%s-AdaBoostResults-param-Tlabel%d", argv[9], Tlabel);
  ofstream AdaBoostParamFile(fileName, ios::out);
  AdaBoostParamFile << " " << dilateR << " " << DX << " " << DY << " " << DZ << endl;
  AdaBoostParamFile.close();
  sprintf(fileName, "%s-AdaBoostResults-Tlabel%d", argv[9], Tlabel);
  AdaBoostTrain(&X[0], &Y[0], totalSample, NFeature, iterN, fileName);
  return 0;
}
