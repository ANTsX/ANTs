/*
Learns segmentation bias produced by one host segmentation method using AdaBoost and context features.
Hongzhi Wang, 07/18/2010
*/
#include "util.h"
#include "AdaBoost.h"
using namespace std;

int usage()
{
  cout << "Corrective Learning: Learning classifiers for correcting segmentation errors for a host segmentation method. " << endl;
  cout << "usage: " << endl;
  cout << " bl dim [options] AdaBoostOutputPrefix" << endl;
  cout << endl;
  cout << "required options:" << endl;
  cout << "  dim                                 Image dimension (2 or 3)" << endl;
  cout << "  -ms sub1_mseg.nii ... subN_mseg.nii Manual segmentations" << endl;
  cout << "  -as sub1_aseg.nii ... subN_aseg.nii Automatic segmentationds produced for the training images by the host segmentation method" << endl;
  cout << "  -tl Tlabel                          Our method learns to correct one label at a time. This parameter specifies"<<endl;
  cout << "                                      the target label that is going to be learned by this learning task."<<endl<<endl;
  cout << "  -rd dilationR                       This parameter specifies the region of interest (ROI) for this learning task."<<endl;
  cout << "                                      A label's working ROI is defined as the region assigned with the Tlabel by the host"<<endl;
  cout << "                                      segmentation plus some dilation (unless a mask is provided, see below)." <<endl;
  cout << "                                      dilationR specifies the dilation radius."<<endl;
  cout << "                                      The dilated ROI should cover most or all voxels that are manually assigned to Tlabel"<<endl;
  cout << "                                      Default: 2x2x2" << endl;
  cout << "  -rf featureRadius                   Patch radius for feature extraction. " << endl;
  cout << "                                      scalar or vector (AxBxC) " << endl;
  cout << "                                      Default: 2x2x2" << endl;
  cout << "  -rate sampleRate                    0<= sampleRate <=1. When the ROI is large or there are too many training images," <<endl;
  cout << "                                      loadng every voxel in ROI as a training sample may be impossible due to the memory limit."<<endl;
  cout << "                                      To address this problem, sampleRate specifies the percentage voxels from ROI that will be"<<endl;
  cout << "                                      problem, sampleRate specifies the percentage voxels from ROI that will be used"<<endl;
  cout << "                                      used. If sampleRate=0.01, 1 percent voxels will be used."<<endl;
  cout << "                                      Default: 1" << endl;
  cout << "  -i iteration                        Training iteration for AdaBoost learning." << endl;
  cout << "                                      Default: 100" << endl;
  cout << " options: " << endl
       << "  -c featureChannel                   Number of feature channels" << endl
       << "  -f sub1_feature1.nii sub1_feature2.nii ... subN_feature1.nii subN_feature2.nii ... "<<endl
       << "                                      Feature images for the training subjects. The number of feautres for each"<<endl
       << "                                      subject should  equal featureChannel." << endl
       << "  -m mask1 ... maskN:                 Specify labels' working ROI for the training images." <<endl
       << "                                      ROI will be derived by performing dilation on this mask then."<<endl;

  return -1;
}

template<unsigned int VDim>
struct BLParam
{
  vector<string> fnManual;
  vector<string> fnAuto;
  vector<string> fnFeature;
  vector<string> fnMask;
  vector<string> fnPosterior;

  string fnOutput;
  
  int featureChannel;
  int targetLabel;
  int DilateR;
  double sampleRate;
  int iteration,DX,DY,DZ;

  BLParam()
    {
    featureChannel = 0;
    DilateR = 1;
    sampleRate = 1;
    iteration = 100;
    DX=2;
    DY=2;
    DZ=2;
    }
  void Print(std::ostream &oss)
  {
    oss << " dilation radius: " << DilateR << endl;
    oss << " sampling rate: " << sampleRate << endl;
    oss << " iteration: " << iteration << endl;
    oss << " feature channels: " << featureChannel << endl;
    for (size_t i=0;i<fnManual.size();i++)
    {
      oss << "manual segmentation: " << fnManual[i] << endl;
      oss << "automatic segmentation: " << fnAuto[i] << endl;
      for (int j=0;j<featureChannel;j++)
      {
        oss << "feature "<<j<<" : " << fnFeature[i*featureChannel+j] << endl;
      }
    }
    oss << "Output prefix: " << fnOutput << endl;
    oss << "Feature Radius: (" << DX <<","<<DY<<","<<DZ<<")" << endl;
  }
};

template <unsigned int VDim>
int lfapp(int argc, char *argv[])
{
  // Parameter vector
  BLParam<VDim> p;

  // Read the parameters from command line
  p.fnOutput = argv[argc-1];
  int argend = argc-2;
  for(int j = 2; j < argc-1; j++)
    {
    string arg = argv[j];

    if(arg == "-ms")
    {
      // Read the following options as manual segmentations
      while(argv[j+1][0] != '-' && j < argend)
      {
        p.fnManual.push_back(argv[++j]);
        cout<<p.fnManual[p.fnManual.size()-1]<<endl;
      }
    }
    else if(arg == "-as")
    {
      // Read the following options as automatic segmentations
      while(argv[j+1][0] != '-' && j < argend)
      {
        p.fnAuto.push_back(argv[++j]);
        cout<<p.fnAuto[p.fnAuto.size()-1]<<endl;
      }
    }

    else if(arg == "-f")
    {
      // Read the following options as feature images
      while(argv[j+1][0] != '-' && j < argend)
      {
        p.fnFeature.push_back(argv[++j]);
        cout<<p.fnFeature[p.fnFeature.size()-1]<<endl;
      }
    }
    else if(arg == "-m")
      {
      // Read the following options as ROI images
      while(argv[j+1][0] != '-' && j < argend)
        p.fnMask.push_back(argv[++j]);
      }
    else if(arg == "-tl")
      {
        p.targetLabel=atoi(argv[++j]);
      }
    else if(arg == "-c")
      {
        p.featureChannel=atoi(argv[++j]);
        cout<<"feature channel: "<<p.featureChannel<<endl;
      }
    else if(arg == "-rd")
      {
      p.DilateR = atoi(argv[++j]);
      }
    else if(arg == "-rate")
      {
      p.sampleRate = atof(argv[++j]);
      }
    else if(arg == "-i")
      {
      p.iteration = atoi(argv[++j]);
      }
    else if(arg == "-rf" && j < argend)
      {
      string rstr(argv[++j]);
      string t;
      int tp =-1,p1;
      tp=rstr.find("x");
      p1=rstr.find("x",tp+1);
      t.assign(rstr,0,tp);
      p.DX=atoi(t.c_str());
      t.assign(rstr,tp+1,p1-tp);
      p.DY=atoi(t.c_str());
      t.assign(rstr,p1+1,rstr.length()-1-p1);
      p.DZ=atoi(t.c_str());
      cout<<tp<<" "<<p1<<" "<<p.DX<<" "<<p.DY<<" "<<p.DZ<<endl;
      }
    else
      {
      cerr << "Unknown option " << arg << endl;
      return -1;
      }
    }
  // We have the parameters now. Check for validity
  if(p.fnFeature.size() != p.fnManual.size()*p.featureChannel)
  {
    cerr << "Number of features and segmentations does not match! "<<p.fnFeature.size()/p.featureChannel
         <<" atlas images "<<p.fnManual.size()<<"segmentatons" << endl;
    return -1;
  }
  if(p.fnAuto.size() != p.fnManual.size())
    {
    cerr << "Number of automatic segmentations and manaul segmentations does not match!"<<p.fnAuto.size()
         <<" automatic segmentations "<<p.fnManual.size()<<" manual segmentatons" << endl;
    return -1;
    }

  // Print parametes
  cout << "segAdapter PARAMETERS:" << endl;
  p.Print(cout);

//  std::vector<ImageType::Pointer> ims;
//  std::vector<IteratorType> imits;
//  std::vector<NeighborhoodIteratorType> imnits;

  ifstream imfile;
  string imstr;
  string manualstr;
  string autostr;
  string maskstr;

  int imFlag=0,maskFlag=0;

  int NFeature;
  NFeature = ((1+p.featureChannel)*(p.DX*2+1)*(p.DY*2+1)*(p.DZ*2+1)+3)*4-3;

  int MD=p.DX;
  if (p.DY>MD) MD=p.DY;
  if (p.DZ>MD) MD=p.DZ;

  NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(MD);
  int fileInd=0;

  vector<double> X;
  vector<int> Y;
  int j, k, l, r, c, d,  iNFeature=0, iNFeature1=0;
  int totalSample=0;
  for (size_t mi=0;mi<p.fnManual.size();mi++)
  {
    std::vector<ImageType::Pointer> ims;
    std::vector<IteratorType> imits;
    std::vector<NeighborhoodIteratorType> imnits;
    manualstr=p.fnManual[mi];
    cout<<manualstr<<endl;
    autostr=p.fnAuto[mi];
    cout<<autostr<<endl;
    if (p.fnMask.size()>0)
    {
      maskstr=p.fnMask[mi];
      cout<<maskstr<<endl;
    }

    if (manualstr.length() && autostr.length())
    {
      cout<<++fileInd<<endl;
      ReaderType::Pointer manualSeg = ReaderType::New();
      ReaderType::Pointer autoSeg = ReaderType::New();
      ReaderType::Pointer autoPosterior = ReaderType::New();
      ReaderType::Pointer ROIMask = ReaderType::New();

      cout<<"loading training subject: "<<fileInd<<" ... "<<endl;
      if (p.fnFeature.size())
      {
        for (int  it=0;it<p.featureChannel;it++)
        {
          imstr=p.fnFeature[mi*p.featureChannel+it].c_str();
          cout<<imstr<<endl;
          cout<<"loading "<<imstr<<endl;
          ReaderType::Pointer im1 = ReaderType::New();

          im1->SetFileName(imstr);
          try
          {
            im1->Update();
            cout<<im1->GetOutput()->GetRequestedRegion().GetSize()[0]<<":"
                <<im1->GetOutput()->GetRequestedRegion().GetSize()[1]<<":"
                <<im1->GetOutput()->GetRequestedRegion().GetSize()[2]<<endl;
          }
          catch ( itk::ExceptionObject &err)
          {
            cout << "ExceptionObject caught !" << std::endl;
            cout << err << std::endl;
            return -1;
          }
          ims.push_back(im1->GetOutput());
        }
      }
      if (p.fnMask.size())
      {
        cout<<"loading "<<maskstr<<endl;
        ROIMask->SetFileName(maskstr);
        try
        {
          ROIMask->Update();
        }
        catch ( itk::ExceptionObject &err)
        {
          cout << "ExceptionObject caught !" << std::endl;
          cout << err << std::endl;
          return -1;
        }
      }

      cout<<"loading "<<manualstr<<endl;
      manualSeg->SetFileName(manualstr);
      try
      {
        manualSeg->Update();
        r=manualSeg->GetOutput()->GetRequestedRegion().GetSize()[0];
        c=manualSeg->GetOutput()->GetRequestedRegion().GetSize()[1];
        d=manualSeg->GetOutput()->GetRequestedRegion().GetSize()[2];
        cout<<manualSeg->GetOutput()->GetRequestedRegion().GetSize()[0]<<":"
            <<manualSeg->GetOutput()->GetRequestedRegion().GetSize()[1]<<":"
            <<manualSeg->GetOutput()->GetRequestedRegion().GetSize()[2]<<endl;
      }
      catch ( itk::ExceptionObject &err)
      {
        cout << "ExceptionObject caught !" << std::endl;
        cout << err << std::endl;
        return -1;
      }

      autoSeg->SetFileName(autostr);
      try
      {
        autoSeg->Update();
        cout<<autoSeg->GetOutput()->GetRequestedRegion().GetSize()[0]<<":"
            <<autoSeg->GetOutput()->GetRequestedRegion().GetSize()[1]<<":"
            <<autoSeg->GetOutput()->GetRequestedRegion().GetSize()[2]<<endl;
      }
      catch ( itk::ExceptionObject &err)
      {
        cout << "ExceptionObject caught !" << std::endl;
        cout << err << std::endl;
        return -1;
      }

      ImageType::Pointer mask = ImageType::New();
      ImageType::Pointer segmask = ImageType::New();

      segmask->SetRegions(autoSeg->GetOutput()->GetRequestedRegion());
      segmask->SetSpacing( autoSeg->GetOutput()->GetSpacing() );
      segmask->SetOrigin( autoSeg->GetOutput()->GetOrigin() );
      segmask->SetDirection(autoSeg->GetOutput()->GetDirection());
      segmask->Allocate();

      mask->SetRegions(autoSeg->GetOutput()->GetRequestedRegion());
      mask->SetSpacing( autoSeg->GetOutput()->GetSpacing() );
      mask->SetOrigin( autoSeg->GetOutput()->GetOrigin() );
      mask->SetDirection(autoSeg->GetOutput()->GetDirection());
      mask->Allocate();

      cout<<"preparing the data..."<<endl;
      // compute ROI,  mask, using the dilation radius
      // The dilated ROI for the whole non-background segmentation, i.e. regions with non-zero labels.
      IndexIteratorType segmaskit(segmask, segmask->GetRequestedRegion());
      IndexIteratorType imaskit(mask, mask->GetRequestedRegion());
      IndexIteratorType autosegit(autoSeg->GetOutput(), autoSeg->GetOutput()->GetRequestedRegion());
      if (maskFlag)
      {
        myDilate(ROIMask->GetOutput(),segmask,-1,p.DilateR);
        myDilate(autoSeg->GetOutput(),mask,-1,p.DilateR);
      }
      else
      {
        cout<<"DilateR: "<<p.DilateR<<endl;
        // The dilated ROI for the specific label. Note that this label can be the background label, which is 0.
        myDilate(autoSeg->GetOutput(),segmask,p.targetLabel,p.DilateR);
        myDilate(autoSeg->GetOutput(),mask,-1,p.DilateR);
      }

      IteratorType maskit0(mask, mask->GetRequestedRegion() );
      IteratorType segmaskit0(segmask, segmask->GetRequestedRegion() );

      if (p.fnFeature.size())
      {
        for (size_t it=0;it<ims.size();it++)
        {
          double mim=0;
          int tc=0;
          IteratorType imit0(ims[it], ims[it]->GetRequestedRegion() );
          NeighborhoodIteratorType imnit(radius, ims[it], ims[it]->GetRequestedRegion());
          imits.push_back(imit0);
          imnits.push_back(imnit);
          //normalize image intensity by the mean intensity of the whole figure region
          for (imit0.GoToBegin(), maskit0.GoToBegin(); !imit0.IsAtEnd(); ++imit0, ++maskit0)
          {
            if (maskit0.Value()>0)
            {
              mim += imit0.Value();
              tc++;
            }
          }
          mim /= tc;
          cout<<"mean intensity : sample #  "<<mim<<":"<<tc<<endl;
          for (imit0.GoToBegin(); !imit0.IsAtEnd(); ++imit0)
          {
            imit0.Set(imit0.Value()-mim);
          }
        }
      }

      //compute the spatial center of the specific label
      double mx=0,my=0,mz=0,totalval=0;
      int tc=0;
//      int sampleFlag=0;
      for (segmaskit.GoToBegin(),autosegit.GoToBegin(); !segmaskit.IsAtEnd(); ++segmaskit,++autosegit)
      {
        ImageType::IndexType idx = segmaskit.GetIndex();
        //excluding voxels which are too close to the edge of the image from consideration
        if (idx[0]<p.DX || idx[0]>=r-p.DX || idx[1]<p.DY || idx[1]>=c-p.DY || idx[2]<p.DZ || idx[2]>=d-p.DZ)
          continue;
        if (segmask->GetPixel(idx)>0 && (mask->GetPixel(idx)>0))
        {
          mx+=idx[0];
          my+=idx[1];
          mz+=idx[2];
          totalval+=1;
          tc++;
        }
      }
      if (tc<1)
        continue;
      mx=mx/totalval;
      my=my/totalval;
      mz=mz/totalval;
      cout<<"center location: ("<<mx<<","<<my<<","<<mz<<") feature #:"<<NFeature<<endl;

      NeighborhoodIteratorType imnit;
      NeighborhoodIteratorType autosegnit(radius, autoSeg->GetOutput(), autoSeg->GetOutput()->GetRequestedRegion());
      NeighborhoodIteratorType::OffsetType offset;
      cout<<"computing features..."<<endl;
      tc=0;
      for (segmaskit.GoToBegin();!segmaskit.IsAtEnd(); ++segmaskit)
      {
        ImageType::IndexType idx = segmaskit.GetIndex();
        //excluding voxels which are too close to the edge of the image from consideration
        if (idx[0]<p.DX || idx[0]>=r-p.DX || idx[1]<p.DY || idx[1]>=c-p.DY || idx[2]<p.DZ || idx[2]>=d-p.DZ)
          continue;
        if (mask->GetPixel(idx)>0 && segmask->GetPixel(idx)>0)
        {
          double randtest=((double) rand() / (RAND_MAX));
          if (randtest>p.sampleRate)
            continue;
          tc++;
          //class label  1 : voxels are labeled with the same label by human
          //            -1 : voxels are not labeled with the same label by human
          if (manualSeg->GetOutput()->GetPixel(idx)==p.targetLabel)
            Y.push_back(1);
          else
            Y.push_back(-1);
          autosegnit.SetLocation(idx);
          for (int it=0;it<p.featureChannel;it++)
          {
            imnits[it].SetLocation(idx);
          }

          iNFeature=X.size();
          // spatial features: xyz coordinates normalized by the center
          X.push_back(idx[0]-mx);
          X.push_back(idx[1]-my);
          X.push_back(idx[2]-mz);
          // apearance and segmentation features
          for (j=-p.DX;j<p.DX+1;j++)
          {
            for (k=-p.DY;k<p.DY+1;k++)
            {
              for (l=-p.DZ;l<p.DZ+1;l++)
              {
                offset[0]=j;
                offset[1]=k;
                offset[2]=l;
                X.push_back(autosegnit.GetPixel(offset));
                for (int it=0;it<p.featureChannel;it++)
                {
                  X.push_back(imnits[it].GetPixel(offset));
                }
              }
            }
          }

          // enhance spatial correlations by multiple spatial features with appearnace and segmentation features
          iNFeature1=X.size();
          for (j=iNFeature;j<iNFeature+3;j++)
          {
            for (k=j;k<iNFeature1;k++)
            {
              X.push_back(X[j]*X[k]);
            }
          }
        }
      }
      totalSample+=tc;
      cout<<X.size()<<endl;
      cout<< "Total samples:"<<tc<<"  feature #:"<<NFeature<<"  # of training data:"<<totalSample<<" feature:"<<imFlag<<endl;
      if (totalSample>5000000)
      {
         cout<<"Too many training samples...";
         return -1;
      }
    }
  }
  if (totalSample<50)
  {
    cout<<"Too few training samples...";
    return -1;
  }
  char fileName[1024];
  sprintf (fileName, "%s-AdaBoostResults-param-Tlabel%d",p.fnOutput.c_str(),p.targetLabel);
  ofstream AdaBoostParamFile(fileName, ios::out);
  AdaBoostParamFile<<" "<<p.DilateR<<" "<<p.DX<<" "<<p.DY<<" "<<p.DZ<<" "<<p.featureChannel<<endl;
  AdaBoostParamFile.close();
  sprintf (fileName, "%s-AdaBoostResults-Tlabel%d",p.fnOutput.c_str(),p.targetLabel);

  AdaBoostTrain(&X[0],&Y[0],totalSample,NFeature,p.iteration,fileName);

  return 0;
}



int main( int argc, char ** argv )
{
  // Parse user input
  if(argc < 5) return usage();

  // Get the first option
  int dim = atoi(argv[1]);

  // Call the templated method
  if(dim == 2)
    return lfapp<2>(argc, argv);
  else if(dim == 3)
    return lfapp<3>(argc, argv);
  else
    {
    cerr << "Dimension " << argv[1] << " is not supported" << endl;
    return -1;
    }
}

