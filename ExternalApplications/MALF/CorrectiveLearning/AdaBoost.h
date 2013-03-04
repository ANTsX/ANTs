/*
Implements AdaBoost learning and classification.

Implemented by Hongzhi Wang 07/18/2010
*/

#include <stdio.h>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <time.h>
using namespace std;

typedef struct
  {
  double x;
  int y;
  double w;
  } node;

bool operator<(const node & n1, const node & n2)
{
  return n1.x < n2.x;
}

bool compare(const node & n1, const node & n2)
{
  return n1.x < n2.x;
}

typedef struct
  {
  int featureID;
  double weightedRate;
  int sign;
  double threshold;
  } weakLearner;

// Construct an optimal linear weak learner using one feature by a single thresholding
// X:  stores response of each training sample at one feature
// wl: stores the optimal weak learner for this feature
// void weakClassifier(const vector<node> X, weakLearner &wl){
void weakClassifier(node *X, int NSample, weakLearner & wl)
{
  double p1 = 0, p2, mp;
  int    i; // , NSample=X.size();

  for( i = 0; i < NSample; i++ )
    {
    if( X[i].y == 1 )
      {
      p1 += X[i].w;
      }
    }
  p2 = 1 - p1;
  if( p1 >= p2 )
    {
    wl.threshold = X[0].x - 0.001;
    wl.sign = 1;
    mp = p1;
    }
  else
    {
    wl.threshold = X[0].x - 0.001;
    wl.sign = -1;
    mp = p2;
    }

  i = 0;
  while( i < NSample )
    {
    while( i < NSample && X[i].x == X[i + 1].x )
      {
      if( X[i].y == 1 )
        {
        p1 -= X[i].w;
        }
      else
        {
        p1 += X[i].w;
        }
      i++;
      }

    if( X[i].y == 1 )
      {
      p1 -= X[i].w;
      }
    else
      {
      p1 += X[i].w;
      }
    p2 = 1 - p1;
    if( p1 > mp )
      {
      mp = p1;
      wl.threshold = X[i].x;
      wl.sign = 1;
      }
    if( p2 > mp )
      {
      mp = p2;
      wl.threshold = X[i].x;
      wl.sign = -1;
      }
    i++;
    }

  wl.weightedRate = mp;
}

// Implement AdaBoost learning
// X: input training samples. An 1D array of size NSamplexNFeature,
// Y: ground truth labels of the training samples. A array of size NSample with value 1 of -1.
// NSample: #of train samples.
// NFeature: #of features used to describe each sample.
// iterN: #of AdaBoost learning iterations.
// fileName: the file names of the learning output.
int AdaBoostTrain(double* X, int* Y, int NSample, int NFeature, int iterN, char* fileName)
{
  double * W = new double[NSample], tw = 1.0 / NSample, R;
  double * CH = new double[NSample], *F = new double[NSample];
  double * H = new double[NSample];
  node *   tX = new node[NSample];
  double   weightedError, alpha, totalW;
  double   trueError;
  int      i, j, CC = 0;
  int      LP = 0, LN = 0;

  for( i = 0; i < NSample; i++ )
    {
    if( Y[i] == 1 )
      {
      LP++;
      }
    W[i] = tw;
    H[i] = 0;
    }
  LN = NSample - LP;
  if( LP > LN )
    {
    R = double(LN) / NSample;
    }
  else
    {
    R = double(LP) / NSample;
    }

  if( R == 0 )
    {
    cout << "Errors in AdaBoostTrain!" << endl
         << "All training samples are from one class!" << endl
         << "This error was produced either because the initial segmentation produced by the" << endl
         << "host segmentation algorithm does not contain this label or because the initial" << endl
         << "are too far away from the manually labeled regions." << endl
         << "For the first case, our current algorithm can not recover errors for this target." << endl
         << "But we still should be able to recover errors for other labeles." << endl
         << "For the second case, use a larger dilation radius should fix this problem." << endl;
    delete [] tX;
    delete [] H;
    delete [] F;
    delete [] CH;
    delete [] W;
    return -1;
    }

  weakLearner wl, bwl;
  bwl.weightedRate = 0;
  bwl.sign = 0;
  bwl.threshold = 0;
  bwl.featureID = -1;
  ofstream AdaBoostTrainFile(fileName, ios::out);

  cout << NSample << " " << NFeature << " " << iterN << " " << fileName << endl;
  time_t second1 = time(NULL), second2;
  while( CC < iterN )
    {
    CC++;
    bwl.weightedRate = -1;
    for( i = 0; i < NFeature; i++ )
      {
      for( j = 0; j < NSample; j++ )
        {
        tX[j].x = X[j * NFeature + i];
        tX[j].y = Y[j];
        tX[j].w = W[j];
        }
      sort(tX, tX + NSample);

      weakClassifier(tX, NSample, wl);
      if( wl.weightedRate > bwl.weightedRate )
        {
        bwl.sign = wl.sign;
        bwl.weightedRate = wl.weightedRate;
        bwl.threshold = wl.threshold;
        bwl.featureID = i;
        }
      }
    alpha = 0.5 * log(bwl.weightedRate / (1 - bwl.weightedRate) );
    weightedError = 1;
    totalW = 0;
    trueError = 0;
    for( i = 0; i < NSample; i++ )
      {
      if( bwl.sign == 1 )
        {
        CH[i] = 2 * double(X[i * NFeature + bwl.featureID] > bwl.threshold) - 1;
        }
      else
        {
        CH[i] = 2 * double(X[i * NFeature + bwl.featureID] <= bwl.threshold) - 1;
        }
      F[i] = double(CH[i] == Y[i]);
      H[i] += alpha * CH[i];
      if( H[i] * Y[i] < 0 )
        {
        trueError++;
        }
      weightedError -= F[i] * W[i];
      W[i] *= exp(-alpha * (CH[i] * Y[i]) );
      totalW += W[i];
      }
    cout << "weighted error: " << 1 - bwl.weightedRate << ":" << weightedError << endl;
    for( i = 0; i < NSample; i++ )
      {
      W[i] /= totalW;
      }
    second2 = time(NULL);
    cout << "iter" << CC << "  trainingError:" << trueError << "  errorRatio:" << trueError / NSample
         << "  error(randomGuess):" << R << "  featureID:" << bwl.featureID << "  alpha:" << alpha
         << "  sign:" << bwl.sign << "  threshold:" << bwl.threshold
         << "  time(s): " << (second2 - second1) << endl;

    AdaBoostTrainFile << CC << " " << alpha << " " << weightedError << " " << bwl.featureID << " " << bwl.sign << " "
                      << bwl.threshold << " " << trueError / NSample << endl;
    second1 = second2;
    }

  AdaBoostTrainFile.close();
  delete [] tX;
  delete [] W;
  delete [] H;
  delete [] F;
  delete [] CH;
  return 0;
}

// implement AdaBoost classification
// X: input testing sampleis. An 1D array of size NSamplexNFeature,
// NSample: #of the testing samples
// featureID, alpha, sign, threshold are parameters specifying AdaBoost weaklearners.
// H: the discretized response of the AdaBoost classiffier, 1 or -1.
// cH: AdaBoost classifier's responses transferred to probabilities.
void AdaBoostClassify(double *X, int NSample, int NFeature, int LC, int *featureID, double *alpha, int *sign,
                      double *threshold, double* H, double* cH)
{
  int i, j;

  for( i = 0; i < NSample; i++ )
    {
    cH[i] = 0;
    H[i] = 0;
    }
  for( i = 0; i < LC; i++ )
    {
    if( sign[i] == 1 )
      {
      for( j = 0; j < NSample; j++ )
        {
        if( X[j * NFeature + featureID[i]] > threshold[i] )
          {
          cH[j] += alpha[i];
          }
        else
          {
          cH[j] -= alpha[i];
          }
        }
      }
    else
      {
      for( j = 0; j < NSample; j++ )
        {
        if( X[j * NFeature + featureID[i]] <= threshold[i] )
          {
          cH[j] += alpha[i];
          }
        else
          {
          cH[j] -= alpha[i];
          }
        }
      }
    }
  for( i = 0; i < NSample; i++ )
    {
    double t = exp(cH[i]);
    cH[i] = t / (t + 1 / t);
    if( cH[i] > 0.5 )
      {
      H[i] = 1;
      }
    else
      {
      H[i] = -1;
      }
    }
}
