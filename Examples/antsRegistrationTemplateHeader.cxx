#include "antsRegistrationTemplateHeader.h"

namespace ants
{
const char *
RegTypeToFileName(const std::string & type, bool & writeInverse, bool & writeVelocityField, bool minc)
{
  std::string str(type);

  ConvertToLowerCase(str);
  if (str == "syn" || str == "symmetricnormalization" || str == "bsplinesyn" ||
      str == "timevaryingbsplinevelocityfield" || str == "tvdmffd" || str == "timevaryingvelocityfield" ||
      str == "tvf" || str == "exponential" || str == "bsplineexponential")
  {
    writeInverse = true;
  }
  else
  {
    writeInverse = false;
  }

  if (str == "timevaryingvelocityfield" || str == "tvf" || str == "exp" || str == "exponential" ||
      str == "bsplineexponential")
  {
    writeVelocityField = true;
  }
  else
  {
    writeVelocityField = false;
  }

  if (str == "rigid")
  {
    if (minc)
      return "_Rigid.xfm";
    else
      return "Rigid.mat";
  }
  else if (str == "affine" || str == "compositeaffine" || str == "compaff")
  {
    if (minc)
      return "_Affine.xfm";
    else
      return "Affine.mat";
  }
  else if (str == "similarity")
  {
    if (minc)
      return "_Similarity.xfm";
    else
      return "Similarity.mat";
  }
  else if (str == "translation")
  {
    if (minc)
      return "_Translation.xfm";
    else
      return "Translation.mat";
  }
  else if (str == "bspline" || str == "ffd")
  {
    if (minc)
      return "_BSpline.txt";
    else
      return "BSpline.txt";
  }
  else if (str == "genericaffine")
  {
    if (minc)
      return "_GenericAffine.xfm";
    else
      return "GenericAffine.mat";
  }
  else if (str == "gaussiandisplacementfield" || str == "gdf" || str == "bsplinedisplacementfield" || str == "dmffd" ||
           str == "syn" || str == "symmetricnormalization" || str == "bsplinesyn" || str == "exp" ||
           str == "exponential" || str == "bsplineexponential")
  {
    if (minc)
      return "_NL.xfm";
    else
      return "Warp.nii.gz";
  }
  else if (str == "timevaryingvelocityfield" || str == "tvf" || str == "timevaryingbsplinevelocityfield" ||
           str == "tvdmffd")
  {
    if (minc)
      return "_Warp.mnc";
    else
      return "Warp.nii.gz";
  }
  return "BOGUS.XXXX";
}
} // end namespace ants
