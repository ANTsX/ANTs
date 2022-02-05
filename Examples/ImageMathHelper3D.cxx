#include "ImageMath_Templates.hxx"

namespace ants
{

int
ImageMathHelper3D(int argc, char ** argv)
{
  int returnval = ImageMathHelperAll<3>(argc, argv);
  if (returnval == EXIT_FAILURE)
  {
    returnval = ImageMathHelper3DOnly<3>(argc, argv);
  }
  if (returnval == EXIT_FAILURE)
  {
    returnval = ImageMathHelper2DOr3D<3>(argc, argv);
  }
  if (returnval == EXIT_FAILURE)
  {
    returnval = ImageMathHelper3DOr4D<3>(argc, argv);
  }
  return returnval;
}

} // namespace ants
