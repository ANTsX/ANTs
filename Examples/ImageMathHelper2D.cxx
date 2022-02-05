#include "ImageMath_Templates.hxx"

namespace ants
{

int
ImageMathHelper2D(int argc, char ** argv)
{
  int returnval = ImageMathHelperAll<2>(argc, argv);
  if (returnval == EXIT_FAILURE)
  {
    returnval = ImageMathHelper2DOnly<2>(argc, argv);
  }
  if (returnval == EXIT_FAILURE)
  {
    returnval = ImageMathHelper2DOr3D<2>(argc, argv);
  }
  return returnval;
}

} // namespace ants
