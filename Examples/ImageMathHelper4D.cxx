#include "ImageMath_Templates.hxx"

namespace ants
{

int
ImageMathHelper4D(int argc, char ** argv)
{
  int returnval = ImageMathHelperAll<4>(argc, argv);
  if (returnval == EXIT_FAILURE)
  {
    returnval = ImageMathHelper4DOnly<4>(argc, argv);
  }
  if (returnval == EXIT_FAILURE)
  {
    returnval = ImageMathHelper3DOr4D<4>(argc, argv);
  }
  return returnval;
}


} // namespace ants
