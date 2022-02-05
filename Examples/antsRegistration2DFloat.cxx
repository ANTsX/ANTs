#include "antsRegistrationTemplateHeader.h"

namespace ants
{

// Instantiate the 2DFloat version
int
antsRegistration2DFloat(ParserType::Pointer & parser)
{
  return DoRegistration<float, 2>(parser);
}

} // end namespace ants
