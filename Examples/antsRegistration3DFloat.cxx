#include "antsRegistrationTemplateHeader.h"

namespace ants
{

// Instantiate the 3DFloat version
int
antsRegistration3DFloat(ParserType::Pointer & parser)
{
  return DoRegistration<float, 3>(parser);
}

} // end namespace ants
