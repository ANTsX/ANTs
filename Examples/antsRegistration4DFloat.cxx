#include "antsRegistrationTemplateHeader.h"

namespace ants
{

// Instantiate the 4DFloat version
int
antsRegistration4DFloat(ParserType::Pointer & parser)
{
  return DoRegistration<float, 4>(parser);
}

} // end namespace ants
