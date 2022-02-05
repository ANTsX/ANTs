#include "antsRegistrationTemplateHeader.h"

namespace ants
{

// Instantiate the 4DDouble version
int
antsRegistration4DDouble(ParserType::Pointer & parser)
{
  return DoRegistration<double, 4>(parser);
}

} // end namespace ants
