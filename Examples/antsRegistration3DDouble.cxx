#include "antsRegistrationTemplateHeader.h"

namespace ants {

//Instantiate the 3DDouble version
int antsRegistration3DDouble(ParserType::Pointer & parser)
{
    return  DoRegistration<3>( parser );
}

} //end namespace ants
