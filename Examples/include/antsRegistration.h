
#ifndef ANTSREGISTRATION_H
#define ANTSREGISTRATION_H

namespace ants
{
int antsRegistration( std::vector<std::string>, // equivalent to argv of command line parameters to main()
                      std::ostream* out_stream  // [optional] output stream to write
                      );
} // namespace ants

#endif // ANTSREGISTRATION_H
