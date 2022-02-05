#ifndef SIMPLESYNREGISTRATION_H
#define SIMPLESYNREGISTRATION_H

namespace ants
{
extern int
simpleSynRegistration(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                      std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // SIMPLESYNREGISTRATION_H
