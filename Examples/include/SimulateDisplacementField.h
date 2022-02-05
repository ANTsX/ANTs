#ifndef SIMULATEDISPLACEMENTFIELD_H
#define SIMULATEDISPLACEMENTFIELD_H

namespace ants
{
extern int
SimulateDisplacementField(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                          std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // SIMULATEIMAGE_H
