
#ifndef ANTSINTEGRATEVELOCITYFIELD_H
#define ANTSINTEGRATEVELOCITYFIELD_H

namespace ants
{
extern int
ANTSIntegrateVelocityField(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                           std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // ANTSINTEGRATEVELOCITYFIELD_H
