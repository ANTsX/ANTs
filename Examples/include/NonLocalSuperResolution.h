
#ifndef NONLOCALSUPERRESOLUTION_H
#define NONLOCALSUPERRESOLUTION_H

namespace ants
{
extern int
NonLocalSuperResolution(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                        std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // NONLOCALSUPERRESOLUTION_H
