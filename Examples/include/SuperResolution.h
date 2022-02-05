
#ifndef SUPERRESOLUTION_H
#define SUPERRESOLUTION_H

namespace ants
{
extern int
SuperResolution(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // SUPERRESOLUTION_H
