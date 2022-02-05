
#ifndef DENOISEIMAGE_H
#define DENOISEIMAGE_H

namespace ants
{
extern int
DenoiseImage(std::vector<std::string>, // equivalent to argv of command line parameters to main()
             std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // DENOISEIMAGE_H
