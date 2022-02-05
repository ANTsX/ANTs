
#ifndef EXTRACTREGIONFROMIMAGEBYMASK_H
#define EXTRACTREGIONFROMIMAGEBYMASK_H

namespace ants
{
extern int
ExtractRegionFromImageByMask(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                             std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // EXTRACTREGIONFROMIMAGEBYMASK_H
