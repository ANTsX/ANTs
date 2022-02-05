#ifndef IMAGEINTENSITYSTATISTICS_H
#define IMAGEINTENSITYSTATISTICS_H

namespace ants
{
extern int
ImageIntensityStatistics(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                         std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // IMAGEINTENSITYSTATISTICS_H
