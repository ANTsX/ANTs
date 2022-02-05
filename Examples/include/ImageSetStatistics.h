
#ifndef IMAGESETSTATISTICS_H
#define IMAGESETSTATISTICS_H

namespace ants
{
extern int
ImageSetStatistics(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                   std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // IMAGESETSTATISTICS_H
