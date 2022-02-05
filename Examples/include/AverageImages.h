
#ifndef AVERAGEIMAGES_H
#define AVERAGEIMAGES_H

namespace ants
{
extern int
AverageImages(std::vector<std::string>, // equivalent to argv of command line parameters to main()
              std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // AVERAGEIMAGES_H
