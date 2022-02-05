
#ifndef AVERAGETENSORIMAGES_H
#define AVERAGETENSORIMAGES_H

namespace ants
{
extern int
AverageTensorImages(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                    std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // AVERAGETENSORIMAGES_H
