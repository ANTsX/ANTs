
#ifndef THRESHOLDIMAGE_H
#define THRESHOLDIMAGE_H

namespace ants
{
extern int
ThresholdImage(std::vector<std::string>, // equivalent to argv of command line parameters to main()
               std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // THRESHOLDIMAGE_H
