
#ifndef MEASUREMINMAXMEAN_H
#define MEASUREMINMAXMEAN_H

namespace ants
{
extern int
MeasureMinMaxMean(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                  std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // MEASUREMINMAXMEAN_H
