
#ifndef TIMESCCAN_H
#define TIMESCCAN_H

namespace ants
{
extern int
TimeSCCAN(std::vector<std::string>, // equivalent to argv of command line parameters to main()
          std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // TIMESCCAN_H
