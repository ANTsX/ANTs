
#ifndef SCCAN_H
#define SCCAN_H

namespace ants
{
extern int
sccan(std::vector<std::string>, // equivalent to argv of command line parameters to main()
      std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // SCCAN_H
