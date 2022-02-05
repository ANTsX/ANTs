
#ifndef CHECKTOPOLOGY_H
#define CHECKTOPOLOGY_H

namespace ants
{
extern int
CheckTopology(std::vector<std::string>, // equivalent to argv of command line parameters to main()
              std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // CHECKTOPOLOGY_H
