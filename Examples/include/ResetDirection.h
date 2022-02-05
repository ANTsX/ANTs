
#ifndef RESETDIRECTION_H
#define RESETDIRECTION_H

namespace ants
{
extern int
ResetDirection(std::vector<std::string>, // equivalent to argv of command line parameters to main()
               std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // RESETDIRECTION_H
