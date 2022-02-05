
#ifndef SETSPACING_H
#define SETSPACING_H

namespace ants
{
extern int
SetSpacing(std::vector<std::string>, // equivalent to argv of command line parameters to main()
           std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // SETSPACING_H
