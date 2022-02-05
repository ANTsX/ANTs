
#ifndef SETORIGIN_H
#define SETORIGIN_H

namespace ants
{
extern int
SetOrigin(std::vector<std::string>, // equivalent to argv of command line parameters to main()
          std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // SETORIGIN_H
