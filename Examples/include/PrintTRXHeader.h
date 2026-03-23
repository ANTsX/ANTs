
#ifndef PRINTTRXHEADER_H
#define PRINTTRXHEADER_H

namespace ants
{
extern int
PrintTRXHeader(std::vector<std::string>, // equivalent to argv of command line parameters to main()
               std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // PRINTTRXHEADER_H
