
#ifndef PRINTHEADER_H
#define PRINTHEADER_H

namespace ants
{
extern int
PrintHeader(std::vector<std::string>, // equivalent to argv of command line parameters to main()
            std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // PRINTHEADER_H
