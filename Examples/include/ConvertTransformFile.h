
#ifndef CONVERTTRANSFORMFILE_H
#define CONVERTTRANSFORMFILE_H

namespace ants
{
extern int
ConvertTransformFile(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                     std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // PRINTHEADER_H
