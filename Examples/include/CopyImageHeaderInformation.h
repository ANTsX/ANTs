
#ifndef COPYIMAGEHEADERINFORMATION_H
#define COPYIMAGEHEADERINFORMATION_H

namespace ants
{
extern int
CopyImageHeaderInformation(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                           std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // COPYIMAGEHEADERINFORMATION_H
