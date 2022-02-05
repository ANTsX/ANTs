
#ifndef CONVERTIMAGEPIXELTYPE_H
#define CONVERTIMAGEPIXELTYPE_H

namespace ants
{
extern int
ConvertImagePixelType(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                      std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // CONVERTIMAGEPIXELTYPE_H
