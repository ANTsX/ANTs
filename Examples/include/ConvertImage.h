
#ifndef CONVERTIMAGE_H
#define CONVERTIMAGE_H

namespace ants
{
extern int
ConvertImage(std::vector<std::string>, // equivalent to argv of command line parameters to main()
             std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // TILEIMAGES_H
