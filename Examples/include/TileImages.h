
#ifndef TILEIMAGES_H
#define TILEIMAGES_H

namespace ants
{
extern int
TileImages(std::vector<std::string>, // equivalent to argv of command line parameters to main()
           std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // TILEIMAGES_H
