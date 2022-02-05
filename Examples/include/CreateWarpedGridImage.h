
#ifndef CREATEWARPEDGRIDIMAGE_H
#define CREATEWARPEDGRIDIMAGE_H

namespace ants
{
extern int
CreateWarpedGridImage(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                      std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // CREATEWARPEDGRIDIMAGE_H
