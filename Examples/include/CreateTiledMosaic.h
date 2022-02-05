
#ifndef CREATETILEDMOSAIC_H
#define CREATETILEDMOSAIC_H

namespace ants
{
extern int
CreateTiledMosaic(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                  std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // CREATETILEDMOSAIC_H
