
#ifndef WARPTIMESERIESIMAGEMULTITRANSFORM_H
#define WARPTIMESERIESIMAGEMULTITRANSFORM_H

namespace ants
{
extern int
WarpTimeSeriesImageMultiTransform(std::vector<std::string>, // equivalent to argv of command line parameters to
                                                            // main()
                                  std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // WARPTIMESERIESIMAGEMULTITRANSFORM_H
