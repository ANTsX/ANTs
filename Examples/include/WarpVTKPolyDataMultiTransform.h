
#ifndef WARPVTKPOLYDATAMULTITRANSFORM_H
#define WARPVTKPOLYDATAMULTITRANSFORM_H

namespace ants
{
extern int
WarpVTKPolyDataMultiTransform(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                              std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // WARPVTKPOLYDATAMULTITRANSFORM_H
