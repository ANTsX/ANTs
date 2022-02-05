
#ifndef WARPTENSORIMAGEMULTITRANSFORM_H
#define WARPTENSORIMAGEMULTITRANSFORM_H

namespace ants
{
extern int
WarpTensorImageMultiTransform(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                              std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // WARPTENSORIMAGEMULTITRANSFORM_H
