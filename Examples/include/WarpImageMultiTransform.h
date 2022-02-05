
#ifndef WARPIMAGEMULTITRANSFORM_H
#define WARPIMAGEMULTITRANSFORM_H

namespace ants
{
extern int
WarpImageMultiTransform(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                        std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // WARPIMAGEMULTITRANSFORM_H
