
#ifndef COMPOSEMULTITRANSFORM_H
#define COMPOSEMULTITRANSFORM_H

namespace ants
{
extern int
ComposeMultiTransform(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                      std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // COMPOSEMULTITRANSFORM_H
