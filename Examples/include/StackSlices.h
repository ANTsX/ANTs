
#ifndef STACKSLICES_H
#define STACKSLICES_H

namespace ants
{
extern int
StackSlices(std::vector<std::string>, // equivalent to argv of command line parameters to main()
            std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // STACKSLICES_H
