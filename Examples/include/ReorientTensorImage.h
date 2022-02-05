
#ifndef REORIENTTENSORIMAGE_H
#define REORIENTTENSORIMAGE_H

namespace ants
{
extern int
ReorientTensorImage(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                    std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // REORIENTTENSORIMAGE_H
