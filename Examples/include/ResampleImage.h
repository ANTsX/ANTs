
#ifndef RESAMPLEIMAGE_H
#define RESAMPLEIMAGE_H

namespace ants
{
extern int
ResampleImage(std::vector<std::string>, // equivalent to argv of command line parameters to main()
              std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // RESAMPLEIMAGE_H
