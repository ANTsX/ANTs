
#ifndef CREATEIMAGE_H
#define CREATEIMAGE_H

namespace ants
{
extern int
CreateImage(std::vector<std::string>, // equivalent to argv of command line parameters to main()
            std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // CREATEIMAGE_H
