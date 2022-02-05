
#ifndef MULTIPLYIMAGES_H
#define MULTIPLYIMAGES_H

namespace ants
{
extern int
MultiplyImages(std::vector<std::string>, // equivalent to argv of command line parameters to main()
               std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // MULTIPLYIMAGES_H
