
#ifndef DENRRD_H
#define DENRRD_H

namespace ants
{
extern int
DeNrrd(std::vector<std::string>, // equivalent to argv of command line parameters to main()
       std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // DENRRD_H
