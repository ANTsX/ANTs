
#ifndef REBASETENSORIMAGE_H
#define REBASETENSORIMAGE_H

namespace ants
{
extern int
RebaseTensorImage(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                  std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // REBASETENSORIMAGE_H
