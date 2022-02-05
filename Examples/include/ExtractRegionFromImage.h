
#ifndef EXTRACTREGIONFROMIMAGE_H
#define EXTRACTREGIONFROMIMAGE_H

namespace ants
{
extern int
ExtractRegionFromImage(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                       std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // EXTRACTREGIONFROMIMAGE_H
