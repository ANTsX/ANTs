
#ifndef EXTRACTSLICEFROMIMAGE_H
#define EXTRACTSLICEFROMIMAGE_H

namespace ants
{
extern int
ExtractSliceFromImage(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                      std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // EXTRACTSLICEFROMIMAGE_H
