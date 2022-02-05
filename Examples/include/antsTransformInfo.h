
#ifndef ANTSTRANSFORMINFO_H
#define ANTSTRANSFORMINFO_H

namespace ants
{
extern int
antsTransformInfo(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                  std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // ANTSMOTIONCORRSTATS_H
