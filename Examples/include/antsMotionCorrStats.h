
#ifndef ANTSMOTIONCORRSTATS_H
#define ANTSMOTIONCORRSTATS_H

namespace ants
{
extern int
antsMotionCorrStats(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                    std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // ANTSMOTIONCORRSTATS_H
