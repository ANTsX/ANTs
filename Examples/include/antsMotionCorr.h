
#ifndef ANTSMOTIONCORR_H
#define ANTSMOTIONCORR_H

namespace ants
{
extern int
antsMotionCorr(std::vector<std::string>, // equivalent to argv of command line parameters to main()
               std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // ANTSMOTIONCORR_H
