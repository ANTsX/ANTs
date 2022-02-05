
#ifndef N4BIASFIELDCORRECTION_H
#define N4BIASFIELDCORRECTION_H

namespace ants
{
extern int
N4BiasFieldCorrection(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                      std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // N4BIASFIELDCORRECTION_H
