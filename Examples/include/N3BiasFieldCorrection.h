
#ifndef N3BIASFIELDCORRECTION_H
#define N3BIASFIELDCORRECTION_H

namespace ants
{
extern int
N3BiasFieldCorrection(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                      std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // N3BIASFIELDCORRECTION_H
