
#ifndef AVERAGEAFFINETRANSFORMNORIGID_H
#define AVERAGEAFFINETRANSFORMNORIGID_H

namespace ants
{
extern int
AverageAffineTransformNoRigid(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                              std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // AVERAGEAFFINETRANSFORMNORIGID_H
