
#ifndef AVERAGEAFFINETRANSFORM_H
#define AVERAGEAFFINETRANSFORM_H

namespace ants
{
extern int
AverageAffineTransform(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                       std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // AVERAGEAFFINETRANSFORM_H
