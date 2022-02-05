
#ifndef FITBSPLINECURVETOPOINTS_H
#define FITBSPLINECURVETOPOINTS_H

namespace ants
{
extern int
FitBSplineToPoints(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                   std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // FITBSPLINECURVETOPOINTS_H
