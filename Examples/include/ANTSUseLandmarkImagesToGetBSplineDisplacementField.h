
#ifndef ANTSUSELANDMARKIMAGESTOGETBSPLINEDISPLACEMENTFIELD_H
#define ANTSUSELANDMARKIMAGESTOGETBSPLINEDISPLACEMENTFIELD_H

namespace ants
{
extern int
ANTSUseLandmarkImagesToGetBSplineDisplacementField(std::vector<std::string>, // equivalent to argv of command line
                                                                             // parameters to main()
                                                   std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // ANTSUSELANDMARKIMAGESTOGETBSPLINEDISPLACEMENTFIELD_H
