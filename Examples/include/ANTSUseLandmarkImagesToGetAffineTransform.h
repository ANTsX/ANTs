
#ifndef ANTSUSELANDMARKIMAGESTOGETAFFINETRANSFORM_H
#define ANTSUSELANDMARKIMAGESTOGETAFFINETRANSFORM_H

namespace ants
{
extern int
ANTSUseLandmarkImagesToGetAffineTransform(std::vector<std::string>, // equivalent to argv of command line
                                                                    // parameters to main()
                                          std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // ANTSUSELANDMARKIMAGESTOGETAFFINETRANSFORM_H
