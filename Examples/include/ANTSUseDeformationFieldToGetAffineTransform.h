
#ifndef ANTSUSEDEFORMATIONFIELDTOGETAFFINETRANSFORM_H
#define ANTSUSEDEFORMATIONFIELDTOGETAFFINETRANSFORM_H

namespace ants
{
extern int
ANTSUseDeformationFieldToGetAffineTransform(std::vector<std::string>, // equivalent to argv of command line
                                                                      // parameters to main()
                                            std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // ANTSUSEDEFORMATIONFIELDTOGETAFFINETRANSFORM_H
