
#ifndef PERMUTEFLIPIMAGEORIENTATIONAXES_H
#define PERMUTEFLIPIMAGEORIENTATIONAXES_H

namespace ants
{
extern int
PermuteFlipImageOrientationAxes(std::vector<std::string>, // equivalent to argv of command line parameters to
                                                          // main()
                                std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // PERMUTEFLIPIMAGEORIENTATIONAXES_H
