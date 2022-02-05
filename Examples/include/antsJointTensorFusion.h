#ifndef ANTSJOINTTENSORFUSION_H
#define ANTSJOINTTENSORFUSION_H

namespace ants
{
extern int
antsJointTensorFusion(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                      std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // ANTSJOINTTENSORFUSION_H
