#ifndef ANTSJOINTFUSION_H
#define ANTSJOINTFUSION_H

namespace ants
{
extern int
antsJointFusion(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // ANTSJOINTFUSION_H
