#ifndef ANTSSURF_H
#define ANTSSURF_H

namespace ants
{
extern int
antsSurf(std::vector<std::string>, // equivalent to argv of command line parameters to main()
         std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // ANTSSURF_H
