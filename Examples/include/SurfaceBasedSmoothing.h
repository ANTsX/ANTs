
#ifndef SURFACEBASEDSMOOTHING_H
#define SURFACEBASEDSMOOTHING_H

namespace ants
{
extern int
SurfaceBasedSmoothing(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                      std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // SURFACEBASEDSMOOTHING_H
