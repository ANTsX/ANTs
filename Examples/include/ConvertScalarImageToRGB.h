
#ifndef CONVERTSCALARIMAGETORGB_H
#define CONVERTSCALARIMAGETORGB_H

namespace ants
{
extern int
ConvertScalarImageToRGB(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                        std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // CONVERTSCALARIMAGETORGB_H
