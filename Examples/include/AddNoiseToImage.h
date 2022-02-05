
#ifndef ADDNOISETOIMAGE_H
#define ADDNOISETOIMAGE_H

namespace ants
{
extern int
AddNoiseToImage(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // ADDNOISETOIMAGE_H
