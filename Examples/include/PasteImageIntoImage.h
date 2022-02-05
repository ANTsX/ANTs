
#ifndef PASTEIMAGEINTOIMAGE_H
#define PASTEIMAGEINTOIMAGE_H

namespace ants
{
extern int
PasteImageIntoImage(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                    std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // PASTEIMAGEINTOIMAGE_H
