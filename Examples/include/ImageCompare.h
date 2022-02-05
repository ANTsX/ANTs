
#ifndef IMAGECOMPARE_H
#define IMAGECOMPARE_H

namespace ants
{
extern int
ImageCompare(std::vector<std::string>, // equivalent to argv of command line parameters to main()
             std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // IMAGECOMPARE_H
