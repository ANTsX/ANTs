
#ifndef IMAGEMATH_H
#define IMAGEMATH_H

namespace ants
{
extern int
ImageMath(std::vector<std::string>, // equivalent to argv of command line parameters to main()
          std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // IMAGEMATH_H
