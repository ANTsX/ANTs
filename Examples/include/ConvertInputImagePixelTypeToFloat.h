
#ifndef CONVERTINPUTIMAGEPIXELTYPETOFLOAT_H
#define CONVERTINPUTIMAGEPIXELTYPETOFLOAT_H

namespace ants
{
extern int
ConvertInputImagePixelTypeToFloat(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                                  std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // CONVERTINPUTIMAGEPIXELTYPETOFLOAT_H
