
#ifndef SETDIRECTIONBYMATRIX_H
#define SETDIRECTIONBYMATRIX_H

namespace ants
{
extern int
SetDirectionByMatrix(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                     std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // SETDIRECTIONBYMATRIX_H
