
#ifndef CREATEJACOBIANDETERMINANTIMAGE_H
#define CREATEJACOBIANDETERMINANTIMAGE_H

namespace ants
{
extern int
CreateJacobianDeterminantImage(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                               std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // CREATEJACOBIANDETERMINANTIMAGE_H
