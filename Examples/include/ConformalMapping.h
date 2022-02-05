
#ifndef CONFORMALMAPPING_H
#define CONFORMALMAPPING_H

namespace ants
{
extern int
ConformalMapping(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                 std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // CONFORMALMAPPING_H
