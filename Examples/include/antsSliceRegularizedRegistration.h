
#ifndef ANTSLICEREGULARIZEDREGISTRATION_H
#define ANTSLICEREGULARIZEDREGISTRATION_H

namespace ants
{
extern int
antsSliceRegularizedRegistration(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                                 std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // ANTSLICEREGULARIZEDREGISTRATION_H
