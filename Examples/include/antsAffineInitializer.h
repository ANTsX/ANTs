
#ifndef ANTSAFFINEINITIALIZER_H
#define ANTSAFFINEINITIALIZER_H

namespace ants
{
extern int
antsAffineInitializer(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                      std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // ANTSAFFINEINITIALIZER_H
