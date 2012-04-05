
#ifndef ANTSINTEGRATEVECTORFIELD_H
#define ANTSINTEGRATEVECTORFIELD_H

namespace ants
{
int ANTSIntegrateVectorField( std::vector<std::string>, // equivalent to argv of command line parameters to main()
                              std::ostream* out_stream  // [optional] output stream to write
                              );
} // namespace ants

#endif // ANTSINTEGRATEVECTORFIELD_H
