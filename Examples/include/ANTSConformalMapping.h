
#ifndef ANTSCONFORMALMAPPING_H
#define ANTSCONFORMALMAPPING_H

namespace ants
{
int ANTSConformalMapping( std::vector<std::string>, // equivalent to argv of command line parameters to main()
                          std::ostream* out_stream  // [optional] output stream to write
                          );
} // namespace ants

#endif // ANTSCONFORMALMAPPING_H
