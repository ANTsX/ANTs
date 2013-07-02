
#ifndef ANTSJACOBIAN_H
#define ANTSJACOBIAN_H

namespace ants
{
extern int ANTSJacobian( std::vector<std::string>,  // equivalent to argv of command line parameters to main()
                  std::ostream* out_stream   // [optional] output stream to write
                  );
} // namespace ants

#endif // ANTSJACOBIAN_H
