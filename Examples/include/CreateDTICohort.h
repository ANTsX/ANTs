
#ifndef CREATEDTICOHORT_H
#define CREATEDTICOHORT_H

namespace ants
{
int CreateDTICohort( std::vector<std::string>, // equivalent to argv of command line parameters to main()
                     std::ostream* out_stream  // [optional] output stream to write
                     );
} // namespace ants

#endif // CREATEDTICOHORT_H
