#ifndef ANTSAI_H
#define ANTSAI_H

namespace ants
{
extern int
antsAI(std::vector<std::string>, // equivalent to argv of command line parameters to main()
       std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // ANTSAI_H
