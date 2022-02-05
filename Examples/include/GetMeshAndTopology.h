
#ifndef GETMESHANDTOPOLOGY_H
#define GETMESHANDTOPOLOGY_H

namespace ants
{
extern int
GetMeshAndTopology(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                   std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // GETMESHANDTOPOLOGY_H
