
#ifndef ANTS_MOCO_H
#define ANTS_MOCO_H

namespace ants
{
extern int
ants_moco(std::vector<std::string>, // equivalent to argv of command line parameters to main()
          std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // ANTS_MOCO_H
