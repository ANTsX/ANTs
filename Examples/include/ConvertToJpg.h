
#ifndef CONVERTTOJPG_H
#define CONVERTTOJPG_H

namespace ants
{
extern int
ConvertToJpg(std::vector<std::string>, // equivalent to argv of command line parameters to main()
             std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // CONVERTTOJPG_H
