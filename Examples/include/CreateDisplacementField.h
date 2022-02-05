
#ifndef CREATEDISPLACEMENTFIELD_H
#define CREATEDISPLACEMENTFIELD_H

namespace ants
{
extern int
CreateDisplacementField(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                        std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // CREATEDISPLACEMENTFIELD_H
