#ifndef ANTSVOL_H
#define ANTSVOL_H

namespace ants
{
extern int
antsVol(std::vector<std::string>, // equivalent to argv of command line parameters to main()
        std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // ANTVOL_H
