
#ifndef ANTSALIGNORIGIN_H
#define ANTSALIGNORIGIN_H

namespace ants
{
extern int
antsAlignOrigin(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // ANTSLIGNORIGIN_H
