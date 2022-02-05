
#ifndef ANTSAPPLYTRANSFORMS_H
#define ANTSAPPLYTRANSFORMS_H

namespace ants
{
extern int
antsApplyTransforms(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                    std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // ANTSAPPLYTRANSFORMS_H
