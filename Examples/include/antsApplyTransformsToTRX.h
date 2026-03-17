#ifndef ANTSAPPLYTRANSFORMSTOTRIX_H
#define ANTSAPPLYTRANSFORMSTOTRIX_H

namespace ants
{
extern int
antsApplyTransformsToTRX(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                          std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // ANTSAPPLYTRANSFORMSTOTRIX_H
