
#ifndef ANTSAPPLYTRANSFORMSTOGIFTI_H
#define ANTSAPPLYTRANSFORMSTOGIFTI_H

namespace ants
{
extern int
antsApplyTransformsToGifti(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                            std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // ANTSAPPLYTRANSFORMSTOGIFTI_H
