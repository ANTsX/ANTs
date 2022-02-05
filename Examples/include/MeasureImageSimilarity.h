
#ifndef MEASUREIMAGESIMILARITY_H
#define MEASUREIMAGESIMILARITY_H

namespace ants
{
extern int
MeasureImageSimilarity(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                       std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // MEASUREIMAGESIMILARITY_H
