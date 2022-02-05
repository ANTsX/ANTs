
#ifndef COMPUTESIMILARITYMETRIC_H
#define COMPUTESIMILARITYMETRIC_H

namespace ants
{
extern int
ComputeSimilarityMetric(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                        std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // COMPUTESIMILARITYMETRIC_H
