
#ifndef GETCONNECTEDCOMPONENTSFEATUREIMAGES_H
#define GETCONNECTEDCOMPONENTSFEATUREIMAGES_H

namespace ants
{
extern int
GetConnectedComponentsFeatureImages(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                                    std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // GETCONNECTEDCOMPONENTSFEATUREIMAGES_H
