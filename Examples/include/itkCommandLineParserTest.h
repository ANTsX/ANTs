
#ifndef ITKCOMMANDLINEPARSERTEST_H
#define ITKCOMMANDLINEPARSERTEST_H

namespace ants
{
extern int
itkCommandLineParserTest(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                         std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // ITKCOMMANDLINEPARSERTEST_H
