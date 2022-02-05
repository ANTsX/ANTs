
#ifndef MEMORYTEST_H
#define MEMORYTEST_H

namespace ants
{
extern int
MemoryTest(std::vector<std::string>, // equivalent to argv of command line parameters to main()
           std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // MEMORYTEST_H
