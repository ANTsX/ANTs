
#ifndef TENSORDERIVEDIMAGE_H
#define TENSORDERIVEDIMAGE_H

namespace ants
{
int TensorDerivedImage( std::vector<std::string>, // equivalent to argv of command line parameters to main()
                        std::ostream* out_stream  // [optional] output stream to write
                        );
} // namespace ants

#endif // TENSORDERIVEDIMAGE_H
