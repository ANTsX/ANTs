
#ifndef ANTSIMAGEWRITE_H
#define ANTSIMAGEWRITE_H

#include "itkImage.h"

namespace ants
{
extern int
antsImageWrite(itk::Image<float, 3>::Pointer, // image to write
               std::string                    // filename of the target file
);
} // namespace ants

#endif // ANTSIMAGEWRITE_H
