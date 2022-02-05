
#ifndef ANTSIMAGEREAD_H
#define ANTSIMAGEREAD_H

#include "itkImage.h"

namespace ants
{
extern itk::Image<float, 3>::Pointer antsImageRead(std::string // filename of the image to be read
);
} // namespace ants

#endif // ANTSIMAGEREAD_H
