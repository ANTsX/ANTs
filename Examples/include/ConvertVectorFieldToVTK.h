
#ifndef CONVERTVECTORFIELDTOVTK_H
#define CONVERTVECTORFIELDTOVTK_H

namespace ants
{
extern int
ConvertVectorFieldToVTK(std::vector<std::string>, // equivalent to argv of command line parameters to main()
                        std::ostream * out_stream // [optional] output stream to write
);
} // namespace ants

#endif // CONVERTVECTORFIELDTOVTK_H
