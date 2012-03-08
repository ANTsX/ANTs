
find_package(SlicerExecutionModel REQUIRED GenerateCLP)
include(${GenerateCLP_USE_FILE})

set(antsRegistrationCLP_SOURCE antsRegistrationCLP.cxx)

generateclp(antsRegistrationCLP_SOURCE antsREGISTRATIONCLP.xml)

add_executable(antsRegistrationCLP ${antsRegistrationCLP_SOURCE})

target_link_libraries(antsRegistrationCLP ${ITK_LIBRARIES})
