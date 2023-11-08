set(DEPS tbb tbbmalloc tbbmalloc_proxy Eigen3::Eigen igl::core igl::opengl_glfw igl::opengl_glfw_imgui igl::png gmp mpfr mkl_intel_lp64 mkl_sequential mkl_core)

file(GLOB HEADERS "*.h" "include/*.hpp" )
file(GLOB SOURCES "*.c" "src/*.cpp" )

add_executable(ProgrammableDigitalWeaves ${HEADERS} ${SOURCES})
target_include_directories(ProgrammableDigitalWeaves PUBLIC ../../Libs/spectra/include)
target_link_libraries(ProgrammableDigitalWeaves PUBLIC ${DEPS})
