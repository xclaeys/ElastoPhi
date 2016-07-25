
# Eigen package
find_package (Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR}) 

add_executable(main cluster.hpp hmatrix.hpp  lrmat.hpp  matrix.hpp  point.hpp  prog.cpp export.hpp loading.hpp) 

target_link_libraries (main  ${EIGEN3_LIBRARIES})
