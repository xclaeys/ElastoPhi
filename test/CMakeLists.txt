
# Eigen package
find_package (Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(Tests test_algebra.cpp)

target_link_libraries (Tests  ${EIGEN3_LIBRARIES})
