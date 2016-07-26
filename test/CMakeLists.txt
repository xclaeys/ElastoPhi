
# Eigen package
find_package (Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(Test_algebra test_algebra.cpp)
add_executable(Test_ACA test_ACA.cpp)

target_link_libraries (Test_algebra  ${EIGEN3_LIBRARIES})
target_link_libraries (Test_ACA  ${EIGEN3_LIBRARIES})
