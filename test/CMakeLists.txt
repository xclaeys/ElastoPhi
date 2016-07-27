
add_executable(Test_algebra test_algebra.cpp)
add_executable(Test_lrmatrix test_lrmatrix.cpp)

target_link_libraries (Test_algebra  ${EIGEN3_LIBRARIES})
target_link_libraries (Test_lrmatrix  ${EIGEN3_LIBRARIES})
