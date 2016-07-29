
# Creation du test Test_algebra a partir de test_algebra
add_executable(Test_algebra test_algebra.cpp)
target_link_libraries (Test_algebra  ${EIGEN3_LIBRARIES})

# Creation du test Test_cluster a partir de test_cluster
add_executable(Test_cluster test_cluster.cpp)
target_link_libraries (Test_cluster  ${EIGEN3_LIBRARIES})

# Creation du test Test_lrmatrix a partir de test_lrmatrix
add_executable(Test_lrmatrix test_lrmatrix.cpp)
target_link_libraries (Test_lrmatrix  ${EIGEN3_LIBRARIES})

# Creation du test Test_hmatrix a partir de test_hmatrix
add_executable(Test_hmatrix test_hmatrix.cpp)
target_link_libraries (Test_hmatrix  ${EIGEN3_LIBRARIES})
