
# Creation de l executable main a partir de prog
add_executable(main prog.cpp ${MY_HEADERS})
target_link_libraries (main  ${EIGEN3_LIBRARIES})

# Creation de l executable CompaSparse a partir de comparison_sparse
add_executable(CompaSparse comparison_sparse.cpp)
target_link_libraries (CompaSparse ${EIGEN3_LIBRARIES})

# Creation de l executable Compression a partir de compress
add_executable(Compression compress.cpp)
target_link_libraries (Compression ${EIGEN3_LIBRARIES})

# Creation de l executable MultiCompression a partir de multi_compress.cpp
add_executable(MultiCompression multi_compress.cpp)
target_link_libraries (MultiCompression ${EIGEN3_LIBRARIES})

# Creation de l executable Test_ACA a partir de test_ACA
add_executable(Test_ACA test_ACA.cpp)
target_link_libraries (Test_ACA  ${EIGEN3_LIBRARIES})

# Creation de l executable ErrDecrease a partir de err_decrease
add_executable(ErrDecrease err_decrease.cpp)
target_link_libraries (ErrDecrease ${EIGEN3_LIBRARIES})

# Creation de l executable VisuCluster a partir de visucluster
add_executable(VisuCluster visucluster.cpp)
target_link_libraries (VisuCluster ${EIGEN3_LIBRARIES})

# Creation de l executable VisuMatrix a partir de visumatrix
add_executable(VisuMatrix visumatrix.cpp)
target_link_libraries (VisuMatrix ${EIGEN3_LIBRARIES})

# Creation de l executable VisuMesh a partir de visumesh
add_executable(VisuMesh visumesh.cpp)
target_link_libraries (VisuMesh ${EIGEN3_LIBRARIES})
