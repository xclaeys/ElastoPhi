
# Creation de l executable main a partir de prog
add_executable(main prog.cpp ${MY_HEADERS}) 
target_link_libraries (main  ${EIGEN3_LIBRARIES})

# Creation de l executable Visu a partir de visualization
add_executable(Visu visualization.cpp)
target_link_libraries (Visu ${EIGEN3_LIBRARIES})

# Creation de l executable Compression a partir de compress
add_executable(Compression compress.cpp)
target_link_libraries (Compression ${EIGEN3_LIBRARIES})

# Creation de l executable MultiCompression a partir de multi_compress.cpp
add_executable(MultiCompression multi_compress.cpp)
target_link_libraries (MultiCompression ${EIGEN3_LIBRARIES})

# Creation de l executable Test_ACA a partir de test_ACA
add_executable(Test_ACA test_ACA.cpp) 
target_link_libraries (Test_ACA  ${EIGEN3_LIBRARIES})
