
# Creation de l executable main a partir de prog
add_executable(main prog.cpp ${MY_HEADERS}) 
target_link_libraries (main  ${EIGEN3_LIBRARIES})

# Creation de l executable compression a partir de compress
add_executable(Compression compress.cpp)
target_link_libraries (Compression ${EIGEN3_LIBRARIES})

# Creation de l executable Test_ACA a partir de test_ACA
add_executable(Test_ACA test_ACA.cpp) 
target_link_libraries (Test_ACA  ${EIGEN3_LIBRARIES})
