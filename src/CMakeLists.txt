
# Eigen package
find_package (Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR}) 

add_executable(main prog.cpp ${MY_HEADERS}) 

target_link_libraries (main  ${EIGEN3_LIBRARIES})
