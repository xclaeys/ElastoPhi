cmake_minimum_required(VERSION 2.6)
Project(ElastoPhi)

#Add the headers into the project
include_directories(include)

#Add sources file
add_subdirectory(src)

#Add test files
add_subdirectory(test)


# Testing project
enable_testing()
add_test(NAME Algebra_test COMMAND Test_algebra)
add_test(NAME ACA_test COMMAND Test_ACA)
