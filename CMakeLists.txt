cmake_minimum_required(VERSION 2.6)
Project(ElastoPhi)

# List the headers and add them into the project
FILE(GLOB_RECURSE MY_HEADERS "include/*.hpp")
include_directories(include)
message("${MY_HEADERS}")
# Add sources file
add_subdirectory(src)

# Add test files
add_subdirectory(test)


# Testing project
enable_testing()
add_test(NAME Algebra_test COMMAND Test_algebra)
add_test(NAME ACA_test COMMAND Test_ACA)


# Add documentation
add_subdirectory(doc)