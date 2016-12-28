cmake_minimum_required(VERSION 3.0)

add_subdirectory(.. zbessel.dir)

find_library(SLATEC_LIBRARY slatec)
if(NOT SLATEC_LIBRARY)
  message(FATAL_ERROR "the SLATEC library is required for tests")
endif()

add_executable(test-zbessel test.cc)
target_link_libraries(test-zbessel zbessel "${SLATEC_LIBRARY}")
