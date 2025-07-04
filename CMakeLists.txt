cmake_minimum_required(VERSION 3.13) 

set(PROJECT_NAME project)
project(${PROJECT_NAME})

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Configs" FORCE)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)  


foreach(CONFIG ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${CONFIG} CONFIG)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${CONFIG} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endforeach()

set(LIBRARY "${PROJECT_NAME}")
set(TESTS   "test_${PROJECT_NAME}")
set(INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/include")


find_package(MPI REQUIRED)


include_directories("${INCLUDE}" gtest ${MPI_INCLUDE_PATH})


add_subdirectory(src)  
add_subdirectory(samples)
add_subdirectory(gtest)
add_subdirectory(test)

# REPORT
message( STATUS "")
message( STATUS "General configuration for ${PROJECT_NAME}")
message( STATUS "======================================")
message( STATUS "")
message( STATUS "   Configuration: ${CMAKE_BUILD_TYPE}")
message( STATUS "   MPI include path: ${MPI_INCLUDE_PATH}")
message( STATUS "   MPI libraries: ${MPI_LIBRARIES}")
message( STATUS "")