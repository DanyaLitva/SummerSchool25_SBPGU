set(target ${LIBRARY})

file(GLOB hdrs "*.h*" "${INCLUDE}/*.h*")
file(GLOB srcs "*.cpp")

add_library(${target} STATIC ${srcs} ${hdrs})

target_include_directories(${target} PRIVATE ${MPI_INCLUDE_PATH})
target_link_libraries(${target} ${MPI_LIBRARIES} ${LIBRARY_DEPS})

