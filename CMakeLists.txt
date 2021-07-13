cmake_minimum_required(VERSION 3.10.2)
project(nmp)



find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS system thread)
find_package(OpenMP)
find_package(igraph REQUIRED)

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_definitions(${EIGEN3_DEFINITIONS})

include_directories(
    include
    ${EIGEN3_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

# target link libraries's target == target of add_executable

add_executable(test_ 
src/utils.cpp 
src/main.cpp 
src/factor_graph.cpp
src/dbscan.cpp
)


target_link_libraries(test_ ${Boost_LIBRARIES} igraph::igraph)