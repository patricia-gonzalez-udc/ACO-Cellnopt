cmake_minimum_required(VERSION 3.0)

# define a project name
project(aco-cellnopt)

# enforce the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(HIGHFIVE_USE_BOOST OFF)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

add_subdirectory(HighFive)

include_directories(src)
include_directories(aco)
set(SOURCE_FILES src/vector_funcs.hpp 
                 src/data.hpp
                 src/compute_score_t1.hpp
                 src/simulator.hpp
                 src/get_fit.hpp)

set(SOURCE_ACO  aco/ants.cpp
		aco/utilities_aco.cpp
		aco/parallel.cpp
		aco/aco.hpp) 

add_library(aco_lib STATIC ${SOURCE_ACO})
add_executable(aco aco/cellnopt_aco.cpp ${SOURCE_FILES})
target_include_directories(aco PRIVATE ${MPI_C_INCLUDE_PATH})
target_compile_options(aco PRIVATE ${OpenMP_CXX_FLAGS})
target_link_libraries(aco PRIVATE aco_lib HighFive ${OpenMP_CXX_FLAGS}
	${MPI_C_LIBRARIES})
