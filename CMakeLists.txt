cmake_minimum_required(VERSION 3.28)
project(V_obc_NCell_inv_12_6)

set(CMAKE_CXX_STANDARD 17)

add_executable(run_mc run_mc.cpp mc_subroutine/mc_read_load_compute.cpp potentialFunction/potentialFunctions.cpp)

add_executable(main main.cpp potentialFunction/potentialFunctions.cpp mc_subroutine/mc_read_load_compute.cpp)

set_source_files_properties(
        mc_subroutine/mc_read_load_compute.cpp
        PROPERTIES

        COMPILE_OPTIONS "-O3"
)

set_source_files_properties(
        potentialFunction/potentialFunctions.cpp
        PROPERTIES

        COMPILE_OPTIONS "-O3"
)
set(Python3_INCLUDE_DIRS "/usr/include/python3.10")
set(Python3_LIBRARIES "/usr/lib/x86_64-linux-gnu/libpython3.10.so")
include_directories(${Python3_INCLUDE_DIRS})

set(Boost_NO_SYSTEM_PATHS YES)
if (Boost_NO_SYSTEM_PATHS)
    set(Boost_ROOT "/home/polya/cppLocal/boostBuild")
    set(Boost_INCLUDE_DIR "/home/polya/cppLocal/boostBuild/include")
    set(Boost_LIBRARY_DIR "/home/polya/cppLocal/boostBuild/lib")
endif (Boost_NO_SYSTEM_PATHS)


find_package(Boost   REQUIRED COMPONENTS filesystem python${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR})

if(Boost_FOUND)
    message(STATUS "boost root is: ${Boost_ROOT}")
    message(STATUS "boost include path is : ${Boost_INCLUDE_DIR}")
    message(STATUS "boost library path is : ${Boost_LIBRARY_DIR}")
    message(STATUS "boost libraries is : ${Boost_LIBRARIES}")
    include_directories(${Boost_INCLUDE_DIR})
    link_directories(${Boost_LIBRARY_DIR})
else()
    message(WARNING "boost not found.")
endif()

target_link_libraries(run_mc Boost::filesystem)
#target_link_libraries(run_mc Boost::serialization)
target_link_libraries(run_mc Boost::python)
# Find OpenMP
#find_package(OpenMP REQUIRED)
#
## Link the OpenMP libraries
#if(OpenMP_CXX_FOUND)
#    target_link_libraries(run_mc  OpenMP::OpenMP_CXX)
#endif()


#python
target_link_libraries(run_mc ${Python3_LIBRARIES})


