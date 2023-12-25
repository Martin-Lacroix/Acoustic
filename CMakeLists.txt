cmake_minimum_required(VERSION 3.10)
project(Acoustic)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE SOURCE_FILES "source/*.cpp")

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

find_package(Gmsh REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(include)
include_directories(${GMSH_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIR})

target_link_libraries(${PROJECT_NAME} ${GMSH_LIBRARIES})
