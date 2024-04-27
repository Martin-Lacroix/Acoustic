cmake_minimum_required(VERSION 3.10)
project(Acoustic)

# Make CXX 11 required

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the source files in the project

file(GLOB_RECURSE SOURCE_FILES "source/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCE_FILES})

# Find Eigen and GMSH

find_package(Eigen3 REQUIRED)
find_library(GMSH_LIBRARIES "gmsh" HINTS $ENV{GMSH_LIBRARIES} NO_CACHE)
find_path(GMSH_INCLUDE_DIRS "gmsh.h" HINTS $ENV{GMSH_INCLUDE_DIRS} NO_CACHE)

# Add the include durectories

include_directories(include)
include_directories(${GMSH_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIR})

# Link the libraries to the project

target_link_libraries(${PROJECT_NAME} ${GMSH_LIBRARIES})
