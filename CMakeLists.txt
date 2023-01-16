cmake_minimum_required(VERSION 3.16)
project(MoebiusDeform)

set(CMAKE_CXX_STANDARD 14)

message(STATUS "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")

### conditionally compile certain modules depending on libraries found on the system
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

### libIGL options:
option(LIBIGL_EMBREE           "Build target igl::embree"           ON)
option(LIBIGL_GLFW             "Build target igl::glfw"             ON)
option(LIBIGL_IMGUI            "Build target igl::imgui"            ON)
option(LIBIGL_OPENGL           "Build target igl::opengl"           ON)
option(LIBIGL_PNG              "Build target igl::png"              ON)
### set(LIBIGL_WITH_PNG CACHE BOOL ON)

### Adding libIGL and Directional: choose the path to your local copy
include(libigl)
include(libhedra)
find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

### Output directories
if(MSVC)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR})
else()
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../")
endif()



message("libhedra includes: ${LIBHEDRA_INCLUDE_DIRS}")


# Add your project files
FILE(GLOB SRCFILES *.cpp *.h)
add_executable(${PROJECT_NAME}_bin ${SRCFILES})
target_link_libraries(${PROJECT_NAME}_bin PUBLIC igl::core igl::glfw igl::opengl ${CERES_LIBRARIES})
