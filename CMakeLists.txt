set(CMAKE_CXX_STANDARD 11)
cmake_minimum_required(VERSION 3.3)
project(Animation)

INCLUDE(FindPkgConfig)
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
PKG_SEARCH_MODULE(SDL2 REQUIRED sdl2)
PKG_SEARCH_MODULE(SDL2IMAGE REQUIRED SDL2_image>=2.0.0)

include_directories(include)
file( GLOB gKit_SOURCES gKit/*.cpp )
file( GLOB SOURCES src/*.cpp )
set(SOURCE_FILES ${gKit_SOURCES} ${SOURCES})

include_directories(gKit ${SDL2_INCLUDE_DIRS} ${SDL2IMAGE_INCLUDE_DIRS} gKit)
add_executable(Animation ${SOURCE_FILES})
target_link_libraries(Animation ${OPENGL_LIBRARIES} ${SDL2_LIBRARIES} ${SDL2IMAGE_LIBRARIES} ${GLEW_LIBRARY})