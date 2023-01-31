cmake_minimum_required(VERSION 3.24)
project(vgis10)

set(CMAKE_CXX_STANDARD 20)

include_directories(include)
include_directories(external/stb_image)

set(SRC_FILES
        src/Main.cpp
        src/FrameClass.cpp
        src/Core.cpp
        src/Logger.cpp
        )

add_executable(vgis10 ${SRC_FILES})
