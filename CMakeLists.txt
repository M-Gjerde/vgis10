cmake_minimum_required(VERSION 3.24)
project(vgis10)

set(CMAKE_CXX_STANDARD 20)

include_directories(include)
include_directories(external/stb_image)

find_package(OpenCV REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${OpenCV_INCLUDE_DIRS})

include(cmake/Common.cmake)

set(SRC_FILES
        src/Main.cpp
        src/FrameClass.cpp
        src/Core.cpp
        src/Logger.cpp
        src/CameraCalibration.cpp
        src/PointSelection.cpp)

add_executable(vgis10 ${SRC_FILES})
target_link_libraries(vgis10 fmt::fmt ${OpenCV_LIBS} Eigen3::Eigen sophus)

add_custom_command(TARGET vgis10
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/calibration ${CMAKE_BINARY_DIR}/calibration
        )

message(${CMAKE_SOURCE_DIR}/calibration)