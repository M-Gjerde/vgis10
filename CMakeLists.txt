cmake_minimum_required(VERSION 3.24)
project(vgis10)

set(CMAKE_CXX_STANDARD 20)

include(cmake/Common.cmake)

include_directories(include)
include_directories(external/stb_image)

find_package(OpenCV REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${OpenCV_INCLUDE_DIRS})


set(SRC_FILES
        src/Main.cpp
        src/FrameClass.cpp
        src/Core.cpp
        src/Logger.cpp
        src/CameraCalibration.cpp
        src/PointSelection.cpp
        src/Tracker.cpp
        src/Optimized/AccumulatedSCHessian.cpp
        src/Optimized/AccumulatedTopHessian.cpp
        src/Optimized/Energy.cpp
        src/ImmaturePoint.cpp
        src/TrackerOptimize.cpp
        src/CoarseTracker.cpp
        src/Optimized/PointFrameResidual.cpp
        src/PointHessian.cpp
        src/Optimized/EFResidual.cpp
        src/Frame.cpp)

add_executable(vgis10 ${SRC_FILES})
target_link_libraries(vgis10 fmt::fmt ${OpenCV_LIBS} Eigen3::Eigen sophus yaml-cpp::yaml-cpp -lboost_thread)

add_custom_command(TARGET vgis10
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/calibration ${CMAKE_BINARY_DIR}/calibration
        )

message(${CMAKE_SOURCE_DIR}/calibration)