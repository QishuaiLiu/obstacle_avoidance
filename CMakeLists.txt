cmake_minimum_required(VERSION 3.5)
project(optimized_test)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(catkin REQUIRED
        roscpp
        rospy
        std_msgs
        nav_msgs
        circular_map)

catkin_package()

include_directories(include ${catkin_INCLUDE_DIRS})

include_directories(/usr/local/include)
include_directories(include)
find_package(casadi REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories( ${EIGEN3_INCLUDE_DIR} )


add_executable(${PROJECT_NAME} src/main.cpp src/casadiOptim.cpp src/obstacleMap.cpp)
target_link_libraries(${PROJECT_NAME} casadi ${catkin_LIBRARIES}
)
