#include <iostream>
#include <casadi/casadi.hpp>
#include "obstacleMap.h"
#include "ros/ros.h"
#include "nav_msgs/Path.h"
#include "geometry_msgs/PoseStamped.h"
#include "tf/transform_datatypes.h"
#include "polynomialOptim.h"

using namespace casadi;

ros::Publisher optimized_path_pub, occupied_pub, free_pub, distance_pub, voxel_map_pub, obstacle_point_path_pub;
geometry_msgs::PoseStamped getPoint(int segment, int order, double time,
                                    const std::vector<casadi::DM>& optimal_solution) {
    geometry_msgs::PoseStamped pt;
    int start_index = segment * (order + 1);
    for (int j = 0; j <= order; ++j) {

        // pt.pose.position.x += static_cast<double>(optimal_solution[0](start_index + j)) * casadiOptim::powInt(time, j);
        // pt.pose.position.y += static_cast<double>(optimal_solution[1](start_index + j)) * casadiOptim::powInt(time, j);
        // pt.pose.position.z += static_cast<double>(optimal_solution[2](start_index + j)) * casadiOptim::powInt(time, j);
        pt.pose.orientation = tf::createQuaternionMsgFromYaw(0);
    }
    // printf("time: %f, x: %.2f, y: %.2f, z: %.2f\n", time, pt.pose.position.x,
    //        pt.pose.position.y, pt.pose.position.z);
    return pt;
}
void getOptimizedPath(int order, const std::vector<double>& times,
                      const std::vector<casadi::DM>& optimal_solution) {
    nav_msgs::Path optimized_path;
    optimized_path.header.frame_id = "world";
    optimized_path.header.stamp = ros::Time::now();

    optimized_path.poses.clear();
    for (int i = 0; i + 1 < times.size(); ++i) {
        int segment_number = 10;
        double step_span = (times[i + 1] - times[i]) / segment_number;
        for (int time_step = 0; time_step < segment_number; time_step++) {
            double current_time = times[i] + step_span * time_step;
            geometry_msgs::PoseStamped pt = getPoint(i, order, current_time, optimal_solution);
            optimized_path.poses.push_back(pt);
        }
    }
    // for (int i = 0; i < optimized_path.poses.size(); ++i) {
    //     printf("x: %.2f, y: %.2f, z: %.2f\n", optimized_path.poses[i].pose.position.x,
    //            optimized_path.poses[i].pose.position.y,
    //            optimized_path.poses[i].pose.position.z);
    // }
    optimized_path_pub.publish(optimized_path);

}

void getObstacleCheckPath(const std::vector<Eigen::Vector3f>& path) {
    nav_msgs::Path obstacle_check_path;
    obstacle_check_path.header.frame_id = "world";
    obstacle_check_path.header.stamp = ros::Time::now();
    obstacle_check_path.poses.clear();

    for (int i = 0; i < path.size(); ++i) {
        geometry_msgs::PoseStamped pt;
        pt.pose.position.x = path[i].x();
        pt.pose.position.y = path[i].y();
        pt.pose.position.z = path[i].z();
        pt.pose.orientation = tf::createQuaternionMsgFromYaw(0);
        obstacle_check_path.poses.push_back(pt);
    }

    obstacle_point_path_pub.publish(obstacle_check_path);
}

int main(int argc, char** argv) {

    int segment = 4;
    int order = 6;
    int derivative = 4;

    int dimension = (1 + order) * segment;
    Eigen::Vector3d initial_point_pos, initial_point_vel, initial_point_acc;
    Eigen::Vector3d final_point_pos, final_point_vel, final_point_acc;
    // prepare already known info.
    std::vector<double> times = {0, 1, 3.5, 5.5, 10.0};
    Eigen::Vector3d velocity_constraint;
    initial_point_pos << -3.0, -2.1, -3.5;   // pos for initial point
    final_point_pos << 6.1, 3.2, 3.2;    // pos for final point
    initial_point_vel << 0.0, 0.0, 0.0;
    final_point_vel << 0.0, 0.0, 0.0;
    initial_point_acc << 0.0, 0.0, 0.0;
    final_point_acc << 0.0, 0.0, 0.0;

    velocity_constraint << 1.2, 1.2, 1.2;


    std::shared_ptr<optim::polynomialOptim> polynomial_instance = std::make_shared<optim::polynomialOptim>(segment, order, derivative, times);
    polynomial_instance->setInitialBoundaryCondition(initial_point_pos, initial_point_vel, initial_point_acc);
    polynomial_instance->setFinalBoundaryCondition(final_point_pos, final_point_vel, final_point_acc);
    polynomial_instance->setVelocityConstraints(velocity_constraint);

    std::shared_ptr<optim::obstacleMap> obstacleMap = std::make_shared<optim::obstacleMap>();
    obstacleMap->insertPointCloud();
    sensor_msgs::PointCloud occupied_cloud, free_cloud;
    obstacleMap->getPointCloud(occupied_cloud, free_cloud);

    pathplan_msgs::VoxelMap map = obstacleMap->visualizeMap(0.6);
    map.header.frame_id = "world";

    auto optimal_solution = polynomial_instance->optimize();

    ros::init(argc, argv, "optim_test");
    ros::NodeHandle nh;
    optimized_path_pub = nh.advertise<nav_msgs::Path>("/optim_test/optimizedPath", 5, false);
    occupied_pub =
            nh.advertise<sensor_msgs::PointCloud>("/circular_map/occupied", 1, false);
    // getOptimizedPath(order, times, optimal_solution);
    distance_pub =
        nh.advertise<sensor_msgs::PointCloud>("/circular_map/distance", 1, false);
    voxel_map_pub = nh.advertise<pathplan_msgs::VoxelMap>("/circular_map/voxel_map", 1);
    obstacle_point_path_pub = nh.advertise<nav_msgs::Path>("/optim_test/obstacle_check_path", 5, false);



    ros::Rate rate(10);
    while(ros::ok()) {
        // std::cout << "test" << std::endl;
        occupied_pub.publish(occupied_cloud);
        voxel_map_pub.publish(map);
        // getOptimizedPath(order, times, optimal_solution);
        std::cout << "begin get check path" << std::endl;
        // auto pose = casadiInstance.getSegmentCost(three_optimal_solution);
        // getObstacleCheckPath(pose);
        ros::spinOnce();
    }

    return 0;
}