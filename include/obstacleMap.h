//
// Created by qishuai on 2021/8/2.
//

#ifndef CASADI_TEST_OBSTACLECOST_H
#define CASADI_TEST_OBSTACLECOST_H

#include "circular_map/circular_map.h"
#include "common_files/common_def.h"
#include "common_files/time.h"
#include "sensor_msgs/PointCloud.h"
#include <pathplan_msgs/VoxelMap.h>

namespace optim {
    class obstacleMap {
    public:
        obstacleMap();
        void insertPointCloud();
        void getPointCloud(sensor_msgs::PointCloud& occupied_cloud, sensor_msgs::PointCloud& free_cloud);
        pathplan_msgs::VoxelMap visualizeMap(float dis_threshold);
        float getDistance(const circular_map::MapPoint3f& point, float& dis) { dis = cmap_->getDistance(point); };
        void getDistanceAndGrad(const circular_map::MapPoint3f& point, float& dis, Eigen::Vector3f& grad);
        void setGrad(const std::vector<Eigen::Vector3f>& grad) {
            grad_ = grad;
        }
        void getGrad(std::vector<Eigen::Vector3f>& grad) {
            grad = grad_;
        }
    private:
        void generatePointCloud();
    private:
        std::shared_ptr<circular_map::CircularMap> cmap_;
        PointCloudWithInfo outer_cloud_world_;
        PointCloudWithInfo obstacle_cloud_;
        std::vector<Eigen::Vector3f> grad_;
    };
}

#endif //CASADI_TEST_OBSTACLECOST_H
