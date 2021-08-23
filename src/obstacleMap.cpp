//
// Created by qishuai on 2021/8/2.
//

#include "obstacleMap.h"
namespace casadiOptim {
    obstacleMap::obstacleMap() {
        cmap_ = std::make_shared<circular_map::CircularMap>();
        cmap_->initializeMap(6, 6, 0.3, 2.0,
                             SubFlightMode::TRACK_FREE);
        circular_map::MapPoint3f origin(0, 0, 0);
        cmap_->updateOrigin(origin);
        cmap_->resetMap();
    }

    void obstacleMap::generatePointCloud() {

        PointCloudWithInfo cloud_world;
        cloud_world.id = (CamGroup)1;
        float resol = 0.3;
        // cloud_world.timestamp = getCurrentTime();
        int size = 686;
        cloud_world.cloud.resize(4, size);
        int i = 0;
        for (float x = -1.0; x <= 1.0; x += resol) {
            for (float y = -0.5; y <= 1.0; y += resol) {
                for (float z = -2.0; z <= 2.0; z += resol) {
                    // cloud_world.cloud.conservativeResize(4, ++size);
                    // cloud_world.cloud.conservativeResize(4, ++size);
                    if (i >= size)
                        break;
                    cloud_world.cloud.col(i)(0) = x;
                    cloud_world.cloud.col(i)(1) = y;
                    cloud_world.cloud.col(i++)(2) = z;
                    // std::cout << "x: " << x << " y: " << y << " z: " << z << std::endl;
                }
            }
        }


        obstacle_cloud_ = cloud_world;
        PointCloudWithInfo obstacle_cloud;
        obstacle_cloud.id = (CamGroup) 1;
        int outer_size = 100000;
        obstacle_cloud.cloud.resize(4, outer_size);
        int outer_i = 0;

        float min_value = -11;
        float max_value = 11;

        for (float x = min_value; x <= max_value; x += resol) {
            for (float y = min_value; y <= max_value; y += resol) {
                if (outer_i >= outer_size)
                    break;
                float z = min_value;
                obstacle_cloud.cloud.col(outer_i)(0) = x;
                obstacle_cloud.cloud.col(outer_i)(1) = y;
                obstacle_cloud.cloud.col(outer_i++)(2) = z;
                z = max_value;
                obstacle_cloud.cloud.col(outer_i)(0) = x;
                obstacle_cloud.cloud.col(outer_i)(1) = y;
                obstacle_cloud.cloud.col(outer_i++)(2) = z;
            }
        }

        for (float x = min_value; x <= max_value; x += resol) {
            for (float z = min_value; z <= max_value; z += resol) {
                if (outer_i >= outer_size)
                    break;
                float y = min_value;
                obstacle_cloud.cloud.col(outer_i)(0) = x;
                obstacle_cloud.cloud.col(outer_i)(1) = y;
                obstacle_cloud.cloud.col(outer_i++)(2) = z;
                y = max_value;
                obstacle_cloud.cloud.col(outer_i)(0) = x;
                obstacle_cloud.cloud.col(outer_i)(1) = y;
                obstacle_cloud.cloud.col(outer_i++)(2) = z;
            }
        }

        for (float y = min_value; y <= max_value; y += resol) {
            for (float z = min_value; z <= max_value; z += resol) {
                if (outer_i >= outer_size)
                    break;
                float x = min_value;
                obstacle_cloud.cloud.col(outer_i)(0) = x;
                obstacle_cloud.cloud.col(outer_i)(1) = y;
                obstacle_cloud.cloud.col(outer_i++)(2) = z;
                x = max_value;
                obstacle_cloud.cloud.col(outer_i)(0) = x;
                obstacle_cloud.cloud.col(outer_i)(1) = y;
                obstacle_cloud.cloud.col(outer_i++)(2) = z;
            }
        }

        outer_cloud_world_ = obstacle_cloud;




        // circular_map::MapPoint3f origin(0, 0, 0);
        //
        // cmap_->insertPointCloud(cloud_world, origin);
        // cmap_->updateDistance();

    }
    void obstacleMap::insertPointCloud() {
        generatePointCloud();
        circular_map::MapPoint3f origin(0.0, 0.0, 0.0);
        cmap_->insertPointCloud(outer_cloud_world_, origin);
        cmap_->insertPointCloud(obstacle_cloud_, origin);
        cmap_->updateDistance();
    }

    void obstacleMap::getPointCloud(sensor_msgs::PointCloud& occupied_cloud, sensor_msgs::PointCloud& free_cloud) {
        cmap_->getMarkerOccupied(occupied_cloud);
        cmap_->getMarkerFree(free_cloud);
    }

    pathplan_msgs::VoxelMap obstacleMap::visualizeMap(float distance) {
        return cmap_->visualizeMap(distance);
    }

    void obstacleMap::getDistanceAndGrad(const circular_map::MapPoint3f& point,
                                         float& dis, Eigen::Vector3f& grad) {
        auto res = cmap_->getResolution();
        circular_map::MapPoint3f ComputerPoint = point - 0.5 * res * Eigen::Vector3f(1, 1, 1);
        auto idx = cmap_->getIdx(ComputerPoint);
        auto idx_point = cmap_->getPoint(idx);

        Eigen::Vector3f diff = (point - idx_point) / res;

        float values[2][2][2];

        for (int x = 0; x < 2; ++x) {
            for (int y = 0; y < 2; ++y) {
                for (int z = 0; z < 2; ++z) {
                    circular_map::MapPoint3i cur_idx = idx + Eigen::Vector3i(x, y, z);
                    values[x][y][z] = cmap_->getDistance(cur_idx);
                }
            }
        }

        float v00 = (1 - diff[0]) * values[0][0][0] + diff[0] * values[1][0][0];
        float v01 = (1 - diff[0]) * values[0][0][1] + diff[0] * values[1][0][1];
        float v10 = (1 - diff[0]) * values[0][1][0] + diff[0] * values[1][1][0];
        float v11 = (1 - diff[0]) * values[0][1][1] + diff[0] * values[1][1][1];

        float v0 = (1 - diff[1]) * v00 + diff[1] * v10;
        float v1 = (1 - diff[1]) * v01 + diff[1] * v11;

        dis = (1 - diff[2]) * v0 + diff[2] * v1;

        grad[2] = (v1 - v0) / res;
        grad[1] = ((1 - diff[2]) * (v10 - v00) + diff[2] * (v11 - v01)) / res;
        grad[0] = (1 - diff[2]) * (1 - diff[1]) * (values[1][0][0] - values[0][0][0]);
        grad[0] += (1 - diff[2]) * diff[1] * (values[1][1][0] - values[0][1][0]);
        grad[0] += diff[2] * (1 - diff[1]) * (values[1][0][1] - values[0][0][1]);
        grad[0] += diff[2] * diff[1] * (values[1][1][1] - values[0][1][1]);

        grad[0] / res;
    }
}