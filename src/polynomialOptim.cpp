//
// Created by qishuai on 2021/8/25.
//

#include "polynomialOptim.h"
namespace optim {
    polynomialOptim::polynomialOptim(int segment, int order, int derivative, std::vector<double> times):
    order_(order), segment_(segment), derivative_(derivative), time_(times) {
        dimension_ = (order + 1) * segment;

    }

    void polynomialOptim::setInitialBoundaryCondition(Eigen::Vector3d &initial_pose,
                                                  Eigen::Vector3d &initial_velocity,
                                                  Eigen::Vector3d &initial_accelerate) {
        initial_pose_ = initial_pose;
        initial_velocity_ = initial_velocity;
        initial_accelerate_ = initial_accelerate;
    }

    void polynomialOptim::setFinalBoundaryCondition(Eigen::Vector3d &final_pose,
                                                Eigen::Vector3d &final_velocity,
                                                Eigen::Vector3d &final_accelerate) {
        final_pose_ = final_pose;
        final_velocity_ = final_velocity;
        final_accelerate_ = final_accelerate;
    }


    void polynomialOptim::setVelocityConstraints(const Eigen::Vector3d &velocity_constraints) {
        velocity_constraints_ = velocity_constraints;
    }

}
