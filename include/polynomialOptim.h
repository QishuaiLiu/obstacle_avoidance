//
// Created by qishuai on 2021/8/25.
//

#ifndef SRC_POLYNOMIALOPTIM_H
#define SRC_POLYNOMIALOPTIM_H
#include <vector>
#include <Eigen/Eigen>
#include <nlopt.hpp>
#include <memory>

namespace optim {

    class polynomialOptim {
    public:
        polynomialOptim(int segment, int order, int derivative, std::vector<double> times);
        void setInitialBoundaryCondition(Eigen::Vector3d& initial_pose,
                                         Eigen::Vector3d& initial_velocity,
                                         Eigen::Vector3d& initial_accelerate);
        void setFinalBoundaryCondition(Eigen::Vector3d& final_pose,
                                       Eigen::Vector3d& final_velocity,
                                       Eigen::Vector3d& final_accelerate);
        void setVelocityConstraints(const Eigen::Vector3d& velocity_constraints);
        void setOptimizer();
    private:
        const int order_;
        const int segment_;
        const int derivative_;
        int dimension_;
        std::vector<double> time_;
        std::shared_ptr<nlopt::opt> optimizer;


        std::vector<double> lower_bound_, upper_bound_;

        Eigen::Vector3d initial_pose_;
        Eigen::Vector3d initial_velocity_;
        Eigen::Vector3d initial_accelerate_;
        Eigen::Vector3d final_pose_;
        Eigen::Vector3d final_velocity_;
        Eigen::Vector3d final_accelerate_;

        Eigen::Vector3d velocity_constraints_;
    };

} // namespace optim
#endif //SRC_POLYNOMIALOPTIM_H
