//
// Created by qishuai on 2021/7/21.
//

#ifndef CAADI_TEST_CASADIOPTIM_H
#define CAADI_TEST_CASADIOPTIM_H

#include <utility>
#include <vector>
#include <Eigen/Eigen>
#include <casadi/casadi.hpp>
#include "obstacleMap.h"

namespace casadiOptim {
    template <typename T>
    T powInt(T x, unsigned int n) {
        if (n == 0) return T{1};

        auto y = T{1};
        while (n > 1) {
            if (n % 2 == 1) y *= x;
            x *= x;
            n /= 2;
        }
        return x * y;
    }

    template <typename T>
    T fact(T num) {
        if (num == 1 || num == 0) {
            return 1;
        } else {
            return num * fact(num - 1);
        }
    }


    class casadiOptim {
    public:
        casadiOptim(int segment, int order, int derivative, std::vector<double> times);
        void setInitialBoundaryCondition(Eigen::Vector3d& initial_pose,
                                         Eigen::Vector3d& initial_velocity,
                                         Eigen::Vector3d& initial_accelerate);
        void setFinalBoundaryCondition(Eigen::Vector3d& final_pose,
                                       Eigen::Vector3d& final_velocity,
                                       Eigen::Vector3d& final_accelerate);
        void setVelocityConstraints(const Eigen::Vector3d& velocity_constraints);
        void setObjectiveMatrix();
        void setConstraintMatrix();
        int getNumOfVariables() { return dimension_; }
        casadi::DMDict solve();

    private:
        casadi::DM setConstraintMatrixOneDim();
        void setBound();
        std::vector<double> setBound(int dim);

    private:
        const int order_;
        const int segment_;
        const int derivative_;
        int dimension_;

        casadi::MX x;  // variables
        casadi::DM Q;  // Hessian matrix
        casadi::DM A;  // constraints matrix

        std::vector<double> lower_bound_, upper_bound_;

        Eigen::Vector3d initial_pose_;
        Eigen::Vector3d initial_velocity_;
        Eigen::Vector3d initial_accelerate_;
        Eigen::Vector3d final_pose_;
        Eigen::Vector3d final_velocity_;
        Eigen::Vector3d final_accelerate_;

        std::vector<double> time_;
        Eigen::Vector3d velocity_constraints_;

    };
} // namespace casadiOptim



#endif //CAADI_TEST_CASADIOPTIM_H
