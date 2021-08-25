//
// Created by qishuai on 2021/8/25.
//

#include "polynomialOptim.h"
#include <iostream>

namespace optim {
    polynomialOptim::polynomialOptim(int segment, int order, int derivative, std::vector<double> times):
    order_(order), segment_(segment), derivative_(derivative), time_(times) {
        dimension_ = (order + 1) * segment;
        coeff_ = Eigen::MatrixXd(order_ + 1, order_ + 1);
        quadratic_coefficients_ = std::vector<Eigen::MatrixXd>(derivative_);
        setQuadraticCoeff();
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

    void polynomialOptim::setOptimizer() {
        optimizer = std::make_shared<nlopt::opt>(nlopt::LD_MMA, 3 * dimension_);
        std::vector<double> lb(3 * dimension_, -1e8);
        optimizer->set_lower_bounds(lb);
        optimizer->set_min_objective(costWarp, this);
        optimizer->set_xtol_rel(1e-4);

    }

    double polynomialOptim::smooth_objective(const std::vector<double>&x, std::vector<double>& grad) {

    }

    void polynomialOptim::setQuadraticCoeff() {
        coeff_.setZero();
        coeff_.row(0).setOnes();
        int order = order_;
        for (int i = 1; i <= order_; ++i) {
            for (int j = order_ - order; j <= order_; ++j) {
                coeff_(i, j) = (order - order_ + j) * coeff_(i - 1, j);
            }
            order--;
        }

        int degree = order_ + 1;
        for (int derive = 1; derive <= derivative_; ++derive) {
            quadratic_coefficients_[derive - 1] = Eigen::MatrixXd(order_ + 1, order_ + 1);
            quadratic_coefficients_[derive - 1].setZero();
            for (int col = 0; col < degree - derive; col++) {
                for (int row = 0; row < degree - derive; row++) {
                    double exp = (order_ - derive) * 2 + 1 - row - col;
                    quadratic_coefficients_[derive - 1](order_ - row, order_ - col) =
                            coeff_(derive, order_ - row) *
                            coeff_(derive, order_ - col) * 2 / exp;
                }
            }
        }
    }


}
