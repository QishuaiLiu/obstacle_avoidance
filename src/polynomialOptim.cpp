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
        smooth_objective_matrix_.resize(segment);
        setQuadraticCoeff();
        computeObjMatrix();
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
        std::vector<double> smooth_grad(x.size());

        double res = 0;

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < segment_; ++j) {
                std::vector<double> segment_parameter = std::vector<double>(x.begin() + i * dimension_ + j * (order_ + 1),
                                                                            x.begin() + i * dimension_ + (j + 1) * (order_ + 1));
                Eigen::VectorXd param(order_ + 1);
                param = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(segment_parameter.data(), segment_parameter.size());
                res += param.transpose() * smooth_objective_matrix_[i] * param;
                Eigen::VectorXd segment_grad = smooth_objective_matrix_[i] * param;
                for (int k = 0; k <= order_; ++k) {
                    grad[i * dimension_ + j * segment_ + k] += segment_grad(k);
                }
            }
        }
        return res;

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

    Eigen::Vector3d polynomialOptim::evaluate(double t, const std::vector<double> &x, int derivative) {
        Eigen::Vector3d res;
        int segment_index = 0;
        while (t > time_[segment_index]) {
            segment_index++;
        }
        segment_index--;
        for (int i = order_; i >= derivative; --i) {
            int index = segment_index * (order_ + 1);
            res(0) = res(0) * t + coeff_(derivative, i) * x[0 * dimension_ + index + i];
            res(1) = res(1) * t + coeff_(derivative, i) * x[1 * dimension_ + index + i];
            res(2) = res(2) * t + coeff_(derivative, i) * x[2 * dimension_ + index + i];
        }

        return res;
    }

    void polynomialOptim::computeObjMatrix() {

        int degree = order_ + 1;
        int derive = derivative_;
        ///////////////////////////////////////////////////
        auto time_exp = [this]( int exp) {
            std::vector<std::vector<double>> time_segment_exp(time_.size());

            for (int i = 0; i < time_.size(); ++i) {
                time_segment_exp[i] = std::vector<double>(exp, 1);
                for (int j = 1; j < exp; ++j) {
                    time_segment_exp[i][j] = time_segment_exp[i][j - 1] * time_[i];
                }
            }
            return time_segment_exp;
        };
        std::vector<std::vector<double>> time_point_exp = time_exp((order_ - derive) * 2 - 1);
        //////////////////////////////////////////////////

        for (int i = 0; i < segment_; ++i) {
            smooth_objective_matrix_[i] = Eigen::MatrixXd(order_ + 1, order_ + 1);
            smooth_objective_matrix_[i].setZero();
            for (int col = 0; col < degree - derive; col++) {
                for (int row = 0; row < degree - derive; row++) {
                    int exp = (order_ - derive) * 2 + 1 - row - col;
                    smooth_objective_matrix_[i](order_ - row, order_ - col) =
                            quadratic_coefficients_[derive - 1](order_ - row, order_ - col)
                    * (time_point_exp[i + 1][exp] - time_point_exp[i][exp]);
                }
            }
        }

        for (int i = 0; i < segment_; ++i) {
            std::cout << smooth_objective_matrix_[i] << std::endl;
            std::cout << "=========" << std::endl;
        }
    }

}
