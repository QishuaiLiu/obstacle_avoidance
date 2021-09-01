//
// Created by qishuai on 2021/8/25.
//

#include "polynomialOptim.h"
#include <iostream>

namespace optim {

    polynomialOptim::polynomialOptim(int segment, int order, int derivative, std::vector<double> times):
    order_(order), segment_(segment), derivative_(derivative), time_(times) {
        dimension_ = (order + 1) * segment;
        // setQuadraticCoeff();
        computeObjMatrix();
        computeConstraintMatrix();
        setOptimizer();
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
        optimizer_ = std::make_shared<nlopt::opt>(nlopt::LD_SLSQP, 3 * dimension_); //LD_SLSQP
        std::vector<double> lb(3 * dimension_, -1e8);
        std::vector<double> ub(3 * dimension_, 1e8);
        optimizer_->set_lower_bounds(lb);
        optimizer_->set_upper_bounds(ub);
        optimizer_->set_min_objective(costWarp, this);
        optimizer_->set_xtol_rel(1e-4);
        std::vector<double> tol_constraint(5, 1e-4);
        // optimizer_->add_equality_mconstraint(equalConstraintWarp, this, tol_constraint);
    }

    void polynomialOptim::optimize() {

        std::vector<double> x(3 * dimension_);
        for (int i = 0; i < x.size(); ++i) {
            x[i] = 1.0;
        }
        double minf;
        nlopt::result result = optimizer_->optimize(x, minf);
        std::cout << "result is: " << result << " minf: " << minf << std::endl;
    }

    double polynomialOptim::smooth_objective(const std::vector<double>&x, std::vector<double>& grad) {
        std::vector<double> smooth_grad(x.size());

        double res = 0;
        static int num = 0;

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < segment_; ++j) {
                std::vector<double> segment_parameter = std::vector<double>(x.begin() + i * dimension_ + j * (order_ + 1),
                                                                            x.begin() + i * dimension_ + (j + 1) * (order_ + 1));
                Eigen::VectorXd param(order_ + 1);
                param = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(segment_parameter.data(), segment_parameter.size());
                res += param.transpose() * smooth_objective_matrix_[j] * param;
                Eigen::VectorXd segment_grad = smooth_objective_matrix_[j] * param;
                for (int k = 0; k <= order_; ++k) {
                    grad[i * dimension_ + j * (order_ + 1) + k] += segment_grad(k);
                }
            }
        }
        // for (int i = 0; i < grad.size(); ++i) {
        //     if (i % 28 == 0)
        //         std::cout << std::endl;
        //     std::cout << grad[i] << " ";
        // }
        // std::cout << std::endl;
        //
        // std::cout << "time is: " << num++ << " res is: " << res << " size: " <<  grad.size()<< std::endl;
        return res;

    }

    void polynomialOptim::setQuadraticCoeff() {
        coeff_ = Eigen::MatrixXd(order_ + 1, order_ + 1);
        coeff_.setZero();
        coeff_.row(0).setOnes();
        int order = order_;
        for (int i = 1; i <= order_; ++i) {
            for (int j = order_ - order; j <= order_; ++j) {
                coeff_(i, j) = (order - order_ + j) * coeff_(i - 1, j);
            }
            order--;
        }

        quadratic_coefficients_ = std::vector<Eigen::MatrixXd>(derivative_);
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
        setQuadraticCoeff();
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
        time_point_exp_ = time_exp(order_ + 1);
        //////////////////////////////////////////////////
        smooth_objective_matrix_.resize(segment_);
        for (int i = 0; i < segment_; ++i) {
            smooth_objective_matrix_[i] = Eigen::MatrixXd(order_ + 1, order_ + 1);
            smooth_objective_matrix_[i].setZero();
            for (int col = 0; col < degree - derive; col++) {
                for (int row = 0; row < degree - derive; row++) {
                    int exp = (order_ - derive) * 2 + 1 - row - col;
                    smooth_objective_matrix_[i](order_ - row, order_ - col) =
                            quadratic_coefficients_[derive - 1](order_ - row, order_ - col)
                    * (time_point_exp_[i + 1][exp] - time_point_exp_[i][exp]);
                }
            }
        }
    }



    void polynomialOptim::totalEqualConstraint(unsigned m, double *result, unsigned n, const double* x, double* grad) {
        // Eigen::VectorXd coeff_with_time_first, coeff_with_time_last;
        //
        // int dim = 3 * dimension_;
        // for (int i = 0; i < 3; ++i) {
            // getCoeffWithTime(coeff_with_time_first, i, time_[0]);
            // getCoeffWithTime(coeff_with_time_last, i, time_.back());
            // for (int j = 0; j < order_ + 1; ++j) {
            //     grad[i * dim + j] = coeff_with_time_first[j];
            //     result[i] += coeff_with_time_first[j] * x[j] - 0;
            //
            //     grad[(3 + i) * dim  + 2 * dimension_ + (segment_ - 1) * (order_ + 1) + j] =
            //             coeff_with_time_last[j];
            //     result[3 + i] += coeff_with_time_last[j] * x[2 * dimension_ + (segment_ - 1) * (order_ + 1) + j] - 0;
            // }
        //     for (int j = 0; j < 128; ++j) {
        //         grad[j  + i] = 2;
        //         result[j] = 4;
        //     }
        // }

        Eigen::VectorXd vals(dimension_);
        for (int i = 0; i < 3; ++i) {
            for (int k = 0; k < dimension_; ++k) {
                vals[k] = x[i * dimension_ + k];
            }

            for (int j = 0; j < 3 * segment_ + 3; ++j) {
                result[i * dimension_ + j] = constraint_matrix_.row(j) * vals;
            }
        }


        for (int i = 0; i < n; ++i) {
            for (int k = 0; k < dimension_; ++k) {
                vals[k] = x[i * dimension_ + k];
            }
            for (int j = 0; j < m; ++j) {
                grad[j * n + i] = 1;
                result[j] += grad[j * n + i] * 0.1;
            }
        }


    }

    void polynomialOptim::computeConstraintMatrix() {
        constraint_matrix_ = Eigen::MatrixXd(3 * segment_ + 3, dimension_);
        for (int i = 0; i < segment_ - 1; ++i) {
            for (int j = 0; j <= order_; ++j) {
                constraint_matrix_(i * 3, i * (order_ + 1) + j) = time_point_exp_[i + 1][j];
                if (j >= 1) {
                    constraint_matrix_(i * 3 + 1, i * (order_ + 1) + j) = j * time_point_exp_[i + 1][j - 1];
                }
                if (j >= 2) {
                    constraint_matrix_(i * 3 + 2, i * (order_ + 1) + j) = j * (j - 1) * time_point_exp_[i + 1][j - 2];
                }
            }

            for (int j = 0; j <= order_; ++j) {
                constraint_matrix_(i * 3, (i + 1) * (order_ + 1) + j) = -time_point_exp_[i + 1][j];
                if (j >= 1) {
                    constraint_matrix_(i * 3 + 1, (i + 1) * (order_ + 1) + j) = -j * time_point_exp_[i + 1][j - 1];
                }
                if (j >= 2) {
                    constraint_matrix_(i * 3 + 2, (i + 1) * (order_ + 1) + j) = -j * (j - 1) * time_point_exp_[i + 1][j - 2];
                }
            }
        }

        for (int i = 0; i <= order_; ++i) {
            constraint_matrix_(3 * (segment_ - 1), i) = time_point_exp_[0][i];
            constraint_matrix_(3 * segment_, (segment_ - 1) * (order_ + 1) + i) = time_point_exp_[time_.size() - 1][i];

            if (i >= 1) {
                constraint_matrix_(3 * (segment_ - 1) + 1, i) = i * time_point_exp_[0][i - 1];
                constraint_matrix_(3 * segment_ + 1, (segment_ - 1) * (order_ + 1) + i) = i * time_point_exp_[time_.size() - 1][i - 1];
            }
            if (i >= 2) {
                constraint_matrix_(3 * (segment_ - 1) + 2, i) = i * (i - 1) * time_point_exp_[0][i - 2];
                constraint_matrix_(3 * segment_ + 2, (segment_ - 1) * (order_ + 1) + i) =
                        i * (i - 1) * time_point_exp_[time_.size() - 1][i - 2];
            }
        }

    }

    void polynomialOptim::getCoeffWithTime(Eigen::VectorXd& coeff_with_time, int derivative, double t) {
        int size = coeff_.cols();
        coeff_with_time = Eigen::VectorXd(size);
        coeff_with_time.setZero();
        coeff_with_time[derivative] = coeff_(derivative, derivative);

       // int segment = 0;
       // while(t > time_[segment]) {
       //     segment++;
       // }
       // segment--;

        double ti = t;
        for (int i = derivative + 1; i < size; ++i) {
            coeff_with_time[i] = coeff_(derivative, i) * ti;
            ti *= t;
        }

    }

}
