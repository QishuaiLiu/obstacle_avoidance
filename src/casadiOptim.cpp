//
// Created by qishuai on 2021/7/21.
//

#include "casadiOptim.h"
namespace casadiOptim {

    casadiOptim::casadiOptim(int segment, int order, int derivative, std::vector<double> times) :
    segment_(segment), order_(order), derivative_(derivative), time_(std::move(times)) {
        dimension_ = (order + 1) * segment;
        int three_dimension = dimension_ * 3;
        x = casadi::MX::sym("x", three_dimension, 1);
        Q = casadi::DM(three_dimension, three_dimension);
        A = casadi::DM((3 * segment_ + 6) * 3, three_dimension);
        lower_bound_ = std::vector<double>((3 * segment_ + 6) * 3);
        upper_bound_ = std::vector<double>((3 * segment_ + 6) * 3);
    }

    void casadiOptim::setInitialBoundaryCondition(Eigen::Vector3d &initial_pose,
                                                  Eigen::Vector3d &initial_velocity,
                                                  Eigen::Vector3d &initial_accelerate) {
        initial_pose_ = initial_pose;
        initial_velocity_ = initial_velocity;
        initial_accelerate_ = initial_accelerate;
    }

    void casadiOptim::setFinalBoundaryCondition(Eigen::Vector3d &final_pose,
                                                Eigen::Vector3d &final_velocity,
                                                Eigen::Vector3d &final_accelerate) {
        final_pose_ = final_pose;
        final_velocity_ = final_velocity;
        final_accelerate_ = final_accelerate;
    }

    void casadiOptim::setObjectiveMatrix() {
        for (int dimension = 0; dimension < 3; ++dimension) {
            for (int i = 0; i < segment_; ++i) {
                for (int r = derivative_; r <= order_; ++r) {
                    for (int c = derivative_; c <= order_; ++c) {
                        double time_intergrate = powInt(time_[i + 1], r + c - 7) - powInt(time_[i], r + c - 7);
                        double numerator = fact(r) * fact(c) * time_intergrate;
                        double denominator = fact(r - 4) * fact(c - 4) * (r + c - 7);
                        Q(i * (order_ + 1) + r + dimension * dimension_,
                          i * (order_ + 1) + c + dimension * dimension_) = numerator / denominator;
                    }
                }
            }
        }

    }

    casadi::DM casadiOptim::setConstraintMatrixOneDim() {
        casadi::DM A = casadi::DM(3 * segment_ + 6, dimension_);
        for (int i = 0; i < segment_ - 1; ++i) {
            for (int j = 0; j <= order_; ++j) {
                A(i * 3, i * (order_ + 1) + j) = powInt(time_[i + 1], j);
                if (j >= 1) {
                    A(i * 3 + 1, i * (order_ + 1) + j) = j * powInt(time_[i + 1], j - 1);
                }
                if (j >= 2) {
                    A(i * 3 + 2, i * (order_ + 1) + j) = static_cast<double>(fact(j)) / fact(j - 2) * powInt(time_[i + 1], j - 2);
                }
            }

            for (int j = 0; j <= order_; ++j) {
                A(i * 3, (i + 1) * (order_ + 1) + j) = -powInt(time_[i + 1], j);
                if (j >= 1) {
                    A(i * 3 + 1, (i + 1) * (order_ + 1) + j) = -j * powInt(time_[i + 1], j - 1);
                }
                if (j >= 2) {
                    A(i * 3 + 2, (i + 1) * (order_ + 1) + j) = -static_cast<double>(fact(j)) / fact(j - 2) * powInt(time_[i + 1], j - 2);
                }
            }
        }

        for (int i = 0; i <= order_; ++i) {
            A(3 * (segment_ - 1), i) = powInt(time_[0], i);
            A(3 * segment_, (segment_ - 1) * (order_ + 1) + i) = powInt(time_.back(), i);

            if (i >= 1) {
                A(3 * (segment_ - 1) + 1, i) = i * powInt(time_[0], i - 1);
                A(3 * segment_ + 1, (segment_ - 1) * (order_ + 1) + i) = i * powInt(time_.back(), i - 1);
            }

            if (i >= 2) {
                A(3 * (segment_ - 1) + 2, i) = static_cast<double>(fact(i)) / fact(i - 2) * powInt(time_[0], i - 2);
                A(3 * segment_ + 2, (segment_ - 1) * (order_ + 1) + i) = static_cast<double>(fact(i)) / fact(i - 2) * powInt(time_.back(), i - 2);
            }
        }

        for (int i = 0; i < segment_ - 1; ++i) {
            for (int j = 1; j <= order_; ++j) {
                A(3 * segment_ + 3 + i, (i + 1) * (order_ + 1) + j) = j * powInt(time_[i + 1], j - 1);
            }
        }

        return A;
    }

    void casadiOptim::setConstraintMatrix() {
        int row_num = 3 * segment_ + 6;
        auto temp_A = setConstraintMatrixOneDim();
        for (int dim = 0; dim < 3; ++dim) {
            A(casadi::Slice(dim * row_num, (dim + 1) * row_num), casadi::Slice(dim * dimension_, (dim + 1) * dimension_))
             = temp_A;
        }
    }

    std::vector<double> casadiOptim::setBound(int dim) {
        std::vector<double> equal_constraints(3 * segment_ + 6);
        for (int i = 0; i < 3 * (segment_ - 1); ++i) {
            equal_constraints[i] = 0.0;
        }

        equal_constraints[3 * (segment_ - 1)] = initial_pose_(dim);
        equal_constraints[3 * (segment_ - 1) + 1] = initial_velocity_(dim);
        equal_constraints[3 * (segment_ - 1) + 2] = initial_accelerate_(dim);
        equal_constraints[3 * segment_] = final_pose_(dim);
        equal_constraints[3 * segment_ + 1] = final_velocity_(dim);
        equal_constraints[3 * segment_ + 2] = final_accelerate_(dim);

        return equal_constraints;

        // lower_bound_ = equal_constraints;
        // upper_bound_ = equal_constraints;
        // for (int i = 0; i < (segment_ - 1); ++i) {   // intermediate segment velocity constraints
        //     lower_bound_[lower_bound_.size() - 1 - i] = -velocity_constraints_(dim);
        //     upper_bound_[upper_bound_.size() - 1 - i] = velocity_constraints_(dim);
        // }
    }

    void casadiOptim::setBound() {
        int constraint_num = (3 * segment_ + 6);
        for (int dim = 0; dim < 3; ++dim) {
            auto equal_constraint = setBound(dim);
            int start_index = dim * constraint_num;
            for (int i = 0; i < constraint_num; ++i) {
                lower_bound_[start_index + i] = equal_constraint[i];
                upper_bound_[start_index + i] = equal_constraint[i];
            }

            for (int i = 0; i < (segment_ - 1); ++i) {
                lower_bound_[(dim + 1) * constraint_num - 1 - i] = -velocity_constraints_(dim);
                upper_bound_[(dim + 1) * constraint_num - 1 - i] = velocity_constraints_(dim);
            }
        }
    }

    void casadiOptim::setVelocityConstraints(const Eigen::Vector3d &velocity_constraints) {
        velocity_constraints_ = velocity_constraints;
    }

    casadi::DMDict casadiOptim::solve() {
        setObjectiveMatrix();
        setConstraintMatrix();
        std::vector<double> initial_value(dimension_ * 3);

        casadi::Dict callopt;
        callopt["enable_fd"] = true;
        // MyCallback callback("f", obstacle_map_, time_, dimension_, Q, callopt);

        // casadi::MX smooth_cost = mtimes(mtimes(x.T(), Q), x) + callback(x).at(0);
        // casadi::MX smooth_cost = callback(x).at(0);

        casadi::MX f = 4;
        casadi::MX g = mtimes(A, x);
        casadi::MXDict nlp;
        nlp["x"] = x;
        nlp["f"] = f;
        nlp["g"] = g;
        casadi::Dict opts_ipopt;
        opts_ipopt["ipopt"] = casadi::Dict{{"max_iter", 100}, {"hessian_approximation", "limited-memory"}
                                           };
        // opts_ipopt["ipopt"] = casadi::Dict{{"max_iter", 20},
        //                                    {"print_level", 0},
        //                                    {"acceptable_tol", 1e-8},
        //                                    {"acceptable_obj_change_tol", 1e-6},
        //                                    {"hessian_approximation", "limited-memory"}};

        // std::vector<casadi::DMDict> optimal_solution(3);

        ////////////////solver
        setBound();
        for (int i = 0; i < dimension_ * 3; ++i) {
            initial_value[i] = 1.0;
        }

        std::cout << "initial value size: " << initial_value.size() <<
        "lower bound size: " << lower_bound_.size() <<
        "upper bound size: " << upper_bound_.size() << std::endl;


        casadi::Function solver = nlpsol("F","ipopt",nlp, opts_ipopt);
        // casadi::Function solver = nlpsol("F","ipopt",nlp);
        casadi::DMDict optimal_solution = solver(casadi::DMDict
                   {{"x0",casadi::DM(initial_value)},{"ubg",upper_bound_},
                   {"lbg",lower_bound_}});
        return optimal_solution;
    }


} // namespace casadiOptim