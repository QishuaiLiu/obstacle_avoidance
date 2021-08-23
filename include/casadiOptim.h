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
        void setMap(std::shared_ptr<obstacleMap> obstacle_map) { obstacle_map_ = std::move(obstacle_map); }
        std::vector<Eigen::Vector3f> getSegmentCost(const casadi::DM& optimal_x);

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
        std::shared_ptr<obstacleMap> obstacle_map_;

    };
} // namespace casadiOptim

class callbackJacobian : public casadi::Callback {
private:
    std::shared_ptr<casadiOptim::obstacleMap> cmap_;
    casadi::DM Q_;
public:
    // Constructor.
    callbackJacobian(const std::string& name, std::shared_ptr<casadiOptim::obstacleMap>  map, const casadi::DM& Q,
                     const casadi::Dict &opts = casadi::Dict()): cmap_(std::move(map)), Q_(Q) {
        construct(name, opts);
    }

    // Destructor.
    ~callbackJacobian() override { }

    // Initialize callback.
    void init() override {
        std::cout << "initializing callback jacobian" << std::endl;
    }

    // Number of inputs.
    casadi_int get_n_in() override { return 2; }

    // Number of outputs.
    casadi_int get_n_out() override { return 1; }

    // Set sparsity of input.
    casadi::Sparsity get_sparsity_in(casadi_int i) override {
        if (i == 0)
        return casadi::Sparsity::dense(84, 1);
    }

    // Set sparsity of output.
    casadi::Sparsity get_sparsity_out(casadi_int i) override {
        return casadi::Sparsity::dense(1, 84);
    }

    // Evaluate numerically.
    std::vector<casadi::DM> eval(const std::vector<casadi::DM>& arg) const override {
        std::cout << "evaluating jacobian callback" << std::endl;
        // std::vector<Eigen::Vector3f> grad;
        // cmap_->getGrad(grad);
        casadi::DM out(1, 84);
        casadi::DM smooth_jacobian = mtimes(arg[0].T(), Q_);

        // for (int j = 0; j < grad.size(); ++j) {
        //     out(0, 0 * grad.size() + j) = grad[j][0];
        //     out(0, 1 * grad.size() + j) = grad[j][1];
        //     out(0, 2 * grad.size() + j) = grad[j][2];
        // }

        return { smooth_jacobian };
    }
};



class MyCallback : public casadi::Callback {
private:
    std::shared_ptr<casadiOptim::obstacleMap> cmap_;
    std::vector<double> time_;
    int dimension_;
    casadi::DM Q_;
public:
    MyCallback(const std::string& name, std::shared_ptr<casadiOptim::obstacleMap>  map, std::vector<double> time,
               int dimension, const casadi::DM& Q, const casadi::Dict& opts = casadi::Dict()):
            cmap_(std::move(map)), time_(std::move(time)), dimension_(dimension), Q_(Q){
        construct(name, opts);
    }
    ~MyCallback() override {}

    casadi_int get_n_in() override {return 1; }
    casadi_int get_n_out () override {return 1; }

    casadi::Sparsity get_sparsity_in(casadi_int i) override {
        return casadi::Sparsity::dense(84, 1);
    }

    casadi::Sparsity get_sparsity_out(casadi_int i) override {
        return casadi::Sparsity::dense(1, 1);
    }

    void init() override {
        std::cout << "initializing_object! " << std::endl;
    }

    std::vector<casadi::DM> eval(const std::vector<casadi::DM>& arg) const override {
        int order = 7;
        int segment = 4;
        std::vector<std::vector<float>> coef(segment, std::vector<float>(order));
        float truncate_distance = 1.97; // pow((0.5 / 0.3), 2);
        std::vector<std::vector<std::vector<float>>> three_coef(3, coef);

        // for (int i = 0; i < segment; ++i) {
        //     for (int j = 0; j < order; ++j) {
        //         coef[i][j] = arg.at(0)(21 + j, 0).scalar();
        //     }
        // }
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < segment; ++j) {
                for (int k = 0; k < order; ++k) {
                    three_coef[i][j][k] = arg.at(0)(i * dimension_ + j * (order) + k).scalar();
                }
            }
        }

        casadi::DM f = 0;
        float eta = 12.2;
        int num = 0;
        std::vector<Eigen::Vector3f> grad;


        int time_num = 8;

        for (int i = 0; i < segment; ++i) {
            for (int k = 0; k < time_num; ++k) {
                double delta_time = (time_[i + 1] - time_[i]) / time_num;
                std::vector<float> pos(3), vel(3);
                for (int dim = 0; dim < 3; ++dim) {
                    const std::vector<float>& coef_one_dim = three_coef[dim][i];
                    for (int j = 0; j <= order; ++j) {
                        pos[dim] += (coef_one_dim[j]) * casadiOptim::powInt(time_[i] + delta_time * k, j);
                        // if (j > 0) {
                        //     vel[dim] += (coef_one_dim[j]) * j * casadiOptim::powInt(time_[i] + delta_time * k, j - 1);
                        // }
                    }
                }
                // printf("pos x: %.2f, y: .2f, z: .2f", pos[0], pos[1], pos[2]);
                circular_map::MapPoint3f point(pos[0], pos[1], pos[2]);
                float distance;
                Eigen::Vector3f point_grad;
                cmap_->getDistanceAndGrad(point, distance, point_grad);
                if (distance > truncate_distance) {
                    grad.emplace_back(Eigen::Vector3f(0, 0, 0));
                } else {
                    float temp_f = eta * pow((distance - truncate_distance), 2);
                    f += temp_f;
                    grad.emplace_back(point_grad);
                    // printf("f is: %.2f \n", temp_f);
                }
            }
        }
        cmap_->setGrad(grad);
        return {mtimes(mtimes(arg[0].T(), Q_), arg[0])};
    }

    // bool has_reverse(casadi_int nadj) const override {return true;}

    bool has_jacobian() const override {return true; }

    casadi::Function get_jacobian(const std::string&name, const std::vector<std::string>& inames, const std::vector<std::string> &onames,
                         const casadi::Dict &opts) const override {
        callbackJacobian cbJacobian(name, cmap_, Q_);
        casadi::Function jacobianFunction = cbJacobian;
        return jacobianFunction;
    }
};


#endif //CAADI_TEST_CASADIOPTIM_H
