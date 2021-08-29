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
        // EIGEN_MAKE_ALIGNED_OPERATOR_NEW
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
        static double costWarp(const std::vector<double>&x, std::vector<double>& grad, void *data) {
            return reinterpret_cast<polynomialOptim*>(data)->smooth_objective(x, grad);
        };
        void equalConstraintWarp(unsigned m, double *result, unsigned n, const double* x, double* grad, void* data) {
            reinterpret_cast<polynomialOptim*>(data)->totalEqualConstraint(x, result, grad);
        };
        double smooth_objective(const std::vector<double>&x, std::vector<double>& grad);

        void totalEqualConstraint(const double* x, double* result, double* grad);

        void setQuadraticCoeff();

        Eigen::Vector3d evaluate(double t, const std::vector<double> &x, int derivative);

        void getCoeffWithTime(Eigen::VectorXd& coeff_with_time, int derivative, double t);

        void computeObjMatrix();


    private:
        const int order_;
        const int segment_;
        const int derivative_;
        int dimension_;
        std::vector<double> time_;
        std::shared_ptr<nlopt::opt> optimizer_;


        std::vector<double> lower_bound_, upper_bound_;

        Eigen::Vector3d initial_pose_;
        Eigen::Vector3d initial_velocity_;
        Eigen::Vector3d initial_accelerate_;
        Eigen::Vector3d final_pose_;
        Eigen::Vector3d final_velocity_;
        Eigen::Vector3d final_accelerate_;

        Eigen::Vector3d velocity_constraints_;

        Eigen::MatrixXd coeff_;
        std::vector<Eigen::MatrixXd> quadratic_coefficients_;
        std::vector<Eigen::MatrixXd> smooth_objective_matrix_;

    };

} // namespace optim
#endif //SRC_POLYNOMIALOPTIM_H
