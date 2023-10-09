//
// Created by hongbo on 08/10/23.
//

#ifndef METRICCURVATURE_SMOOTHER_H
#define METRICCURVATURE_SMOOTHER_H
#include <Eigen/Core>
#include <set>
#include <memory>
namespace IGLUtils {
    class Smoother {
    public:
        static std::shared_ptr<Eigen::MatrixXd> SmoothVec(unsigned iter,
                                                      const std::vector<std::set<int>> &adjacencyList,
                                                      std::shared_ptr<Eigen::MatrixXd> ori_mat);
        static std::shared_ptr<Eigen::VectorXd> SmoothScalar(unsigned iter,
                                                         const std::vector<std::set<int>> &adjacencyList,
                                                         std::shared_ptr<Eigen::VectorXd> ori_vec);

        static Eigen::Matrix3d LogEU_tensor_interpolation(std::vector<Eigen::Matrix3d> etens);

        static Eigen::Matrix3d linear_average(std::vector<Eigen::Matrix3d> etens, double smooth_coeff = 1.0f);
        static Eigen::VectorXd linear_average(std::vector<Eigen::VectorXd> vec, double smooth_coeff = 1.0f);


    };
}

#endif //METRICCURVATURE_SMOOTHER_H
