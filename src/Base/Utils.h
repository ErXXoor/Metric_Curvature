//
// Created by hongbo on 04/10/23.
//
#include <Eigen/Core>

#ifndef METRICCURVATURE_UTILS_H
#define METRICCURVATURE_UTILS_H
namespace IGLUtils{
    Eigen::Matrix3d LogEU_tensor_interpolation(std::vector<Eigen::Matrix3d> etens)
    {
        Eigen::Matrix3d linear_ten = Eigen::Matrix3d::Zero();
        for (auto i = 0; i < etens.size(); i++)
        {
            Eigen::Matrix3d log_t =etens[i].array().log() * 1.0 / etens.size() ;
            linear_ten += log_t;
        }
        Eigen::Matrix3d ten_out = linear_ten.array().exp();
        return ten_out;
    }
    Eigen::Matrix3d linear_average(std::vector<Eigen::Matrix3d> etens,double smooth_coeff = 1.0f)
    {
        Eigen::Matrix3d linear_ten = Eigen::Matrix3d::Zero();
        for(const auto & eten : etens)
        {
            linear_ten+=eten;
        }
        linear_ten/=etens.size();
        Eigen::Matrix3d result = smooth_coeff*linear_ten+(1-smooth_coeff)*etens[0];
        return result;
    }

}


#endif //METRICCURVATURE_UTILS_H
