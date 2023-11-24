//
// Created by hongbo on 08/10/23.
//
#include "Base/Smoother.h"

namespace IGLUtils{
        std::shared_ptr<Eigen::MatrixXd> Smoother::SmoothVec(unsigned iter,
                                                      const std::vector<std::set<int>> &adjacency_list,
                                                      std::shared_ptr<Eigen::MatrixXd> ori_mat) {
            std::shared_ptr<Eigen::MatrixXd> result = ori_mat;
            for (auto i = 0; i < iter; i++) {
                std::vector<Eigen::VectorXd> vec_smt_list;

                for (const auto &v_set: adjacency_list) {
                    std::vector<Eigen::VectorXd> nghb_tensor;
                    for (auto id: v_set) {
                        nghb_tensor.emplace_back(result->row(id));
                    }
                    auto vec_smt = Smoother::linear_average(nghb_tensor);
                    vec_smt.normalize();
                    vec_smt_list.emplace_back(vec_smt);
                }

                result = std::make_shared<Eigen::MatrixXd>();
                result->resize(ori_mat->rows(), ori_mat->cols());
                for(auto v_i=0;v_i<vec_smt_list.size();v_i++)
                {
                    result->row(v_i) = vec_smt_list[v_i];
                }
            }
            return result;
        }

        std::shared_ptr<Eigen::VectorXd> Smoother::SmoothScalar(unsigned int iter,
                                                                const std::vector<std::set<int>> &adjacencyList,
                                                                std::shared_ptr<Eigen::VectorXd> ori_vec) {
            std::shared_ptr<Eigen::VectorXd> result = ori_vec;
            for (auto i = 0; i < iter; i++) {
                std::vector<Eigen::VectorXd> vec_smt_list;

                for (const auto &v_set: adjacencyList) {
                    std::vector<Eigen::VectorXd> nghb_tensor;
                    for (auto id: v_set) {
                        nghb_tensor.emplace_back(result->segment(id,1));
                    }
                    auto vec_smt = Smoother::linear_average(nghb_tensor);
                    vec_smt_list.emplace_back(vec_smt);
                }

                result = std::make_shared<Eigen::VectorXd>();
                result->resize(ori_vec->rows(), ori_vec->cols());
                auto aaa = vec_smt_list.size();
                for(auto v_i=0;v_i<vec_smt_list.size();v_i++)
                {
                    result->row(v_i) = vec_smt_list[v_i];
                }
            }
            return result;
        }

        Eigen::Matrix3d Smoother::LogEU_tensor_interpolation(std::vector<Eigen::Matrix3d> etens) {
            Eigen::Matrix3d linear_ten = Eigen::Matrix3d::Zero();
            for (auto i = 0; i < etens.size(); i++) {
                Eigen::Matrix3d log_t = etens[i].array().log() * 1.0 / etens.size();
                linear_ten += log_t;
            }
            Eigen::Matrix3d ten_out = linear_ten.array().exp();
            return ten_out;
        }

        Eigen::Matrix3d Smoother::linear_average(std::vector<Eigen::Matrix3d> etens, double smooth_coeff) {
            Eigen::Matrix3d linear_ten = Eigen::Matrix3d::Zero();
            for (const auto &eten: etens) {
                linear_ten += eten;
            }
            linear_ten /= etens.size();
            Eigen::Matrix3d result = smooth_coeff * linear_ten + (1 - smooth_coeff) * etens[0];
            return result;
        }

        Eigen::VectorXd Smoother::linear_average(std::vector<Eigen::VectorXd> vec, double smooth_coeff) {
            Eigen::VectorXd linear_vec = Eigen::VectorXd::Zero(vec[0].size());
            for (const auto &v: vec) {
                linear_vec += v;
            }
            linear_vec /= double(vec.size());
            Eigen::VectorXd result = smooth_coeff * linear_vec + (1 - smooth_coeff) * vec[0];
            return result;
        }
}