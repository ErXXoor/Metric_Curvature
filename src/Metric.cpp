//
// Created by hongbo on 04/10/23.
//
#include"Base/Metric.h"
#include "Base/Utils.h"
#include<algorithm>
namespace IGLUtils{
    void Metric::ComposeMetric(std::shared_ptr<Eigen::MatrixXd> m_n,
                                std::shared_ptr<Eigen::MatrixXd> m_min_pd,
                                std::shared_ptr<Eigen::MatrixXd> m_max_pd,
                                std::shared_ptr<Eigen::VectorXd> m_min_pv,
                                std::shared_ptr<Eigen::VectorXd> m_max_pv)
                                {
        for(auto i=0;i<m_n->rows();i++) {
            auto k_max = m_max_pv->coeff(i);
            auto k_min = m_min_pv->coeff(i);
            k_max = std::clamp(abs(k_max),1.0,400.0);
            k_min = std::clamp(abs(k_min),1.0,400.0);

            auto k_max_dir = m_max_pd->row(i);
            auto k_min_dir = m_min_pd->row(i);
            auto n = m_n->row(i);
            Eigen::Matrix3d r = Eigen::Matrix3d::Zero();
            r.block(0,0,1,3) = k_min_dir;
            r.block(1,0,1,3) = k_max_dir;
            r.block(2,0,1,3) = n;

            Eigen::Matrix3d s = Eigen::Matrix3d::Zero();
            s(0,0) = k_min;
            s(1,1) = k_max;
            s(2,2) = 0.0;
            auto tensor = r.transpose()*s*r;
            m_metric.emplace_back(tensor);
        }
    }

    void Metric::SmoothMetric(unsigned int iter,
                              const std::vector<std::set<int>> &adjacency_list)
    {
        for(auto i=0;i<iter;i++)
        {
            std::vector<Eigen::Matrix3d> metric_smt_list;

            for(const auto& v_set:adjacency_list)
            {
                std::vector<Eigen::Matrix3d> nghb_tensor;
                for(auto id:v_set)
                {
                    nghb_tensor.emplace_back(m_metric[id]);
                }
                auto metric_smt = linear_average(nghb_tensor);
                metric_smt_list.emplace_back(metric_smt);
            }
            m_metric = metric_smt_list;
        }
    }

}