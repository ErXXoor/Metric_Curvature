//
// Created by hongbo on 04/10/23.
//
#include "Base/Metric.h"
#include <fstream>
#include <algorithm>

namespace IGLUtils {

    void Metric::ComposeMetric(std::shared_ptr<Eigen::MatrixXd> m_n,
                               std::shared_ptr<Eigen::MatrixXd> m_min_pd,
                               std::shared_ptr<Eigen::MatrixXd> m_max_pd,
                               std::shared_ptr<Eigen::VectorXd> m_min_pv,
                               std::shared_ptr<Eigen::VectorXd> m_max_pv) {
        for (auto i = 0; i < m_n->rows(); i++) {
            auto k_max = m_max_pv->coeff(i);
            auto k_min = m_min_pv->coeff(i);
            k_max = std::clamp(abs(k_max), 1.0, 400.0);
            k_min = std::clamp(abs(k_min), 1.0, 400.0);

            auto k_max_dir = m_max_pd->row(i);
            auto k_min_dir = m_min_pd->row(i);
            auto n = m_n->row(i);
            Eigen::Matrix3d r = Eigen::Matrix3d::Zero();
            r.block(0, 0, 1, 3) = k_min_dir;
            r.block(1, 0, 1, 3) = k_max_dir;
            r.block(2, 0, 1, 3) = n;

            Eigen::Matrix3d s = Eigen::Matrix3d::Zero();
            s(0, 0) = k_max/k_min;
            s(1, 1) = 1;
            s(2, 2) = 0.0;

            m_s.emplace_back(s);
            m_r.emplace_back(r);
            auto tensor = s * r;
            m_metric.emplace_back(tensor);
        }
    }

    void Metric::SaveMetric(const std::string &filepath) {
        std::ofstream out(filepath);
        for (auto &i: m_metric) {
            out << i(0, 0) << "," << i(0, 1) << "," << i(0, 2) << ",";
            out << i(1, 0) << "," << i(1, 1) << "," << i(1, 2) << ",";
            out << i(2, 0) << "," << i(2, 1) << "," << i(2, 2) << std::endl;
        }
        out.close();
    }

    void Metric::SaveSR(const std::string& s_filepath, const std::string& r_filepath){
        std::ofstream out(s_filepath);
        for(auto &i:m_s){
            out << i(0, 0) << "," << i(0, 1) << "," << i(0, 2) << ",";
            out << i(1, 0) << "," << i(1, 1) << "," << i(1, 2) << ",";
            out << i(2, 0) << "," << i(2, 1) << "," << i(2, 2) << std::endl;
        }
        out.close();

        std::ofstream out2(r_filepath);
        for(auto &i:m_r){
            out2 << i(0, 0) << "," << i(0, 1) << "," << i(0, 2) << ",";
            out2 << i(1, 0) << "," << i(1, 1) << "," << i(1, 2) << ",";
            out2 << i(2, 0) << "," << i(2, 1) << "," << i(2, 2) << std::endl;
        }
    }

}