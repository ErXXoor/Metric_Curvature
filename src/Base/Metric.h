//
// Created by hongbo on 04/10/23.
//

#ifndef METRICCURVATURE_METRIC_H
#define METRICCURVATURE_METRIC_H

#include <Eigen/Core>
#include <memory>
#include<set>

namespace IGLUtils {
    class Metric {
    public:
        Metric() = default;

        ~Metric() = default;

        void ComposeMetric(std::shared_ptr<Eigen::MatrixXd> m_n,
                           std::shared_ptr<Eigen::MatrixXd> m_min_pd,
                           std::shared_ptr<Eigen::MatrixXd> m_max_pd,
                           std::shared_ptr<Eigen::VectorXd> m_min_pv,
                           std::shared_ptr<Eigen::VectorXd> m_max_pv);

        void SaveMetric(const std::string &filepath);
        void SaveSR(const std::string &s_filepath, const std::string &r_filepath);

    private:
        std::vector<Eigen::Matrix3d> m_metric;
        std::vector<Eigen::Matrix3d> m_s;
        std::vector<Eigen::Matrix3d> m_r;
    };

}
#endif //METRICCURVATURE_METRIC_H
