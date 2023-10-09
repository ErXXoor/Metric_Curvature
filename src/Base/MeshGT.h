//
// Created by hongbo on 13/09/23.
//

#ifndef METRIC_CURVATURE_MESHIO_H
#define METRIC_CURVATURE_MESHIO_H

#include <string>
#include <Eigen/Core>
#include <memory>
#include "Base/Metric.h"

namespace IGLUtils {
    class MeshGT {
    public:
        MeshGT();

        ~MeshGT() = default;

        bool LoadMesh(const std::string &filename);

        bool CalculateCurvature();

        void ProcessMetric(unsigned int smooth_ring, unsigned int smooth_iter);

        void ViewCurvature();

        void SaveCurvature(const std::string &filepath);

        void SaveMetric(const std::string &filepath);

        void SaveMeshInfo(const std::string &filepath);

    private:
        std::string m_filename;
        std::shared_ptr<Eigen::MatrixXd> m_v;
        std::shared_ptr<Eigen::MatrixXi> m_f;
        std::shared_ptr<Eigen::MatrixXd> m_n;
        std::shared_ptr<Eigen::MatrixXd> m_min_pd;
        std::shared_ptr<Eigen::MatrixXd> m_max_pd;
        std::shared_ptr<Eigen::VectorXd> m_max_pv;
        std::shared_ptr<Eigen::VectorXd> m_min_pv;
        std::shared_ptr<Metric> m_metric;
    };
}
#endif //METRIC_CURVATURE_MESHIO_H
