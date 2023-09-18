//
// Created by hongbo on 13/09/23.
//

#ifndef METRIC_CURVATURE_MESHIO_H
#define METRIC_CURVATURE_MESHIO_H
#include <string>
#include <Eigen/Core>
#include <memory>

namespace IGLUtils {
    class MeshGT {
    public:
        MeshGT();

        bool LoadMesh(const std::string &filename);

        void CalculateCurvature();

        void ComposeSRMatVert();

        void ComposeQMatFace();

        void ViewCurvature();

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
        std::vector<Eigen::Matrix3d> m_S_mat;
        std::vector<Eigen::Matrix3d> m_R_mat;
        std::vector<Eigen::Matrix3d> m_QF_mat;
    };
}
#endif //METRIC_CURVATURE_MESHIO_H
