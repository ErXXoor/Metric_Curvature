//
// Created by hongbo on 13/09/23.
//

#include "MeshIO/MeshGT.h"
#include <igl/read_triangle_mesh.h>
#include <igl/principal_curvature.h>
#include <igl/opengl/glfw/Viewer.h>
#include <ostream>
#include <filesystem>

namespace IGLUtils {
    namespace fs = std::filesystem;

    MeshGT::MeshGT() {
        m_v = std::make_shared<Eigen::MatrixXd>();
        m_f = std::make_shared<Eigen::MatrixXi>();
        m_n = std::make_shared<Eigen::MatrixXd>();
    }

    bool MeshGT::LoadMesh(const std::string &filename) {
        m_filename = fs::path(filename).stem().string();
        auto success = igl::read_triangle_mesh(filename, *m_v, *m_f);
        if (!success) {
            std::cerr << "Failed to load mesh from " << filename << std::endl;
            return false;
        }
        igl::per_vertex_normals(*m_v, *m_f, *m_n);
        return success;
    }

    void MeshGT::CalculateCurvature() {
        m_min_pd = std::make_shared<Eigen::MatrixXd>();
        m_max_pd = std::make_shared<Eigen::MatrixXd>();
        m_max_pv = std::make_shared<Eigen::VectorXd>();
        m_min_pv = std::make_shared<Eigen::VectorXd>();
        igl::principal_curvature(*m_v, *m_f, *m_max_pd, *m_min_pd, *m_max_pv, *m_min_pv);
    }

    void MeshGT::ComposeSRMatVert() {
        if (m_min_pd == nullptr || m_max_pd == nullptr) {
            CalculateCurvature();
        }

        for (int i = 0; i < m_n->rows(); ++i) {
            Eigen::Matrix3d mat_r;
            mat_r << m_min_pd->row(i), m_max_pd->row(i), m_n->row(i);
            m_R_mat.emplace_back(mat_r);

            Eigen::DiagonalMatrix<double, 3> diag;
            diag.diagonal() << std::sqrt(abs(m_min_pv->coeff(i))), std::sqrt(abs(m_max_pv->coeff(i))), 0.0;
            m_S_mat.emplace_back(diag);
        }
    }

    void MeshGT::ComposeQMatFace() {
        if (m_S_mat.empty() || m_R_mat.empty()) {
            ComposeSRMatVert();
        }
        for (auto i = 0; i < m_f->rows(); i++) {
            auto v_ids = m_f->row(i);
            Eigen::Matrix3d FS;
            FS.setZero();
            for (auto j = 0; j < v_ids.size(); j++) {
                FS += m_S_mat[v_ids[j]];
            }
            FS /= double(v_ids.size());
            auto FR = m_R_mat[v_ids[0]];
            auto Q_mat = FS * FR;
            m_QF_mat.emplace_back(Q_mat);
        }
    }


    void MeshGT::SaveMetric(const std::string &filepath) {
        std::filesystem::path path(filepath);
        path = path / (m_filename + "_cur.csv");
        std::ofstream out(path.string());
        // vertex id,max,maxcurvature,min,mincurvature,normal
        //csv reader has problem with the first line so no header
//        std::vector<std::string> col_name = {"vertex", "max_curvature","max_direction_x","max_direction_y","max_direction_z",
//                                             "min_curvature","min_direction_x","min_direction_y","min_direction_z",
//                                             "normal_x","normal_y","normal_z"};
//        for(auto &name:col_name)
//        {
//            out<<name<<",";
//        }
//        out<<std::endl;

        for(int i=0;i<m_v->rows();i++)
        {
            auto max_pd = m_max_pd->row(i);
            auto min_pd = m_min_pd->row(i);
            auto normal = m_n->row(i);
            out<<i<<",";
            out<<m_max_pv->coeff(i)<<",";
            out<<max_pd(0)<<","<<max_pd(1)<<","<<max_pd(2)<<",";
            out<<m_min_pv->coeff(i)<<",";
            out<<min_pd(0)<<","<<min_pd(1)<<","<<min_pd(2)<<",";
            out<<normal(0)<<","<<normal(1)<<","<<normal(2)<<std::endl;
        }
        out.close();
    }

    void MeshGT::SaveMeshInfo(const std::string &filepath) {
        //Save vertex
        std::filesystem::path path(filepath);
        auto path_s = path / (m_filename + "_vert.csv");
        std::ofstream out_s(path_s.string());
        // vertex x,y,z
        for(int i=0;i<m_v->rows();i++)
        {
            auto v = m_v->row(i);
            out_s<<v(0)<<","<<v(1)<<","<<v(2)<<std::endl;
        }
        out_s.close();

        //Save face
        auto path_f = path / (m_filename + "_face.csv");
        std::ofstream out_f(path_f.string());
        // face v1,v2,v3
        for(int i=0;i<m_f->rows();i++)
        {
            auto f = m_f->row(i);
            out_f<<f(0)<<","<<f(1)<<","<<f(2)<<std::endl;
        }
        out_f.close();
    }

    void MeshGT::ViewCurvature() {
        igl::opengl::glfw::Viewer viewer;
        viewer.data().set_mesh(*m_v, *m_f);
        // Average edge length for sizing
        const double avg = igl::avg_edge_length(*m_v, *m_f);

        // Draw a red segment parallel to the maximal curvature direction
        const Eigen::RowVector3d red(0.8, 0.2, 0.2), blue(0.2, 0.2, 0.8);
        viewer.data().add_edges(*m_v + *m_min_pd * avg, *m_v - *m_min_pd * avg, red);

        // Draw a blue segment parallel to the minimal curvature direction
        viewer.data().add_edges(*m_v + *m_max_pd * avg, *m_v - *m_max_pd * avg, blue);

        // Hide wireframe
        viewer.data().show_lines = false;

        viewer.launch();
    }

}