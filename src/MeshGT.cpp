//
// Created by hongbo on 13/09/23.
//

#include "Base/MeshGT.h"
#include <igl/read_triangle_mesh.h>
#include <igl/principal_curvature.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/adjacency_list.h>

#include <ostream>
#include <filesystem>
#include "Base/Smoother.h"
#include <igl/copyleft/comiso/nrosy.h>
#include <igl/uniformly_sample_two_manifold.h>

namespace IGLUtils {
    namespace fs = std::__fs::filesystem;

    MeshGT::MeshGT() {
        m_v = std::make_shared<Eigen::MatrixXd>();
        m_f = std::make_shared<Eigen::MatrixXi>();
        m_n = std::make_shared<Eigen::MatrixXd>();
        m_metric = std::make_shared<Metric>();
        m_min_pd = std::make_shared<Eigen::MatrixXd>();
        m_max_pd = std::make_shared<Eigen::MatrixXd>();
        m_max_pv = std::make_shared<Eigen::VectorXd>();
        m_min_pv = std::make_shared<Eigen::VectorXd>();
        tmp_pd = std::make_shared<Eigen::MatrixXd>();
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

    bool MeshGT::NormalizeMesh() {
        Eigen::RowVector3d min_corner, max_corner;
        min_corner = m_v->colwise().minCoeff();
        max_corner = m_v->colwise().maxCoeff();
        auto diag_size = (max_corner - min_corner).norm();
        if(diag_size>1.2||diag_size<0.8)
        {
            auto mid_point = (max_corner + min_corner) / 2;
            auto scale = 1.0 / diag_size;
            m_v->rowwise() -= mid_point;
            (*m_v) *= scale;
        }
        return true;
    }


    bool MeshGT::CalculateCurvature() {
        try {
            igl::principal_curvature(*m_v, *m_f, *m_max_pd, *m_min_pd, *m_max_pv, *m_min_pv);
        }
        catch (const std::exception &e) {
            // Catch and handle exceptions
            std::cerr << "An exception occurred: " << e.what() << std::endl;
            return false;
        }
        return true;
    }

    void MeshGT::ProcessMetric(unsigned int smooth_ring, unsigned int smooth_iter) {
        std::vector<std::vector<int>> adj;
        igl::adjacency_list(*m_f, adj);

        std::vector<std::set<int>> ring_neighbor;
        for (auto &i: adj) {
            std::set<int> ring;
            for (auto &j: i) {
                ring.insert(j);
            }
            ring_neighbor.emplace_back(ring);
        }
        if (smooth_ring > 1) {
            for (auto i = 0; i < smooth_ring - 1; i++) {
                std::vector<std::set<int>> new_ring_neighbor;
                std::copy(ring_neighbor.begin(), ring_neighbor.end(), std::back_inserter(new_ring_neighbor));
                for (auto set_i = 0; set_i < ring_neighbor.size(); set_i++) {
                    for (auto v_id: ring_neighbor[set_i]) {
                        new_ring_neighbor[set_i].insert(ring_neighbor[v_id].begin(), ring_neighbor[v_id].end());
                    }
                }
                ring_neighbor = new_ring_neighbor;
            }
        } else {
            for (int i = 0; i < ring_neighbor.size(); i++) {
                ring_neighbor[i].insert(i);
            }

        }

        smt_n = Smoother::SmoothVec(smooth_iter, ring_neighbor, m_n);
        smt_max_pd = Smoother::SmoothVec(smooth_iter, ring_neighbor, m_max_pd);
        smt_min_pd = Smoother::SmoothVec(smooth_iter, ring_neighbor, m_min_pd);

        Eigen::VectorXd tmp_max_pv = m_max_pv->array().abs();
        tmp_max_pv = tmp_max_pv.cwiseMin(25.0).cwiseMax(1.0);
        std::shared_ptr<Eigen::VectorXd> max_ptr = std::make_shared<Eigen::VectorXd>(tmp_max_pv);

        Eigen::VectorXd tmp_min_pv = m_min_pv->array().abs();
        tmp_min_pv = tmp_min_pv.cwiseMin(25.0).cwiseMax(1.0);
        std::shared_ptr<Eigen::VectorXd> min_ptr = std::make_shared<Eigen::VectorXd>(tmp_min_pv);

        auto smt_max_pv = Smoother::SmoothScalar(smooth_iter, ring_neighbor, max_ptr);
        auto smt_min_pv = Smoother::SmoothScalar(smooth_iter, ring_neighbor, min_ptr);

        m_metric->ComposeMetric(smt_n,
                                smt_min_pd,
                                smt_max_pd,
                                smt_min_pv,
                                smt_max_pv);
    }

    void MeshGT::SaveTmpPD(const std::string &filepath) {
        std::ofstream out(filepath);
        for (auto i = 0; i < tmp_pd->rows(); i++) {
            out << (*tmp_pd)(i, 0) << "," << (*tmp_pd)(i, 1) << "," << (*tmp_pd)(i, 2) << std::endl;
        }
    }

    void MeshGT::SaveCurvature(const std::string &filepath) {
        fs::path path(filepath);
        path = path / (m_filename + "_cur.csv");
        std::ofstream out(path.string());

        //Header
        out << "min_cur" << "," << "min_pd_x" << "," << "min_pd_y" << "," << "min_pd_z" << ",";
        out << "max_cur" << "," << "max_pd_x" << "," << "max_pd_y" << "," << "max_pd_z" << ",";
        out << "normal_x" << "," << "normal_y" << "," << "normal_z" << std::endl;
        for (auto i = 0; i < m_min_pd->rows(); i++) {
            out << (*m_min_pv)(i) << "," << (*m_min_pd)(i, 0) << "," << (*m_min_pd)(i, 1) << "," << (*m_min_pd)(i, 2)
                << ",";
            out << (*m_max_pv)(i) << "," << (*m_max_pd)(i, 0) << "," << (*m_max_pd)(i, 1) << "," << (*m_max_pd)(i, 2)
                << ",";
            out << (*m_n)(i, 0) << "," << (*m_n)(i, 1) << "," << (*m_n)(i, 2) << std::endl;

        }

        out.close();
    }

    void MeshGT::SaveMetric(const std::string &filepath, bool save_sr) {
        fs::path path(filepath);
        auto path_m = path / (m_filename + "_m.csv");
        m_metric->SaveMetric(path_m.string());

        if (save_sr) {
            auto path_s = path / (m_filename + "_s.csv");
            auto path_r = path / (m_filename + "_r.csv");
            m_metric->SaveSR(path_s.string(), path_r.string());
        }
    }

    void MeshGT::SaveMeshInfo(const std::string &filepath) {
        //Save vertex
        fs::path path(filepath);
        auto path_s = path / (m_filename + "_vert.csv");
        std::ofstream out_s(path_s.string());
        // vertex x,y,z
        for (int i = 0; i < m_v->rows(); i++) {
            auto v = m_v->row(i);
            out_s << v(0) << "," << v(1) << "," << v(2) << std::endl;
        }
        out_s.close();

        //Save face
        auto path_f = path / (m_filename + "_face.csv");
        std::ofstream out_f(path_f.string());
        // face v1,v2,v3
        for (int i = 0; i < m_f->rows(); i++) {
            auto f = m_f->row(i);
            out_f << f(0) << "," << f(1) << "," << f(2) << std::endl;
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
        viewer.data().add_edges(*m_v + *smt_min_pd * avg, *m_v, red);

        // Draw a blue segment parallel to the minimal curvature direction
        viewer.data().add_edges(*m_v + *smt_max_pd * avg, *m_v, blue);

        // Hide wireframe
        viewer.data().show_lines = false;

        viewer.launch();
    }

    void MeshGT::SaveMesh(const std::string &filepath) {
        std::filesystem::path path(filepath);
        auto path_s = path / (m_filename + "_norm.obj");
        igl::writeOBJ(path_s.string(), *m_v, *m_f);
    }

}