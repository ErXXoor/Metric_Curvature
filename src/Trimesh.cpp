//
// Created by Hongbo on 4/5/24.
//
#include "Base/TriMesh.h"
#include "igl/massmatrix.h"
// This way of computing it tends to be faster than using %
#define NEXT_MOD3(i) ((i) < 2 ? (i) + 1 : (i) - 2)
#define PREV_MOD3(i) ((i) > 0 ? (i) - 1 : (i) + 2)
namespace IGLUtils {

    void rot_coord_sys(const Eigen::Vector3d &old_u,
                       const Eigen::Vector3d &old_v,
                       const Eigen::Vector3d &new_norm,
                       Eigen::Vector3d &new_u, Eigen::Vector3d &new_v) {
        new_u = old_u;
        new_v = old_v;
        Eigen::Vector3d old_norm = old_u.cross(old_v);
        float ndot = old_norm.dot(new_norm);
        if (ndot <= -1.0f) {
            new_u = -new_u;
            new_v = -new_v;
            return;
        }

        Eigen::Vector3d perp_old = new_norm - ndot * old_norm;
        Eigen::Vector3d dperp = 1.0f / (1 + ndot) * (old_norm + new_norm);
        new_u -= dperp * (new_u.dot(perp_old));
        new_v -= dperp * (new_v.dot(perp_old));
    }

    void proj_curv(const Eigen::Vector3d &old_u,
                   const Eigen::Vector3d &old_v,
                   float old_ku, float old_kuv, float old_kv,
                   const Eigen::Vector3d &new_u, const Eigen::Vector3d &new_v,
                   float &new_ku, float &new_kuv, float &new_kv) {
        Eigen::Vector3d r_new_u, r_new_v;
        rot_coord_sys(new_u, new_v, old_u.cross(old_v), r_new_u, r_new_v);

        float u1 = r_new_u.dot(old_u);
        float v1 = r_new_u.dot(old_v);
        float u2 = r_new_v.dot(old_u);
        float v2 = r_new_v.dot(old_v);
        new_ku = old_ku * u1 * u1 + old_kuv * (2.0f * u1 * v1) + old_kv * v1 * v1;
        new_kuv = old_ku * u1 * u2 + old_kuv * (u1 * v2 + u2 * v1) + old_kv * v1 * v2;
        new_kv = old_ku * u2 * u2 + old_kuv * (2.0f * u2 * v2) + old_kv * v2 * v2;
    }

    void diagonalize_cur(const Eigen::Vector3d &old_u, const Eigen::Vector3d &old_v,
                         float ku, float kuv, float kv,
                         const Eigen::Vector3d &new_n,
                         Eigen::Vector3d &new_pd1, Eigen::Vector3d &new_pd2,
                         float &k1, float &k2) {
        Eigen::Vector3d r_old_u, r_old_v;
        rot_coord_sys(old_u, old_v, new_n, r_old_u, r_old_v);

        double c = 1., s = 0., tt = 0.;
        if (kuv != 0.0) {
            auto h = 0.5 * (kv - ku) / kuv;
            tt = (h < 0.0) ? 1.0 / (h - std::sqrt(1.0 + h * h)) : 1.0 / (h + sqrt(1.0 + h * h));
            c = 1.0 / sqrt(1.0 + tt * tt);
            s = tt * c;
        }

        k1 = ku - tt * kuv;
        k2 = kv + tt * kuv;

        if (std::fabs(k1) >= std::fabs(k2)) {
            new_pd1 = c * r_old_u - s * r_old_v;
        } else {
            std::swap(k1, k2);
            new_pd1 = s * r_old_u + c * r_old_v;
        }
        new_pd2 = new_n.cross(new_pd1);
    }

    void TriMesh::PointArea(MeshGT *mesh) {
        Eigen::MatrixXd v = *(mesh->m_v);
        Eigen::MatrixXi f = *(mesh->m_f);

        m_point_areas = Eigen::VectorXd::Zero(v.rows());
        m_corner_areas = Eigen::MatrixXd::Zero(f.rows(), 3);

#pragma omp parallel for
        for (int i = 0; i < f.rows(); i++) {
            std::vector<Eigen::Vector3d> e = {v.row(f(i, 2)) - v.row(f(i, 1)),
                                              v.row(f(i, 0)) - v.row(f(i, 2)),
                                              v.row(f(i, 1)) - v.row(f(i, 0))};

            float area = 0.5f * e[0].cross(e[1]).norm();
            std::vector<double> l2 = {e[0].norm(), e[1].norm(), e[2].norm()};

            std::vector<double> bcw = {l2[0] * (l2[1] + l2[2] - l2[0]),
                                       l2[1] * (l2[2] + l2[0] - l2[1]),
                                       l2[2] * (l2[0] + l2[1] - l2[2])};

            if (bcw[0] <= 0.0f) {
                m_corner_areas(i, 1) = -0.25f * l2[2] * area / (e[0].dot(e[2]));
                m_corner_areas(i, 2) = -0.25f * l2[1] * area / (e[0].dot(e[1]));
                m_corner_areas(i, 0) = area - m_corner_areas(i, 1) - m_corner_areas(i, 2);
            } else if (bcw[1] <= 0.0f) {
                m_corner_areas(i, 2) = -0.25f * l2[0] * area / (e[1].dot(e[0]));
                m_corner_areas(i, 0) = -0.25f * l2[2] * area / (e[1].dot(e[2]));
                m_corner_areas(i, 1) = area - m_corner_areas(i, 2) - m_corner_areas(i, 0);
            } else if (bcw[2] <= 0.0f) {
                m_corner_areas(i, 0) = -0.25f * l2[1] * area / (e[2].dot(e[1]));
                m_corner_areas(i, 1) = -0.25f * l2[0] * area / (e[2].dot(e[0]));
                m_corner_areas(i, 2) = area - m_corner_areas(i, 0) - m_corner_areas(i, 1);
            } else {
                float scale = 0.5f * area / (bcw[0] + bcw[1] + bcw[2]);
                for (int j = 0; j < 3; j++)
                    m_corner_areas(i, j) = scale * (bcw[NEXT_MOD3(j)] + bcw[PREV_MOD3(j)]);
            }

#pragma omp atomic
            m_point_areas(f(i, 0)) += m_corner_areas(i, 0);
#pragma omp atomic
            m_point_areas(f(i, 1)) += m_corner_areas(i, 1);
#pragma omp atomic
            m_point_areas(f(i, 2)) += m_corner_areas(i, 2);
        }
    }


    void TriMesh::CalculateCurvature(MeshGT *mesh) {
        Eigen::MatrixXd v = *(mesh->m_v);
        Eigen::MatrixXi f = *(mesh->m_f);
        Eigen::MatrixXd n = *(mesh->m_n);
        Eigen::SparseMatrix<double> M;

        Eigen::MatrixXd pdir1 = Eigen::MatrixXd::Zero(v.rows(), 3);
        Eigen::MatrixXd pdir2 = Eigen::MatrixXd::Zero(v.rows(), 3);
        for (auto i = 0; i < mesh->m_f->rows(); i++) {
            pdir1.row(f(i, 0)) = v.row(f(i, 1)) - v.row(f(i, 0));
            pdir2.row(f(i, 0)) = v.row(f(i, 2)) - v.row(f(i, 1));
            pdir1.row(f(i, 1)) = v.row(f(i, 0)) - v.row(f(i, 2));
        }
        Eigen::VectorXd k1 = Eigen::VectorXd::Zero(v.rows());
        Eigen::VectorXd k2 = Eigen::VectorXd::Zero(v.rows());
        Eigen::VectorXd k12 = Eigen::VectorXd::Zero(v.rows());

        PointArea(mesh);

#pragma omp parallel for
        for (int i = 0; i < v.rows(); i++) {
            Eigen::Vector3d tmp_dir = pdir1.row(i);
            Eigen::Vector3d tmp_n = n.row(i);
            pdir1.row(i) = tmp_dir.cross(tmp_n);
            pdir1.row(i).normalize();
            pdir2.row(i) = tmp_n.cross(tmp_dir);
        }

#pragma omp parallel for
        for (int i = 0; i < f.rows(); i++) {
            Eigen::Vector3d e0 = v.row(f(i, 2)) - v.row(f(i, 1));
            Eigen::Vector3d e1 = v.row(f(i, 0)) - v.row(f(i, 2));
            Eigen::Vector3d e2 = v.row(f(i, 1)) - v.row(f(i, 0));

            Eigen::Vector3d t = e0;
            t.normalize();
            Eigen::Vector3d tmp_n = e0.cross(e1);
            Eigen::Vector3d b = tmp_n.cross(t);
            b.normalize();

            Eigen::Vector3d m = Eigen::Vector3d::Zero();
            Eigen::Matrix3d w = Eigen::Matrix3d::Zero();
            Eigen::Matrix3d edge;
            edge << e0, e1, e2;
            for (auto j = 0; j < 3; j++) {
                auto u = edge.col(j).dot(t);
                auto v = edge.col(j).dot(b);
                w(0, 0) += u * u;
                w(0, 1) += u * v;
                w(1, 0) += u * v;

                Eigen::Vector3d dn = n.row(f(i, PREV_MOD3(j))) - n.row(f(i, NEXT_MOD3(j)));
                double dnu = dn.dot(t);
                double dnv = dn.dot(b);
                m(0) += dnu * u;
                m(1) += dnu * v + dnv * u;
                m(2) += dnv * v;
            }
            w(1, 1) = w(0, 0) + w(2, 2);
            w(1, 2) = w(0, 1);

            m = w.ldlt().solve(m);

            for (int j = 0; j < 3; j++) {
                int vj = f(i, j);
                float c1, c12, c2;
                proj_curv(t, b, m(0), m(1), m(2), pdir1.row(vj), pdir2.row(vj), c1, c12, c2);
                auto wt = m_corner_areas(i, j) / m_point_areas(vj);
#pragma omp atomic
                k1(vj) += wt * c1;
#pragma omp atomic
                k12(vj) += wt * c12;
#pragma omp atomic
                k2(vj) += wt * c2;
            }
        }

        // Compute principal directions and curvatures at each vertex
#pragma omp parallel for
        for (int i = 0; i < v.rows(); i++) {
            Eigen::Vector3d new_pd1, new_pd2;
            float new_k1, new_k2;
            diagonalize_cur(pdir1.row(i), pdir2.row(i),
                            k1(i), k12(i), k2(i),
                            n.row(i), new_pd1, new_pd2,
                            new_k1, new_k2);
            pdir1.row(i) = new_pd1;
            pdir2.row(i) = new_pd2;
            k1(i) = new_k1;
            k2(i) = new_k2;
        }
        mesh->m_min_pd = std::make_shared<Eigen::MatrixXd>(pdir2);
        mesh->m_max_pd = std::make_shared<Eigen::MatrixXd>(pdir1);
        mesh->m_min_pv = std::make_shared<Eigen::VectorXd>(k2);
        mesh->m_max_pv = std::make_shared<Eigen::VectorXd>(k1);
    }
}