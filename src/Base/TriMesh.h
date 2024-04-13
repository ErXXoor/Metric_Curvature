//
// Created by Hongbo on 4/5/24.
//

#ifndef METRICCURVATURE_TRIMESH_H
#define METRICCURVATURE_TRIMESH_H

#include "Base/MeshGT.h"
#include "MeshGT.h"
#include "Base/Metric.h"

namespace IGLUtils {
    class TriMesh {
    public:
        TriMesh() = default;

        ~TriMesh() = default;

        void PointArea(MeshGT *mesh);

        void CalculateCurvature(MeshGT *mesh);

        Eigen::VectorXd m_point_areas;
        Eigen::MatrixXd m_corner_areas;
    };
}
#endif //METRICCURVATURE_TRIMESH_H
