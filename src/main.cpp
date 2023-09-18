//
// Created by hongbo on 13/09/23.
//
#include <igl/opengl/glfw/Viewer.h>
#include "MeshIO/MeshGT.h"
#include <string>

int main(int argc, char *argv[]) {
    std::string filename = "/home/hongbo/Desktop/code/datasets/kitten_remesh.obj";

    IGLUtils::MeshGT mesh;
    if(mesh.LoadMesh(filename))
    {
        mesh.CalculateCurvature();
//    mesh.ViewCurvature();
        mesh.SaveMeshInfo("../../tmp/");
    mesh.SaveMetric("../../tmp/");
    }
    else
    {
        std::cerr<<"Mesh load failed!"<<std::endl;
    }
}