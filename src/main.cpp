//
// Created by hongbo on 13/09/23.
//
#include <igl/opengl/glfw/Viewer.h>
#include "Base/MeshGT.h"
#include <string>
#include <fstream>
#include <filesystem>
#include <thread>
#include <chrono>

namespace fs = std::filesystem;

bool RunPipeline(const std::string &filename, const std::string &outpath) {
    auto file_idx = fs::path(filename).parent_path().filename().string();
    IGLUtils::MeshGT mesh;
    if (!mesh.LoadMesh(filename)) {
        std::cout << "file_idx: " << file_idx << " " << "load mesh failed" << std::endl;
        return false;
    }
    if(!mesh.NormalizeMesh()){
        std::cout<<"file_idx: "<<file_idx<<" "<<"normalize mesh failed"<<std::endl;
        return false;
    }
    if (!mesh.CalculateCurvature()) {
        std::cout << "file_idx: " << file_idx << " " << "calculate curvature failed" << std::endl;
        return false;
    }

    try {
        mesh.ProcessMetric(0, 0);
//        mesh.SaveMetric(outpath);
        mesh.SaveCurvature(outpath);
        mesh.SaveMeshInfo(outpath);
//        mesh.SaveMesh(outpath);
    }
    catch (const std::exception &e) {
        std::cout << "file_idx: " << file_idx << " " << e.what() << " " << "save failed" << std::endl;
        return false;
    }
    std::cout << "file_idx: " << file_idx << " " << "success" << std::endl;
    return true;
}

int main(int argc, char *argv[]) {
    auto start = std::chrono::steady_clock::now();

    std::string list_path = "/home/hongbo/Desktop/code/Metric_Curvature/res/thing10k_list.txt";
    std::filesystem::path base_path = "/media/hongbo/45ad552c-e83b-4f01-9864-7d87cfa1377e/hongbo/Thing10k_tetwild/10k_surface/";

    std::vector<std::string> folder_list;
    //read txt file to folder_list
    std::ifstream fin(list_path);
    std::string line;
    while(std::getline(fin,line))
    {
        folder_list.push_back(line);
    }
    fin.close();

    const int num_thread = 16;
    std::vector<std::thread> threads;

    for (auto i=0;i<folder_list.size();i+=num_thread)
    {
        for(auto j=0;j<num_thread;j++)
        {
            if(i+j>=folder_list.size())
                break;
            auto obj_name = folder_list[i+j]+"_norm_rxn.obj";
            auto obj_path = base_path / folder_list[i+j]/obj_name;
            auto output_path = base_path / folder_list[i+j];
            threads.emplace_back(RunPipeline, obj_path, output_path);
        }
        if(!threads.empty())
        {
            for(auto &t:threads)
                t.join();
        }
        threads.clear();
        std::cout<<"iteration "<<i<<" finished"<<std::endl;

    }

//    auto obj_path = "/home/hongbo/Desktop/code/Metric_Curvature/res/111240_sf.obj";
//    auto output_path = "../../tmp";
//    RunPipeline(obj_path, output_path);

    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "final duration: " << duration.count() << std::endl;
}