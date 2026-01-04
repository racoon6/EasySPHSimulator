//
// Created by Wanghuanxiong on 2026/1/4.
//
#include <iostream>
#include <string>
#include <fstream>  // 用于文件写入
#include "include/SPHSimulator.h"

void saveParticlesToPLY(const std::span<const SPHParticle>& particles, const std::string& filename)
{
    // 打开文件（输出模式，覆盖已有文件）
    std::ofstream plyFile(filename);
    if (!plyFile.is_open()) {
        std::cerr << "Error: 无法打开文件 " << filename << " 进行写入！" << std::endl;
        return;
    }

    // 1. 写入PLY文件头部（ASCII格式，点云类型）
    plyFile << "ply" << std::endl;
    plyFile << "format ascii 1.0" << std::endl;
    plyFile << "element vertex " << particles.size() << std::endl;  // 顶点数量=粒子数量
    plyFile << "property float x" << std::endl;                     // 位置x
    plyFile << "property float y" << std::endl;                     // 位置y
    plyFile << "property float z" << std::endl;                     // 位置z
    plyFile << "end_header" << std::endl;

    // 2. 写入每个粒子的顶点数据
    for (const auto& p : particles) {
        plyFile << p.pos.x << " "
            << p.pos.y << " "
            << p.pos.z << " "
            << std::endl;
    }

    // 关闭文件
    plyFile.close();
    std::cout << "save" << filename << std::endl;
}
int main()
{
    SPHSimulator simulator(300);
    std::cout << "SPH Fluid Simulation Start..." << std::endl;
    // std::cout << "printing particles..." << std::endl;
    // std::span<const SPHParticle> particles = simulator.getParticles();;
    // for (auto const &particle : particles)
    // {
    //     std::cout << particle.pos.x <<" "<< particle.pos.y <<" " <<particle.pos.z <<" " << std::endl;
    // }
    int j = 0;
    for (int i = 0; i < 10000; i++) {
        simulator.step();
        // 输出第100步的粒子位置
        if (i % 100 == 0) {
            const auto& particles = simulator.getParticles();
            // 拼接文件名（带步数，方便序列识别）
            std::string filename = "sph_step_" + std::to_string(j++) + ".ply";
            saveParticlesToPLY(particles, filename);
            // std::cout << "Step " << i << ", Particle 0 Pos: ("
            //           << particles[0].pos.x << ", "
            //           << particles[0].pos.y << ", "
            //           << particles[0].pos.z << ")" << std::endl;
        }
    }
    std::cout << "Simulation End." << std::endl;
    return 0;
}
