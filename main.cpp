//
// Created by Wanghuanxiong on 2026/1/4.
//
#include <iostream>
#include <string>
#include <fstream>  // 用于文件写入
#include "include/SPHSimulator.h"

void saveParticlesToPLY(const std::span<const SPHParticle>& particles, const std::string& filename)
{
    std::cout << filename<<std::endl;
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

