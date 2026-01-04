//
// Created by Wanghuanxiong on 2026/1/4.
//
#include <iostream>

#include "include/SPHSimulator.h"

int main()
{
    SPHSimulator simulator(100);
    std::cout << "SPH Fluid Simulation Start..." << std::endl;
    // std::cout << "printing particles..." << std::endl;
    // std::span<const SPHParticle> particles = simulator.getParticles();;
    // for (auto const &particle : particles)
    // {
    //     std::cout << particle.pos.x <<" "<< particle.pos.y <<" " <<particle.pos.z <<" " << std::endl;
    // }

    for (int i = 0; i < 1000; i++) {
        simulator.step();
        // 输出第100步的粒子位置
        if (i % 100 == 0) {
            const auto& particles = simulator.getParticles();
            std::cout << "Step " << i << ", Particle 0 Pos: ("
                      << particles[0].pos.x << ", "
                      << particles[0].pos.y << ", "
                      << particles[0].pos.z << ")" << std::endl;
        }
    }
    return 0;
}
