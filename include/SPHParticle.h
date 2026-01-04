//
// Created by Wanghuanxiong on 2025/12/17.
//

#ifndef SPHPARTICLE_H
#define SPHPARTICLE_H
#include "glm/glm.hpp"

// SPH粒子类（WCSPH：弱可压缩SPH）
struct SPHParticle
{
    glm::vec3 pos;//粒子位置
    glm::vec3 vel;//粒子速度
    glm::vec3 acc;//加速度
    float density;//密度
    float pressure;//压力  该粒子所在位置的 “局部压强值”（标量）用于代表粒子周围流体的压缩 / 膨胀状态
    float mass = 8.0f;//质量（取固定值）

    // 重置加速度（每次仿真步长前）
    constexpr void  resetAcc() noexcept
    {
        acc = glm::vec3(0.0f);
    }

    constexpr void integrate(float dt) noexcept
    {
        vel += dt * acc;
        pos += dt * vel;
        // 可选：速度阻尼（防止数值爆炸）
        vel *= 0.999f;
    }
};
#endif //SPHPARTICLE_H
