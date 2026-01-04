//
// Created by Wanghuanxiong on 2026/1/3.
//

#ifndef SPHSIMULATOR_H
#define SPHSIMULATOR_H
#include <span>
#include <vector>
#include "SPHParticle.h"


class SPHSimulator {
public:
    static constexpr float REST_DENSITY = 1000.0f; //静止密度（水）
    static constexpr float GAS_CONSTANT = 40000.0f; //气体常数（控制压缩性）
    static constexpr float VISCOSITY = 0.01f;//粘性系数
    static constexpr float KERNEL_RADIUS = 0.3f;//核函数半径（邻居搜索范围）
    static constexpr float DT = 0.0005f;//仿真步长

    // SPHSimulator();
    explicit SPHSimulator(int particleCount);

    //仿真主循环（单步）
    void step();

    //获取粒子列表（只读）
    std::span<const SPHParticle> getParticles() const noexcept
    {
        return _particles;
    }
private:
    std::vector<SPHParticle> _particles;//粒子数组
    //网格哈希（邻居搜索优化：按空间网格划分粒子）
    std::unordered_map<int, std::vector<int>> _grid;

    //核函数：Poly（密度计算）
    float kernelPoly6(const glm::vec3 &r,float h) const noexcept;
    //核函数梯度：Spiky（压力力计算）
    glm::vec3 kernelSpikyGrad(const glm::vec3 &r,float h) const noexcept;
    // 核函数拉普拉斯：Viscosity（粘性力计算）
    float kernelViscosityLaplacian(const glm::vec3 &r, float h) const noexcept;

    //1.构建空间网格（邻居搜索的前置）
    void buildSpatialGrid();
    //2.计算所有粒子的密度/压力
    void computeDensityPressure();
    // 3. 计算所有粒子的受力（压力力+粘性力+重力）
    void computeForces();
    // 4. 积分更新粒子状态
    void integrateParticles();
    // 5. 边界处理（简单弹性边界）
    void handleBoundary();

};



#endif //SPHSIMULATOR_H
