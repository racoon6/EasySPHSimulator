//
// Created by Wanghuanxiong on 2026/1/4.
//
#include "../include/SPHSimulator.h"
#include <numbers>
#include <random>

//构造函数：初始化粒子（随机分布在立方体区域）
SPHSimulator::SPHSimulator(int particleCount)
{
    _particles.resize(particleCount);

    // 1. 确定立方体的粒子分布（比如300个粒子：x=6, y=10, z=5 → 6*10*5=300）
    const int nx = 6;   // x方向粒子数
    const int ny = 10;  // y方向粒子数（决定水块高度）
    const int nz = 5;   // z方向粒子数
    const float step = 0.2f;  // 粒子间距（0.2足够大，Houdini里不会挤）

    // 2. 按网格排列粒子，组成立方体水块
    int idx = 0;
    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                if (idx >= particleCount) break;
                auto& p = _particles[idx];
                // 立方体位置：x/z从0开始，y从2.0开始（让水块在高空，自由落体到地面）
                p.pos = glm::vec3(
                    x * step,
                    y * step + 2.0f,  // y轴偏移，保证水块初始在地面上方
                    z * step
                );
                p.vel = glm::vec3(0.0f);       // 初始静止
                p.density = REST_DENSITY;      // 初始密度设为静止密度
                p.pressure = 0.0f;
                idx++;
            }
        }
    }
}

void SPHSimulator::step()
{
    buildSpatialGrid();       // 1. 构建空间网格
    computeDensityPressure(); // 2. 计算密度/压力
    computeForces();          // 3. 计算受力
    integrateParticles();     // 4. 积分更新
    handleBoundary();         // 5. 边界处理
}

//Poly6核函数（密度计算）
float SPHSimulator::kernelPoly6(const glm::vec3 &r, const float h) const noexcept
{
    const float r2 = glm::dot(r, r);//距离的平法，减少开根计算
    const float h2 = h * h;
    if (r2 > h2) return 0.0f;

    const float h3 = h2 * h;
    const float h9 = h3 * h3 * h3;
    const float x = (h2 - r2);

    return (315.0f * x * x * x)/ (64.0f * std::numbers::pi_v<float> * h9);
}

// Spiky核函数梯度（压力力）
glm::vec3 SPHSimulator::kernelSpikyGrad(const glm::vec3& r, float h) const noexcept
{
    const float rLen = glm::length(r);
    if (rLen > h || rLen < 1e-6f) return glm::vec3(0.0f);
    const float h3 = h * h * h;
    const float val = h - rLen;
    return -45.0f / (std::numbers::pi_v<float> * h3) * val * val * glm::normalize(r);
}

// Viscosity核函数拉普拉斯（粘性力）
float SPHSimulator::kernelViscosityLaplacian(const glm::vec3& r, float h) const noexcept
{
    const float rLen = glm::length(r);
    if (rLen > h) return 0.0f;
    return 45.0f / (std::numbers::pi_v<float> * h * h * h) * (h - rLen);
}

// 构建空间网格（哈希键：网格坐标转整数）
void SPHSimulator::buildSpatialGrid()
{
    _grid.clear();
    for (int i = 0; i< _particles.size(); i++)
    {
        const auto& particle = _particles[i];
        // 计算粒子所在网格坐标
        int gridX = static_cast<int>(particle.pos.x / KERNEL_RADIUS);
        int gridY = static_cast<int>(particle.pos.y / KERNEL_RADIUS);
        int gridZ = static_cast<int>(particle.pos.z / KERNEL_RADIUS);
        // 哈希键（简单拼接）
        int key = (gridX << 20) | (gridY << 10) | gridZ;
        _grid[key].push_back(i);
    }
}

void SPHSimulator::computeDensityPressure()
{
    for (int i = 0; i< _particles.size(); i++)
    {
        auto& particle = _particles[i];
        particle.density = 0.0f;

        // 找当前粒子所在网格及相邻网格的粒子（邻居）
        int gridX = static_cast<int>(particle.pos.x / KERNEL_RADIUS);
        int gridY = static_cast<int>(particle.pos.y / KERNEL_RADIUS);
        int gridZ = static_cast<int>(particle.pos.z / KERNEL_RADIUS);

        // 遍历3x3x3相邻网格
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {
                    int key = ((gridX + dx) << 20) | ((gridY + dy) << 10) | (gridZ + dz);
                    if (!_grid.count(key)) continue;
                    // 遍历网格内的邻居粒子
                    for (int j : _grid[key]) {
                        if (i == j) continue;
                        const auto& pj = _particles[j];
                        glm::vec3 r = particle.pos - pj.pos;
                        particle.density += pj.mass * kernelPoly6(r, KERNEL_RADIUS);
                    }
                }
            }
        }

        // 压力计算：P = k*(ρ - ρ0)（WCSPH状态方程）
        particle.pressure = GAS_CONSTANT * (particle.density - REST_DENSITY);
    }
}

// 计算受力（压力力+粘性力+重力）
void SPHSimulator::computeForces()
{
    const glm::vec3 gravity(0.0f, -9.81f, 0.0f);  // 重力加速度
    for (int i = 0; i< _particles.size(); i++)
    {
        auto& particle = _particles[i];
        particle.resetAcc();
        particle.acc += gravity;  // 先加重力

        // 找邻居计算压力力/粘性力
        int gridX = static_cast<int>(particle.pos.x / KERNEL_RADIUS);
        int gridY = static_cast<int>(particle.pos.y / KERNEL_RADIUS);
        int gridZ = static_cast<int>(particle.pos.z / KERNEL_RADIUS);
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {
                    int key = ((gridX + dx) << 20) | ((gridY + dy) << 10) | (gridZ + dz);
                    if (!_grid.count(key)) continue;
                    for (int j : _grid[key]) {
                        if (i == j) continue;
                        const auto& pj = _particles[j];
                        glm::vec3 r = particle.pos - pj.pos;
                        float rLen = glm::length(r);
                        if (rLen > KERNEL_RADIUS) continue;

                        // 1. 压力力（排斥力）
                        glm::vec3 pressureForce = -pj.mass * (particle.pressure + pj.pressure) / (2.0f * pj.density)
                                                  * kernelSpikyGrad(r, KERNEL_RADIUS);
                        // 2. 粘性力（阻尼）
                        glm::vec3 viscosityForce = VISCOSITY * pj.mass * (pj.vel - particle.vel) / pj.density
                                                   * kernelViscosityLaplacian(r, KERNEL_RADIUS);
                        // 累加受力（F=ma → a=F/m）
                        particle.acc += (pressureForce + viscosityForce) / particle.mass;
                    }
                }
            }
        }
    }
}

void SPHSimulator::handleBoundary()
{
    const float boundary = 0.1f;    // 边界厚度
    const float bounce = 0.8f;      // 反弹系数
    for (auto& p : _particles) {
        // x轴边界
        if (p.pos.x < boundary) {
            p.pos.x = boundary;
            p.vel.x *= -bounce;
        } else if (p.pos.x > 3.0f - boundary) {
            p.pos.x = 3.0f - boundary;
            p.vel.x *= -bounce;
        }
        // y轴边界（地面）
        if (p.pos.y < boundary) {
            p.pos.y = boundary;
            p.vel.y *= -bounce;
        }else if (p.pos.y > 3.0f - boundary) {
            p.pos.y = 3.0f - boundary;
            p.vel.y *= -bounce;
        }
        // z轴边界
        if (p.pos.z < boundary) {
            p.pos.z = boundary;
            p.vel.z *= -bounce;
        } else if (p.pos.z > 3.0f - boundary) {
            p.pos.z = 3.0f - boundary;
            p.vel.z *= -bounce;
        }
    }
}

void SPHSimulator::integrateParticles()
{
    for (auto& p : _particles) {
        p.integrate(DT);
    }
}






