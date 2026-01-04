cmake_minimum_required(VERSION 3.25)
project(EasySPHSimulator)

set(CMAKE_CXX_STANDARD 20)
# 强制要求该C++标准，若编译器不支持则报错（确保生效）
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
add_executable(EasySPHSimulator main.cpp
        include/SPHParticle.h
        include/SPHSimulator.h
        src/SPHSimulator.cpp
)
