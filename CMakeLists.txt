cmake_minimum_required(VERSION 3.20)
project(EasySPHSimulator)

# 关键：强制MSVC编译器使用UTF-8编码，彻底解决C4819警告
if (MSVC)
    add_compile_options("$<$<C_COMPILER_ID:MSVC>:/utf-8>")
    add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/utf-8>")
endif()

set(CMAKE_CXX_STANDARD 20)
# 强制要求该C++标准，若编译器不支持则报错（确保生效）
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
add_executable(EasySPHSimulator main.cpp
        include/SPHParticle.h
        include/SPHSimulator.h
        src/SPHSimulator.cpp
)