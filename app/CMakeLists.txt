cmake_minimum_required(VERSION 3.10.0)

include_directories(include ../lib/include)

add_executable(pll-smc main.cpp)

target_link_libraries(pll-smc LINK_PUBLIC pll-smc-lib)
