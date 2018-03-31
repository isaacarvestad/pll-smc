cmake_minimum_required(VERSION 3.10.0)

include_directories(include libpll)

file(GLOB SOURCES "src/*.cpp")

add_library(pll-smc-lib STATIC ${SOURCES})

target_link_libraries(pll-smc-lib libpll.a)
