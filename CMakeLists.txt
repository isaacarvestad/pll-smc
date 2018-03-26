cmake_minimum_required(VERSION 3.10.0)

project(pll-smc-project)

# Make sure we are searching /usr/local for libpll headers and library.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -I/usr/local/include -L/usr/local/lib")

add_subdirectory(lib)
add_subdirectory(app)
