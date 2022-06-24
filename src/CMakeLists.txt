cmake_minimum_required(VERSION 3.10)


project(SqueezePIR VERSION 0.0 LANGUAGES CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../bin)

add_executable(SqueezePIR main.cpp)

find_package(SEAL 3.7 REQUIRED)

target_link_libraries(SqueezePIR SEAL::seal)