cmake_minimum_required(VERSION 3.0)

project(zbessel VERSION 1.0 LANGUAGES CXX)
add_library(zbessel STATIC zbessel.cc)
target_include_directories(zbessel PUBLIC "${PROJECT_SOURCE_DIR}")
