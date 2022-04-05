cmake_minimum_required(VERSION 2.4)

PROJECT(qs_test)
include_directories(${PROJ_SOURCE_DIR})
add_executable(qsort qs_teste.cc quicksort.cc)