project(eigen)

cmake_minimum_required(VERSION 2.6)

set (CMAKE_CXX_STANDARD 11)

find_package(VISP REQUIRED)

set(cpp_
  main.cpp )


 
foreach(cpp ${cpp_})
  visp_add_target(${cpp})
endforeach()
