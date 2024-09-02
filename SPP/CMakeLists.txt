cmake_minimum_required(VERSION 3.25)
project(SPP C)

set(CMAKE_C_STANDARD 99)


add_executable(SPP main.c spp.c spp.h UNB3m.h UNB3m.c bdgim.c bdgim.h)
