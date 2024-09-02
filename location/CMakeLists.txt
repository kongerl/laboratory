cmake_minimum_required(VERSION 3.25)
project(location C)

set(CMAKE_C_STANDARD 11)

add_executable(location main.c all.h satellite.c Spp.c bdgim.c UNB3m.c)
