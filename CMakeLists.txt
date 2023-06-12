cmake_minimum_required(VERSION 3.20)
project(Hadrons VERSION 1.5.0 LANGUAGES CXX)

include(GNUInstallDirs)
include(FetchContent)
set(FETCHCONTENT_BASE_DIR ${CMAKE_BINARY_DIR}/deps)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

find_library(GRID NAMES Grid PATH_SUFFIXES lib REQUIRED)