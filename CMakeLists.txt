## Copyright 2009-2021 Intel Corporation
## SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.13)
FIND_PACKAGE(embree 3.0 REQUIRED)

#add_library(GMSH ../gmsh-4.8.4-Linux64-sdk/include)
ADD_EXECUTABLE(VF ~/embree/kernels/embree.rc ViewFactor.cpp functions.h)
add_subdirectory(MshIO)
target_link_libraries(VF mshio::mshio)
TARGET_LINK_LIBRARIES(VF embree)
#TARGET_LINK_LIBRARIES(HEmbree ${PROJECT_SOURCE_DIR}/../gmsh-4.8.4-Linux64-sdk/include)
TARGET_INCLUDE_DIRECTORIES(VF PRIVATE ${PROJECT_SOURCE_DIR})
## SET_PROPERTY(TARGET minimal PROPERTY FOLDER tutorials/single)
## SET_PROPERTY(TARGET minimal APPEND PROPERTY COMPILE_FLAGS " ${FLAGS_LOWEST}")
## INSTALL(TARGETS minimal DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT examples)
## SIGN_TARGET(minimal)
