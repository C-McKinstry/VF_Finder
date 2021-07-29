# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.13.4)

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(sanitizer-populate NONE)

include(ExternalProject)
ExternalProject_Add(sanitizer-populate
                     "UPDATE_DISCONNECTED" "False" "GIT_REPOSITORY" "https://github.com/arsenm/sanitizers-cmake.git" "GIT_TAG" "99e159ec9bc8dd362b08d18436bd40ff0648417b"
                    SOURCE_DIR          "/home/cmckinstry/projecta/_deps/sanitizer-src"
                    BINARY_DIR          "/home/cmckinstry/projecta/_deps/sanitizer-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
)
