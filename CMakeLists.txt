include(FetchContent)

set(FETCHCONTENT_QUIET OFF)
set(EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external")
set(FETCHCONTENT_UPDATES_DISCONNECTED TRUE)

project(BoundaryDeformationMaps)
cmake_minimum_required(VERSION 3.19)


set (CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++17")
#set (CMAKE_BUILD_TYPE Release)

#set(CHOLMOD_INCLUDE_DIR /opt/local/include/)
#include(${CMAKE_SOURCE_DIR}/cmake/FindSUITESPARSE.cmake)

#message(CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH})
#set(CMAKE_PREFIX_PATH "/opt/local/include/")
#message(CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH})

#NOTE: if you get CHOLMOD issues again, do this:
# export OpenMP_ROOT=$(brew --prefix)/opt/libomp

#find_package (SUITESPARSE REQUIRED MODULE)
find_package (CGAL REQUIRED NO_MODULE)

if(NOT TARGET Eigen3::Eigen)
    FetchContent_Declare(eigen
            GIT_REPOSITORY https://gitlab.com/libeigen/eigen
            GIT_TAG 5.0.0
            #GIT_REPOSITORY https://gitlab.com/mheistermann/eigen
            #GIT_TAG fix/spqr-index-vs-storageindex
            SOURCE_DIR "${EXTERNAL_DIR}/eigen"
    )
    FetchContent_Populate(eigen)
    message("Downloaded Eigen3 to ${eigen_SOURCE_DIR}")
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen3::Eigen INTERFACE "$<BUILD_INTERFACE:${eigen_SOURCE_DIR}>")
    #target_compile_definitions(Eigen3::Eigen INTERFACE -DEIGEN_HAS_STD_RESULT_OF=0)
endif()




#if(NOT TARGET OpenVolumeMesh::OpenVolumeMesh)
if(NOT TARGET OpenVolumeMesh::OpenVolumeMesh)

    FetchContent_Declare(openvolumemesh
            GIT_REPOSITORY  git@gitlab.vci.rwth-aachen.de:OpenVolumeMesh/OpenVolumeMesh.git
            GIT_TAG common-laplacians
            SOURCE_DIR "${EXTERNAL_DIR}/OpenVolumeMesh"
    )
    FetchContent_MakeAvailable(openvolumemesh)
endif()

#if(NOT TARGET OpenMesh)
#    FetchContent_Declare(openmesh
#            GIT_REPOSITORY  https://gitlab.vci.rwth-aachen.de:9000/OpenMesh/OpenMesh.git
#            GIT_TAG OpenMesh-11.0
#            SOURCE_DIR "${EXTERNAL_DIR}/OpenMesh"
#            )
#    FetchContent_MakeAvailable(openmesh)
#    include_directories("${EXTERNAL_DIR}/OpenMesh/src/")
#endif()

if(NOT DEFINED TINYAD_INCLUDE_DIR)
    set(TINYAD_INCLUDE_DIR ${EXTERNAL_DIR}/TinyAD/include/)
    message(" set TINYAD_INCLUDE_DIR to ${TINYAD_INCLUDE_DIR}")
else()
    message(" TINYAD_INCLUDE_DIR already set to ${TINYAD_INCLUDE_DIR}")
endif()

if(NOT TARGET CoMISo)
    FetchContent_Declare(comiso
            GIT_REPOSITORY  git@gitlab.vci.rwth-aachen.de:CoMISo/CoMISo.git
            GIT_TAG master
            SOURCE_DIR "${EXTERNAL_DIR}/CoMISo"
            )
    FetchContent_MakeAvailable(comiso)
endif()


if(NOT TARGET TinyAD)
    FetchContent_Declare(tinyad
            GIT_REPOSITORY  https://github.com/patr-schm/TinyAD.git
            GIT_TAG main
            SOURCE_DIR "${EXTERNAL_DIR}/TinyAD"
            )
    FetchContent_MakeAvailable(tinyad)
    set(TINYAD_INCLUDE_DIR "${EXTERNAL_DIR}/TinyAD/include")
endif()




if(NOT TARGET nlohmann_json)
    FetchContent_Declare(json
            GIT_REPOSITORY  https://github.com/nlohmann/json.git
            SOURCE_DIR "${EXTERNAL_DIR}/nlohmann_json"
    )
    FetchContent_MakeAvailable(json)
endif()

if(NOT TARGET argparse)
    FetchContent_Declare(argparse
            GIT_REPOSITORY  https://github.com/p-ranav/argparse.git
            SOURCE_DIR "${EXTERNAL_DIR}/argparse"
    )
    FetchContent_MakeAvailable(argparse)
endif()

if(NOT TARGET CDT)
    FetchContent_Declare(CDT
            GIT_REPOSITORY  https://github.com/artem-ogre/CDT.git
            SOURCE_DIR "${EXTERNAL_DIR}/CDT"
    )
    FetchContent_MakeAvailable(CDT)
endif()

if(NOT TARGET CGAL_OVM_wrapper)

    FetchContent_Declare(CGAL_OVM_wrapper
            GIT_REPOSITORY git@gitlab.inf.unibe.ch:cgg/misc/cgal-ovm-wrapper.git
            GIT_TAG main
            SOURCE_DIR "${EXTERNAL_DIR}/CGAL_OVM_wrapper"
    )
    FetchContent_MakeAvailable(CGAL_OVM_wrapper)
    #set(CGAL_OVM_wrapper 1)
endif()

if(NOT TARGET TetSweep)

    FetchContent_Declare(TetSweep
            GIT_REPOSITORY git@github.com:vnigolian/TetSweep.git
            GIT_TAG main
            SOURCE_DIR "${EXTERNAL_DIR}/TetSweep"
    )
    FetchContent_MakeAvailable(TetSweep)

    include_directories("${EXTERNAL_DIR}/TetSweep/include/")
endif()




#message("Copying parameters file to build folder...")
#configure_file("SAE-parameters-brutus.json" ${CMAKE_CURRENT_BINARY_DIR}/SAE-parameters.json COPYONLY)

FILE(GLOB SOURCE_FILES
        src/*.cc
        #../OpenVolumeMesh/src/OpenVolumeMesh/Core/detail/TopologicalFaceSet.cc
        )

include_directories(src/)
#include_directories(/opt/local/include/)
include_directories(/usr/local/include/suitesparse/)
include_directories(/opt/homebrew/Cellar/suite-sparse/7.10.1/include/suitesparse/)
include_directories(/usr/local/Cellar/suite-sparse/7.10.1/include/suitesparse/)
#/opt/homebrew/Cellar/suite-sparse/7.9.0
add_subdirectory("${EXTERNAL_DIR}/CDT/CDT")
include_directories("${EXTERNAL_DIR}")



add_library (BoundaryDeformationMaps ${SOURCE_FILES})
target_link_libraries (BoundaryDeformationMaps PUBLIC
        OpenVolumeMesh::OpenVolumeMesh
        #    OpenMeshCore
        CGAL::CGAL
        Eigen3::Eigen
        nlohmann_json
        argparse
        CoMISo
        CDT
        CGAL_OVM_wrapper)


if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(" enabling sanitizers")
    add_compile_options(-fsanitize=address -fsanitize=undefined)
    add_link_options(-fsanitize=address -fsanitize=undefined)

    target_compile_options(BoundaryDeformationMaps PRIVATE -fsanitize=address -fsanitize=undefined)
    target_link_options(BoundaryDeformationMaps PRIVATE -fsanitize=address -fsanitize=undefined)
endif()

execute_process(COMMAND ln -s ${CMAKE_SOURCE_DIR}/opti-parameters.json ${CMAKE_BINARY_DIR}/opti-parameters.json)

add_subdirectory(cli)


