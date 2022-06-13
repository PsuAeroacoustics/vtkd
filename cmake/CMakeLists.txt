cmake_minimum_required(VERSION 3.14)

# Project version
project(vtk_shim)

if(${BUILD_TYPE} STREQUAL "")
	set(CMAKE_BUILD_TYPE "Debug")
elseif(${BUILD_TYPE} MATCHES "debug")
	set(CMAKE_BUILD_TYPE "Debug")
elseif(${BUILD_TYPE} MATCHES "release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

find_package(VTK COMPONENTS vtkCommonCore vtkCommonDataModel vtkIOXML)
add_library(vtk_shim STATIC ../source/vtkd/vtk.cpp)

set_property(TARGET vtk_shim PROPERTY CXX_STANDARD 20)

target_link_libraries(vtk_shim VTK)

if(${VTK_MAJOR_VERSION} GREATER_EQUAL "9")
	target_include_directories(vtk_shim PRIVATE "${VTK_PREFIX_PATH}/include/vtk-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
else()
	message(NOTICE "Old way")
	target_include_directories(vtk_shim PRIVATE ${VTK_INCLUDE_DIRS})
endif()

target_compile_options(vtk_shim PRIVATE -fPIE)
