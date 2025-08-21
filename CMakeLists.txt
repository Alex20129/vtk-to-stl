CMAKE_MINIMUM_REQUIRED(VERSION 3.8)

PROJECT(vtk-to-stl)

SET(CMAKE_CXX_FLAGS "-O3 -Wall")

FIND_PACKAGE(VTK COMPONENTS CommonCore FiltersCore IOGeometry ImagingCore ImagingGeneral REQUIRED)

LINK_LIBRARIES(VTK::CommonCore VTK::FiltersCore VTK::IOGeometry VTK::ImagingCore VTK::ImagingGeneral)

vtk_module_autoinit(TARGETS ${CMAKE_PROJECT_NAME} MODULES VTK::CommonCore VTK::FiltersCore VTK::IOGeometry VTK::ImagingCore VTK::ImagingGeneral)

ADD_EXECUTABLE(${CMAKE_PROJECT_NAME} main.cpp)
