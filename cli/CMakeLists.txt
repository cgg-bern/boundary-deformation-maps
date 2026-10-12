

cmake_minimum_required(VERSION 3.4)


add_executable (boundary_deformation_maps boundary_deformation_maps.cc)
target_link_libraries(boundary_deformation_maps  BoundaryDeformationMaps)
