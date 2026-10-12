#pragma once


#include <cstdio>
#include <iostream>
#include <queue>
#include <random>
#include <stack>
 	

#include <OpenVolumeMesh/Core/Properties/PropertyPtr.hh>
#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>
#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <OpenVolumeMesh/Core/detail/TopologicalLinkT_impl.hh>


//#include <OpenVolumeMesh/Core/detail/TopologicalLink.hh>

//#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
//#include <OpenMesh/Core/IO/MeshIO.hh>

#include <Eigen/Dense>

//#include <CGAL/Tetrahedron_3.h>

#include <CoMISo/Config/config.hh>
#include <CoMISo/Utils/StopWatch.hh>
#include <CoMISo/NSolver/NPDerivativeChecker.hh>
#include <CoMISo/NSolver/NPTiming.hh>
#include <CoMISo/NSolver/NewtonSolver.hh>
#include <CoMISo/NSolver/TruncatedNewtonPCG.hh>
#include <CoMISo/NSolver/ConstraintTools.hh>
#include <CoMISo/NSolver/LinearConstraint.hh>
#include <CoMISo/NSolver/FiniteElementProblem.hh>

#include "CDT.h"

#include <CGAL_OVM_wrapper/CGAL_OVM_wrapper.hh>

#include <argparse/argparse.hpp>

#include <nlohmann/json.hpp>


typedef OpenVolumeMesh::Vec2d Vec2d;
typedef OpenVolumeMesh::Vec3d Vec3d;

typedef OpenVolumeMesh::TetrahedralGeometryKernel<OpenVolumeMesh::Vec3d, OpenVolumeMesh::TetrahedralMeshTopologyKernel> TetrahedralMesh;


struct TimeOutException : public std::exception{
    const char * what () const throw ()
    {
        return "timeout time reached";
    }
};


struct VirtualVertexSetupException : public std::exception{
    const char * what () const throw ()
    {
        return "failed to setup virtual vertex";
    }
};


