//
// Created by Valentin Nigolian on 14.12.2024.
//

#pragma once

#include "includes.hh"
#include "OptimizationElements2D.hh"
#include "OptimizationElements3D.hh"
#include "EntityPriorityQueue.hh"
#include "QualityVector.hh"
#include "Logger.hh"

//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <CGAL/Cartesian.h>

#warning TODO: make ref/domain mesh naming consistent


//designed so the higher is better
enum QUALITY_METRIC {UNKNOWN_METRIC=-1,
        DIHEDRAL_ANGLE,  // penalizes small dihedral angles more
        NEGATIVE_DIHEDRAL_ANGLE, // penalizes large dihedral angles more
        DIHEDRAL_ANGLE_DIFF, // difference between the dihedral angle what is closest between 0 and 2π. Penalizes small and large angles equally
        VOLUME_LENGTH, // volume / RMS(edge lengths)^3
        WICKE, // don't recall right now but the metric coming from the Wicke paper on elasto-plastic deformations
        CONSTANT_METRIC, // can be used to remesh without conditions
        METRIC_COUNT};

QUALITY_METRIC metric_name_to_enum(const std::string& metric_name);




/* splits are "safe" if they strictly decrease the max angle in the new mid-vertex' neighborhood */
void split_boundary_edges_opposite_to_high_angles(TetrahedralMesh& domain_mesh,
                                                  TetrahedralMesh& codomain_mesh,
                                                  double angle_threshold_deg,
                                                  int max_iterations,
                                                  bool only_allow_safe_splits);


void split_edges_to_minimize_max_boundary_angle_with_terminal_edge_propagation(TetrahedralMesh& domain_mesh,
                                                                               TetrahedralMesh& codomain_mesh,
                                                                               double angle_threshold_deg,
                                                                               int max_iterations);


OpenVolumeMesh::EdgeHandle find_terminal_boundary_edge(TetrahedralMesh& mesh,
                                                       const OpenVolumeMesh::EdgeHandle& eh);


double edge_length(const TetrahedralMesh& mesh,
                   const OpenVolumeMesh::EdgeHandle& eh);

double compute_average_incident_edge_length(const TetrahedralMesh &mesh,
                                            const OpenVolumeMesh::VertexHandle &vh,
                                            bool only_boundary);


void max_face_angle(TetrahedralMesh& mesh,
                    OpenVolumeMesh::FaceHandle fh,
                    double& max_angle,
                    OpenVolumeMesh::EdgeHandle& op_eh);

void invert_mesh(TetrahedralMesh& mesh,
                 bool move_interior_to_centroid);


/** replaces the mesh with each other, except the OUTPUT domain mesh will have a single interior vertex.
 * This assumes that the INPUT codomain is star-shaped.
 * @return: 0 -> success, 1 -> codomain is not star-shaped*/
int invert_problem(TetrahedralMesh& domain_mesh,
                    TetrahedralMesh& codomain_mesh);


/**
 * translates and rotates the mesh, so that the first vertex of the halfface is at p,
 * the halfface's normal is n, and the edge vector (v0->v1) is aligned with e1 */
int translate_and_rotate_to_target(TetrahedralMesh& mesh,
                                    const OpenVolumeMesh::HalfFaceHandle& hfh,
                                    const Vec3d& p,
                                    const Vec3d& n,
                                    const Vec3d& e1);

void translate_to_pos(TetrahedralMesh &mesh,
                      const OpenVolumeMesh::VertexHandle &anchor_vh,
                      const Vec3d& target_pos);

void scale_to_volume(TetrahedralMesh& mesh,
                    const OpenVolumeMesh::VertexHandle& anchor_vh,
                    const double& target_vol);


void scale_to_boundary_area(TetrahedralMesh &mesh,
                            const OpenVolumeMesh::VertexHandle &anchor_vh,
                            const double& target_boundary_area);


void print_boundary_quality_stats(const TetrahedralMesh& mesh);

int count_fully_boundary_elements(const TetrahedralMesh& mesh);

int split_fully_boundary_elements(TetrahedralMesh& mesh);

int split_cord_edges_and_faces(TetrahedralMesh& mesh);

bool collapse_would_create_cord(const TetrahedralMesh& mesh, const OpenVolumeMesh::EdgeHandle& eh);

/** @return the number of boundary edges whose tet fan covers more than a fulll turn (i.e. sum of dihedral angles ≥ 360° */
int check_boundary_holonomy(const TetrahedralMesh& mesh);

bool valid_holonomy(const TetrahedralMesh &mesh,
                            const OpenVolumeMesh::EdgeHandle& eh);

void check_boundary_SDE(const TetrahedralMesh& domain_mesh,
                        const TetrahedralMesh& codomain_mesh,
                        double& min,
                        double& avg,
                        double& max);

void check_full_SDE(const TetrahedralMesh& mesh,
                    const TetrahedralMesh& ref_mesh,
                    double& min_SDE,
                    double& avg_SDE,
                    double& max_SDE);

void check_full_conformal_distortion(const TetrahedralMesh& mesh,
                                     const TetrahedralMesh& ref_mesh,
                                     double& min_SDE,
                                     double& avg_SDE,
                                     double& max_SDE);

double halfface_SDE(const TetrahedralMesh& domain_mesh,
                    const TetrahedralMesh& codomain_mesh,
                    const OpenVolumeMesh::HalfFaceHandle& hfh,
                    const double ref_scaling = 1.0);


double cell_SDE(const TetrahedralMesh& mesh,
                const OpenVolumeMesh::CellHandle& ch,
                const double ref_scaling = 1.0);

double cell_SDE(const TetrahedralMesh& mesh,
                const TetrahedralMesh& ref_mesh,
                const OpenVolumeMesh::CellHandle& ch);

double cell_conformal_distortion(const TetrahedralMesh& mesh,
                                 const TetrahedralMesh& ref_mesh,
                                 const OpenVolumeMesh::CellHandle& ch);




int euler_characteristic(const TetrahedralMesh& mesh);

/** IMPORTANT NOTE: this assumes that @arg mesh is a SURFACE mesh connectivity-wise. i.e. you should delete its interior vertices first
 * Additionally, the surface is assumed to be 2-manifold (with or without boundary)
 * The reason I'm using a tet mesh here is to avoid having to deal with correspondences between the original tet mesh and a surface mesh
 *
 * returns -1 if the mesh is not a 2-manifold */
int is_surface_mesh_boundary(const TetrahedralMesh& boundary_mesh,
                                    const OpenVolumeMesh::VertexHandle& vh);

template<typename VecT>
VecT solve_Chebyshev(const std::vector<VecT>& normals,
                     const std::vector<VecT>& points,
                     double& radius,
                     bool print_debug){

    VecT new_position;
    //warning: this assumes the size is known at compile time (kind of)
    const int vec_size = new_position.size();

    using DataType = CGAL::Gmpq;


    using Program  = CGAL::Quadratic_program<DataType>;
    using Solution = CGAL::Quadratic_program_solution<DataType>;

    int max_precision(0);
    int n_constraints(0);

    CGAL::Quadratic_program_options options;
    //options.set_pricing_strategy(CGAL::QP_BLAND);     // Bland's rule
    //options.set_pricing_strategy(CGAL::QP_FILTERED_DANTZIG);


    //create solver
    // by default, we have an LP with Ax <= b
    Program lp(CGAL::SMALLER, false, 0, false, 0);

    // now set the non-default entries
    //const int X(0), Y(1), R(2);
    const int R(0);
    //R, x,y,z,...
    std::vector<int> vars(vec_size + 1);


    for(int i(0); i< (int)normals.size(); i++){

        const VecT& n = normals[i];
        const VecT& p = points[i];
        if(print_debug){
            std::cout<<" - centroid = "<<p.transpose()<<std::endl;
            std::cout<<" -   normal = "<<n.transpose()<<std::endl;
            std::cout<<" ---------"<<std::endl;
            //std::cout<<" points.push_back({"<<p[0]<<", "<<p[1]<<"});"<<std::endl;
            //std::cout<<" normals.push_back({"<<n[0]<<", "<<n[1]<<"});"<<std::endl;
        }

        DataType b = n.dot(p);

        //Ai^T . x + |Ai| . r <= Bi
        for(int i(0); i<n.size(); i++) {
            lp.set_a(i+1, n_constraints, (n[i]));
            //lp.set_a(Y, n_constraints, (n[1]));
            //std::cout<<" - A_"<<(i+1)<<" = "<<n[i]<<std::endl;
        }

        lp.set_a(R, n_constraints, (n.norm()));

        lp.set_b(n_constraints, (b));

        //std::cout<<" - a_r = "<<n.norm()<<std::endl;
        //std::cout<<" - b = "<<b<<std::endl;

        n_constraints++;
    }

    //r >= 0
    lp.set_l(R, true, 0.0);

    // Objective function: -r = (0x + 0y - 1r)
    for(int i(0); i<vec_size; i++){
        lp.set_c(i, 0);
    }
    //lp.set_c(X,  0);
    //lp.set_c(Y,  0);
    lp.set_c(R, -1);


    //solve LP

    //std::cout<<" - about to solve system..."<<std::endl;
    Solution s = CGAL::solve_linear_program(lp, DataType(), options);

    radius = 0;
    int i(0);
    for(auto val_it = s.variable_values_begin(); val_it != s.variable_values_end(); val_it++){

        DataType sol_value = (val_it->numerator()/ val_it->denominator());

        if(i == 0){
            radius = CGAL::to_double(sol_value);
        }else{
            new_position[i-1] = CGAL::to_double(sol_value);
        }
        i++;
    }

    if(print_debug){
        std::cout<<" --> infeasible: "<<s.is_infeasible()<<std::endl;
        std::cout<<" -->  unbounded: "<<s.is_unbounded()<<std::endl;
        std::cout<<" -->      valid: "<<s.is_valid()<<std::endl;
    }

    if(s.is_infeasible() || s.is_unbounded() || !s.is_valid()){
        if(print_debug){
            std::cout<<" --> infeasible: "<<s.is_infeasible()<<std::endl;
            std::cout<<" -->  unbounded: "<<s.is_unbounded()<<std::endl;
            std::cout<<" -->      valid: "<<s.is_valid()<<std::endl;
        }
        VecT zero_pos;
        for(int i(0); i<vec_size; i++){
            zero_pos[i]=0;
        }
        radius = 0.0;
        return zero_pos;
    }

    /*if(!s.solves_quadratic_program(lp)){
        std::cout<<" ---> solution does not solve LP"<<std::endl;
        std::cout<<" error: "<<s.get_error()<<std::endl;
    }*/

    if(print_debug){
        std::cout<<" radius = "<<radius<<std::endl;
    }

    return new_position;
}


Vec3d compute_visible_point(const TetrahedralMesh &mesh,
                            const OpenVolumeMesh::VertexHandle &vh,
                            const std::vector<OpenVolumeMesh::HalfFaceHandle>& incident_hfs);

/* For a given boundary edge, computes a point on the (local) outside of the mesh that is angle-wise halfway between the two incident boundary faces.
 * The norm is the average between the distance between the edge midpoint and the vertices opposite to the edge (on the incident faces) */
Vec3d compute_rotational_mid_point(const TetrahedralMesh& mesh,
                                   const OpenVolumeMesh::EdgeHandle& eh,
                                   const std::pair<OpenVolumeMesh::HalfFaceHandle, OpenVolumeMesh::HalfFaceHandle>& incident_hfs);


/* For a given boundary edge, computes the point that minimizes conformal distortion (from a regular tet)
 * for the two virtual tets incident to the two given halffaces */
Vec3d compute_optimal_virtual_vertex_pos(const TetrahedralMesh &mesh,
                                         const OpenVolumeMesh::EdgeHandle &eh,
                                         const std::vector<OpenVolumeMesh::HalfFaceHandle>& incident_hfs);


double compute_min_IVL_energy(const TetrahedralMesh& mesh,
                              const OpenVolumeMesh::CellHandle ch);

void compute_centroid_and_normal(const TetrahedralMesh& mesh,
                                  const OpenVolumeMesh::HalfFaceHandle fh,
                                  Vec3d& centroid,
                                  Vec3d& normal);

double compute_dihedral_angle_sum_deg(const TetrahedralMesh &mesh,
                                      const OpenVolumeMesh::EdgeHandle& eh);

double compute_dihedral_angle_deg(const TetrahedralMesh& mesh,
                                  const OpenVolumeMesh::CellHandle ch,
                                  const OpenVolumeMesh::EdgeHandle& eh);

double compute_boundary_dihedral_angle_deg(const TetrahedralMesh &mesh,
                                           const OpenVolumeMesh::EdgeHandle eh);

double compute_interior_angle_deg(const TetrahedralMesh& mesh,
                                  const OpenVolumeMesh::FaceHandle fh,
                                  const OpenVolumeMesh::VertexHandle vh);

double compute_volume(const TetrahedralMesh& mesh,
                      const OpenVolumeMesh::CellHandle ch);

double compute_area(const TetrahedralMesh& mesh,
                    const OpenVolumeMesh::HalfFaceHandle  fh);

double compute_total_volume(const TetrahedralMesh& mesh);

double compute_min_edge_length(const TetrahedralMesh& mesh);


double compute_total_edge_length(const TetrahedralMesh& mesh);

int count_boundary_faces(const TetrahedralMesh& mesh);

int count_boundary_edges(const TetrahedralMesh& mesh);

int count_boundary_vertices(const TetrahedralMesh& mesh);

int count_boundary_vertices(const TetrahedralMesh& mesh,
                            const OpenVolumeMesh::CellHandle& ch);

double compute_boundary_surface_area(const TetrahedralMesh& mesh);

double compute_average_boundary_surface_area(const TetrahedralMesh& mesh);

void compute_boundary_edge_lengths(const TetrahedralMesh& mesh,
                                   double& min_length,
                                   double& max_length,
                                   double& avg_length);


// current heuristic: total volume divided by the volume of a regular tet with the same edge length as
// a regular triangle of average boundary surface area
double compute_optimal_number_of_elements(const TetrahedralMesh& mesh);

Eigen::Vector3d compute_barycenter(const TetrahedralMesh& mesh);

void translate_barycenter_to_origin(TetrahedralMesh& mesh);

//NOTE: this is a quick-and-dirty version that supposed to be called AFTER translate_barycenter_to_origin
void scale_to_unit_bounding_box_diagonal(TetrahedralMesh& mesh);

void compute_bounding_box(const TetrahedralMesh &mesh,
                          Eigen::Vector3d &min_corner,
                          Eigen::Vector3d &max_corner);

Eigen::Vector3d compute_element_bounding_box_diagonal(const TetrahedralMesh& mesh,
                                                      const OpenVolumeMesh::CellHandle ch);

Eigen::Vector3d compute_global_bounding_box_diagonal(const TetrahedralMesh& mesh);


Eigen::Vector3d vec2vec(const OpenVolumeMesh::Vec3d& vec);

Eigen::Vector2d vec2vec(const OpenVolumeMesh::Vec2d &vec);

OpenVolumeMesh::Vec3d vec2vec(const Eigen::Vector3d& vec);

OpenVolumeMesh::Vec2d vec2vec(const Eigen::Vector2d& vec);

Eigen::Vector3d eigen_conversion(const Eigen::Vector2d& vec);
Eigen::Vector3d eigen_conversion(const Eigen::Vector3d& vec);



/** IMPORTANT NOTE: this assumes that the codomain_mesh boundary vertex positions are already set to the prescribed/target positions.
 * If it's NOT the case in your situation you can temporarily set them to the prescribed boundary position, since that's all that's used from the codomain mesh */
void setup_boundary_adapted_ref_positions(const TetrahedralMesh& domain_mesh,
                                          const TetrahedralMesh& codomain_mesh,
                                          const OpenVolumeMesh::CellHandle& ch,
                                          double global_scaling_factor,
                                          std::vector<Eigen::Vector3d>& P_ref);


void setup_constraints_adapted_ref_positions(const TetrahedralMesh& domain_mesh,
                                             const TetrahedralMesh& codomain_mesh,
                                             const OpenVolumeMesh::CellHandle& ch,
                                             double global_scaling_factor,
                                             const std::array<bool, 4>& constraints_mask,
                                             std::vector<Eigen::Vector3d>& P_ref);

int setup_mesh_as_one_ring_submesh(const TetrahedralMesh& mesh,
                                    const OpenVolumeMesh::EdgeHandle& he,
                                    TetrahedralMesh& submesh,
                                    OpenVolumeMesh::EdgeHandle& submesh_he);

double deg_to_rad(const double& theta_deg);
double rad_to_deg(const double& theta_rad);



/** NOTE: area variation for elements with three vertices on the boundary, and length variation for two vertices on the boundary.
 * returns 0 otherwise */
double boundary_size_variation(const TetrahedralMesh& ref_mesh,
                               const TetrahedralMesh& mesh,
                               const OpenVolumeMesh::CellHandle& ch,
                               double global_scaling_factor);


bool Tutte_embedding(TetrahedralMesh& mesh);


void export_elements_with_extreme_angles(const TetrahedralMesh& mesh,
                                         const std::string& output_path,
                                         double angle_thresh_deg);

void export_boundary_as_triangle_mesh(const TetrahedralMesh& mesh,
                                      const std::string& output_path);

void export_one_ring_neighborhood(const TetrahedralMesh& mesh,
                                  const OpenVolumeMesh::VertexHandle& vh,
                                  const std::string& output_path);


/** axis: 0->x, 1->y, 2->z */
void fill_layer_grouping_prop(TetrahedralMesh& ref_mesh,
                            TetrahedralMesh& mesh,
                            int axis,
                            int cycle_count,
                            int steps_per_cycle);


std::string edge_to_str(const TetrahedralMesh& mesh,
                            const OpenVolumeMesh::EdgeHandle& eh);

std::string halfface_to_str(const TetrahedralMesh& mesh,
                            const OpenVolumeMesh::HalfFaceHandle& hfh);

std::string cell_to_str(const TetrahedralMesh& mesh,
                        const OpenVolumeMesh::CellHandle& hfh);

void manual_break(std::string msg);
void manual_stop(std::string msg);

/**
 * Return t such that n dot ((p + t * d) - c) = 0
 * This can be used, for instance, to know how much a vertex at p can move in direction d,
 * before inverting an element whose opposite edge/face has centroid c and normal n
 * NOTE: d should be normalized
 * */
template<typename Vec>
double compute_t_factor(const Vec& c,
                        const Vec& n,
                        const Vec& p,
                        const Vec& d){

    // parallel case
    if (d.dot(n) == 0.0) {
        return std::numeric_limits<double>::infinity();
    }

    return (c-p).dot(n) / d.dot(n);
}

/** Checks whether vertex could move "around" the face or if it would cross it, if it moved in direction d */
bool face_is_blocking(const TetrahedralMesh& mesh,
                      const OpenVolumeMesh::HalfFaceHandle& hfh,
                      const OpenVolumeMesh::VertexHandle& vh,
                      const Vec3d& d);


/** If vertex vh would go AROUND halfface hfh by following direction d,
 * then it would pass through either the inside or the outside of the halfspaces defined by each edge (in the face's plane).
 * If the vertex would go through the outside of an edge, we consider it as blocking.
 * That means we can have either one, two or zero blocking edges.
 * The latter case is when the vertex would cross the face (inside all halfspaces)
 *
 * NOTE: This assumes that the vertex actually crosses the face's plane.
 * i.e. there's no check on whether the direction is parallel to the face's plane */
std::vector<OpenVolumeMesh::EdgeHandle> find_edges_blocking_movement(const TetrahedralMesh &mesh,
                                                                     const OpenVolumeMesh::HalfFaceHandle &hfh,
                                                                     const OpenVolumeMesh::VertexHandle &vh,
                                                                     const Vec3d &d,
                                                                     double epsilon = 1e-7);

//std::ostream& operator<<(std::ostream& os, const std::vector<OpenVolumeMesh::VertexHandle>& verts);

OpenVolumeMesh::EdgeHandle opposite_edge_in_cell(const TetrahedralMesh &mesh,
                                         const OpenVolumeMesh::CellHandle &ch,
                                         const OpenVolumeMesh::EdgeHandle &eh);

std::string resize_str(const std::string str, const int size);


template<unsigned long S>
std::ostream& operator<<(std::ostream& os, const std::array<OpenVolumeMesh::VertexHandle, S>& verts){
    for(auto v: verts){
        os<<" "<<v;
    }
    return os;
}