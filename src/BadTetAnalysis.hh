#pragma once

#include "Helpers.hh"
#include "NonConformingMeshesT.hh"


struct PerElementStats{
    bool domain_or_codomain_deformation = false;
    int n_boundary_vertices = 0;

    bool converged = false;
    int opti_iters = 0;
    double initial_SDE = 0;
    double final_SDE = 0;
    double ref_SDE = 0; //SDE for ref element (with ref of ref being a regular tet of the same volume)

    // the dihedral angles are not set for triangle meshes (obviously)
    double min_dihedral_angle_deg = std::numeric_limits<double>::max();;
    double avg_dihedral_angle_deg = 0;
    double max_dihedral_angle_deg = std::numeric_limits<double>::lowest();;

    double min_interior_angle_deg = std::numeric_limits<double>::max();;
    double avg_interior_angle_deg = 0;
    double max_interior_angle_deg = std::numeric_limits<double>::lowest();;

    double volume = 0;
    double relative_volume = 0;
    double bounding_box_diagonal_norm = 0;
    double boundary_size_variation = 0;
    double relative_height = 0;
    double volume_length = 0;

    double param1 = 0;
    double param2 = 0;
    double param3 = 0;
};

using TetOptiExperimentFunction = void (*)(bool, std::vector<PerElementStats>&);


void set_mesh_as_regular_tet(TetrahedralMesh& mesh);

void run_cap_tet_experiments(bool change_domain_or_codomain,
                             std::vector<PerElementStats>& stats);


void run_sliver_tet_experiments(bool change_domain_or_codomain,
                                std::vector<PerElementStats>& stats);


void run_spade_tet_experiments(bool change_domain_or_codomain,
                               std::vector<PerElementStats>& stats);

//not sure about the name
void run_triangle_tet_experiments(bool change_domain_or_codomain,
                                  std::vector<PerElementStats>& stats);

void run_sliver_tet_rotation_experiments(bool change_domain_or_codomain,
                                         std::vector<PerElementStats>& stats);

void run_flat_tet_experiments(bool change_domain_or_codomain,
                              std::vector<PerElementStats>& stats);

void run_flat_tet_exp_stretch_experiments(bool change_domain_or_codomain,
                                          std::vector<PerElementStats>& stats);

void swap_meshes(TetrahedralMesh& domain_mesh,
                 TetrahedralMesh& codomain_mesh);


void optimize_tet_and_fill_stats(const TetrahedralMesh& domain_mesh,
                                 const TetrahedralMesh& codomain_mesh,
                                 OpenVolumeMesh::CellHandle ch,
                                 bool change_domain_or_codomain,
                                 bool constrain_boundary,
                                 double scale_factor,
                                 std::vector<PerElementStats>& stats);

void optimize_tet(const TetrahedralMesh& domain_mesh,
                  const TetrahedralMesh& codomain_mesh,
                  OpenVolumeMesh::CellHandle ch,
                  bool constrain_boundary,
                  double scale_factor,
                  PerElementStats& stat);

void optimize_tet(const std::array<Vec3d, 4>& domain_positions,
                  const std::array<Vec3d, 4>& codomain_positions,
                  const std::array<bool, 4>& constraints_mask,
                  PerElementStats& stat);

void export_stats_to_csv(const std::vector<PerElementStats>& stats,
                         const std::string& csv_path);

void write_experiment_results(TetOptiExperimentFunction test_function, std::string test_name);



void analyze_all_elements(const NonConformingTetrahedralMesh &ref_mesh,
                                         const NonConformingTetrahedralMesh &mesh,
                                         double ref_to_mesh_scale_factor,
                                         bool use_ref_mesh_as_reference,
                                         double reg_tet_target_vol,
                                         VERTEX_POSITION_SOURCE pos_source,
                                         std::string csv_file_path);

