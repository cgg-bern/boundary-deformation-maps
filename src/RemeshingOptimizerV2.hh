#pragma once

#pragma clang diagnostic ignored "-W#warnings"
#pragma clang diagnostic ignored "-Woverloaded-virtual"

#include "includes.hh"
#include "EntityPriorityQueue.hh"
#include "LightWeightStopWatch.hh"
#include "Helpers.hh"
#include "Logger.hh"
//#include "NonConformingMeshesT.hh"
#include "BadTetAnalysis.hh"
#include "NonConformingHelpers.hh"

//#include "SAE_StarShapifyableExpansionCone.hh"



namespace OVM = OpenVolumeMesh;

enum WRAPPING_RESULT {
    DRY_RUN_FAILURE=-3,
        WRAPPING_NUMERICS_ERROR, //-2
    WRAPPING_ERROR, //-1
    WRAPPING_SUCCESS, //0
    WRAPPING_FAILURE, //1
    WRAPPING_ANGLE_FAILURE, //2
    WRAPPING_TIMEOUT, //3
    WRAPPING_HOLONOMY_ERROR, //4
    FINAL_ORIENTATION_CHECK_ERROR, //5
    FINAL_CONSTRAINTS_SATISFACTION_ERROR, //6
    WRAPPING_UNDEFINED_RESULT // 7
};

enum WRAPPING_STAGE {
    ANGLES_STAGE,
    EDGES_STAGE,
    POSITIONS_STAGE,
    WRAPPING_STAGE_COUNT
};

enum REMESHING_METHOD {
    GLOBAL_REMESHING,
    FOCUSED_REMESHING
};

struct EnergyRecord {

    double initial_quality_energy_ = 0;
    double final_quality_energy_ = 0;
    double final_quality_gradient_norm_ = 0;

    EnergyArray initial_virtual_energy_;
    EnergyArray final_virtual_energy_;
    EnergyArray final_virtual_gradient_norm_;

    double final_quality_gradient_max_coeff_ = 0;
    //double final_quality_H_max_coeff_ = 0;

    EnergyArray final_virtual_gradient_max_coeff_;
    //EnergyArray final_virtual_H_max_coeff_;

    double final_total_gradient_max_coeff_ = 0;
    double final_total_Hessian_max_coeff_ = 0;
};


//maybe this would work better with a json record.
struct OptimizerStats {

    int outer_iters_ = 0;
    int inner_iters_ = 0;
    int first_snap_iter_ = -1;

    int total_split_count_ = 0;
    int total_swap_count_ = 0;
    int total_collapse_count_ = 0;

    double total_opti_time_s_ = 0;
    double total_upkeep_time_s_ = 0;
    double total_remeshing_time_s_ = 0;
    double total_unaccounted_for_time_s_ = 0;

    bool backup_mode_ = false;

    std::vector<int> deep_remeshing_iterations_;

    std::vector<double> wrapped_ratio_history_;
    std::vector<double> virtual_weight_scaling_history_;
    std::vector<double> worst_quality_history_;
    std::vector<double> min_relative_volume_history_;
    std::vector<double> min_relative_edge_length_history_;
    std::vector<int> newton_iters_history_;

    std::vector<double> min_quality_before_remeshing_history_;
    std::vector<double> avg_quality_before_remeshing_history_;

    std::vector<double> min_quality_after_remeshing_history_;
    std::vector<double> avg_quality_after_remeshing_history_;

    std::vector<double> unweighted_initial_quality_energy_history_;
    std::vector<double> unweighted_initial_primary_virtual_energy_history_;
    std::vector<double> unweighted_initial_secondary_virtual_energy_history_;

    std::vector<double> unweighted_final_quality_energy_history_;
    std::vector<double> unweighted_final_primary_virtual_energy_history_;
    std::vector<double> unweighted_final_secondary_virtual_energy_history_;

    std::vector<double> unweighted_final_quality_gradient_norm_history_;
    std::vector<double> unweighted_final_primary_virtual_gradient_norm_history_;
    std::vector<double> unweighted_final_secondary_virtual_gradient_norm_history_;

    std::vector<double> final_total_gradient_max_coeff_history_;
    std::vector<double> final_total_Hessian_max_coeff_history_;

    /*static constexpr int quality_histo_min_pow = -7;
    std::vector<double> pre_deep_remeshing_min_quality_history_;
    std::vector<std::array<double, -quality_histo_min_pow + 1>> pre_deep_remeshing_quality_histogram_history_; //going from 1e-7 to 1e0
    std::vector<double> post_deep_remeshing_min_quality_history_;
    std::vector<std::array<double, -quality_histo_min_pow + 1>> post_deep_remeshing_quality_histogram_history_;*/

    std::vector<std::vector<double>> quality_history_; //quite a bit of data but it'll do for now
    std::vector<std::vector<double>> final_virtual_energies_history_; //quite a bit of data but it'll do for now

    std::vector<double> opti_time_s_history_;
    std::vector<double> max_angle_diff_deg_history_;
    std::vector<double> avg_angle_diff_deg_history_;

    double init_worst_quality_;
    double init_min_relative_volume_;
    double init_min_relative_edge_length_;

    double final_worst_quality_;
    double final_min_relative_volume_;
    double final_min_relative_edge_length_;
};


//untangling optimization
template <typename DistortionElemT, typename VirtualElemT>
class BoundaryGlueingRemeshingOptimizerT {
public:


    using EigenVec = Eigen::Vector3d;

    /// Constructor
    BoundaryGlueingRemeshingOptimizerT(NonConformingTetrahedralMesh &_ref_mesh,
                                       NonConformingTetrahedralMesh &_mesh,
                                       nlohmann::json& _json_data,
                                       const std::string& mesh_name);


    WRAPPING_RESULT run(double max_allocated_time_s);

    //this version doesn't use the virtual elements, and instead drives deformation with positional boundary constraints
    WRAPPING_RESULT run_positional(double max_allocated_time_s);

    OptimizerStats opti_stats() const;

    std::vector<std::vector<Vec3d>>& boundary_vertices_trajectories();

    std::vector<std::vector<Vec3d>>& boundary_vertices_target_paths();

private:


    //double ref_to_mesh_volume_scaling_factor() const;

    double compute_ref_to_mesh_volume_scaling_factor() const;

    int translate_and_rotate_mesh_to_codomain(const OVM::HalfFaceHandle& hfh);

    OpenVolumeMesh::HalfFaceHandle find_best_root_halfface(int ring_k) const;

    void setup_codomain_boundary_props();

    void setup_min_IVL_prop();

    void setup_codomain_average_incident_edge_length_prop();

    void setup_virtual_vertices_positions();

    void setup_virtual_vertex_position(const OVM::EdgeHandle& eh);

    int initialize_codomain();


    /** NOTE: "possible" means that
     * 1. we're not already in the snap phase
     * 2. the max boundary angle is below a certain threshold (this will probably be replaced by the virtual energy at some point)*/
    void enter_snap_phase_if_possible(bool stop_at_dihedral_angle_energy_convergence,
                                      bool force_snap_phase,
                                      double max_dihedral_angle_diff_deg,
                                      double max_dihedral_angle_diff_deg_thresh);

    int try_to_snap_next_candidates();


    void compute_element_quality(std::vector<double>& quality_vector, double& min_quality) const;

    // relative to the average
    double compute_min_relative_volume() const;

    void compute_boundary_dihedral_angle_diff(const std::vector<OVM::HalfFaceHandle>& to_wrap_hfs,
                                              const OVM::VertexHandle& to_wrap_vh,
                                              std::vector<double>& angles_diff_deg,
                                              double& max_abs_angle_diff_deg) const;


    void compute_dihedral_diff_progress(double& avg_diff,
                                        double& max_diff,
                                        double& weighted_diff_norm);

    void update_log_IVL_prop();
    void update_log_rel_vol_prop();

    double compute_boundary_edge_vector_rel_diff(double& rel_diff_norm,
                                                 double& max_rel_diff) const;

    OVM::EdgeHandle opposite_edge(const OVM::HalfFaceHandle &hfh,
                                  const OVM::VertexHandle &vh) const;

    void compute_target_angles_deg(const std::vector<OVM::HalfFaceHandle> &to_wrap_hfs,
                                   const OVM::VertexHandle &to_wrap_vh,
                                   std::vector<double> &current_angles_deg,
                                   double &max_abs_current_angle_deg,
                                   std::vector<double> &target_angles_deg,
                                   double& max_abs_target_angle_deg,
                                   int& max_abs_target_angle_deg_idx) const;



    void mark_as_wrapped(const OVM::VertexHandle& vh);


    WRAPPING_RESULT run_safety_checks(const std::string& step_description) const;

    WRAPPING_RESULT run_final_checks(double constr_eps,
                                     LightWeightStopWatch &sw) const;

    /** @return the max constraints violation */
    double snap_satisfied_boundary_constraints(double constr_eps);

    bool boundary_constraints_are_all_satisifed(double constr_eps) const;


    /** NOTE: the parameter is just for debugging
     * @return
     * -1: error
     *  0: success
     *  1: nothing was remeshed */
    int run_remeshing(REMESHING_METHOD method,
                      int opti_step,
                      int remeshing_outer_iters,
                      double candidate_quality_threshold,
                      double target_quality_threshold);

    int run_global_wrapping(int newton_iters,
                            bool& converged,
                            EnergyRecord& energy_record,
                            Eigen::VectorXd& quality_g,
                            std::array<Eigen::VectorXd, 3>& virtual_gs,
                            bool mock_run = false);

    int run_global_positional_wrapping(int newton_iters,
                                       bool &converged,
                                       EnergyRecord &energy_record,
                                       Eigen::VectorXd& quality_g,
                                       double& constraints_violation_inf_norm);

    void setup_growing_front_equality_constraints(COMISO::FiniteElementProblem& fe_problem,
                                                  const OpenVolumeMesh::VertexPropertyT<int>& global_vid_prop,
                                                  std::vector<COMISO::LinearConstraint>& eq_constraints);


    void check_for_timeout() const;

    //void export_histories() const;

    //NOTE: only useful for the hill climbing remeshing stuff, because I already have it implemented for the non-conforming mesh types.
    //However, we'll only be using the conforming positions in this method.
    NonConformingTetrahedralMesh& ref_mesh_;
    NonConformingTetrahedralMesh& mesh_;

    nlohmann::json& json_data_;

    const std::string& mesh_name_;

    const double ref_mesh_total_volume_; //domain volume
    const double target_total_volume_;//current codomain volume
    //const double ref_to_mesh_volume_scaling_factor_;

    const double ref_mesh_total_boundary_area_;
    const double mesh_total_boundary_area_;
    const double mesh_to_ref_boundary_scaling_factor_;

    LightWeightStopWatch main_sw_;

    double max_allocated_time_s_;


    OVM::VertexPropertyT<Vec3d> codomain_pos_prop_;
    OVM::EdgePropertyT<double> boundary_dihedral_angle_deg_prop_;
    OVM::EdgePropertyT<double> boundary_dihedral_angle_diff_deg_prop_;
    OVM::EdgePropertyT<double> boundary_edge_length_prop_;
    OVM::EdgePropertyT<Vec3d> virtual_vertex_pos_prop_;
    OVM::CellPropertyT<double> min_IVL_prop_;
    OVM::CellPropertyT<double> log_IVL_prop_;
    OVM::CellPropertyT<double> log_rel_vol_prop_;
    OVM::VertexPropertyT<double> average_codomain_incident_edge_length_prop_;

    double min_codomain_boundary_edge_length_;
    double avg_codomain_boundary_edge_length_;
    int initial_boundary_edges_count_;
    int initial_cell_count_;

    OVM::VertexPropertyT<bool> wrapped_prop_;

    int wrapped_count_ = 0;
    //OVM::VertexHandle final_vertex_;
    OVM::HalfFacePropertyT<int> root_hfh_prop_;
    OVM::HalfFaceHandle root_hfh_;
    OVM::VertexHandle anchor_vh_;
    bool snap_phase_ = false;

    double quality_weight_scaling_ = 1.0;
    WeightArray virtual_weight_scaling_ = std::vector<double>({1.0, 1.0, 1.0});

    //optional, to draw trajectories
    std::vector<std::vector<Vec3d>> boundary_vertices_trajectories_;

    std::vector<std::vector<Vec3d>> boundary_vertices_target_paths_;


    // --------- stats and stuff

    OptimizerStats opti_stats_;


};

#if !defined(REMESHING_OPTIMIZER_V2_C)
#define REMESHING_OPTIMIZER_V2_TEMPLATES
#include "RemeshingOptimizerV2.cc"
#endif


#if !defined(POSITIONAL_REMESHING_OPTIMIZER_C)
#define POSITIONAL_REMESHING_OPTIMIZER_TEMPLATES
#include "PositionalRemeshingOptimizer.cc"
#endif


