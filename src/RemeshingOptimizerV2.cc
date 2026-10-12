
#define REMESHING_OPTIMIZER_V2_C

#include "RemeshingOptimizerV2.hh"

using namespace OpenVolumeMesh;


template <typename DistortionElemT, typename VirtualElemT>
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
BoundaryGlueingRemeshingOptimizerT(NonConformingTetrahedralMesh &_ref_mesh,
                                   NonConformingTetrahedralMesh &_mesh,
                                   nlohmann::json &_json_data,
                                   const std::string& _mesh_name) :
ref_mesh_(_ref_mesh),
mesh_(_mesh),
json_data_(_json_data),
mesh_name_(_mesh_name),
ref_mesh_total_volume_(compute_total_volume(ref_mesh_)),
target_total_volume_(compute_total_volume(mesh_)),
ref_mesh_total_boundary_area_(compute_boundary_surface_area(ref_mesh_)),
mesh_total_boundary_area_(compute_boundary_surface_area(mesh_)),
mesh_to_ref_boundary_scaling_factor_(std::sqrt(ref_mesh_total_boundary_area_ / mesh_total_boundary_area_)),
wrapped_prop_(mesh_.request_vertex_property<bool>("wrapped_vertex", false)),
codomain_pos_prop_(mesh_.request_vertex_property<Vec3d>("codomain_pos", {0,0,0})),
boundary_dihedral_angle_deg_prop_(mesh_.request_edge_property<double>("boundary_dihedral_angle_deg", 0.0)),
boundary_dihedral_angle_diff_deg_prop_(mesh_.request_edge_property<double>("boundary_dihedral_angle_diff_deg", 0.0)),
boundary_edge_length_prop_(mesh_.request_edge_property<double>("boundary_edge_length_prop", 0.0)),
virtual_vertex_pos_prop_(mesh_.request_edge_property<Vec3d>("virtual_vertex_pos", {0,0,0})),
min_IVL_prop_(mesh_.request_cell_property<double>("min_IVL", std::numeric_limits<double>::max())),
log_IVL_prop_(mesh_.request_cell_property<double>("log_IVL", std::numeric_limits<double>::max())),
log_rel_vol_prop_(mesh_.request_cell_property<double>("log_rel_vol", std::numeric_limits<double>::max())),
average_codomain_incident_edge_length_prop_(mesh_.request_vertex_property<double>("average_incident_codomain_edge_length", std::numeric_limits<double>::max())),
root_hfh_prop_(mesh_.request_halfface_property<int>("root_halfface", 0)){

#if 0
    if (const char* env_p = std::getenv("OMP_NUM_THREADS")){
        int omp_num_threads_env = std::stoi(env_p);
        if(omp_num_threads_env != 1){
            log(LogLevel::ERROR)<<" ERROR - OMP_NUM_THREADS is set to "<<omp_num_threads_env<<" instead of 1"<<std::endl;
            exit(EXIT_FAILURE);
        }
    }else{
        log(LogLevel::ERROR)<<" ERROR - OMP_NUM_THREADS environment variable not set. It should be set to 1"<<std::endl;
        exit(EXIT_FAILURE);
    }
#endif

    log(LogLevel::INFO)<<" -           target mesh volume = "<<target_total_volume_<<std::endl;
    log(LogLevel::INFO)<<" -              ref mesh volume = "<<ref_mesh_total_volume_<<std::endl;
    //log(LogLevel::INFO)<<" -   ref to mesh volume scaling = "<<ref_to_mesh_scaling_factor_<<std::endl;

    mesh_.set_persistent(codomain_pos_prop_);
    mesh_.set_persistent(wrapped_prop_);
    mesh_.set_persistent(boundary_dihedral_angle_deg_prop_);
    mesh_.set_persistent(boundary_dihedral_angle_diff_deg_prop_);
    mesh_.set_persistent(virtual_vertex_pos_prop_);
    mesh_.set_persistent(min_IVL_prop_);
    mesh_.set_persistent(average_codomain_incident_edge_length_prop_);
    mesh_.set_persistent(root_hfh_prop_);
    mesh_.set_persistent(log_IVL_prop_);
    mesh_.set_persistent(log_rel_vol_prop_);

    for (auto vh: mesh_.vertices()) {
        wrapped_prop_[vh] = false;
    }

    //not using a prop because it's optional
    if(json_data_["enable_vertex_trajectory_mesh"]){
        boundary_vertices_trajectories_.resize(mesh_.n_vertices());
        boundary_vertices_target_paths_.resize(mesh_.n_vertices());
    }
}


template <typename DistortionElemT, typename VirtualElemT>
int
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
initialize_codomain(){
    log(LogLevel::IMPORTANT) << " -> initializing codomain mesh..." << std::endl;


    bool warm_start(false);
    for (auto vh: mesh_.vertices()) {
        if (wrapped_prop_[vh]) {
            /*log(LogLevel::IMPORTANT) << " -> found at least one wrapped vertex "<<vh<<", skipping initialization" << std::endl;
            warm_start = true;
            break;*/

            log(LogLevel::ERROR)<<" ERROR - wrapping warm start is not supported anymore"<<std::endl;
            return -1;
        }
    }

    if (true || !warm_start) {

        setup_codomain_boundary_props();
        setup_min_IVL_prop();
        setup_codomain_average_incident_edge_length_prop();

        double max_codomain_boundary_edge_length(0);
        compute_boundary_edge_lengths(mesh_, min_codomain_boundary_edge_length_, max_codomain_boundary_edge_length, avg_codomain_boundary_edge_length_);


        // save the current position and replace the codomain mesh with the domain one
        anchor_vh_ = VertexHandle(-1);
        for (auto vh: mesh_.vertices()) {
            if (mesh_.is_boundary(vh)) {
                codomain_pos_prop_[vh] = mesh_.vertex(vh);
                if (!anchor_vh_.is_valid()) {
                    anchor_vh_ = vh;
                }
            } else {
                codomain_pos_prop_[vh] = {0, 0, 0};
            }
            mesh_.set_vertex(vh, ref_mesh_.vertex(vh));
        }

        translate_to_pos(mesh_, anchor_vh_, codomain_pos_prop_[anchor_vh_]);
        //scale_to_volume(mesh_, anchor_vh_, mesh_total_volume_);
        scale_to_boundary_area(mesh_, anchor_vh_, mesh_total_boundary_area_);
        log(LogLevel::IMPORTANT)<<" ->              target volume: "<<target_total_volume_ <<std::endl;
        log(LogLevel::IMPORTANT)<<" ->        updated mesh volume: "<<compute_total_volume(mesh_) <<std::endl;
        log(LogLevel::IMPORTANT)<<" ->  target mesh boundary area: "<<mesh_total_boundary_area_ <<std::endl;
        log(LogLevel::IMPORTANT)<<" -> updated mesh boundary area: "<<compute_boundary_surface_area(mesh_) <<std::endl;


        //kind of a hack, so the anchor vertex is constrained
        wrapped_prop_[anchor_vh_] = true;

        auto updated_volume = compute_total_volume(mesh_);
        auto updated_boundary_area = compute_boundary_surface_area(mesh_);

        auto volume_rel_diff = std::abs(target_total_volume_ - updated_volume) / target_total_volume_;
        auto boundary_area_rel_diff = std::abs(mesh_total_boundary_area_ - updated_boundary_area) / mesh_total_boundary_area_;

        if (volume_rel_diff > 1e-7) {
            log(LogLevel::WARNING) << " WARNING - volume after scaling is " << updated_volume << " but should be " << target_total_volume_ << " (rel diff = " << volume_rel_diff << ")"
                                 << std::endl;
        }
        if (boundary_area_rel_diff > 1e-7) {
            log(LogLevel::ERROR) << " ERROR - boundary area after scaling is " << updated_boundary_area << " but should be " << mesh_total_boundary_area_ << " (rel diff = " << boundary_area_rel_diff << ")"
                                 << std::endl;
            return -1;
        }

        //do that AFTER initializing the work mesh
        setup_virtual_vertices_positions();

        initial_boundary_edges_count_ = count_boundary_edges(mesh_);

    }else{
        //root_hfh_ = root_hfh_prop_[0];
        log(LogLevel::WARNING)<<" WARNING - warm start isn't implemented yet for this version, re-initializing the codomain"<<std::endl;
        //return WRAPPING_ERROR;
    }

    return 0;
}


template <typename DistortionElemT, typename VirtualElemT>
void
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
enter_snap_phase_if_possible(bool stop_at_dihedral_angle_energy_convergence,
                             bool force_snap_phase,
                             double max_dihedral_angle_diff_deg,
                             double max_dihedral_angle_diff_deg_thresh) {

    if (stop_at_dihedral_angle_energy_convergence) {
        log(LogLevel::WARNING) << " WARNING - skipping snap phase for first stage experiment" << std::endl;
        return;
    }

    if (!snap_phase_) {

        log(LogLevel::IMPORTANT) << " ----------------------------------------------------------------------------" << std::endl;
        log(LogLevel::IMPORTANT) << " - trying to enter snap phase..." << std::endl;

        if (force_snap_phase || max_dihedral_angle_diff_deg < max_dihedral_angle_diff_deg_thresh) {
            log(LogLevel::IMPORTANT) << " --> (unweighted) virtual energy below threshold, or forcing snap phase. Moving on to snap phase" <<
                    std::endl;

            snap_phase_ = true;
            opti_stats_.first_snap_iter_ = opti_stats_.outer_iters_;

            //un-wrap the anchor because it might not be compatible with the root face
            wrapped_prop_[anchor_vh_] = false;
            wrapped_count_ = 0;


            if (json_data_["enable_debugging_meshes_exportation"]) {
                std::string ext = json_data_["debugging_meshes_exportation_extension"];
                mesh_.export_conforming_mesh(mesh_name_+"_codomain_mesh_after_before_rotation" + ext);
            }

            root_hfh_ = find_best_root_halfface(json_data_["global_wrapping_root_hfh_ring_k"]);
            root_hfh_prop_[root_hfh_] = 1;
            log(LogLevel::IMPORTANT) << " -> using halfface " << halfface_to_str(mesh_, root_hfh_) << " as root" <<
                    std::endl;

            //scale_to_volume(mesh_, anchor_vh_, mesh_total_volume_);
            scale_to_boundary_area(mesh_, anchor_vh_, mesh_total_boundary_area_);

            log(LogLevel::IMPORTANT)<<" -       target        volume: "<<target_total_volume_ <<std::endl;
            log(LogLevel::IMPORTANT)<<" - updated mesh        volume: "<<compute_total_volume(mesh_) <<std::endl;
            log(LogLevel::IMPORTANT)<<" -  target mesh boundary area: "<<mesh_total_boundary_area_ <<std::endl;
            log(LogLevel::IMPORTANT)<<" - updated mesh boundary area: "<<compute_boundary_surface_area(mesh_) <<std::endl;

            translate_and_rotate_mesh_to_codomain(root_hfh_);

            //re-compute the virtual vertices after the rigid transformation
            setup_virtual_vertices_positions();

            auto root_verts = mesh_.get_halfface_vertices(root_hfh_);
            wrapped_prop_[root_verts[0]] = true;
            wrapped_count_ = 1;

            if (json_data_["enable_debugging_meshes_exportation"]) {
                std::string ext = json_data_["debugging_meshes_exportation_extension"];
                mesh_.export_conforming_mesh(mesh_name_+"_codomain_mesh_after_before_snap_phase" + ext);
            }
        } else {
            log(LogLevel::IMPORTANT) << " --> (unweighted) virtual energy is still above threshold" << std::endl;
        }
        log(LogLevel::IMPORTANT) << " ----------------------------------------------------------------------------" << std::endl;
    }
}


template <typename DistortionElemT, typename VirtualElemT>
WRAPPING_RESULT
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
run(double max_allocated_time_s) {


    main_sw_.reset();
    max_allocated_time_s_ = max_allocated_time_s;

    WRAPPING_RESULT wrapping_result(WRAPPING_UNDEFINED_RESULT);

    try{

        log(LogLevel::IMPORTANT)<<" ================================================================================================================================ "<<std::endl;
        log(LogLevel::IMPORTANT)<<" ================================================================================================================================ "<<std::endl;
        log(LogLevel::IMPORTANT)<<" -> running global wrapping process..."<<std::endl;
        log(LogLevel::IMPORTANT)<<" ================================================================================================================================ "<<std::endl;


        snap_phase_ = false;

        const bool enable_focused_remeshing = json_data_["global_wrapping_enable_focused_remeshing"];
        const int deep_remeshing_inner_iters = json_data_["global_wrapping_deep_remeshing_inner_iters"];
        const double quick_remeshing_candidate_quality_threshold = json_data_["global_wrapping_quick_remeshing_candidate_quality_threshold"];
        const double quick_remeshing_target_quality_threshold = json_data_["global_wrapping_quick_remeshing_target_quality_threshold"];
        const double deep_remeshing_candidate_quality_threshold = json_data_["global_wrapping_deep_remeshing_candidate_quality_threshold"];
        const double deep_remeshing_target_quality_threshold = json_data_["global_wrapping_deep_remeshing_target_quality_threshold"];


        const int max_outer_iters = json_data_["global_wrapping_max_outer_iters"];
        const int max_inner_iters = json_data_["global_wrapping_max_inner_iters"];
        const int max_total_inner_iters = json_data_["global_wrapping_max_total_inner_iters"];
        const int newton_iters = json_data_["global_wrapping_max_newton_iters"];

        const double final_max_H_coeff_thresh = json_data_["global_wrapping_max_abs_hessian_coeff_thresh"];

        const double virtual_w_scale_base_variation = json_data_["global_wrapping_virtual_elements_w_variation"];

        const double virtual_energy_inner_progress_thresh = json_data_["global_wrapping_virtual_energy_inner_progress_thresh"];
        const double virtual_energy_outer_progress_thresh = json_data_["global_wrapping_virtual_energy_outer_progress_thresh"];
        const double unweighted_virtual_energy_success_thresh = json_data_["global_wrapping_virtual_energy_success_thresh"];
        const double max_dihedral_diff_deg_thresh = json_data_["global_wrapping_max_angle_diff_snap_thresh"];

        const double vol_bound_w = json_data_["global_wrapping_IVL_vol_w"];
        const double min_rel_vol_bound = json_data_["global_wrapping_IVL_rel_vol_bound"];

        const bool stop_once_below_dihedral_angle_diff_threshold = json_data_["global_wrapping_stop_once_below_dihedral_angle_diff_threshold"];
        const bool stop_at_virtual_energy_convergence = json_data_["global_wrapping_stop_at_virtual_energy_convergence"];
        const bool enable_full_history = json_data_["global_wrapping_enable_full_history"];

        //const double start_virtual_w_min_scaling = json_data_["global_wrapping_virtual_energy_start_w_min_scaling"];
        const double start_virtual_w = json_data_["global_wrapping_virtual_energy_start_w"];
        const double start_quality_w = json_data_["global_wrapping_quality_energy_start_w"];
        const double max_start_virtual_w = json_data_["global_wrapping_virtual_energy_max_start_w"];
        const double max_virtual_w = json_data_["global_wrapping_virtual_energy_max_w"];
        const bool enable_adaptive_start_w = json_data_["global_wrapping_enable_adaptive_start_w"];
        const bool enable_relaxation_pass = json_data_["global_wrapping_enable_relaxation_pass"];
        const bool enable_fixed_boundary_relaxation_pass = json_data_["global_wrapping_enable_fixed_boundary_relaxation_pass"];
        const bool enable_backup_DA_energy_mode = json_data_["global_wrapping_enable_backup_DA_energy_mode"];

        //last (4th) entry is for the quality energy
        const std::vector<double> adaptive_start_w_target_g_norm = std::vector<double>(json_data_["global_wrapping_adaptive_start_w_target_g_norm"]);

        const double constr_eps = json_data_["constr_eps"];

        const bool enable_forced_snapping = json_data_["global_wrapping_enable_forced_snapping"];

        int initial_boundary_vertices_count = count_boundary_vertices(mesh_);


        if(initialize_codomain()){
            log(LogLevel::ERROR)<<" ERROR - couldn't initialize codomain"<<std::endl;
            return WRAPPING_ERROR;
        }
        if (auto res = run_safety_checks("codomain initialization")) {
            return res;
        }


        IO::FileManager fm;
        mesh_.export_conforming_mesh(mesh_name_ + "_codomain_mesh_init.ovm");
        export_boundary_as_triangle_mesh(mesh_, mesh_name_+"_codomain_mesh_boundary_init.obj");
        export_boundary_as_triangle_mesh(ref_mesh_, mesh_name_+"_target_boundary.obj");


        log(LogLevel::IMPORTANT) << " ---------------------------------------------------------------------- " << std::endl;
        log(LogLevel::IMPORTANT) << "  -> initial boundary dihedral angles: " << std::endl;

        double avg_abs_dihedral_diff_deg, max_abs_dihedral_diff_deg, dihedral_diff_weighted_norm;
        compute_dihedral_diff_progress(avg_abs_dihedral_diff_deg,
                                       max_abs_dihedral_diff_deg,
                                       dihedral_diff_weighted_norm);

        double min_quality;
        std::vector<double> quality_vector;
        compute_element_quality(quality_vector, min_quality);

        double min_rel_vol = compute_min_relative_volume();
        double min_rel_edge_len = compute_min_edge_length(mesh_) / min_codomain_boundary_edge_length_;

        opti_stats_.wrapped_ratio_history_.push_back((double)wrapped_count_ / (double)initial_boundary_vertices_count);
        opti_stats_.max_angle_diff_deg_history_.push_back(max_abs_dihedral_diff_deg);
        opti_stats_.avg_angle_diff_deg_history_.push_back(avg_abs_dihedral_diff_deg);
        opti_stats_.worst_quality_history_.push_back(min_quality);
        opti_stats_.min_relative_volume_history_.push_back(min_rel_vol);
        opti_stats_.min_relative_edge_length_history_.push_back(min_rel_edge_len);
        opti_stats_.init_worst_quality_ = min_quality;
        opti_stats_.init_min_relative_volume_ = min_rel_vol;
        opti_stats_.init_min_relative_edge_length_ = min_rel_edge_len;
        if (enable_full_history) {
            opti_stats_.quality_history_.push_back(quality_vector);
        }



        log(LogLevel::IMPORTANT) << "  avg abs dihedral angle difference =  " <<avg_abs_dihedral_diff_deg<<"°"<< std::endl;
        log(LogLevel::IMPORTANT) << "  max abs dihedral angle difference =  " <<max_abs_dihedral_diff_deg<<"°"<< std::endl;
        log(LogLevel::IMPORTANT) << "     edge-weighted diff vector norm =  " <<dihedral_diff_weighted_norm<<"°"<< std::endl;
        log(LogLevel::IMPORTANT) << "      initial worst element quality =  " <<min_quality<< std::endl;
        log(LogLevel::IMPORTANT) << "        initial min relative volume =  " <<min_rel_vol<< std::endl;
        log(LogLevel::IMPORTANT) << "   initial min relative edge length =  " <<min_rel_edge_len<< std::endl;
        log(LogLevel::IMPORTANT) << " ---------------------------------------------------------------------- " << std::endl;

        if (enable_fixed_boundary_relaxation_pass && vol_bound_w > 0 && min_rel_vol < min_rel_vol_bound) {
            log(LogLevel::IMPORTANT) << " min relative volume is below threshold, running a fixed-boundary relaxation pass " << std::endl;

            if(optimize_mesh<NonConformingTetrahedralMesh, DistortionElemT>(mesh_, ref_mesh_total_volume_/ mesh_.n_logical_cells(), 1000, true, nullptr, false, json_data_, LogConfig::get_min_level() > LogLevel::INFO)){
                log(LogLevel::ERROR)<<" error while optimizing reference mesh"<<std::endl;
                return WRAPPING_ERROR;
            }
        }


        if(json_data_["enable_vertex_trajectory_mesh"]){
            for (auto vh: mesh_.vertices()) {
                if(mesh_.is_boundary(vh)){
                    boundary_vertices_trajectories_[vh.idx()].push_back(mesh_.vertex(vh));
                }
            }
        }


        bool converged(false);
        EnergyRecord energy_record;

        Eigen::VectorXd quality_g;
        std::array<Eigen::VectorXd, 3> virtual_gs;

        double last_dihedral_diff_weighted_norm = dihedral_diff_weighted_norm;

        quality_weight_scaling_ = start_quality_w;
        virtual_weight_scaling_ = start_virtual_w;

        opti_stats_.backup_mode_ = false;

        EnergyArray last_pass_unweighted_initial_virtual_e(std::numeric_limits<double>::max());
        EnergyArray last_pass_unweighted_final_virtual_e(std::numeric_limits<double>::max());

        int& outer_wrap_it = opti_stats_.outer_iters_;
        while(outer_wrap_it < max_outer_iters && wrapped_count_ != initial_boundary_vertices_count) {

            LightWeightStopWatch pass_sw;


            log(LogLevel::IMPORTANT) << " ================================================================================================================================ " << std::endl;
            log(LogLevel::IMPORTANT) << " ====================================================================================== running global wrap pass " << outer_wrap_it  << std::endl;
            log(LogLevel::IMPORTANT) << " ====================================================================================== last unweighted energies: " << last_pass_unweighted_final_virtual_e << std::endl;
            log(LogLevel::IMPORTANT) << " ======================================================================================            wrapped count: " << wrapped_count_ << " / " << initial_boundary_vertices_count << std::endl;
            log(LogLevel::IMPORTANT) << " ======================================================================================       total elapsed time: " << main_sw_.total_duration() << "s" << std::endl;
            log(LogLevel::IMPORTANT) << " ================================================================================================================================ " << std::endl;



            if(enable_adaptive_start_w){

                if (!enable_relaxation_pass || outer_wrap_it) {
                    log(LogLevel::IMPORTANT) << " ---------------------------------------------------------------------- " << std::endl;
                    log(LogLevel::IMPORTANT) << " - computing optimal starting weight from gradients... " << std::endl;

                    quality_weight_scaling_ = start_quality_w;
                    virtual_weight_scaling_ = start_virtual_w;

                    if (!enable_backup_DA_energy_mode || (enable_backup_DA_energy_mode && !opti_stats_.backup_mode_)) {
                        log(LogLevel::IMPORTANT) << " --> still not in backup mode, setting virtual weights for secondary and tertiary to zero " << std::endl;
                        virtual_weight_scaling_.secondary = 0;
                        virtual_weight_scaling_.tertiary = 0;
                    }


                    if (run_global_wrapping(1, converged, energy_record, quality_g, virtual_gs, true)) {
                        log(LogLevel::ERROR) << " error while running global wrap opti" << std::endl;
                        return WRAPPING_ERROR;
                    }

                    if (!enable_relaxation_pass && !outer_wrap_it) {
                        quality_weight_scaling_ = adaptive_start_w_target_g_norm[3] / quality_g.norm();
                        log(LogLevel::IMPORTANT) << " -         quality g norm: " << quality_g.norm() << std::endl;
                        log(LogLevel::IMPORTANT) << " -          target g norm: " << adaptive_start_w_target_g_norm[3] << std::endl;
                        log(LogLevel::IMPORTANT) << " -   quality start weight: " << quality_weight_scaling_ << std::endl;
                        log(LogLevel::IMPORTANT) << " --------------------------------- " << std::endl;
                    }

                    EnergyArray target_g_norm = adaptive_start_w_target_g_norm;

                    EnergyArray virtual_g_norm = {virtual_gs[0].norm(), virtual_gs[1].norm(), virtual_gs[2].norm()};
                    virtual_weight_scaling_ = target_g_norm / virtual_g_norm;
                    virtual_weight_scaling_.make_finite();

                    log(LogLevel::IMPORTANT) << " -         virtual g norm: "<<virtual_g_norm << std::endl;
                    log(LogLevel::IMPORTANT) << " -          target g norm: "<< target_g_norm<< std::endl;
                    log(LogLevel::IMPORTANT) << " -   virtual start weight: "<<virtual_weight_scaling_<< std::endl;
                    log(LogLevel::IMPORTANT) << " ---------------------------------------------------------------------- " << std::endl;
                }else {

                    log(LogLevel::IMPORTANT) << " ---------------------------------------------------------------------- " << std::endl;
                    log(LogLevel::IMPORTANT) << " - running relaxation pass, setting virtual weights to zero " << std::endl;

                    virtual_weight_scaling_ = 0.0;

                    log(LogLevel::IMPORTANT) << " ---------------------------------------------------------------------- " << std::endl;
                }
            }else{
                //a bit brute force, but it'll do for now

                log(LogLevel::IMPORTANT) << " ---------------------------------------------------------------------- " << std::endl;
                log(LogLevel::IMPORTANT) << " - using default starting weights "<< virtual_weight_scaling_<< std::endl;

                log(LogLevel::IMPORTANT) << " ---------------------------------------------------------------------- " << std::endl;
            }


            //only do this after the first pass
            while (outer_wrap_it) {

                log(LogLevel::IMPORTANT) << " ---------------------------------------------------------------------- " << std::endl;
                log(LogLevel::IMPORTANT) << " - checking virtual energy variation with weights "<<virtual_weight_scaling_ << std::endl;
                if (run_global_wrapping(3, converged, energy_record, quality_g, virtual_gs, true)) {
                    log(LogLevel::ERROR) << " error while running global wrap opti" << std::endl;
                    return WRAPPING_ERROR;
                }

                auto  virtual_energy_progress = (energy_record.initial_virtual_energy_ - energy_record.final_virtual_energy_) / energy_record.initial_virtual_energy_;

                log(LogLevel::IMPORTANT) << " -> initial / final primary virtual energy: "<< energy_record.initial_virtual_energy_<<" / "<<energy_record.final_virtual_energy_<< std::endl;
                log(LogLevel::IMPORTANT) << " ->                      relative progress: "<< virtual_energy_progress<< std::endl;


                if (energy_record.final_virtual_energy_.primary < energy_record.initial_virtual_energy_.primary ||
                    energy_record.final_virtual_energy_.secondary < energy_record.initial_virtual_energy_.secondary ||
                    energy_record.final_virtual_energy_.tertiary < energy_record.initial_virtual_energy_.tertiary ) {

                    log(LogLevel::IMPORTANT) << " ---> either of the virtual energies decreased -> valid starting virtual weights: "<< virtual_weight_scaling_<< std::endl;
                    last_pass_unweighted_initial_virtual_e = energy_record.initial_virtual_energy_ / virtual_weight_scaling_;
                    break;
                }
                log(LogLevel::IMPORTANT) << " ---> energy did not decrease, increasing start virtual weights "<< std::endl;

                //initially only increasing the weight for the primary virtual elements until it becomes dominant
                virtual_weight_scaling_ *= virtual_w_scale_base_variation;
                if (virtual_weight_scaling_.primary > max_start_virtual_w ||
                    virtual_weight_scaling_.secondary > max_start_virtual_w ||
                    virtual_weight_scaling_.tertiary > max_start_virtual_w) {

                    log(LogLevel::IMPORTANT) << " ---> either of the virtual weights are higher than the max value ("<< max_start_virtual_w<<"), stopping here (they'll be capped later) "<< std::endl;
                    break;
                }
                log(LogLevel::IMPORTANT) << " ---------------------------------------------------------------------- " << std::endl;
            }

            virtual_weight_scaling_ = {
                std::min(virtual_weight_scaling_.primary, max_start_virtual_w),
                std::min(virtual_weight_scaling_.secondary, max_start_virtual_w),
                std::min(virtual_weight_scaling_.tertiary, max_start_virtual_w)
            };
            log(LogLevel::IMPORTANT) << " ---> start virtual weights: "<<virtual_weight_scaling_<< std::endl;


            EnergyArray unweighted_final_virtual_e(std::numeric_limits<double>::max());

            bool made_progress_this_pass(false);
            int inner_wrap_it(0);
            while (inner_wrap_it < max_inner_iters && wrapped_count_ != initial_boundary_vertices_count) {

                check_for_timeout();
                LightWeightStopWatch it_sw;

                opti_stats_.virtual_weight_scaling_history_.push_back(virtual_weight_scaling_.primary);

                log(LogLevel::IMPORTANT) << " ================================================================================================================================ " << std::endl;
                log(LogLevel::IMPORTANT) << " ============================================================================   running global wrap iteration " << outer_wrap_it << "." << inner_wrap_it <<" (overall idx:"<<opti_stats_.inner_iters_<<")"<< std::endl;
                log(LogLevel::IMPORTANT) << " ============================================================================ last unweighted virtual energy: " << unweighted_final_virtual_e << std::endl;
                log(LogLevel::IMPORTANT) << " ============================================================================        current virtual weights: " << virtual_weight_scaling_ << std::endl;
                log(LogLevel::IMPORTANT) << " ============================================================================         current quality weight: " << quality_weight_scaling_ << std::endl;
                log(LogLevel::IMPORTANT) << " ============================================================================                  wrapped count: " << wrapped_count_ << " / " << initial_boundary_vertices_count << std::endl;
                log(LogLevel::IMPORTANT) << " ============================================================================             total elapsed time: " << main_sw_.total_duration() << "s" << std::endl;
                log(LogLevel::IMPORTANT) << " ================================================================================================================================ " << std::endl;


                if (run_global_wrapping(newton_iters, converged, energy_record, quality_g, virtual_gs, false)) {
                    log(LogLevel::ERROR) << " error while running global wrap opti" << std::endl;
                    return WRAPPING_ERROR;
                }

                if (json_data_["enable_debugging_meshes_exportation"]) {
                    update_log_IVL_prop();
                    update_log_rel_vol_prop();
                    std::string ext = json_data_["debugging_meshes_exportation_extension"];

                    auto temp_prop = mesh_.request_vertex_property<Vec3d>();
                    //temp
                    {
                        for(auto vh: mesh_.vertices()){
                            temp_prop[vh] = codomain_pos_prop_[vh];
                            if(mesh_.is_boundary(vh)){
                                codomain_pos_prop_[vh] -= mesh_.vertex(vh);
                            }else{
                                codomain_pos_prop_[vh] = {0,0,0};
                            }
                        }
                    }


                    mesh_.export_conforming_mesh(mesh_name_+"_codomain_mesh_after_it_" + std::to_string(outer_wrap_it) + "_" + std::to_string(inner_wrap_it) + ext, false);
                    export_boundary_as_triangle_mesh(mesh_, mesh_name_+"_codomain_mesh_boundary_after_it_" + std::to_string(outer_wrap_it) + "_" + std::to_string(inner_wrap_it) + ".obj");

                    {
                        for(auto vh: mesh_.vertices()){
                            codomain_pos_prop_[vh] = temp_prop[vh];
                        }
                    }
                }

                if (json_data_["global_wrapping_enable_analysis"]) {

                    log(LogLevel::INFO) << " - running post-optimization analysis" << std::endl;

                    analyze_all_elements(ref_mesh_, mesh_, compute_ref_to_mesh_volume_scaling_factor(), true, 1.0, CONFORMING, "");
                    analyze_all_elements(ref_mesh_, ref_mesh_, 1.0, false, ref_mesh_total_volume_ / ref_mesh_.n_logical_cells(), CONFORMING, "");
                }


                if (auto res = run_safety_checks("wrap_it_" + std::to_string(outer_wrap_it) + "_" + std::to_string(inner_wrap_it))) {
                    log(LogLevel::ERROR)<<" safety checks failed, stopping here"<<std::endl;
                    return res;
                }

                log(LogLevel::IMPORTANT) << " ---------------------------------------------------------------------- " << std::endl;
                compute_dihedral_diff_progress(avg_abs_dihedral_diff_deg,
                                               max_abs_dihedral_diff_deg,
                                               dihedral_diff_weighted_norm);

                opti_stats_.max_angle_diff_deg_history_.push_back(max_abs_dihedral_diff_deg);
                opti_stats_.avg_angle_diff_deg_history_.push_back(avg_abs_dihedral_diff_deg);
                double dihedral_diff_weighted_norm_progress = (last_dihedral_diff_weighted_norm - dihedral_diff_weighted_norm) / last_dihedral_diff_weighted_norm;

                double unweighted_initial_quality_e = energy_record.initial_quality_energy_ / quality_weight_scaling_;
                double unweighted_final_quality_e   = energy_record.final_quality_energy_   / quality_weight_scaling_;
                auto      quality_energy_progress = (unweighted_initial_quality_e - unweighted_final_quality_e) / unweighted_initial_quality_e;

                auto unweighted_initial_virtual_e = energy_record.initial_virtual_energy_ / virtual_weight_scaling_;
                unweighted_final_virtual_e      =   energy_record.final_virtual_energy_ / virtual_weight_scaling_;
                auto      virtual_energy_progress = (unweighted_initial_virtual_e - unweighted_final_virtual_e) / unweighted_initial_virtual_e;

                //std::cout<<" diff = "<<(unweighted_initial_virtual_e - unweighted_final_virtual_e)<<std::endl;

                std::vector<double> quality_vector;
                double min_quality;
                compute_element_quality(quality_vector, min_quality);

                min_rel_vol = compute_min_relative_volume();
                min_rel_edge_len = compute_min_edge_length(mesh_) / min_codomain_boundary_edge_length_;

                opti_stats_.unweighted_initial_quality_energy_history_.push_back(unweighted_initial_quality_e);
                opti_stats_.unweighted_initial_primary_virtual_energy_history_.push_back(unweighted_initial_virtual_e.primary);
                opti_stats_.unweighted_initial_secondary_virtual_energy_history_.push_back(unweighted_initial_virtual_e.secondary);

                opti_stats_.unweighted_final_quality_energy_history_.push_back(unweighted_final_quality_e);
                opti_stats_.unweighted_final_primary_virtual_energy_history_.push_back(unweighted_final_virtual_e.primary);
                opti_stats_.unweighted_final_secondary_virtual_energy_history_.push_back(unweighted_final_virtual_e.secondary);

                opti_stats_.unweighted_final_quality_gradient_norm_history_.push_back(energy_record.final_quality_gradient_norm_);
                opti_stats_.unweighted_final_primary_virtual_gradient_norm_history_.push_back(energy_record.final_virtual_gradient_norm_.primary);
                opti_stats_.unweighted_final_secondary_virtual_gradient_norm_history_.push_back(energy_record.final_virtual_gradient_norm_.secondary);

                opti_stats_.final_total_gradient_max_coeff_history_.push_back(energy_record.final_total_gradient_max_coeff_);
                opti_stats_.final_total_Hessian_max_coeff_history_.push_back(energy_record.final_total_Hessian_max_coeff_);

                opti_stats_.worst_quality_history_.push_back(min_quality);
                opti_stats_.min_relative_volume_history_.push_back(min_rel_vol);
                opti_stats_.min_relative_edge_length_history_.push_back(min_rel_edge_len);

                if (enable_full_history) {
                    opti_stats_.quality_history_.push_back(quality_vector);
                }



                log(LogLevel::IMPORTANT) << "  avg abs dihedral angle difference =  " << avg_abs_dihedral_diff_deg << "°" << std::endl;
                log(LogLevel::IMPORTANT) << "  max abs dihedral angle difference =  " << max_abs_dihedral_diff_deg << "°" << std::endl;
                log(LogLevel::IMPORTANT) << "     edge-weighted diff vector norm =  " << dihedral_diff_weighted_norm << "°" << std::endl;
                log(LogLevel::IMPORTANT) << "               vector norm progress =  " << dihedral_diff_weighted_norm_progress << std::endl;
                log(LogLevel::IMPORTANT) << " ------------------------------------------------- " << std::endl;
                log(LogLevel::IMPORTANT) << "             initial quality energy =  " << energy_record.initial_quality_energy_ << std::endl;
                log(LogLevel::IMPORTANT) << "             final   quality energy =  " << energy_record.final_quality_energy_ << std::endl;
                log(LogLevel::IMPORTANT) << "  unweighted initial quality energy =  " << unweighted_initial_quality_e << std::endl;
                log(LogLevel::IMPORTANT) << "  unweighted   final quality energy =  " << unweighted_final_quality_e << std::endl;
                log(LogLevel::IMPORTANT) << "            quality energy progress =  " << quality_energy_progress << std::endl;
                log(LogLevel::IMPORTANT) << " ------------------------------------------------- " << std::endl;
                log(LogLevel::IMPORTANT) << "             initial virtual energy =  " << energy_record.initial_virtual_energy_ << std::endl;
                log(LogLevel::IMPORTANT) << "             final   virtual energy =  " << energy_record.final_virtual_energy_ << std::endl;
                log(LogLevel::IMPORTANT) << "  unweighted initial virtual energy =  " << unweighted_initial_virtual_e << std::endl;
                log(LogLevel::IMPORTANT) << "  unweighted   final virtual energy =  " << unweighted_final_virtual_e << std::endl;
                log(LogLevel::IMPORTANT) << "            virtual energy progress =  " << virtual_energy_progress << std::endl;
                log(LogLevel::IMPORTANT) << " ------------------------------------------------- " << std::endl;
                log(LogLevel::IMPORTANT) << "              worst element quality =  " << min_quality << std::endl;
                log(LogLevel::IMPORTANT) << "                min relative volume =  " << min_rel_vol << std::endl;
                log(LogLevel::IMPORTANT) << "           min relative edge length =  " << min_rel_edge_len << std::endl;
                log(LogLevel::IMPORTANT) << " ------------------------------------------------- " << std::endl;
                log(LogLevel::IMPORTANT) << "   max (abs) quality gradient coeff =  " << energy_record.final_quality_gradient_max_coeff_ << std::endl;
                log(LogLevel::IMPORTANT) << "   max (abs) virtual gradient coeff =  " << energy_record.final_virtual_gradient_max_coeff_ << std::endl;
                log(LogLevel::IMPORTANT) << "   max (abs)   total gradient coeff =  " << energy_record.final_total_gradient_max_coeff_ << std::endl;
                log(LogLevel::IMPORTANT) << "   max (abs)    total Hessian coeff =  " << energy_record.final_total_Hessian_max_coeff_ << std::endl;
                log(LogLevel::IMPORTANT) << " ---------------------------------------------------------------------- " << std::endl;

                last_dihedral_diff_weighted_norm = dihedral_diff_weighted_norm;


                bool remeshed_something = false;
                if (enable_focused_remeshing) {
                    auto hc_result = run_remeshing(FOCUSED_REMESHING, outer_wrap_it, 0, quick_remeshing_candidate_quality_threshold, quick_remeshing_target_quality_threshold);
                    if (hc_result == -1) {
                        log(LogLevel::ERROR) << " error while remeshing" << std::endl;
                        return WRAPPING_ERROR;
                    }
                    remeshed_something = hc_result;
                }

                bool exit_inner_loop(false);
                if (enable_relaxation_pass && !outer_wrap_it) {
                    log(LogLevel::IMPORTANT) << " -> end of relaxation pass, moving on to next pass" << std::endl;
                    exit_inner_loop = true;
                    made_progress_this_pass = true;

                }else if(stop_once_below_dihedral_angle_diff_threshold && max_abs_dihedral_diff_deg < max_dihedral_diff_deg_thresh) {
                    log(LogLevel::IMPORTANT) << " -> final max angle diff is below snapping threshold (" << max_dihedral_diff_deg_thresh << "), stopping here" << std::endl;
                    exit_inner_loop = true;
                    made_progress_this_pass = true;

                    //NOTE: I could override the comparison operator, but I want it to be clear what it means here
                }else if(virtual_energy_progress.primary   > virtual_energy_inner_progress_thresh ||
                         virtual_energy_progress.secondary > virtual_energy_inner_progress_thresh ||
                         virtual_energy_progress.tertiary  > virtual_energy_inner_progress_thresh){

                    auto updated_virtual_weight_scaling = virtual_weight_scaling_ * virtual_w_scale_base_variation;
                    log(LogLevel::IMPORTANT) << " -> progress was made, increasing virtual energy weight from "<<virtual_weight_scaling_<<" to "<< updated_virtual_weight_scaling << std::endl;
                    virtual_weight_scaling_ = updated_virtual_weight_scaling;
                    //actual progress is now measured after the pass was completed
                    //made_progress_this_pass = true;

                }else if(!inner_wrap_it) {
                    auto updated_virtual_weight_scaling = virtual_weight_scaling_ * virtual_w_scale_base_variation;
                    log(LogLevel::IMPORTANT) << " -> no progress made, but first iteration, increasing virtual energy weight from "<<virtual_weight_scaling_<<" to "<< updated_virtual_weight_scaling << std::endl;
                    virtual_weight_scaling_ = updated_virtual_weight_scaling;

                }else if(final_max_H_coeff_thresh > 0.0 && energy_record.final_total_Hessian_max_coeff_ < final_max_H_coeff_thresh && energy_record.final_total_gradient_max_coeff_ < final_max_H_coeff_thresh){
                    auto updated_virtual_weight_scaling = virtual_weight_scaling_ * virtual_w_scale_base_variation;
                    log(LogLevel::IMPORTANT) << " -> no progress made, but max gradient and Hessian coeffs are still below threshold, increasing virtual energy weight from "<<virtual_weight_scaling_<<" to "<< updated_virtual_weight_scaling << std::endl;
                    virtual_weight_scaling_ = updated_virtual_weight_scaling;

                }else if(max_virtual_w > 0 &&
                (virtual_weight_scaling_.primary > max_virtual_w ||
                        virtual_weight_scaling_.secondary > max_virtual_w ||
                        virtual_weight_scaling_.tertiary > max_virtual_w)){
                    log(LogLevel::IMPORTANT) << " -> weight scaling above threshold, moving on to next pass" << std::endl;
                    exit_inner_loop = true;
                }/*else if (remeshed_something) {
                    auto updated_virtual_weight_scaling = virtual_weight_scaling_ * virtual_w_scale_base_variation;
                    log(LogLevel::IMPORTANT) << " -> nothing snapped and no progress made, but something was remeshed, increasing virtual energy weight from "<<virtual_weight_scaling_<<" to "<< updated_virtual_weight_scaling << std::endl;

                }*/else{
                    log(LogLevel::IMPORTANT) << " -> no progress was made, moving on to next pass" << std::endl;
                    exit_inner_loop = true;
                }

                enter_snap_phase_if_possible(stop_once_below_dihedral_angle_diff_threshold, false, max_abs_dihedral_diff_deg, max_dihedral_diff_deg_thresh);

                int snapped_candidates_count = 0;
                if (snap_phase_) {
                    snapped_candidates_count = try_to_snap_next_candidates();
                    log(LogLevel::IMPORTANT) << " --> snapped a total of " << snapped_candidates_count <<" non-wrapped candidates" << std::endl;

                    made_progress_this_pass |= snapped_candidates_count;

                    if (json_data_["enable_debugging_meshes_exportation"]) {
                        std::string ext = json_data_["debugging_meshes_exportation_extension"];
                        mesh_.export_conforming_mesh(mesh_name_+"_codomain_mesh_after_it_" + std::to_string(outer_wrap_it) + "_" + std::to_string(inner_wrap_it) + "_and_snapping" + ext);
                    }
                }

                opti_stats_.wrapped_ratio_history_.push_back((double) wrapped_count_ / (double) initial_boundary_vertices_count);

                if(snap_phase_ && stop_at_virtual_energy_convergence){
                    log(LogLevel::IMPORTANT) << " entered snap phase, stopping here " << std::endl;
                    return WRAPPING_SUCCESS;
                }

                inner_wrap_it++;
                opti_stats_.inner_iters_++;
                log(LogLevel::IMPORTANT) << " ================================================================================================================================ " << std::endl;
                if(exit_inner_loop){
                    break;
                }

                //that's more of a debugging thing
                if(max_total_inner_iters > 0 && opti_stats_.inner_iters_ >= max_total_inner_iters){
                    log(LogLevel::IMPORTANT) << " reached max number of inner iterations, stopping here " << std::endl;
                    return WRAPPING_FAILURE;
                }
            }

            if (enable_relaxation_pass && !outer_wrap_it) {
                log(LogLevel::IMPORTANT) << " -> end of relaxation pass, resetting virtual vertices positions " << std::endl;
                setup_virtual_vertices_positions();
            }

            if(stop_once_below_dihedral_angle_diff_threshold && max_abs_dihedral_diff_deg < max_dihedral_diff_deg_thresh) {
                log(LogLevel::IMPORTANT) << " --> dihedral angle is below snap thresh, counting this as a success" << std::endl;
                return WRAPPING_SUCCESS;
            }

            last_pass_unweighted_final_virtual_e = unweighted_final_virtual_e;

            auto last_pass_virtual_energy_progress = (last_pass_unweighted_initial_virtual_e - last_pass_unweighted_final_virtual_e) / last_pass_unweighted_initial_virtual_e;

            log(LogLevel::IMPORTANT) << " =================================================================  " << std::endl;
            log(LogLevel::IMPORTANT) << "  last pass unweighted initial virtual energy =  " << last_pass_unweighted_initial_virtual_e << std::endl;
            log(LogLevel::IMPORTANT) << "  last pass unweighted   final virtual energy =  " << last_pass_unweighted_final_virtual_e << std::endl;
            log(LogLevel::IMPORTANT) << "            last pass virtual energy progress =  " << last_pass_virtual_energy_progress << std::endl;
            log(LogLevel::IMPORTANT) << " =================================================================  " << std::endl;

            if(last_pass_virtual_energy_progress.primary > virtual_energy_outer_progress_thresh ||
               last_pass_virtual_energy_progress.secondary > virtual_energy_outer_progress_thresh ||
               last_pass_virtual_energy_progress.tertiary > virtual_energy_outer_progress_thresh) {
                auto updated_virtual_weight_scaling = virtual_weight_scaling_ * virtual_w_scale_base_variation;
                log(LogLevel::IMPORTANT) << " -> progress was made during the last pass! Moving on " << std::endl;
                made_progress_this_pass = true;
            }

            if(outer_wrap_it && !made_progress_this_pass) {
                if (enable_backup_DA_energy_mode && !opti_stats_.backup_mode_) {
                    opti_stats_.backup_mode_ = true;
                    log(LogLevel::IMPORTANT) << " --> no progress during this pass, entering backup mode" << std::endl;
                }else {
                    log(LogLevel::IMPORTANT) << " --> no progress during this pass, stopping here" << std::endl;
                    wrapping_result = WRAPPING_FAILURE;
                    break;
                }
            }

            if(!snap_phase_ && enable_forced_snapping) {
                log(LogLevel::IMPORTANT) << " --> forcing snap phase after first pass" << std::endl;

                enter_snap_phase_if_possible(stop_once_below_dihedral_angle_diff_threshold,
                                             true,
                                             max_abs_dihedral_diff_deg,
                                             max_dihedral_diff_deg_thresh);

                int snapped_candidates_count = try_to_snap_next_candidates();
                log(LogLevel::IMPORTANT) << " --> snapped a total of " << snapped_candidates_count <<
                        " non-wrapped candidates" << std::endl;

                if (json_data_["enable_debugging_meshes_exportation"]) {
                    std::string ext = json_data_["debugging_meshes_exportation_extension"];
                    mesh_.export_conforming_mesh(
                        mesh_name_ + "_codomain_mesh_after_it_" + std::to_string(outer_wrap_it) + "_" +
                        std::to_string(inner_wrap_it) + "_and_snapping" + ext);
                }
            }

            if (wrapped_count_ != initial_boundary_vertices_count) {

                log(LogLevel::IMPORTANT) << " ---> still some vertices to snap, remeshing and moving on to next pass " << std::endl;

                compute_element_quality(quality_vector, min_quality);
                double avg_quality(0);
                for(auto& q: quality_vector){ avg_quality += q;}
                opti_stats_.avg_quality_before_remeshing_history_.push_back(avg_quality / quality_vector.size());
                opti_stats_.min_quality_before_remeshing_history_.push_back(min_quality);

                double avg_quality_before = opti_stats_.avg_quality_before_remeshing_history_.back();
                double min_quality_before = opti_stats_.min_quality_before_remeshing_history_.back();;
                double default_quality_bound = json_data_["global_wrapping_IVL_bound"];

                //last-minute hack, not very safe
                //json_data_["global_wrapping_IVL_bound"] = 1.0 / min_quality_before;
                //log(LogLevel::INFO)<<" -> replaced global IVL bound "<<default_quality_bound<<" with min-quality-adapted one: "<<(double)json_data_["global_wrapping_IVL_bound"]<<std::endl;

                auto hc_result = run_remeshing(GLOBAL_REMESHING,
                                               outer_wrap_it,
                                               deep_remeshing_inner_iters,
                                               deep_remeshing_candidate_quality_threshold,
                                               deep_remeshing_target_quality_threshold);

                compute_element_quality(quality_vector, min_quality);
                avg_quality = 0.0;
                for(auto& q: quality_vector){ avg_quality += q;}
                opti_stats_.avg_quality_after_remeshing_history_.push_back(avg_quality / quality_vector.size());
                opti_stats_.min_quality_after_remeshing_history_.push_back(min_quality);

                json_data_["global_wrapping_IVL_bound"] = default_quality_bound;

                double IVL_bound_scaling = json_data_["global_wrapping_IVL_bound_scaling"];

                //extra safety factor is kind of a hack because the bound is measure in IVL, but the quality is VL
                if(json_data_["enable_all_safety_checks"]) {
                    if (min_quality < (min_quality_before / IVL_bound_scaling) * 0.9) {
                        log(LogLevel::ERROR) << " ERROR - min quality decreased from " << min_quality_before << " to " << min_quality << ", which is "
                                             << (min_quality_before / min_quality) << " times worse" << std::endl;
                        return WRAPPING_ERROR;
                    }
                    if (avg_quality < (avg_quality_before / IVL_bound_scaling) * 0.9) {
                        log(LogLevel::ERROR) << " ERROR - avg quality decreased from " << avg_quality_before << " to " << avg_quality << ", which is "
                                             << (avg_quality_before / avg_quality) << " times worse" << std::endl;
                        return WRAPPING_ERROR;
                    }
                }

                if (hc_result == -1) {
                    log(LogLevel::ERROR) << " error while remeshing" << std::endl;
                    return WRAPPING_ERROR;
                }
                opti_stats_.deep_remeshing_iterations_.push_back(opti_stats_.inner_iters_);
                if (snap_phase_) {

                    int snapped_candidates_count = try_to_snap_next_candidates();
                    log(LogLevel::IMPORTANT) << " --> snapped a total of " << snapped_candidates_count <<" non-wrapped candidates" << std::endl;

                    if (json_data_["enable_debugging_meshes_exportation"]) {
                        std::string ext = json_data_["debugging_meshes_exportation_extension"];
                        mesh_.export_conforming_mesh(mesh_name_+"_codomain_mesh_after_it_" + std::to_string(outer_wrap_it) + "_" + std::to_string(inner_wrap_it) + "_and_snapping" + ext);
                    }
                }
            }


            log(LogLevel::IMPORTANT) << " -----> done with pass " << outer_wrap_it << " in " << pass_sw.total_duration() << " seconds" << std::endl;
            last_pass_unweighted_initial_virtual_e = last_pass_unweighted_final_virtual_e;

            outer_wrap_it++;
        }

        double max_constraints_violation = snap_satisfied_boundary_constraints(constr_eps);

        opti_stats_.total_unaccounted_for_time_s_ = main_sw_.total_duration() - (opti_stats_.total_opti_time_s_ + opti_stats_.total_upkeep_time_s_ + opti_stats_.total_remeshing_time_s_);

        log(LogLevel::IMPORTANT)<<" ================================================================================================================================ "<<std::endl;
        log(LogLevel::IMPORTANT)<<" --> done with global optimization process after "<<opti_stats_.outer_iters_<<" passes and "<<opti_stats_.inner_iters_<<" iterations in "<<main_sw_.total_duration()<<" seconds. Timing details: "<<std::endl;
        log(LogLevel::IMPORTANT)<<" --------------------------------------------------------------------------- "<<std::endl;
        log(LogLevel::IMPORTANT)<<" ->      opti time: "<<std::setw(5)<<opti_stats_.total_opti_time_s_<<"s"<<std::endl;
        log(LogLevel::IMPORTANT)<<" ->    upkeep time: "<<std::setw(5)<<opti_stats_.total_upkeep_time_s_<<"s"<<std::endl;
        log(LogLevel::IMPORTANT)<<" -> remeshing time: "<<std::setw(5)<<opti_stats_.total_remeshing_time_s_<<"s"<<std::endl;
        log(LogLevel::IMPORTANT)<<" ->          other: "<<std::setw(5)<<opti_stats_.total_unaccounted_for_time_s_<<"s"<<std::endl;
        log(LogLevel::IMPORTANT)<<" --------------------------------------------------------------------------- "<<std::endl;
        log(LogLevel::IMPORTANT)<<"    Remeshing details: "<<std::endl;
        log(LogLevel::IMPORTANT)<<" --------------------------------------------------------------------------- "<<std::endl;
        log(LogLevel::IMPORTANT)<<" - edge     swaps: "<<opti_stats_.total_swap_count_<<std::endl;
        log(LogLevel::IMPORTANT)<<" - edge collapses: "<<opti_stats_.total_collapse_count_<<std::endl;
        log(LogLevel::IMPORTANT)<<" - edge    splits: "<<opti_stats_.total_split_count_<<std::endl;
        log(LogLevel::IMPORTANT)<<" ================================================================================================================================ "<<std::endl;

        if(wrapping_result == WRAPPING_UNDEFINED_RESULT) {
            wrapping_result = run_final_checks(constr_eps, main_sw_);
        }

    }catch(TimeOutException e) {
        log(LogLevel::INFO)<<"--------------------"<<std::endl;
        log(LogLevel::INFO)<<"||||||||||||||||||||"<<std::endl;
        log(LogLevel::INFO)<<e.what()<<std::endl;
        log(LogLevel::INFO)<<"||||||||||||||||||||"<<std::endl;
        log(LogLevel::INFO)<<"--------------------"<<std::endl;
        wrapping_result = WRAPPING_TIMEOUT;

    }catch(VirtualVertexSetupException e){
        log(LogLevel::INFO)<<"--------------------"<<std::endl;
        log(LogLevel::INFO)<<"||||||||||||||||||||"<<std::endl;
        log(LogLevel::INFO)<<e.what()<<std::endl;
        log(LogLevel::INFO)<<"||||||||||||||||||||"<<std::endl;
        log(LogLevel::INFO)<<"--------------------"<<std::endl;
        wrapping_result = WRAPPING_NUMERICS_ERROR;
    }

    std::vector<double> quality_vector;
    double min_quality, min_rel_vol, min_rel_edge_len;
    compute_element_quality(quality_vector, min_quality);

    min_rel_vol = compute_min_relative_volume();
    min_rel_edge_len = compute_min_edge_length(mesh_) / min_codomain_boundary_edge_length_;

    opti_stats_.final_worst_quality_ = min_quality;
    opti_stats_.final_min_relative_volume_ = min_rel_vol;
    opti_stats_.final_min_relative_edge_length_ = min_rel_edge_len;


    return wrapping_result;
}





template <typename DistortionElemT, typename VirtualElemT>
double
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
compute_ref_to_mesh_volume_scaling_factor() const{

    auto mesh_total_volume = compute_total_volume(mesh_);

    return std::pow(mesh_total_volume / ref_mesh_total_volume_, 1.0 / 3.0);

}



template <typename DistortionElemT, typename VirtualElemT>
int
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
run_remeshing(REMESHING_METHOD remeshing_method, int opti_step, int remeshing_outer_iters, double candidate_quality_threshold, double target_quality_threshold) {



    check_for_timeout();

    //to skip the analysis, output, etc. if there's nothing to remesh
    if (remeshing_method == FOCUSED_REMESHING) {

        bool found_low_qual_tet(false);
        for (auto ch: mesh_.cells()) {
            auto qual = mesh_.compute_quality_vector(ch, VOLUME_LENGTH, CONFORMING)[0];
            if (qual < target_quality_threshold) {
                found_low_qual_tet = true;
                break;
            }
        }
        if (!found_low_qual_tet) {
            log(LogLevel::IMPORTANT)<<" --> no element with quality below target threshold "<<target_quality_threshold<<", skipping remeshing "<<std::endl;
            return 0;
        }
    }

    opti_stats_.total_upkeep_time_s_ += main_sw_.lap_duration();

    log(LogLevel::IMPORTANT)<<" ================================================================================================================================ "<<std::endl;
    log(LogLevel::IMPORTANT)<<" -> entering remeshing step... "<<std::endl;


    if(json_data_["wrapping_enable_remeshing_element_analysis"]) {
        log(LogLevel::INFO)<<" -> running pre-remeshing analysis..."<<std::endl;

        analyze_all_elements(ref_mesh_, mesh_, compute_ref_to_mesh_volume_scaling_factor(), true, compute_total_volume(mesh_) / mesh_.n_logical_cells(), CONFORMING,"");
        analyze_all_elements(ref_mesh_, ref_mesh_, 1.0, false, ref_mesh_total_volume_ / ref_mesh_.n_logical_cells(), CONFORMING,"");
    }


    //kind of a hack because of legacy approach
    json_data_["ref_mesh_pre_remeshing_enable_smoothing"] = true;
    json_data_["ref_mesh_pre_remeshing_max_iterations"] = remeshing_outer_iters;
    json_data_["ref_mesh_pre_remeshing_candidate_quality_threshold"] = candidate_quality_threshold;
    json_data_["ref_mesh_pre_remeshing_target_quality_threshold"] = target_quality_threshold;
    json_data_["ref_mesh_pre_remeshing_relative_volume_threshold"] = json_data_["global_wrapping_IVL_rel_vol_bound"];
    json_data_["ref_mesh_pre_remeshing_enable_early_stop"] = json_data_["global_wrapping_deep_remeshing_enable_early_stop"];


    int new_edges_start_idx(mesh_.n_edges());

    RemeshingStats remeshing_stats;
    if (remeshing_method == GLOBAL_REMESHING) {
        remeshing_stats = hill_climbing_ref_mesh_remeshing_with_phases(mesh_, ref_mesh_, json_data_, true);
    } else if (remeshing_method == FOCUSED_REMESHING) {
        remeshing_stats = focused_remeshing(mesh_, ref_mesh_, json_data_["global_wrapping_enable_focused_splitting"], json_data_, json_data_["global_wrapping_enable_focused_remeshing_refill"]);
    }else{
        log(LogLevel::ERROR)<<" --> unhandled remeshing mode "<<remeshing_method<<std::endl;
        return WRAPPING_ERROR;
    }
    check_for_timeout();

    opti_stats_.total_split_count_    += remeshing_stats.split_count;
    opti_stats_.total_swap_count_     += remeshing_stats.swap_count;
    opti_stats_.total_collapse_count_ += remeshing_stats.collapse_count;


    std::vector<EdgeHandle> new_boundary_edges;
    for(int idx(new_edges_start_idx); idx < mesh_.n_edges(); idx++){
        if(mesh_.is_boundary(EdgeHandle(idx))){
            new_boundary_edges.push_back(EdgeHandle(idx));
        }
    }
    if(!new_boundary_edges.empty()){
        log(LogLevel::ERROR)<<" ERROR - found "<<new_boundary_edges.size()<<" new boundary edges"<<std::endl;
        return WRAPPING_ERROR;
    }

    /*if(true || json_data_["enable_post_remeshing_ref_mesh_opti"]) {
        double ref_mesh_target_volume = ref_mesh_total_volume_ / ref_mesh_.n_logical_elements();

        log(LogLevel::INFO)<<" --------------------------------------------------------- "<<std::endl;
        log(LogLevel::INFO)<<" -> optimizing reference mesh, using target volume = "<<ref_mesh_target_volume<<std::endl;
        if(optimize_mesh<NonConformingTetrahedralMesh, DistortionElemT>(ref_mesh_, ref_mesh_target_volume, 1000, true, nullptr, false, json_data_, LogConfig::get_min_level()> LogLevel::INFO)){
            log(LogLevel::ERROR)<<" error while optimizing reference mesh"<<std::endl;
            return WRAPPING_ERROR;
        }
        log(LogLevel::INFO)<<" --------------------------------------------------------- "<<std::endl;
    }*/

    if (remeshing_stats.remeshed_something < 0) {
        log(LogLevel::ERROR)<<" error during remeshing "<<std::endl;
        return WRAPPING_ERROR;
    }

    if (!remeshing_stats.remeshed_something) {
        log(LogLevel::IMPORTANT)<<" --> found nothing to remesh "<<std::endl;
        return 1;
    }

    if (auto res = run_safety_checks("global remeshing")) {
        return res;
    }

    if(json_data_["wrapping_enable_remeshing_element_analysis"]) {
        log(LogLevel::INFO)<<" -> running post-remeshing analysis..."<<std::endl;

        analyze_all_elements(ref_mesh_, mesh_, compute_ref_to_mesh_volume_scaling_factor(), true, compute_total_volume(mesh_) / mesh_.n_logical_cells(), CONFORMING,"");
        analyze_all_elements(ref_mesh_, ref_mesh_, 1.0, false, ref_mesh_total_volume_ / ref_mesh_.n_logical_cells(), CONFORMING,"");
    }
    log(LogLevel::IMPORTANT)<<" ================================================================================================================================ "<<std::endl;


    opti_stats_.total_remeshing_time_s_ += main_sw_.lap_duration();

    return 0;

}




template <typename DistortionElemT, typename VirtualElemT>
int
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
run_global_wrapping(int newton_iters,
                    bool &converged,
                    EnergyRecord& energy_record,
                    Eigen::VectorXd &quality_g,
                    std::array<Eigen::VectorXd, 3> &virtual_gs,
                    bool mock_run) {

    log(LogLevel::IMPORTANT)<<" ================================================================================================================================ "<<std::endl;
    log(LogLevel::IMPORTANT)<<" - entering optimization step..."<<std::endl;
    //std::cout<<" - scaling target increase: "<<virtual_w_scaling_target_increase<<std::endl;

    if (snap_phase_ && !root_hfh_.is_valid()) {
        log(LogLevel::ERROR)<<" ERROR - in snap phase but the root halfface is not set"<<std::endl;
        return WRAPPING_ERROR;
    }

    using EigenVec = Eigen::Vector3d;
    const int geo_dim = 3;

    std::vector<EigenVec> P_ref;
    //setup_regular_ref_element<MeshT>(target_vol_per_element, P_ref);

    int n_vars(mesh_.n_logical_vertices() * geo_dim);
    int virtual_vertex_idx_start = n_vars;

    //double fof_epsilon(0);
    //std::string ref_opti_element_type = json_data["ref_opti_element_type"];


    auto global_vid_prop = mesh_.request_vertex_property<int>();
    int global_vid(0);
    for(auto v: mesh_.vertices()){
        global_vid_prop[v] = global_vid;
        global_vid++;
    }

    COMISO::FiniteElementProblem fe_problem(n_vars);
    log(LogLevel::IMPORTANT)<<" -> initial n_vars: "<<n_vars<<std::endl;

    const double constr_eps = json_data_["global_wrapping_constr_eps"];

    const bool enable_growing_front_eq_constraints = json_data_["global_wrapping_enable_growing_front_equality_constraints"];
    const bool enable_full_history = json_data_["global_wrapping_enable_full_history"];

    const double edge_length_bound_elements_w = json_data_["global_wrapping_edge_length_bound_elements_w"];
    const double edge_length_bound_elements_boundary_rel_bound= json_data_["global_wrapping_edge_length_bound_elements_boundary_rel_bound"];
    const double edge_length_bound_elements_bound_scaling = json_data_["global_wrapping_edge_length_bound_elements_bound_scaling"];
    const bool enable_edge_length_bound_elements_for_boundary = json_data_["global_wrapping_enable_edge_length_bound_elements_for_boundary"];

    double dihedral_angle_energy_diff_bound_deg = json_data_["global_wrapping_virtual_DA_diff_bound_deg"];

    if (!opti_stats_.inner_iters_ && dihedral_angle_energy_diff_bound_deg > 0) {
        log(LogLevel::WARNING)<<" -> WARNING HACK STUFF JUST FOR QUICK TEST"<<std::endl;
        dihedral_angle_energy_diff_bound_deg = 0;
        log(LogLevel::WARNING)<<" -> first pass, disabling angle energy for now"<<std::endl;
        log(LogLevel::WARNING)<<" -> END OF HACKY STUFF"<<std::endl;
    }


    using VElemT_A = VirtualElemT;
    using VElemT_B = AUG_ANGLE_3D_PH_AD;
    using VElemT_C = AUG_FOF_NO_EPS_3D_AD;

    COMISO::FiniteElementSet<VElemT_A> fe_virtual_A("virtual A (" + VElemT_A::element_name()+")");
    COMISO::FiniteElementSet<VElemT_B> fe_virtual_B("virtual B (" + VElemT_B::element_name()+")");
    COMISO::FiniteElementSet<VElemT_C> fe_virtual_C("virtual C (" + VElemT_C::element_name()+")");

    //std::vector<double> unweighted_virtual_energies;
    NonConformingTetrahedralMesh virtual_elements_mesh, ref_virtual_elements;

    auto virtual_setup_res = setup_holonomy_virtual_elements(mesh_,
                                                             codomain_pos_prop_,
                                                             global_vid_prop,
                                                             virtual_weight_scaling_,
                                                             json_data_,
                                                             virtual_vertex_pos_prop_,
                                                             boundary_dihedral_angle_deg_prop_,
                                                             fe_problem,
                                                             fe_virtual_A,
                                                             fe_virtual_B,
                                                             fe_virtual_C,
                                                             dihedral_angle_energy_diff_bound_deg,
                                                             &virtual_elements_mesh,
                                                             &ref_virtual_elements);

    if(virtual_setup_res == -2){
        log(LogLevel::ERROR) << " error while computing virtual vertices positions during virtual elements setup" << std::endl;
        throw VirtualVertexSetupException();
    }else if (virtual_setup_res) {
        log(LogLevel::ERROR) << " error while setting up holonomy virtual elements" << std::endl;
        return -1;
    }



    if (json_data_["enable_debugging_meshes_exportation"]) {
        Eigen::Vector3d min_corner, max_corner;
        compute_bounding_box(mesh_, min_corner, max_corner);
        double scaling = 1.0 / (max_corner - min_corner).norm();
        auto center = compute_barycenter(mesh_);
        for(auto vh: virtual_elements_mesh.vertices()){
            auto pos = virtual_elements_mesh.vertex(vh);
            virtual_elements_mesh.set_vertex(vh, 1.1 * scaling * (pos - vec2vec(center)));
        }

        //virtual_elements_mesh.export_conforming_mesh("virtual_elements_before_it_" + std::to_string(opti_stats_.outer_iters_) + "_" + std::to_string(opti_stats_.inner_iters_) + ".ovm", false);
        //ref_virtual_elements.export_conforming_mesh("ref_virtual_elements_before_it_" + std::to_string(opti_stats_.outer_iters_) + "_" + std::to_string(opti_stats_.inner_iters_) + ".ovm", true);
    }

    n_vars = fe_problem.n_unknowns();
    log(LogLevel::IMPORTANT)<<" -> updated n_vars to "<<n_vars<<std::endl;


    std::vector<COMISO::LinearConstraint> eq_constraints;
    eq_constraints.reserve(n_vars);
    //COMISO::FiniteElementSet <BaseElemT> fe_barrier("barrier elements");
    COMISO::FiniteElementSet <DistortionElemT> fe_barrier("distortion (" + DistortionElemT::element_name()+")");
    COMISO::FiniteElementSet <ELEN_B_3D_AD> fe_edge_length("edge length (" + ELEN_B_3D_AD::element_name()+")");


    double max_barrier_e(0), max_unw_barrier_e(0), max_w(0);
    double w_sum(0);

    int backup_bound_use_count(0), backup_vol_bound_use_count(0);
    double max_bound(0);

    double w = quality_weight_scaling_ / mesh_.n_logical_elements();

    /*std::vector<double> parameters = {w * Dirichlet_w, Dirichlet_barrier_w, Dirichlet_barrier_min};
    if(ref_opti_element_type == "fof"){
        parameters = {fof_epsilon, conf_w, vol_w, Dirichlet_barrier_w, Dirichlet_barrier_min};
    }else if (ref_opti_element_type == "IVL") {
        parameters = {json_data["wrapping_relaxation_IVL_w"], 0.0, 1.0};
        //log(LogLevel::INFO)<<" - IVL weight = "<<parameters[0]<<std::endl;
    }*/

    //double avg_vol = compute_total_volume(mesh_) / mesh_.n_logical_cells();
    double avg_vol = target_total_volume_ / mesh_.n_logical_cells();

    std::vector<double> parameters;
    log(LogLevel::IMPORTANT)<<" -> average volume (used as optimum for volume bound): "<<avg_vol<<std::endl;
    double vol_bound = double(json_data_["global_wrapping_IVL_rel_vol_bound"]) * avg_vol;
    log(LogLevel::IMPORTANT)<<" -> volume bound: "<<vol_bound<<std::endl;

    //NOTE: those are the "base" weights. They are scaled in a per-element basis during elements setup
    if(DistortionElemT::element_code() == "IVL_B") {
        parameters = {
            json_data_["global_wrapping_IVL_w"],
            json_data_["global_wrapping_IVL_bound_w"],
            json_data_["global_wrapping_IVL_bound"],
            json_data_["global_wrapping_IVL_vol_w"],
            vol_bound,
            avg_vol
        };

    }else{
        log(LogLevel::ERROR)<<" ERROR - global wrapping is not implemented for distortion element "<<DistortionElemT::element_name()<<std::endl;
        return WRAPPING_ERROR;
    }



    //add the barrier elements
    for (auto elh: mesh_.elements()) {

        //log(LogLevel::INFO)<<" --------------------------------------------------- handling element "<<mesh.element_to_str(elh)<<std::endl;

        w_sum += w;
        max_w = std::max(w, max_w);

        typename DistortionElemT::VecI vi_barrier;
        typename DistortionElemT::VecC vc_barrier;
        typename DistortionElemT::VecV vx_barrier;

        auto verts = mesh_.get_element_vertices(elh);

        for (int local_vid(0); local_vid < verts.size(); local_vid++) {
            const auto &vh = verts[local_vid];

            int global_id = geo_dim * global_vid_prop[vh];
            for (int i(0); i < geo_dim; i++) {
                vi_barrier[geo_dim * local_vid + i] = global_id + i;
            }

            for (int i(0); i < geo_dim; i++) {
                fe_problem.x()[global_id + i] = mesh_.vertex(vh)[i];
                vx_barrier[geo_dim * local_vid + i] = fe_problem.x()[global_id + i];
            }
        }

        auto setup_res = DistortionElemT::compute_constants(parameters, P_ref, vc_barrier);
        /*auto setup_res = BaseElemT::compute_constants(w * Dirichlet_w, Dirichlet_barrier_w, Dirichlet_barrier_min,
                                                      P_ref,
                                                      vc_barrier);*/

        DistortionElemT::scale_weights(w, vc_barrier);


        //log(LogLevel::INFO)<<" - constants vector: "<<vc_barrier.transpose()<<std::endl;

        if (setup_res) {
            log(LogLevel::ERROR) << " ERROR: constants couldn't be computed for face: " << elh << std::endl;
            return -1;
        }

        DistortionElemT elem;
        double barrier_e = elem.eval_f(vx_barrier, vc_barrier);

        //std::cout<<" - barrier energy = "<<barrier_e<<", is finite: "<<std::isfinite(barrier_e)<<std::endl;

        double VL_qual = mesh_.compute_quality_vector(elh, VOLUME_LENGTH, CONFORMING)[0];
        double IVL_qual = 1.0 / VL_qual;

        if (!std::isfinite(barrier_e)) {
            log(LogLevel::INFO) << " ----------------------------------------------------------------------------" << std::endl;
            log(LogLevel::INFO) << " distortion energy is infinite for element "<<cell_to_str(mesh_, elh)<<", using a higher bound" << std::endl;
            double bound_w_save = vc_barrier[1];
            double vol_w_save = vc_barrier[3];
            vc_barrier[1] = 0; //disable the bounds to compute the actual value
            vc_barrier[3] = 0;

            double unw_barrier_e = elem.eval_f(vx_barrier, vc_barrier) / vc_barrier[0];

            log(LogLevel::INFO) << " -> unweighted base energy = "<<unw_barrier_e<< std::endl;

            double bound_scale_factor = this->json_data_["global_wrapping_IVL_bound_scaling"];
            double vol_bound_scale_factor = this->json_data_["global_wrapping_IVL_vol_scaling"];

            vc_barrier[1] = bound_w_save;
            vc_barrier[2] = unw_barrier_e * bound_scale_factor;
            vc_barrier[3] = vol_w_save;
            log(LogLevel::INFO) << " -> updated bound to "<<vc_barrier[2]<< std::endl;

            barrier_e = elem.eval_f(vx_barrier, vc_barrier);
            log(LogLevel::INFO) << " -> updated energy = "<<barrier_e<< std::endl;

            if (!std::isfinite(barrier_e)) {

                log(LogLevel::INFO) << " ------------------------------------------------------------" << std::endl;
                log(LogLevel::INFO) << " barrier energy is STILL infinite, using a lower volume bound" << std::endl;
                auto current_vol = compute_volume(mesh_, elh);
                log(LogLevel::INFO) << " -> current volume = "<<current_vol <<" (relative to average = "<<(current_vol / avg_vol)<<"), volume-length = "<<VL_qual<< std::endl;

                vc_barrier[4] = current_vol / vol_bound_scale_factor;
                log(LogLevel::INFO) << " -> updated vol bound to = "<<vc_barrier[4] << std::endl;

                barrier_e = elem.eval_f(vx_barrier, vc_barrier);
                log(LogLevel::INFO) << " -> updated energy = "<<barrier_e<< std::endl;
                backup_vol_bound_use_count++;
            }

            backup_bound_use_count++;
            max_bound = std::max(max_bound, vc_barrier[2]);
        }

        if (!std::isfinite(barrier_e)) {
            log(LogLevel::ERROR) << " barrier energy is still infinite for element "<<cell_to_str(mesh_, elh)<<" after using a higher bound. The volume is most likely too low" << std::endl;
            return -1;
        }

        max_barrier_e = std::max(max_barrier_e, barrier_e);
        max_unw_barrier_e = std::max(max_unw_barrier_e, barrier_e / w);
        //log(LogLevel::INFO)<<" - energy = "<<barrier_e<<std::endl;


        fe_barrier.instances().add_element(vi_barrier,
                                           vc_barrier);

        if(IVL_qual >= parameters[2]){
            log(LogLevel::WARNING)<<" WARNING - IVL energry is "<<(barrier_e / w)<<" but volume-length is "<<VL_qual<<", inverted is "<<IVL_qual<<", which is larger than the bound "<<parameters[2]<<std::endl;
            vc_barrier << 1.0, 0.0, 0.0, 0.0, 0.0;
            double unweighted_barrier_e = elem.eval_f(vx_barrier, vc_barrier);
            log(LogLevel::INFO)<<" unweighted IVL energy = "<<unweighted_barrier_e<<std::endl;
        }
    }

    log(LogLevel::INFO) << " ->            max barrier energy = " << max_barrier_e << std::endl;
    log(LogLevel::INFO) << " -> max unweighted barrier energy = " << max_barrier_e << std::endl;
    log(LogLevel::INFO) << " ->            max barrier weight = " << max_w << std::endl;
    log(LogLevel::INFO) << " ->            barrier weight sum = " << w_sum << std::endl;
    log(LogLevel::INFO) << " ->                     max bound = " << max_bound << std::endl;
    log(LogLevel::INFO) << " -> #elements using a backup energy bound: " << backup_bound_use_count <<" / "<<mesh_.n_logical_cells() << std::endl;
    log(LogLevel::INFO) << " -> #elements using a backup volume bound: " << backup_vol_bound_use_count <<" / "<<mesh_.n_logical_cells() << std::endl;

    fe_problem.add_set(&fe_barrier);


    if (edge_length_bound_elements_w > 0) {
        log(LogLevel::INFO) << " ---------------------------------------------------------" << std::endl;
        log(LogLevel::INFO) << " - setting up edge length bound elements... " << std::endl;

        double el_w_sum(0);
        double max_el_w(0);

        double default_el_target = avg_codomain_boundary_edge_length_;
        double el_bound = edge_length_bound_elements_boundary_rel_bound * min_codomain_boundary_edge_length_;
        double el_w = edge_length_bound_elements_w / mesh_.n_logical_edges();

        log(LogLevel::INFO) << " -                      weight = "<<el_w << std::endl;
        log(LogLevel::INFO) << " -  default target edge length = "<<default_el_target << std::endl;
        log(LogLevel::INFO) << " -    min boundary edge length = "<<min_codomain_boundary_edge_length_ << std::endl;
        log(LogLevel::INFO) << " -  boundary edge length bound = "<<el_bound << std::endl;

        int backup_bound_count(0);
        double min_len = std::numeric_limits<double>::max();
        double max_el_e = 0;

        for (auto eh: mesh_.edges()) {
            if(enable_edge_length_bound_elements_for_boundary || !mesh_.is_boundary(eh)) {
                std::vector<VertexHandle> verts = {mesh_.edge(eh).from_vertex(), mesh_.edge(eh).to_vertex()};

                double el_target = mesh_.is_boundary(eh) ? boundary_edge_length_prop_[eh] : default_el_target;

                el_w_sum += el_w;
                max_w = std::max(el_w, max_el_w);

                ELEN_B_3D_AD::VecI vi_el;
                ELEN_B_3D_AD::VecV vx_el;
                ELEN_B_3D_AD::VecC vc_el;
                vc_el << el_w, el_bound, el_target;

                double len = (mesh_.vertex(mesh_.edge(eh).from_vertex()) - mesh_.vertex(mesh_.edge(eh).to_vertex())).norm();
                if (len < el_bound) {
                    vc_el[1] = edge_length_bound_elements_bound_scaling * len;
                }
                min_len = std::min(min_len, len);

                for (int local_vid(0); local_vid < verts.size(); local_vid++) {
                    const auto &vh = verts[local_vid];

                    int global_id = geo_dim * global_vid_prop[vh];
                    for (int i(0); i < geo_dim; i++) {
                        vi_el[geo_dim * local_vid + i] = global_id + i;
                        fe_problem.x()[global_id + i] = mesh_.vertex(vh)[i];
                        vx_el[geo_dim * local_vid + i] = fe_problem.x()[global_id + i];
                    }
                }

                ELEN_B_3D_AD el_elem;
                auto el_e = el_elem.eval_f(vx_el, vc_el);
                max_el_e = std::max(max_el_e, el_e);
                log(LogLevel::DEBUG) << " - " << (mesh_.is_boundary(eh) ? "boundary" : "interior") << " edge " << edge_to_str(mesh_, eh) << " length: " << len
                                     << " rel length (to min boundary): " << (len / min_codomain_boundary_edge_length_) << ", target: " << el_target << " energy: " << el_e
                                     << std::endl;

                fe_edge_length.instances().add_element(vi_el, vc_el);
            }
        }

        fe_problem.add_set(&fe_edge_length);

        log(LogLevel::INFO) << " - edge length bound elements weight sum = "<<el_w_sum << std::endl;
        log(LogLevel::INFO) << " -                            min length = "<<min_len << std::endl;
        log(LogLevel::INFO) << " -                            max energy = "<<max_el_e << std::endl;
        log(LogLevel::INFO) << " -                 unweighted max energy = "<<(max_el_e / w) << std::endl;
        log(LogLevel::INFO) << " -        #backup edge length bounds used: "<<backup_bound_count << std::endl;
        log(LogLevel::INFO) << " ---------------------------------------------------------" << std::endl;
    }
    
    int constrained_count(0);

    // constrain the root face
    if (snap_phase_) {
        auto root_verts = mesh_.get_halfface_vertices(root_hfh_);
        if (!wrapped_prop_[root_verts[1]] || !wrapped_prop_[root_verts[2]]) {
            auto target_e1 = codomain_pos_prop_[root_verts[1]] - codomain_pos_prop_[root_verts[0]];
            auto target_e2 = codomain_pos_prop_[root_verts[2]] - codomain_pos_prop_[root_verts[0]];
            auto target_n = target_e1.cross(target_e2).normalized();
            log(LogLevel::INFO) << " - target e1: " << target_e1 << std::endl;
            log(LogLevel::INFO) << " - target e2: " << target_e2 << std::endl;
            log(LogLevel::INFO) << " - target  n: " << target_n << std::endl;


            if (!wrapped_prop_[root_verts[1]]) {
                log(LogLevel::INFO) << " --> root vertex " << root_verts[1] <<
                        " is not wrapped yet, constraining it to target e1" << std::endl;
                auto target_n2 = target_e1.cross(target_n).normalized();
                /*log(LogLevel::INFO)<<" - target  n2: "<<target_n2<<std::endl;

                auto p0 = mesh_.vertex(root_verts[0]);
                auto p1 = mesh_.vertex(root_verts[1]);
                auto e1 = p1 - p0;

                log(LogLevel::INFO)<<" - p0 = "<<p0<<std::endl;
                log(LogLevel::INFO)<<" - p1 = "<<p1<<std::endl;
                log(LogLevel::INFO)<<" - n . e1 = "<<target_n.dot(e1)<<std::endl;
                log(LogLevel::INFO)<<" - n . p1 = "<<target_n.dot(p1)<<std::endl;
                log(LogLevel::INFO)<<" - n . p0 = "<<target_n.dot(p0)<<std::endl;
                log(LogLevel::INFO)<<" - n . p1 - n . n0 = "<<(target_n.dot(p1) - target_n.dot(p0))<<std::endl;
                log(LogLevel::INFO)<<" - n2 . e1 = "<<target_n2.dot(e1)<<std::endl;
                log(LogLevel::INFO)<<" - n2 . p1 = "<<target_n2.dot(p1)<<std::endl;
                log(LogLevel::INFO)<<" - n2 . p0 = "<<target_n2.dot(p0)<<std::endl;
                log(LogLevel::INFO)<<" - n2 . p1 - n2 . n0 = "<<(target_n2.dot(p1) - target_n2.dot(p0))<<std::endl;*/

                int var_idx0 = geo_dim * global_vid_prop[root_verts[0]];
                int var_idx1 = geo_dim * global_vid_prop[root_verts[1]];

                COMISO::LinearConstraint::SVectorNC coeffs1, coeffs2;
                coeffs1.resize(n_vars);
                coeffs2.resize(n_vars);

                coeffs1.coeffRef(var_idx1 + 0) = target_n[0];
                coeffs1.coeffRef(var_idx1 + 1) = target_n[1];
                coeffs1.coeffRef(var_idx1 + 2) = target_n[2];

                coeffs1.coeffRef(var_idx0 + 0) = -target_n[0];
                coeffs1.coeffRef(var_idx0 + 1) = -target_n[1];
                coeffs1.coeffRef(var_idx0 + 2) = -target_n[2];

                coeffs2.coeffRef(var_idx1 + 0) = target_n2[0];
                coeffs2.coeffRef(var_idx1 + 1) = target_n2[1];
                coeffs2.coeffRef(var_idx1 + 2) = target_n2[2];

                coeffs2.coeffRef(var_idx0 + 0) = -target_n2[0];
                coeffs2.coeffRef(var_idx0 + 1) = -target_n2[1];
                coeffs2.coeffRef(var_idx0 + 2) = -target_n2[2];


                eq_constraints.push_back(COMISO::LinearConstraint(coeffs1, 0.0));
                eq_constraints.push_back(COMISO::LinearConstraint(coeffs2, 0.0));

                if (COMISO::LinearConstraint(coeffs1, 0.0).eval_constraint(fe_problem.x().data()) > 1e-7 ||
                    COMISO::LinearConstraint(coeffs2, 0.0).eval_constraint(fe_problem.x().data()) > 1e-7) {
                    log(LogLevel::ERROR) << " ERROR - root vertex " << root_verts[1] <<
                            " alignmnent constraints are not satisfied: " << std::endl;
                    std::cout << " -  n.e1 = " << COMISO::LinearConstraint(coeffs1, 0.0).eval_constraint(
                        fe_problem.x().data()) << std::endl;
                    std::cout << " - n2.e1 = " << COMISO::LinearConstraint(coeffs2, 0.0).eval_constraint(
                        fe_problem.x().data()) << std::endl;
                    return WRAPPING_ERROR;
                    }
            }

            if (!wrapped_prop_[root_verts[2]]) {
                log(LogLevel::INFO) << " --> root vertex " << root_verts[2] <<
                        " is not wrapped yet, constraining it to root plane" << std::endl;

                int var_idx0 = geo_dim * global_vid_prop[root_verts[0]];
                int var_idx2 = geo_dim * global_vid_prop[root_verts[2]];

                COMISO::LinearConstraint::SVectorNC coeffs;
                coeffs.resize(n_vars);

                coeffs.coeffRef(var_idx2 + 0) = target_n[0];
                coeffs.coeffRef(var_idx2 + 1) = target_n[1];
                coeffs.coeffRef(var_idx2 + 2) = target_n[2];

                coeffs.coeffRef(var_idx0 + 0) = -target_n[0];
                coeffs.coeffRef(var_idx0 + 1) = -target_n[1];
                coeffs.coeffRef(var_idx0 + 2) = -target_n[2];


                eq_constraints.push_back(COMISO::LinearConstraint(coeffs, 0.0));
                if (COMISO::LinearConstraint(coeffs, 0.0).eval_constraint(fe_problem.x().data()) > constr_eps) {
                    log(LogLevel::ERROR) << " ERROR - root vertex " << root_verts[2] <<
                            " planar constraint is not satisfied: " << std::endl;
                    std::cout << " - n.e2  constr eval = " << COMISO::LinearConstraint(coeffs, 0.0).eval_constraint(
                        fe_problem.x().data()) << std::endl;
                    return WRAPPING_ERROR;
                }
            }
        }
    }

    //constrain the wrapped vertices
    for(auto vh: mesh_.vertices()) {
        if(wrapped_prop_[vh]) {

            auto pos = codomain_pos_prop_[vh];

            int var_idx = geo_dim * global_vid_prop[vh];

            for(int i(0); i<geo_dim; i++) {
                COMISO::LinearConstraint::SVectorNC coeffs;
                coeffs.resize(n_vars);

                coeffs.coeffRef(var_idx + i) = 1;
                eq_constraints.push_back(COMISO::LinearConstraint(coeffs, -pos[i]));
            }
            constrained_count++;
        }

        //x-plane constraints for video and demos
#if 0
        std::cout<<" ERROR - CONSTRAINING TO X PLANE"<<std::endl;
        //constraining to current x plane, temp stuff
        int var_idx = geo_dim * global_vid_prop[vh];
        auto pos = mesh_.vertex(vh);

        for(int i(0); i<1; i++) {
            COMISO::LinearConstraint::SVectorNC coeffs;
            coeffs.resize(n_vars);

            coeffs.coeffRef(var_idx + i) = 1;
            eq_constraints.push_back(COMISO::LinearConstraint(coeffs, -pos[i]));
        }
        constrained_count++;
#endif

    }
    if (!constrained_count) {
            log(LogLevel::ERROR)<<" ERROR - No vertex constrained at all. You should constrain at least one "<<std::endl;
            return -1;
    }

    for (int i(0); i<eq_constraints.size(); i++) {
        double violation = eq_constraints[i].eval_constraint(fe_problem.x().data());
        if (violation > constr_eps) {
            log(LogLevel::ERROR) << "ERROR - violation for already-snapped vertex constraint "<<i<<" is "<<violation<<" > "<<constr_eps << std::endl;
            return WRAPPING_ERROR;
        }
    }


    if(snap_phase_ && enable_growing_front_eq_constraints){
        setup_growing_front_equality_constraints(fe_problem, global_vid_prop, eq_constraints);
    }

    std::vector<COMISO::NConstraintInterface*> constraints_pointers;
    constraints_pointers.reserve(constraints_pointers.size());
    for( auto& c : eq_constraints) {
        constraints_pointers.push_back(&c);
    }

    log(LogLevel::IMPORTANT) << " - total number of constraints pointers: " << constraints_pointers.size() << std::endl;
    /*auto res = COMISO::ConstraintTools::remove_dependent_linear_constraints_only_linear_equality(constraints_pointers);
    log(LogLevel::IMPORTANT) << "#eliminated: " << res.n_constraints_eliminated << ", #infeasible detected: "
                             << res.n_infeasible_detected << std::endl;*/


    log(LogLevel::INFO) << " energies before optimization: " << std::endl;
    if (LogConfig::get_min_level() <= LogLevel::INFO) {
        fe_problem.print_objectives();
    }

    double initial_energy = fe_problem.eval_f(fe_problem.x().data());
    if(initial_energy == std::numeric_limits<double>::infinity()){
        log(LogLevel::INFO)<<" -- infinite initial energy, stopping there."<<std::endl;
        return -1;
    }

    energy_record.initial_quality_energy_ = fe_barrier.eval_f(fe_problem.x().data());
    energy_record.initial_virtual_energy_ = {
        fe_virtual_A.eval_f(fe_problem.x().data()),
        fe_virtual_B.eval_f(fe_problem.x().data()),
        fe_virtual_C.eval_f(fe_problem.x().data())
    };

    bool feasible_solution_found(false);

    COMISO::TruncatedNewtonPCG::Config tn_config;

    tn_config.eps = json_data_["global_wrapping_convergence_eps"];
    tn_config.always_update_preconditioner = true;
    tn_config.max_iters = newton_iters;
    tn_config.max_pcg_iters = json_data_["global_wrapping_max_pcg_iters"];;
    tn_config.allow_warmstart = false;
    tn_config.eps_gdx = json_data_["global_wrapping_tolerance_gdx"];
    tn_config.silent = LogConfig::get_min_level()> LogLevel::INFO;
    tn_config.eps_constraints_violation = constr_eps;
    tn_config.min_acceptable_rel_constraint_violation_decrease = json_data_["global_wrapping_constr_violation_progress_eps"];
    tn_config.min_acceptable_rel_objective_decrease = json_data_["global_wrapping_energy_progress_eps"];


    COMISO::TruncatedNewtonPCG tn(tn_config);


    opti_stats_.total_upkeep_time_s_ += main_sw_.lap_duration();

    tn.solve_projected_normal_equation(&fe_problem, constraints_pointers);
    converged = tn.converged();
    feasible_solution_found = tn.feasible_solution_found();

    auto opti_time = main_sw_.lap_duration();
    opti_stats_.total_opti_time_s_ += opti_time;

    log(LogLevel::INFO) << " energies after optimization: " << std::endl;
    if (LogConfig::get_min_level() <= LogLevel::INFO) {
        fe_problem.print_objectives();
    }

    energy_record.final_quality_energy_ = fe_barrier.eval_f(fe_problem.x().data());
    energy_record.final_virtual_energy_ = {
        fe_virtual_A.eval_f(fe_problem.x().data()),
        fe_virtual_B.eval_f(fe_problem.x().data()),
        fe_virtual_C.eval_f(fe_problem.x().data())
    };

    log(LogLevel::INFO) << " -------------------------------------------" << std::endl;
    log(LogLevel::INFO) << " - computing quality energy gradient..." << std::endl;
    quality_g.resize(fe_problem.x().size());
    quality_g.setZero();
    fe_barrier.accumulate_gradient(fe_problem.x().data(), quality_g.data());
    log(LogLevel::INFO) << " -> quality energy gradient norm = "<<quality_g.norm() << std::endl;
    energy_record.final_quality_gradient_norm_ = quality_g.norm();
    energy_record.final_quality_gradient_max_coeff_ = quality_g.maxCoeff();


    log(LogLevel::INFO) << " -------------------------------------------" << std::endl;
    log(LogLevel::INFO) << " - computing virtual energy gradient..." << std::endl;
    for (int i(0); i < virtual_gs.size(); i++) {
        virtual_gs[i].resize(fe_problem.x().size());
        virtual_gs[i].setZero();
    }
    fe_virtual_A.accumulate_gradient(fe_problem.x().data(), virtual_gs[0].data());
    fe_virtual_B.accumulate_gradient(fe_problem.x().data(), virtual_gs[1].data());
    fe_virtual_C.accumulate_gradient(fe_problem.x().data(), virtual_gs[2].data());
    log(LogLevel::INFO) << " ->   primary virtual energy gradient norm = "<<virtual_gs[0].norm() << std::endl;
    log(LogLevel::INFO) << " -> secondary virtual energy gradient norm = "<<virtual_gs[1].norm() << std::endl;
    log(LogLevel::INFO) << " ->  tertiary virtual energy gradient norm = "<<virtual_gs[2].norm() << std::endl;
    energy_record.final_virtual_gradient_norm_ = {virtual_gs[0].norm(), virtual_gs[1].norm(), virtual_gs[2].norm()};
    energy_record.final_virtual_gradient_max_coeff_ = {virtual_gs[0].maxCoeff(), virtual_gs[1].maxCoeff(), virtual_gs[2].maxCoeff()};

    log(LogLevel::INFO) << " -------------------------------------------" << std::endl;
    log(LogLevel::INFO) << " - computing total virtual energy gradient..." << std::endl;
    Eigen::VectorXd total_g;
    total_g.resize(fe_problem.x().size());
    total_g.setZero();

    fe_barrier.accumulate_gradient(fe_problem.x().data(), total_g.data());
    fe_virtual_A.accumulate_gradient(fe_problem.x().data(), total_g.data());
    fe_virtual_B.accumulate_gradient(fe_problem.x().data(), total_g.data());
    fe_virtual_C.accumulate_gradient(fe_problem.x().data(), total_g.data());
    energy_record.final_total_gradient_max_coeff_ = total_g.maxCoeff();
    log(LogLevel::INFO) << " -> max abs total gradient coeff = "<<energy_record.final_total_gradient_max_coeff_ << std::endl;
    log(LogLevel::INFO) << " ->          total gradient norm = "<<total_g.norm() << std::endl;


    log(LogLevel::INFO) << " -------------------------------------------" << std::endl;
    log(LogLevel::INFO) << " - computing total virtual energy Hessian..." << std::endl;
    energy_record.final_total_Hessian_max_coeff_ = 0.0;
    std::vector<Eigen::Triplet<double>> triplets;
    fe_barrier.accumulate_hessian(fe_problem.x().data(), triplets);
    fe_virtual_A.accumulate_hessian(fe_problem.x().data(), triplets);
    fe_virtual_B.accumulate_hessian(fe_problem.x().data(), triplets);
    fe_virtual_C.accumulate_hessian(fe_problem.x().data(), triplets);
    for(auto& t: triplets){
        energy_record.final_total_Hessian_max_coeff_ = std::max(energy_record.final_total_Hessian_max_coeff_, std::abs(t.value()));
    }
    log(LogLevel::INFO) << " - max abs virtual elements Hessian coeff: "<<energy_record.final_total_Hessian_max_coeff_ << std::endl;
    log(LogLevel::INFO) << " ----------------------------------------------------------" << std::endl;


    if(!mock_run) {
        std::vector<double> unweighted_final_virtual_energies(fe_virtual_A.n_instances());
        fe_virtual_A.eval_f_per_element(fe_problem.x().data(), unweighted_final_virtual_energies.data());
        for(auto& e: unweighted_final_virtual_energies){
            e/= virtual_weight_scaling_.primary;
        }

        if (enable_full_history) {
            opti_stats_.final_virtual_energies_history_.push_back(unweighted_final_virtual_energies);
        }
        opti_stats_.newton_iters_history_.push_back(tn.n_iterations_used());
        opti_stats_.opti_time_s_history_.push_back(opti_time);

        log(LogLevel::INFO) << " --- updating mesh..." << std::endl;
        //update the mesh
        for (auto vh: mesh_.vertices()) {

            if (wrapped_prop_[vh]) {
                mesh_.set_vertex(vh, codomain_pos_prop_[vh]);
            } else {
                int global_id = geo_dim * global_vid_prop[vh];
                Vec3d new_pos;
                for (int i(0); i < 3; i++) {
                    new_pos[i] = fe_problem.x()[global_id + i];
                }

                double diff = (mesh_.vertex(vh) - new_pos).norm();
                mesh_.set_vertex(vh, new_pos);
            }
        }


        if(json_data_["enable_vertex_trajectory_mesh"]){
            for (auto vh: mesh_.vertices()) {

                if(mesh_.is_boundary(vh)){
                    boundary_vertices_trajectories_[vh.idx()].push_back(mesh_.vertex(vh));
                }
            }
        }


        //and update the virtual vertices positions
        int virtual_vertex_idx(0);
        for (auto eh: mesh_.edges()) {
            if (mesh_.is_boundary(eh)) {
                int global_id = virtual_vertex_idx_start + virtual_vertex_idx * 3;

                Vec3d pos;
                for (int i(0); i < 3; i++) {
                    pos[i] = fe_problem.x()[global_id + i];
                }

                virtual_vertex_pos_prop_[eh] = pos;
                //log(LogLevel::INFO) << " - updated virtual vertex prop position to: "<<pos << std::endl;

                virtual_vertex_idx++;
            }
        }
    }else{
        log(LogLevel::INFO) << " --- skipping mesh update..." << std::endl;

    }


    log(LogLevel::IMPORTANT)<<" ================================================================================================================================ "<<std::endl;

    return 0;
}



template <typename DistortionElemT, typename VirtualElemT>
void
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
setup_growing_front_equality_constraints(COMISO::FiniteElementProblem& fe_problem,
                                         const OpenVolumeMesh::VertexPropertyT<int>& global_vid_prop,
                                         std::vector<COMISO::LinearConstraint>& eq_constraints){

    log(LogLevel::INFO)<<" ------------------------------------------------ "<<std::endl;
    log(LogLevel::INFO)<<" adding equality constraints for the non-snapped vertices neighbors to the snapping front... "<<std::endl;

    int dihedral_diff_above_180_count(0);
    const double dihedral_diff_thresh_deg = json_data_["global_wrapping_growing_front_equality_constraints_dihedral_diff_thresh_deg"];
    int constrained_count(0);
    for(auto vh: mesh_.vertices()){
        if(mesh_.is_boundary(vh) && !wrapped_prop_[vh]){
            bool found_snapped_neighbor(false);
            for(auto vv_it = mesh_.vv_iter(vh); vv_it.valid(); vv_it++){
                if(mesh_.is_boundary(*vv_it) && wrapped_prop_[*vv_it]){
                    found_snapped_neighbor = true;
                    break;
                }
            }

            if(found_snapped_neighbor){

                log(LogLevel::DEBUG)<<" ----- checking if possible to add constraints for vertex "<<vh<<" ... "<<std::endl;

                double max_dihedral_diff_deg(0.0);
                for(auto vhf_it = mesh_.vhf_iter(vh); vhf_it.valid(); vhf_it++){

                    if(mesh_.is_boundary(*vhf_it)) {

                        auto op_eh = opposite_edge(*vhf_it, vh);

                        auto diff_deg = std::abs(compute_boundary_dihedral_angle_deg(mesh_, op_eh) - boundary_dihedral_angle_deg_prop_[op_eh]);

                        max_dihedral_diff_deg = std::max(max_dihedral_diff_deg, diff_deg);

                        log(LogLevel::DEBUG)<<" -- dihedral diff for edge "<<edge_to_str(mesh_, op_eh)<<": "<<diff_deg<<std::endl;
                    }
                }
                log(LogLevel::DEBUG)<<" -> max dihedral diff: "<<max_dihedral_diff_deg<<std::endl;



                if(dihedral_diff_thresh_deg >= 360.0 || max_dihedral_diff_deg < dihedral_diff_thresh_deg) {
                    log(LogLevel::DEBUG)<<" --> below threshold ("<<dihedral_diff_thresh_deg<<"), adding constraint"<<std::endl;

                    auto pos = codomain_pos_prop_[vh];

                    int var_idx = 3 * global_vid_prop[vh];

                    for (int i(0); i < 3; i++) {
                        COMISO::LinearConstraint::SVectorNC coeffs;
                        coeffs.resize(fe_problem.n_unknowns());

                        coeffs.coeffRef(var_idx + i) = 1;
                        eq_constraints.push_back(COMISO::LinearConstraint(coeffs, -pos[i]));
                    }
                    constrained_count++;
                }

                if(max_dihedral_diff_deg >= 180.0){
                    dihedral_diff_above_180_count++;
                }
            }
        }
    }

    log(LogLevel::INFO)<<" -> constrained "<<constrained_count<<" boundary vertices neighbor to the snapping front"<<std::endl;
    log(LogLevel::INFO)<<" -> "<<dihedral_diff_above_180_count<<" of them have at least one hinge edge with dihedral diff above 180.0"<<std::endl;
    log(LogLevel::INFO)<<" ------------------------------------------------ "<<std::endl;
}


template <typename DistortionElemT, typename VirtualElemT>
OpenVolumeMesh::HalfFaceHandle
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
find_best_root_halfface(int ring_k) const {

    LightWeightStopWatch sw;
    log(LogLevel::INFO)<<" ------------------------------------------------ "<<std::endl;
    log(LogLevel::INFO)<<" - looking for best root halfface by checking the average dihedral angle diff in their "<<ring_k<<"-ring..."<<std::endl;

    FaceHandle best_fh(-1);
    double min_avg_dihedral_diff(std::numeric_limits<double>::max());

    for(auto fh: mesh_.faces()) {
        if(mesh_.is_boundary(fh)) {
            if(ring_k <= 0){
                log(LogLevel::INFO)<<" - k <= 0, just taking the first boundary halfface instead"<<std::endl;
                best_fh = fh;
                break;
            }
            auto hfh = mesh_.halfface_handle(fh, 0);
            std::set<EdgeHandle> to_check;
            std::set<VertexHandle> current_layer;
            for(auto vh: mesh_.get_halfface_vertices(hfh)){
                current_layer.insert(vh);
            }

            for(int i(0); i<ring_k; i++){
                std::set<VertexHandle> next_layer;
                for(auto vh: current_layer){
                    for(auto out_he: mesh_.outgoing_halfedges(vh)){
                        if(mesh_.is_boundary(mesh_.edge_handle(out_he))){
                            next_layer.insert(mesh_.to_vertex_handle(out_he));

                            to_check.insert(mesh_.edge_handle(out_he));
                        }
                    }
                }

                current_layer = next_layer;
            }

            double dihedral_diff_sum(0);
            for(auto eh: to_check){
                double domain_dihedral_angle = compute_boundary_dihedral_angle_deg(ref_mesh_, eh);
                double codomain_dihedral_angle = boundary_dihedral_angle_deg_prop_[eh];

                auto from_v = mesh_.edge(eh).from_vertex();
                auto   to_v = mesh_.edge(eh).to_vertex();
                auto target_length = (codomain_pos_prop_[from_v] - codomain_pos_prop_[to_v]).norm();

                //target_length = 1.0;

                dihedral_diff_sum += std::abs(domain_dihedral_angle - codomain_dihedral_angle) / target_length;
            }
            dihedral_diff_sum /= to_check.size();

            //log(LogLevel::INFO)<<" -> checked "<<to_check.size()<<" edges in the "<<ring_k<<"-neighborhood of face "<<halfface_to_str(mesh_, hfh)<<" and its average dihedral diff is "<<dihedral_diff_sum<<std::endl;

            if(dihedral_diff_sum < min_avg_dihedral_diff){
                min_avg_dihedral_diff = dihedral_diff_sum;
                best_fh = fh;
                log(LogLevel::INFO)<<" --> updated best face to "<<halfface_to_str(mesh_, hfh)<<" with target-length weighted average dihedral diff: "<<dihedral_diff_sum<<std::endl;
            }
        }
    }

    auto first_hfh = mesh_.halfface_handle(best_fh, 0);

    log(LogLevel::INFO)<<" --> best face is "<<halfface_to_str(mesh_, first_hfh)<<" with average dihedral diff: "<<min_avg_dihedral_diff<<" in "<<sw.total_duration()<<" seconds"<<std::endl;
    log(LogLevel::INFO)<<" ------------------------------------------------ "<<std::endl;

    return mesh_.is_boundary(first_hfh) ? first_hfh : mesh_.opposite_halfface_handle(first_hfh);

}

template <typename DistortionElemT, typename VirtualElemT>
void
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
setup_codomain_boundary_props(){
    log(LogLevel::INFO)<<" ------------------------------------------------ "<<std::endl;
    log(LogLevel::INFO)<<" setting up codomain boundary props (dihedral angles and edge lengths)  "<<std::endl;

    double min_dihedral_angle_deg(std::numeric_limits<double>::max());
    double max_dihedral_angle_deg(0);
    double min_length(std::numeric_limits<double>::max());
    double max_length(0);

    //log(LogLevel::INFO)<<" initializing boundary dihedral angles prop..."<<std::endl;
    for(auto eh: mesh_.edges()){
        if (mesh_.is_boundary(eh)) {

            boundary_dihedral_angle_deg_prop_[eh] = compute_boundary_dihedral_angle_deg(mesh_, eh);
            boundary_edge_length_prop_[eh] = (mesh_.vertex(mesh_.edge(eh).from_vertex()) - mesh_.vertex(mesh_.edge(eh).to_vertex())).norm();
            //log(LogLevel::INFO)<<" - boundary edge "<<eh<<": ("<<mesh_.edge(eh).from_vertex()<<", "<<mesh_.edge(eh).to_vertex()<<") angle = "<<boundary_dihedral_angle_deg_prop_[eh]<<std::endl;

            min_dihedral_angle_deg = std::min(min_dihedral_angle_deg, boundary_dihedral_angle_deg_prop_[eh]);
            max_dihedral_angle_deg = std::max(max_dihedral_angle_deg, boundary_dihedral_angle_deg_prop_[eh]);

            min_length = std::min(min_length, boundary_edge_length_prop_[eh]);
            max_length = std::max(max_length, boundary_edge_length_prop_[eh]);

        }else {
            boundary_dihedral_angle_deg_prop_[eh] = 0.0;
            boundary_edge_length_prop_[eh] = 0.0;
        }
        boundary_dihedral_angle_diff_deg_prop_[eh] = 0.0;
    }

    log(LogLevel::INFO)<<" -> min dihedral angle: "<<min_dihedral_angle_deg<<"°"<<std::endl;
    log(LogLevel::INFO)<<" -> max dihedral angle: "<<max_dihedral_angle_deg<<"°"<<std::endl;
    log(LogLevel::INFO)<<" ->    min edge length: "<<min_length<<std::endl;
    log(LogLevel::INFO)<<" ->    max edge length: "<<max_length<<std::endl;

    log(LogLevel::INFO)<<" ------------------------------------------------ "<<std::endl;

}


template <typename DistortionElemT, typename VirtualElemT>
void
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
setup_min_IVL_prop(){

    double max_min_IVL(0.0);
    double print_thresh(1e-1);
    log(LogLevel::INFO)<<" initializing min IVL energy prop..."<<std::endl;
    log(LogLevel::INFO)<<"(only printing those above "<<print_thresh<<")"<<std::endl;

    for(auto ch: mesh_.cells()){

        int boundary_v_count = count_boundary_vertices(mesh_, ch);

        //NOTE: using the reference mesh so the boundary constraints are already satisfied
        min_IVL_prop_[ch] = compute_min_IVL_energy(ref_mesh_, ch);
        if (min_IVL_prop_[ch] > print_thresh) {
            log(LogLevel::DEBUG)<<" - cell "<<cell_to_str(mesh_, ch)<<" with "<<boundary_v_count<<" boundary vertices min IVL energy = "<<min_IVL_prop_[ch]<<std::endl;
        }

        if(boundary_v_count <= 1 && min_IVL_prop_[ch] != 0.0){
            log(LogLevel::ERROR)<<" ERROR - cell has "<<boundary_v_count<<" boundary vertices but min IVL energy is "<<min_IVL_prop_[ch]<<std::endl;
            exit(EXIT_FAILURE);
        }

        max_min_IVL = std::max(max_min_IVL, min_IVL_prop_[ch]);
    }
    log(LogLevel::INFO)<<" --> max min IVL energy = "<<max_min_IVL<<std::endl;

}


template <typename DistortionElemT, typename VirtualElemT>
void
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
setup_codomain_average_incident_edge_length_prop(){

    double max_min_IVL(0.0);
    double print_thresh(1e-1);

    log(LogLevel::INFO)<<" --------------------------------------------------------------- "<<std::endl;
    log(LogLevel::INFO)<<" initializing average codomain incident edge length prop..."<<std::endl;
    for(auto vh: mesh_.vertices()){
        if(mesh_.is_boundary(vh)){
            average_codomain_incident_edge_length_prop_[vh] = compute_average_incident_edge_length(mesh_, vh, true);
            //log(LogLevel::INFO)<<" - v"<<vh<<": "<<average_codomain_incident_edge_length_prop_[vh]<<std::endl;
        }else{
            average_codomain_incident_edge_length_prop_[vh] = 0.0;
        }
    }
    log(LogLevel::INFO)<<" --------------------------------------------------------------- "<<std::endl;

}



template <typename DistortionElemT, typename VirtualElemT>
void
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
setup_virtual_vertices_positions() {

    for (auto eh: mesh_.edges()) {
        setup_virtual_vertex_position(eh);
    }
}


template <typename DistortionElemT, typename VirtualElemT>
void
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
setup_virtual_vertex_position(const OVM::EdgeHandle& eh){
    if (mesh_.is_boundary(eh)) {

        std::vector<OVM::HalfFaceHandle> incident_boundary_halffaces;
        for (auto ehf_it = mesh_.ehf_iter(eh); ehf_it.valid(); ehf_it++) {
            if (mesh_.is_boundary(*ehf_it)) {
                incident_boundary_halffaces.push_back(*ehf_it);
            }
        }

        auto &left_hf = incident_boundary_halffaces[0];
        auto &right_hf = incident_boundary_halffaces[1];

        auto x = compute_optimal_virtual_vertex_pos(mesh_, eh, {left_hf, right_hf});
        //auto x = compute_rotational_mid_point(mesh_, eh, {left_hf, right_hf});
        //log(LogLevel::INFO) << " - rotational mid-point for edge "<<edge_to_str(mesh_, eh)<<": "<<x << std::endl;


        virtual_vertex_pos_prop_[eh] = x;

    }else {
        virtual_vertex_pos_prop_[eh] = {0,0,0};
    }
}


template <typename DistortionElemT, typename VirtualElemT>
int
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
try_to_snap_next_candidates() {


    log(LogLevel::INFO)<<" ------------------------------------------------------------ "<<std::endl;
    log(LogLevel::INFO)<<" - trying to snap non-wrapped candidates... "<<std::endl;

    const bool enable_min_vol_check = json_data_["global_wrapping_IVL_vol_w"] > 0;
    const double min_rel_vol_bound = json_data_["global_wrapping_IVL_rel_vol_bound"];
    //const double avg_vol = enable_min_vol_check ? 0.0 : compute_total_volume(mesh_) / mesh_.n_logical_cells();
    const double avg_vol = enable_min_vol_check ? 0.0 : target_total_volume_ / mesh_.n_logical_cells();
    const double min_vol_bound = min_rel_vol_bound * avg_vol;


    const bool enable_min_el_check = json_data_["global_wrapping_edge_length_bound_elements_w"] > 0;
    const double min_rel_el_bound = json_data_["global_wrapping_edge_length_bound_elements_boundary_rel_bound"];
    const double min_el_bound = min_rel_el_bound * min_codomain_boundary_edge_length_;

    int total_snapped_count(0);
    int snapped_this_pass_count(0);
    int holonomy_fail_count(0);
    int orientation_fail_count(0);
    int quality_bound_fail_count(0);
    int volume_bound_fail_count(0);
    int edge_length_fail_count(0);

    double quality_bound = json_data_["global_wrapping_IVL_bound"];
    double bound_safety_factor = 0.99;
    log(LogLevel::INFO)<<" - snapping bound: "<<quality_bound<<std::endl;

    DistortionElemT elem;
    typename DistortionElemT::VecC vc;
    typename DistortionElemT::VecV vx;
    vc<<1,0,0;

    if(DistortionElemT::element_code() != "IVL_B"){
        log(LogLevel::ERROR)<<" ERROR - snapping is not implemented for distortion element "<<DistortionElemT::element_name()<<std::endl;
        return WRAPPING_ERROR;
    }

    int new_edges_start_idx(mesh_.n_edges());
    std::set<EdgeHandle> impacted_edges;
    int pass_it(0);
    do {
        snapped_this_pass_count = 0;
        std::vector<VertexHandle> to_wrap;
        auto tried_candidate_prop = mesh_.request_vertex_property<bool>("", false);
        for (auto vh: mesh_.vertices()) {
            if (mesh_.is_boundary(vh) && wrapped_prop_[vh]) {
                for (auto out_he: mesh_.outgoing_halfedges(vh)) {
                    if (mesh_.is_boundary(mesh_.edge_handle(out_he))) {
                        auto neighbor = mesh_.to_vertex_handle(out_he);

                        if (mesh_.is_boundary(neighbor) &&
                            !wrapped_prop_[neighbor] &&
                            !tried_candidate_prop[neighbor]) {
                            //log(LogLevel::INFO)<<" ---- trying to snap candidate neighbor "<<neighbor<<" from "<<mesh_.vertex(neighbor)<<" to "<<codomain_pos_prop_[neighbor]<<" (snapping distance = "<<(mesh_.vertex(neighbor) - codomain_pos_prop_[neighbor]).norm()<<")"<<std::endl;

                            auto current_pos = mesh_.vertex(neighbor);
                            mesh_.set_vertex(neighbor, codomain_pos_prop_[neighbor]);

                            bool found_bad_tet(false);
                            bool found_bad_quality_tet(false);
                            bool found_bad_holonomy(false);
                            bool found_too_small_tet(false);
                            bool found_too_short_edge(false);

                            for (auto vc_it = mesh_.vc_iter(neighbor); vc_it.valid(); vc_it++) {
                                if (CGAL_wrapper::is_degenerate(mesh_, *vc_it) ||
                                    CGAL_wrapper::is_flipped(mesh_, *vc_it)) {
                                    found_bad_tet = true;
                                    //log(LogLevel::INFO) <<" ----->    snapping would create tet "<<cell_to_str(mesh_, *vc_it)<<" with volume "<<CGAL_wrapper::exact_volume(mesh_, *vc_it).to_double()<<std::endl;
                                    orientation_fail_count++;
                                    break;
                                }

                                if (enable_min_vol_check && compute_volume(mesh_, *vc_it) <= min_vol_bound) {
                                    found_too_small_tet = true;
                                    volume_bound_fail_count++;
                                    break;
                                }
                            }

                            if (!found_bad_tet) {
                                for (auto vc_it = mesh_.vc_iter(neighbor); vc_it.valid(); vc_it++) {

                                    auto verts = mesh_.get_cell_vertices(*vc_it);
                                    for (int i(0); i<verts.size(); i++) {
                                        auto pos = mesh_.vertex(verts[i]);
                                        for (int j(0); j<3; j++) {
                                            vx[3 * i + j] = pos[j];
                                        }
                                    }

                                    auto e = elem.eval_f(vx, vc);
                                    if (e >= (quality_bound * bound_safety_factor)) {
                                        found_bad_quality_tet = true;
                                        log(LogLevel::INFO) <<" ----->    snapping would create tet "<<cell_to_str(mesh_, *vc_it)<<" with energy = "<<e<<std::endl;
                                        quality_bound_fail_count++;
                                    }
                                }

                                if (!found_bad_quality_tet) {
                                    for (auto vf_it = mesh_.vf_iter(neighbor); vf_it.valid(); vf_it++) {
                                        if (mesh_.is_boundary(*vf_it)) {
                                            for (auto fe_it = mesh_.fe_iter(*vf_it); fe_it.valid(); fe_it++) {
                                                if (mesh_.is_boundary(*fe_it)) {

                                                    //log(LogLevel::INFO) << " -    checking holonomy of impacted edge "<< edge_to_str(mesh_, *fe_it) << std::endl;
                                                    if (!valid_holonomy(mesh_, *fe_it)) {
                                                        found_bad_holonomy = true;
                                                        //log(LogLevel::INFO) << " ----->    snapping would break holonomy of edge "<< edge_to_str(mesh_, *fe_it) << std::endl;
                                                        holonomy_fail_count++;
                                                        break;
                                                    }
                                                }
                                            }
                                        }


                                        if(!found_bad_holonomy) {
                                            for (auto neighbor_out_he_it: mesh_.outgoing_halfedges(neighbor)) {
                                                auto neighbor2 = mesh_.to_vertex_handle(neighbor_out_he_it);
                                                double len = (mesh_.vertex(neighbor2) - mesh_.vertex(neighbor)).norm();

                                                if (len < min_el_bound) {
                                                    found_too_short_edge = true;
                                                    edge_length_fail_count++;
                                                    break;
                                                }
                                            }
                                        }
                                    }


#warning TODO: check angle energy
                                    //setup the energy
                                    /*std::set<EdgeHandle> local_impacted_edges;
                                    for(auto vf_it = mesh_.vf_iter(neighbor); vf_it.valid(); vf_it++){
                                        for(auto fe_it = mesh_.fe_iter(*vf_it); fe_it.valid(); fe_it++){
                                            if(mesh_.is_boundary(*fe_it)) {
                                                local_impacted_edges.insert();
                                            }
                                        }
                                    }*/
                                }
                            }


                            if (!found_bad_tet && !found_bad_holonomy && !found_bad_quality_tet && !found_too_small_tet && !found_too_short_edge) {
                                log(LogLevel::DEBUG) << " -> snapped vertex " << (neighbor) << " at " << current_pos <<
                                        " to its target position " << codomain_pos_prop_[neighbor]
                                        << " (snapping distance = " << (current_pos - codomain_pos_prop_[neighbor]).
                                        norm() << std::endl;
                                to_wrap.push_back(neighbor);
                            }else{
                                mesh_.set_vertex(neighbor, current_pos);
                            }

                            /*if(check_boundary_holonomy(mesh_)){
                                log(LogLevel::ERROR)<<" ERROR - holonomy issue after snapping vertex "<<neighbor<<std::endl;
                                return WRAPPING_HOLONOMY_ERROR;
                            }*/

                            tried_candidate_prop[neighbor] = true;
                        }
                    }
                }
            } /*else if(wrapped_prop_[vh]){
                mesh_.set_vertex(vh, codomain_pos_prop_[vh]);
            }*/
        }

        for(auto vh: to_wrap){
            mark_as_wrapped(vh);
            snapped_this_pass_count++;
            total_snapped_count++;

            for(auto vf_it = mesh_.vf_iter(vh); vf_it.valid(); vf_it++){
                for(auto fe_it = mesh_.fe_iter(*vf_it); fe_it.valid(); fe_it++){
                    if(mesh_.is_boundary(*fe_it)) {
                        impacted_edges.insert(*fe_it);
                    }
                }
            }
        }
        log(LogLevel::DEBUG)<<" -> snapped at pass "<<std::setw(3)<<pass_it<<": "<<std::setw(5)<<snapped_this_pass_count<<std::endl;
        log(LogLevel::DEBUG)<<" ------------------------------------------ "<<std::endl;
        pass_it++;

    }while (snapped_this_pass_count);

    double max_diff(0);
    int recomputed_count(0);
    for (auto eh: impacted_edges) {
        if (!mesh_.is_deleted(eh)) {

            std::vector<OVM::HalfFaceHandle> incident_boundary_halffaces;
            for (auto ehf_it = mesh_.ehf_iter(eh); ehf_it.valid(); ehf_it++) {
                if (mesh_.is_boundary(*ehf_it)) {
                    incident_boundary_halffaces.push_back(*ehf_it);
                }
            }

            auto backup_x = compute_optimal_virtual_vertex_pos(mesh_, eh, {incident_boundary_halffaces[0], incident_boundary_halffaces[1]});
            auto explicit_x = compute_rotational_mid_point(mesh_, eh, {incident_boundary_halffaces[0], incident_boundary_halffaces[1]});
            virtual_vertex_pos_prop_[eh] = backup_x;

            max_diff = std::max(max_diff, (backup_x - explicit_x).norm());
            log(LogLevel::DEBUG)<<" - recomputed virtual vertex position for edge "<<edge_to_str(mesh_, eh)
            <<" incident to faces "<<halfface_to_str(mesh_, incident_boundary_halffaces[0])
            <<" and "<<halfface_to_str(mesh_, incident_boundary_halffaces[1])
            <<" as "<<backup_x<<" (diff to fof version = "<<(explicit_x - backup_x).norm()<<")"<<std::endl;

            recomputed_count++;
        }
    }

    log(LogLevel::INFO)<<" -> done with snapping! Snapped "<<total_snapped_count<<" vertices in total in "<<pass_it<<" passes "<<std::endl;
    log(LogLevel::INFO)<<" ->            #orientation snapping failures: "<<orientation_fail_count<<std::endl;
    log(LogLevel::INFO)<<" ->               #holonomy snapping failures: "<<holonomy_fail_count<<std::endl;
    log(LogLevel::INFO)<<" ->          #quality bound snapping failures: "<<quality_bound_fail_count<<std::endl;
    log(LogLevel::INFO)<<" ->           #volume bound snapping failures: "<<volume_bound_fail_count<<std::endl;
    log(LogLevel::INFO)<<" ->      #edge length bound snapping failures: "<<edge_length_fail_count<<std::endl;
    log(LogLevel::INFO)<<" ->            #recomputed virtual vertex pos: "<<recomputed_count<<std::endl;
    log(LogLevel::INFO)<<" ->  max diff between new and old updated pos: "<<max_diff<<std::endl;
    log(LogLevel::INFO)<<" ->                     updated wrapped count: "<<wrapped_count_<<std::endl;
    log(LogLevel::INFO)<<" ------------------------------------------------------------ "<<std::endl;


    if (run_safety_checks("snapping")) {
        return -1;
    }

    return total_snapped_count;
}


template <typename DistortionElemT, typename VirtualElemT>
int
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
translate_and_rotate_mesh_to_codomain(const OVM::HalfFaceHandle& hfh) {

    log(LogLevel::IMPORTANT)<<" ----------------------------------------------------------------------------------- "<<std::endl;
    log(LogLevel::IMPORTANT)<<" translating and rotating mesh to codomain... "<<std::endl;

    auto verts = mesh_.get_halfface_vertices(hfh);
    /*std::vector<Vec3d> original_pos = {
        mesh_.vertex(verts[0]),
        mesh_.vertex(verts[1]),
        mesh_.vertex(verts[2])
    };
    mesh_.set_vertex(verts[0], codomain_pos_prop_[verts[0]]);
    mesh_.set_vertex(verts[1], codomain_pos_prop_[verts[1]]);
    mesh_.set_vertex(verts[2], codomain_pos_prop_[verts[2]]);*/

    auto p0 = codomain_pos_prop_[verts[0]];
    auto p1 = codomain_pos_prop_[verts[1]];
    auto p2 = codomain_pos_prop_[verts[2]];

    auto n = (p1 - p0).cross(p2 - p0).normalized();

    auto e1 = (p1 - p0).normalized();

    if (translate_and_rotate_to_target(mesh_, root_hfh_, p0, n, e1)) {
        log(LogLevel::ERROR) << " ERROR - couldn't translate and rotate mesh properly" << std::endl;
        return -1;
    }

    log(LogLevel::IMPORTANT)<<" ----------------------------------------------------------------------------------- "<<std::endl;
    return 0;
}



template <typename DistortionElemT, typename VirtualElemT>
void
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
compute_element_quality(std::vector<double>& quality_vector, double& min_quality) const {

    quality_vector.clear();
    min_quality = std::numeric_limits<double>::max();
    for (auto ch: mesh_.cells()) {
        auto qual = mesh_.compute_quality_vector(ch, VOLUME_LENGTH, CONFORMING)[0];
        min_quality = std::min(min_quality, qual);
        quality_vector.push_back(qual);
    }
}


template <typename DistortionElemT, typename VirtualElemT>
double
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
compute_min_relative_volume() const {

    double min_rel_vol(std::numeric_limits<double>::max());
    double avg_vol = target_total_volume_ / mesh_.n_logical_cells();

    //double avg_vol = compute_total_volume(mesh_) / mesh_.n_logical_cells();
    //std::cout<<" - average vol = "<<avg_vol<<std::endl;

    for (auto ch: mesh_.cells()) {
        double rel_vol = compute_volume(mesh_, ch) / avg_vol;
        min_rel_vol = std::min(rel_vol, min_rel_vol);
        //std::cout<<" -- rel vol for cell "<<cell_to_str(mesh_, ch)<<" = "<<rel_vol<<std::endl;

    }

    return min_rel_vol;
}

template <typename DistortionElemT, typename VirtualElemT>
void
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
compute_boundary_dihedral_angle_diff(const std::vector<OVM::HalfFaceHandle>& to_wrap_hfs,
                                     const OVM::VertexHandle& to_wrap_vh,
                                     std::vector<double>& angles_diff_deg,
                                     double& max_abs_angle_diff_deg) const{

    max_abs_angle_diff_deg = 0;

    log(LogLevel::INFO)<<" ------------ "<<std::endl;
    log(LogLevel::INFO)<<" computing boundary dihedral angle diff... "<<std::endl;

    for(auto hf: to_wrap_hfs){
        log(LogLevel::INFO)<<" - handling halfface "<<halfface_to_str(mesh_, hf)<<std::endl;

        auto op_eh = opposite_edge(hf, to_wrap_vh);
        log(LogLevel::INFO)<<" - opposite (hinge) edge is ("<<mesh_.edge(op_eh).from_vertex()<<", "<<mesh_.edge(op_eh).to_vertex()<<")"<<std::endl;

        double dihedral_sum_deg = compute_dihedral_angle_sum_deg(mesh_, op_eh);
        log(LogLevel::INFO)<<" - dihedral sum = "<<dihedral_sum_deg<<"°"<<std::endl;

        double current_boundary_dihedral_angle_deg = 360.0 - dihedral_sum_deg;
        log(LogLevel::INFO)<<" - current boundary dihedral angle = "<<current_boundary_dihedral_angle_deg<<"°"<<std::endl;

        angles_diff_deg.push_back(current_boundary_dihedral_angle_deg - boundary_dihedral_angle_deg_prop_[op_eh]);
        log(LogLevel::INFO)<<" --> which is = "<<angles_diff_deg.back()<<"° off the codomain boundary dihedral angle of "<<boundary_dihedral_angle_deg_prop_[op_eh]<<"°"<<std::endl;

        //update max
        max_abs_angle_diff_deg = std::max(max_abs_angle_diff_deg,std::abs(angles_diff_deg.back()));
    }

    log(LogLevel::INFO)<<" ------------ "<<std::endl;
}


template <typename DistortionElemT, typename VirtualElemT>
void
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
compute_dihedral_diff_progress(double &avg_diff,
                               double &max_diff,
                               double& weighted_diff_norm) {

    avg_diff = 0;
    max_diff = 0;
    weighted_diff_norm = 0;

    double total_boundary_edge_length = initial_boundary_edges_count_ * avg_codomain_boundary_edge_length_;

    double w_sum(0);
    int boundary_e_count(0);
    for(auto eh: mesh_.edges()){
        if(mesh_.is_boundary(eh)) {
            //log(LogLevel::INFO) << " ------------ checking edge " << eh << ": " << edge_to_str(mesh_, eh) << std::endl;

            auto from_v = mesh_.edge(eh).from_vertex();
            auto   to_v = mesh_.edge(eh).to_vertex();

            double current_angle_deg = compute_boundary_dihedral_angle_deg(mesh_, eh);

            double diff_deg = std::abs(boundary_dihedral_angle_deg_prop_[eh] - current_angle_deg);
            boundary_dihedral_angle_diff_deg_prop_[eh] = diff_deg;

            double codomain_length = (codomain_pos_prop_[from_v] - codomain_pos_prop_[to_v]).norm();
            double w = codomain_length / total_boundary_edge_length;
            w_sum += w;

            weighted_diff_norm += w * diff_deg * diff_deg;

            max_diff = std::max(max_diff, diff_deg);
            avg_diff += diff_deg;

            boundary_e_count++;
        }
    }

    avg_diff     /= boundary_e_count;
    //weighted_diff_norm = std::sqrt(weighted_diff_norm);
}



template <typename DistortionElemT, typename VirtualElemT>
void
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
update_log_IVL_prop() {
    for (auto ch: mesh_.cells()) {
        double VL_qual = mesh_.compute_quality_vector(ch, VOLUME_LENGTH, CONFORMING)[0];
        double IVL_qual = 1.0 / VL_qual;
        log_IVL_prop_[ch] = std::log10(IVL_qual);
        //std::cout<<" - log IVL for cell "<<cell_to_str(mesh_, ch)<<": "<<log_IVL_prop_[ch] <<std::endl;
    }
}


template <typename DistortionElemT, typename VirtualElemT>
void
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
update_log_rel_vol_prop() {
    //auto avg_vol = compute_total_volume(mesh_) / mesh_.n_logical_cells();
    auto avg_vol = target_total_volume_ / mesh_.n_logical_cells();

    for (auto ch: mesh_.cells()) {
        log_rel_vol_prop_[ch] = std::log10(compute_volume(mesh_, ch) / avg_vol);
    }
}


template <typename DistortionElemT, typename VirtualElemT>
double
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
compute_boundary_edge_vector_rel_diff(double& rel_diff_norm,
                                      double& max_rel_diff) const {

    max_rel_diff = 0;

    double rel_diff_acc = 0;
    //log(LogLevel::INFO)<<" ------------ "<<std::endl;
    //log(LogLevel::INFO)<<" computing max boundary edge vector relative diff... "<<std::endl;

    for(auto eh: mesh_.edges()) {
        if(mesh_.is_boundary(eh)) {

            auto from_v = mesh_.edge(eh).from_vertex();
            auto   to_v = mesh_.edge(eh).to_vertex();

            auto current_edge_vec = mesh_.vertex(to_v) - mesh_.vertex(from_v);
            auto target_edge_vec = codomain_pos_prop_[to_v] - codomain_pos_prop_[from_v];

            auto target_edge_vec_norm = target_edge_vec.norm();
            if(target_edge_vec_norm < DBL_EPSILON){
                log(LogLevel::WARNING)<<" WARNING - target edge vector for edge "<<edge_to_str(mesh_, eh)<<" is "<<target_edge_vec<<" with norm = "<<target_edge_vec.norm()<<std::endl;
                target_edge_vec_norm = DBL_EPSILON;
            }

            double rel_diff = (target_edge_vec - current_edge_vec).norm() / target_edge_vec_norm;

            max_rel_diff = std::max(max_rel_diff, rel_diff);
            rel_diff_acc += rel_diff * rel_diff;
        }
    }

    return std::sqrt(rel_diff_acc);
}





template <typename DistortionElemT, typename VirtualElemT>
EdgeHandle
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
opposite_edge(const OVM::HalfFaceHandle &hfh,
              const OVM::VertexHandle &vh) const {

    auto verts = mesh_.get_halfface_vertices(hfh);
    if (verts[0] == vh) {
        return mesh_.edge_handle(mesh_.find_halfedge(verts[1], verts[2]));
    } else if (verts[1] == vh) {
        return mesh_.edge_handle(mesh_.find_halfedge(verts[0], verts[2]));
    } else {
        return mesh_.edge_handle(mesh_.find_halfedge(verts[0], verts[1]));
    }
}



template <typename DistortionElemT, typename VirtualElemT>
void
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
compute_target_angles_deg(const std::vector<OVM::HalfFaceHandle> &to_wrap_hfs,
                          const OVM::VertexHandle &to_wrap_vh,
                          std::vector<double> &current_angles_deg,
                          double &max_abs_current_angle_deg,
                          std::vector<double> &target_angles_deg,
                          double & max_abs_target_angle_deg,
                          int& max_abs_target_angle_deg_idx) const{


    compute_boundary_dihedral_angle_diff(to_wrap_hfs, to_wrap_vh, current_angles_deg, max_abs_current_angle_deg);

    max_abs_target_angle_deg = 0;
    max_abs_target_angle_deg_idx = 0;

    for(int i(0); i < current_angles_deg.size(); i++){
        auto angle_deg = current_angles_deg[i];
        target_angles_deg.push_back(angle_deg < 90.0 ? std::min(0.0, angle_deg + 90.0) : std::max(0.0, angle_deg - 90.0));
        //max_abs_target_angle_deg = std::max(max_abs_target_angle_deg, std::abs(target_angles_deg.back()));

        if(std::abs(target_angles_deg.back()) > max_abs_target_angle_deg){
            max_abs_target_angle_deg = std::abs(target_angles_deg.back());
            max_abs_target_angle_deg_idx = i;
            log(LogLevel::IMPORTANT)<<" - updated halfface with largest abs target angle as n°"<<i<<": "<<halfface_to_str(mesh_, to_wrap_hfs[i])<<" with target angle "<<target_angles_deg.back()<<"°"<<std::endl;
        }
    }
}


template <typename DistortionElemT, typename VirtualElemT>
void
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
mark_as_wrapped(const VertexHandle& vh) {
    wrapped_prop_[vh] = true;
    wrapped_count_++;
}


template <typename DistortionElemT, typename VirtualElemT>
WRAPPING_RESULT
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
run_safety_checks(const std::string& step_description) const {
    if (json_data_["enable_all_safety_checks"]) {

        if (check_boundary_holonomy(mesh_)) {
            log(LogLevel::ERROR) << " ERROR - holonomy issue after " << step_description << std::endl;
            return WRAPPING_HOLONOMY_ERROR;
        }

        //if(check_wrapped_dihedral_angles)

        if (ref_mesh_.contains_degenerate_or_flipped_elements(CONFORMING)) {
            log(LogLevel::ERROR) << " ERROR - input ref mesh contains degenerate or flipped elements after " << step_description << std::endl;

            for (auto elh: ref_mesh_.elements()) {
                if (ref_mesh_.is_flipped_or_degenerate(elh, CONFORMING)) {
                    log(LogLevel::ERROR) << " --> element " << ref_mesh_.element_to_str(elh) << " with volume "
                                         << CGAL_wrapper::exact_volume(ref_mesh_, elh).to_double() << std::endl;
                }
            }
            return WRAPPING_ERROR;
        }

        if(mesh_.contains_degenerate_or_flipped_elements(CONFORMING)){
            log(LogLevel::ERROR)<<" ERROR - mesh contains degenerate or flipped elements after "<<step_description<<std::endl;

            for(auto elh: mesh_.elements()){
                if(mesh_.is_flipped_or_degenerate(elh, CONFORMING)){
                    log(LogLevel::ERROR)<<" --> element "<<mesh_.element_to_str(elh)<<" with volume "<<CGAL_wrapper::exact_volume(mesh_, elh).to_double()<<std::endl;
                }
            }
            return WRAPPING_ERROR;
        }

    }
    return WRAPPING_SUCCESS;
}


template <typename DistortionElemT, typename VirtualElemT>
WRAPPING_RESULT
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
run_final_checks(double constr_eps,
                 LightWeightStopWatch& sw) const {

    if(!boundary_constraints_are_all_satisifed(constr_eps)){
        log(LogLevel::IMPORTANT)<<" +----------------------------------------------+"<<std::endl;
        log(LogLevel::IMPORTANT)<<" | +------------------------------------------+ |"<<std::endl;
        log(LogLevel::IMPORTANT)<<" | |  Boundary constraints validation failed  | |"<<std::endl;
        log(LogLevel::IMPORTANT)<<" | +------------------------------------------+ |"<<std::endl;
        log(LogLevel::IMPORTANT)<<" +----------------------------------------------+"<<std::endl;
        return FINAL_CONSTRAINTS_SATISFACTION_ERROR;
    }

    auto bad_tets = CGAL_wrapper::find_degenerate_or_flipped_tets(mesh_);
    if(!bad_tets.degenerate.empty() || !bad_tets.flipped.empty()){

        log(LogLevel::IMPORTANT)<<" +----------------------------------------------+"<<std::endl;
        log(LogLevel::IMPORTANT)<<" | +------------------------------------------+ |"<<std::endl;
        log(LogLevel::IMPORTANT)<<" | |  Elements orientation validation failed  | |"<<std::endl;
        log(LogLevel::IMPORTANT)<<" | +------------------------------------------+ |"<<std::endl;
        log(LogLevel::IMPORTANT)<<" +----------------------------------------------+"<<std::endl;

        log(LogLevel::IMPORTANT)<<"  degenerate element: "<<std::endl;
        for (auto ch: bad_tets.degenerate) {
            log(LogLevel::IMPORTANT)<<" - "<<cell_to_str(mesh_, ch)<<" with volume "<<CGAL_wrapper::exact_volume(mesh_, ch).to_double()<<std::endl;
        }
        log(LogLevel::IMPORTANT)<<" +----------------------------------------------+"<<std::endl;

        log(LogLevel::IMPORTANT)<<"  flipped element: "<<std::endl;
        for (auto ch: bad_tets.flipped) {
            log(LogLevel::IMPORTANT)<<" - "<<cell_to_str(mesh_, ch)<<" with volume "<<CGAL_wrapper::exact_volume(mesh_, ch).to_double()<<std::endl;
        }

        return FINAL_ORIENTATION_CHECK_ERROR;
    }

    log(LogLevel::IMPORTANT)<<" ********************************************************************************************************************************** "<<std::endl;
    log(LogLevel::IMPORTANT)<<" ----> constraints are satisfied and the map is bijective!"<<std::endl;
    log(LogLevel::IMPORTANT)<<" ----> full process took "<<sw.total_duration()<<" seconds"<<std::endl;
    log(LogLevel::IMPORTANT)<<" ********************************************************************************************************************************** "<<std::endl;

    return WRAPPING_SUCCESS;
}


template <typename DistortionElemT, typename VirtualElemT>
double
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
snap_satisfied_boundary_constraints(double constr_eps) {

    double max_violation = 0.0;
    for(auto vh: mesh_.vertices()){
        if(mesh_.is_boundary(vh)){

            auto diff = (mesh_.vertex(vh) - codomain_pos_prop_[vh]).norm();
            if(diff < constr_eps){
                mesh_.set_vertex(vh, codomain_pos_prop_[vh]);
                //log(LogLevel::INFO)<<" - snapped boundary vertex "<<vh<<" to its target position"<<std::endl;
            }
            if (diff >= max_violation) {
                max_violation = diff;
            }
        }
    }
    return max_violation;
}


template <typename DistortionElemT, typename VirtualElemT>
bool
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
boundary_constraints_are_all_satisifed(double constr_eps) const {

    for(auto vh: mesh_.vertices()){
        if(mesh_.is_boundary(vh)){
            if(mesh_.vertex(vh) != codomain_pos_prop_[vh]){
                log(LogLevel::INFO) << " -> boundary vertex " << vh << " is at " << mesh_.vertex(vh) << " but its prescribed position is " << codomain_pos_prop_[vh] << " (diff = " << (mesh_.vertex(vh) - codomain_pos_prop_[vh]).norm() << std::endl;
                return false;
            }
        }
    }

    return true;
}



template <typename DistortionElemT, typename VirtualElemT>
void
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
check_for_timeout() const{

    if(max_allocated_time_s_ < 0){
        return;
    }

    if(main_sw_.total_duration() > max_allocated_time_s_){
        throw TimeOutException();
    }

}

template <typename DistortionElemT, typename VirtualElemT>
OptimizerStats
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
opti_stats() const {
    return opti_stats_;
}


template <typename DistortionElemT, typename VirtualElemT>
std::vector<std::vector<Vec3d>>&
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
boundary_vertices_trajectories(){
    return boundary_vertices_trajectories_;
}

template <typename DistortionElemT, typename VirtualElemT>
std::vector<std::vector<Vec3d>>&
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
boundary_vertices_target_paths(){
    return boundary_vertices_target_paths_;
}

