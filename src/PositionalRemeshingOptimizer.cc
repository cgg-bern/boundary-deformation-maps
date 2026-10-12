#define POSITIONAL_REMESHING_OPTIMIZER_C


#include "RemeshingOptimizerV2.hh"



//NOTE: having these two member functions of a different class in a separate file is a bit unconventional,
// but I keep making modifications to this version when I want to make them on the actual method, and lose precious time because of it.
// SO it will have to do for now

template <typename DistortionElemT, typename VirtualElemT>
WRAPPING_RESULT
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
run_positional(double max_allocated_time_s) {
    main_sw_.reset();
    max_allocated_time_s_ = max_allocated_time_s;

    log(LogLevel::IMPORTANT)<<" ================================================================================================================================ "<<std::endl;
    log(LogLevel::IMPORTANT)<<" ================================================================================================================================ "<<std::endl;
    log(LogLevel::IMPORTANT)<<" -> running global wrapping process with positional constraints (and without virtual elements)..."<<std::endl;
    log(LogLevel::IMPORTANT)<<" ================================================================================================================================ "<<std::endl;


    try {

        const bool enable_focused_remeshing = json_data_["global_wrapping_enable_focused_remeshing"];
        const int deep_remeshing_inner_iters = json_data_["global_wrapping_deep_remeshing_inner_iters"];
        const double quick_remeshing_candidate_quality_threshold = json_data_["global_wrapping_quick_remeshing_candidate_quality_threshold"];
        const double quick_remeshing_target_quality_threshold = json_data_["global_wrapping_quick_remeshing_target_quality_threshold"];
        const double deep_remeshing_candidate_quality_threshold = json_data_["global_wrapping_deep_remeshing_candidate_quality_threshold"];
        const double deep_remeshing_target_quality_threshold = json_data_["global_wrapping_deep_remeshing_target_quality_threshold"];

        const double constraints_progress_thresh = json_data_["positional_wrapping_constraints_violation_threshold"];
        const int max_outer_iters = json_data_["global_wrapping_max_outer_iters"];
        const int newton_iters = json_data_["global_wrapping_max_newton_iters"];

        //const bool enable_full_history = json_data_["global_wrapping_enable_full_history"];

        const double start_quality_w = json_data_["global_wrapping_quality_energy_start_w"];
        const bool enable_adaptive_start_w = json_data_["global_wrapping_enable_adaptive_start_w"];

        //last (4th) entry is for the quality energy
        const std::vector<double> adaptive_start_w_target_g_norm = std::vector<double>(json_data_["global_wrapping_adaptive_start_w_target_g_norm"]);

        const double constr_eps = json_data_["constr_eps"];

        int initial_boundary_vertices_count = count_boundary_vertices(mesh_);


        if(initialize_codomain()){
            log(LogLevel::ERROR)<<" ERROR - couldn't initialize codomain"<<std::endl;
            return WRAPPING_ERROR;
        }

        double max_abs_dihedral_diff_deg, max_dihedral_diff_deg_thresh;
        enter_snap_phase_if_possible(false, true, max_abs_dihedral_diff_deg, max_dihedral_diff_deg_thresh);


        if (auto res = run_safety_checks("codomain initialization")) {
            return res;
        }

        IO::FileManager fm;
        fm.writeFile(mesh_name_ + "_codomain_mesh_init.ovm", mesh_);
        export_boundary_as_triangle_mesh(mesh_, mesh_name_+"_codomain_mesh_boundary_init.obj");


        double min_quality;
        std::vector<double> quality_vector;
        compute_element_quality(quality_vector, min_quality);

        double min_rel_vol = compute_min_relative_volume();
        double min_rel_edge_len = compute_min_edge_length(mesh_) / min_codomain_boundary_edge_length_;

        opti_stats_.wrapped_ratio_history_.push_back((double)wrapped_count_ / (double)initial_boundary_vertices_count);
        opti_stats_.worst_quality_history_.push_back(min_quality);
        opti_stats_.min_relative_volume_history_.push_back(min_rel_vol);
        opti_stats_.min_relative_edge_length_history_.push_back(min_rel_edge_len);
        opti_stats_.init_worst_quality_ = min_quality;
        opti_stats_.init_min_relative_volume_ = min_rel_vol;
        opti_stats_.init_min_relative_edge_length_ = min_rel_edge_len;
        /*if (enable_full_history) {
            opti_stats_.quality_history_.push_back(quality_vector);
        }*/


        log(LogLevel::IMPORTANT) << "      initial worst element quality =  " <<min_quality<< std::endl;
        log(LogLevel::IMPORTANT) << "        initial min relative volume =  " <<min_rel_vol<< std::endl;
        log(LogLevel::IMPORTANT) << "   initial min relative edge length =  " <<min_rel_edge_len<< std::endl;
        log(LogLevel::IMPORTANT) << " ---------------------------------------------------------------------- " << std::endl;

        bool converged(false);
        EnergyRecord energy_record;

        Eigen::VectorXd quality_g;
        std::array<Eigen::VectorXd, 3> virtual_gs;

        quality_weight_scaling_ = start_quality_w;

        opti_stats_.backup_mode_ = false;

        //EnergyArray last_pass_unweighted_initial_virtual_e(std::numeric_limits<double>::max());
        //EnergyArray last_pass_unweighted_final_virtual_e(std::numeric_limits<double>::max());

        double last_constraint_violation_inf_norm = std::numeric_limits<double>::max();


        if(json_data_["enable_vertex_trajectory_mesh"]){
            for (auto vh: mesh_.vertices()) {
                if(mesh_.is_boundary(vh)){
                    boundary_vertices_trajectories_[vh.idx()].push_back(mesh_.vertex(vh));
                }
            }
        }


        int& outer_wrap_it = opti_stats_.outer_iters_;
        while(outer_wrap_it < max_outer_iters && wrapped_count_ != initial_boundary_vertices_count) {

            LightWeightStopWatch pass_sw;


            log(LogLevel::IMPORTANT) << " ================================================================================================================================ " << std::endl;
            log(LogLevel::IMPORTANT) << " ======================================================================================   running global wrap pass " << outer_wrap_it  << std::endl;
            log(LogLevel::IMPORTANT) << " ====================================================================================== last constraints violation: " << last_constraint_violation_inf_norm << std::endl;
            log(LogLevel::IMPORTANT) << " ======================================================================================              wrapped count: " << wrapped_count_ << " / " << initial_boundary_vertices_count << std::endl;
            log(LogLevel::IMPORTANT) << " ======================================================================================         total elapsed time: " << main_sw_.total_duration() << "s" << std::endl;
            log(LogLevel::IMPORTANT) << " ================================================================================================================================ " << std::endl;



            check_for_timeout();
            LightWeightStopWatch it_sw;

            log(LogLevel::IMPORTANT) << " ================================================================================================================================ " << std::endl;
            log(LogLevel::IMPORTANT) << " ============================================================================   running global wrap iteration " << outer_wrap_it << std::endl;
            log(LogLevel::IMPORTANT) << " ============================================================================                  wrapped count: " << wrapped_count_ << " / " << initial_boundary_vertices_count << std::endl;
            log(LogLevel::IMPORTANT) << " ============================================================================             total elapsed time: " << main_sw_.total_duration() << "s" << std::endl;
            log(LogLevel::IMPORTANT) << " ================================================================================================================================ " << std::endl;

            double constraint_violation_inf_norm;
            if (run_global_positional_wrapping(newton_iters, converged, energy_record, quality_g,
                                               constraint_violation_inf_norm)) {
                log(LogLevel::ERROR) << " error while running global wrap opti" << std::endl;
                return WRAPPING_ERROR;
            }

            if (json_data_["enable_debugging_meshes_exportation"]) {
                update_log_IVL_prop();
                update_log_rel_vol_prop();
                std::string ext = json_data_["debugging_meshes_exportation_extension"];
                mesh_.export_conforming_mesh(mesh_name_+"_codomain_mesh_after_it_" + std::to_string(outer_wrap_it)  + ext);
                export_boundary_as_triangle_mesh(mesh_, mesh_name_+"_codomain_mesh_boundary_after_it_" + std::to_string(outer_wrap_it)  + ".obj");
            }

            if (json_data_["global_wrapping_enable_analysis"]) {

                log(LogLevel::INFO) << " - running post-optimization analysis" << std::endl;

                analyze_all_elements(ref_mesh_, mesh_, compute_ref_to_mesh_volume_scaling_factor(), true, 1.0, CONFORMING, "");
                analyze_all_elements(ref_mesh_, ref_mesh_, 1.0, false, ref_mesh_total_volume_ / ref_mesh_.n_logical_cells(), CONFORMING, "");
            }


            if (auto res = run_safety_checks("wrap_it_" + std::to_string(outer_wrap_it))) {
                log(LogLevel::ERROR)<<" safety checks failed, stopping here"<<std::endl;
                return res;
            }

            log(LogLevel::IMPORTANT) << " ---------------------------------------------------------------------- " << std::endl;

            double unweighted_initial_quality_e = energy_record.initial_quality_energy_ / quality_weight_scaling_;
            double unweighted_final_quality_e   = energy_record.final_quality_energy_   / quality_weight_scaling_;
            auto      quality_energy_progress = (unweighted_initial_quality_e - unweighted_final_quality_e) / unweighted_initial_quality_e;

            std::vector<double> quality_vector;
            double min_quality;
            compute_element_quality(quality_vector, min_quality);

            min_rel_vol = compute_min_relative_volume();
            min_rel_edge_len = compute_min_edge_length(mesh_) / min_codomain_boundary_edge_length_;

            opti_stats_.unweighted_initial_quality_energy_history_.push_back(unweighted_initial_quality_e);

            opti_stats_.unweighted_final_quality_energy_history_.push_back(unweighted_final_quality_e);

            opti_stats_.unweighted_final_quality_gradient_norm_history_.push_back(energy_record.final_quality_gradient_norm_);
            opti_stats_.unweighted_final_primary_virtual_gradient_norm_history_.push_back(energy_record.final_virtual_gradient_norm_.primary);
            opti_stats_.unweighted_final_secondary_virtual_gradient_norm_history_.push_back(energy_record.final_virtual_gradient_norm_.secondary);

            opti_stats_.final_total_gradient_max_coeff_history_.push_back(energy_record.final_total_gradient_max_coeff_);
            opti_stats_.final_total_Hessian_max_coeff_history_.push_back(energy_record.final_total_Hessian_max_coeff_);

            opti_stats_.worst_quality_history_.push_back(min_quality);
            opti_stats_.min_relative_volume_history_.push_back(min_rel_vol);
            opti_stats_.min_relative_edge_length_history_.push_back(min_rel_edge_len);

            /*if (enable_full_history) {
                opti_stats_.quality_history_.push_back(quality_vector);
            }*/

            double constraint_violation_progress = (last_constraint_violation_inf_norm - constraint_violation_inf_norm) / last_constraint_violation_inf_norm;


            log(LogLevel::IMPORTANT) << " ------------------------------------------------- " << std::endl;
            log(LogLevel::IMPORTANT) << "     last constraints violation inf norm = " << last_constraint_violation_inf_norm << std::endl;
            log(LogLevel::IMPORTANT) << "  updated constraints violation inf norm = " << constraint_violation_inf_norm << std::endl;
            log(LogLevel::IMPORTANT) << "          constraints violation progress = " << constraint_violation_progress << std::endl;

            log(LogLevel::IMPORTANT) << " ------------------------------------------------- " << std::endl;
            log(LogLevel::IMPORTANT) << "             initial quality energy =  " << energy_record.initial_quality_energy_ << std::endl;
            log(LogLevel::IMPORTANT) << "             final   quality energy =  " << energy_record.final_quality_energy_ << std::endl;
            log(LogLevel::IMPORTANT) << "  unweighted initial quality energy =  " << unweighted_initial_quality_e << std::endl;
            log(LogLevel::IMPORTANT) << "  unweighted   final quality energy =  " << unweighted_final_quality_e << std::endl;
            log(LogLevel::IMPORTANT) << "            quality energy progress =  " << quality_energy_progress << std::endl;
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

            last_constraint_violation_inf_norm = constraint_violation_inf_norm;


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

            int snapped_candidates_count = 0;
            if (json_data_["positional_wrapping_enable_snapping"]) {
                snapped_candidates_count = try_to_snap_next_candidates();
                log(LogLevel::IMPORTANT) << " --> snapped a total of " << snapped_candidates_count <<
                                         " non-wrapped candidates" << std::endl;

                if (json_data_["enable_debugging_meshes_exportation"]) {
                    std::string ext = json_data_["debugging_meshes_exportation_extension"];
                    mesh_.export_conforming_mesh(
                            mesh_name_ + "_codomain_mesh_after_it_" + std::to_string(outer_wrap_it) + "_and_snapping" + ext);
                }
            }


            opti_stats_.wrapped_ratio_history_.push_back((double) wrapped_count_ / (double) initial_boundary_vertices_count);


            if (constraint_violation_inf_norm > constr_eps && wrapped_count_ != initial_boundary_vertices_count) {

                log(LogLevel::IMPORTANT) << " ---> still some vertices to snap, remeshing and moving on to next pass " << std::endl;
                auto hc_result = run_remeshing(GLOBAL_REMESHING,
                                               outer_wrap_it,
                                               deep_remeshing_inner_iters,
                                               deep_remeshing_candidate_quality_threshold,
                                               deep_remeshing_target_quality_threshold);
                if (hc_result == -1) {
                    log(LogLevel::ERROR) << " error while remeshing" << std::endl;
                    return WRAPPING_ERROR;
                }
                opti_stats_.deep_remeshing_iterations_.push_back(opti_stats_.inner_iters_);
                if (json_data_["positional_wrapping_enable_snapping"]) {

                    snapped_candidates_count = try_to_snap_next_candidates();
                    log(LogLevel::IMPORTANT) << " --> snapped a total of " << snapped_candidates_count <<" non-wrapped candidates" << std::endl;

                    if (json_data_["enable_debugging_meshes_exportation"]) {
                        std::string ext = json_data_["debugging_meshes_exportation_extension"];
                        mesh_.export_conforming_mesh(mesh_name_+"_codomain_mesh_after_it_" + std::to_string(outer_wrap_it) + "_and_snapping" + ext);
                    }
                }
            }

            log(LogLevel::IMPORTANT) << " -----> done with pass " << outer_wrap_it << " in " << pass_sw.total_duration() << " seconds" << std::endl;

            //set default value for potential new edges (otherwise it messes with the file loading in OF)
            for (auto eh: mesh_.edges()) {
                virtual_vertex_pos_prop_[eh] = {0,0,0};
            }

            if (constraint_violation_inf_norm < constr_eps) {
                log(LogLevel::IMPORTANT) << " -----> max constraints violation is "<<constraint_violation_inf_norm<<" < "<<constr_eps<<", stopping here" << std::endl;
                break;
            }
            if (!snapped_candidates_count && constraint_violation_progress < constraints_progress_thresh) {
                log(LogLevel::IMPORTANT) << " -----> no progress was made, stopping here" << std::endl;
                break;
            }


            outer_wrap_it++;
        }

        int snapped_candidates_count = try_to_snap_next_candidates();
        log(LogLevel::IMPORTANT) << " --> snapped a total of " << snapped_candidates_count <<
                                 " non-wrapped candidates" << std::endl;

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

        compute_element_quality(quality_vector, min_quality);

        min_rel_vol = compute_min_relative_volume();
        min_rel_edge_len = compute_min_edge_length(mesh_) / min_codomain_boundary_edge_length_;

        opti_stats_.final_worst_quality_ = min_quality;
        opti_stats_.final_min_relative_volume_ = min_rel_vol;
        opti_stats_.final_min_relative_edge_length_ = min_rel_edge_len;


        if(json_data_["enable_vertex_trajectory_mesh"]){
            for (auto vh: mesh_.vertices()) {
                if(mesh_.is_boundary(vh)){
                    boundary_vertices_target_paths_[vh.idx()].push_back(mesh_.vertex(vh));
                    boundary_vertices_target_paths_[vh.idx()].push_back(codomain_pos_prop_[vh]);
                }
            }
        }

        return run_final_checks(constr_eps, main_sw_);

    }catch(TimeOutException e) {
        log(LogLevel::INFO)<<"--------------------"<<std::endl;
        log(LogLevel::INFO)<<"||||||||||||||||||||"<<std::endl;
        log(LogLevel::INFO)<<e.what()<<std::endl;
        log(LogLevel::INFO)<<"||||||||||||||||||||"<<std::endl;
        log(LogLevel::INFO)<<"--------------------"<<std::endl;
        return WRAPPING_TIMEOUT;
    }catch(VirtualVertexSetupException e){
        log(LogLevel::INFO)<<"--------------------"<<std::endl;
        log(LogLevel::INFO)<<"||||||||||||||||||||"<<std::endl;
        log(LogLevel::INFO)<<e.what()<<std::endl;
        log(LogLevel::INFO)<<"||||||||||||||||||||"<<std::endl;
        log(LogLevel::INFO)<<"--------------------"<<std::endl;
        return WRAPPING_NUMERICS_ERROR;
    }

    return WRAPPING_SUCCESS;
}



template <typename DistortionElemT, typename VirtualElemT>
int
BoundaryGlueingRemeshingOptimizerT<DistortionElemT, VirtualElemT>::
run_global_positional_wrapping(int newton_iters,
                               bool &converged,
                               EnergyRecord &energy_record,
                               Eigen::VectorXd &quality_g,
                               double &constraints_violation_inf_norm) {

    log(LogLevel::IMPORTANT) <<
                             " ================================================================================================================================ "<<std::endl;
    log(LogLevel::IMPORTANT)<<" - entering positional optimization step..."<<std::endl;


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

    const bool enable_full_history = json_data_["global_wrapping_enable_full_history"];

    const double edge_length_bound_elements_w = json_data_["global_wrapping_edge_length_bound_elements_w"];
    const double edge_length_bound_elements_boundary_rel_bound= json_data_["global_wrapping_edge_length_bound_elements_boundary_rel_bound"];
    const double edge_length_bound_elements_bound_scaling = json_data_["global_wrapping_edge_length_bound_elements_bound_scaling"];


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

    double avg_vol = compute_total_volume(mesh_) / mesh_.n_logical_cells();
    std::vector<double> parameters;
    log(LogLevel::IMPORTANT)<<" -> average volume (used as optimum for volume bound): "<<avg_vol<<std::endl;

    //NOTE: those are the "base" weights. They are scaled in a per-element basis during elements setup
    if(DistortionElemT::element_code() == "IVL_B") {
        parameters = {
                json_data_["global_wrapping_IVL_w"],
                json_data_["global_wrapping_IVL_bound_w"],
                json_data_["global_wrapping_IVL_bound"],
                json_data_["global_wrapping_IVL_vol_w"],
                double(json_data_["global_wrapping_IVL_rel_vol_bound"]) * avg_vol,
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

            std::vector<VertexHandle> verts = {mesh_.edge(eh).from_vertex(), mesh_.edge(eh).to_vertex()};

            double el_target = mesh_.is_boundary(eh) ? boundary_edge_length_prop_[eh] : default_el_target;

            el_w_sum += el_w;
            max_w = std::max(el_w, max_el_w);

            ELEN_B_3D_AD::VecI vi_el;
            ELEN_B_3D_AD::VecV vx_el;
            ELEN_B_3D_AD::VecC vc_el;
            vc_el << el_w, el_bound, el_target;

            double len = (mesh_.vertex(mesh_.edge(eh).from_vertex()) - mesh_.vertex(mesh_.edge(eh).to_vertex())).norm();
            if(len < el_bound){
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
            log(LogLevel::DEBUG) << " - "<<(mesh_.is_boundary(eh) ? "boundary" : "interior")<<" edge "<<edge_to_str(mesh_, eh)<<" length: "<<len<<" rel length (to min boundary): "<<(len / min_codomain_boundary_edge_length_)<<", target: "<<el_target<<" energy: " <<el_e<< std::endl;

            fe_edge_length.instances().add_element(vi_el, vc_el);
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

    //constrain all boundary vertices
    for(auto vh: mesh_.vertices()) {
        if(mesh_.is_boundary(vh)) {

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
    }
    if (!constrained_count) {
        log(LogLevel::ERROR)<<" ERROR - No vertex constrained at all. You should constrain at least one "<<std::endl;
        return -1;
    }

    /*for (int i(0); i<eq_constraints.size(); i++) {
        double violation = eq_constraints[i].eval_constraint(fe_problem.x().data());
        if (violation > constr_eps) {
            log(LogLevel::ERROR) << "ERROR - violation for already-snapped vertex constraint "<<i<<" is "<<violation<<" > "<<constr_eps << std::endl;
            return WRAPPING_ERROR;
        }
    }*/

    std::vector<COMISO::NConstraintInterface*> constraints_pointers;
    constraints_pointers.reserve(constraints_pointers.size());
    for( auto& c : eq_constraints) {
        constraints_pointers.push_back(&c);
    }

    log(LogLevel::IMPORTANT) << " - total number of constraints pointers: " << constraints_pointers.size() << std::endl;

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

    constraints_violation_inf_norm = tn.status().constraint_violation_inf_norm;


    auto opti_time = main_sw_.lap_duration();
    opti_stats_.total_opti_time_s_ += opti_time;

    log(LogLevel::INFO) << " energies after optimization: " << std::endl;
    if (LogConfig::get_min_level() <= LogLevel::INFO) {
        fe_problem.print_objectives();
    }

    energy_record.final_quality_energy_ = fe_barrier.eval_f(fe_problem.x().data());


    log(LogLevel::INFO) << " -------------------------------------------" << std::endl;
    log(LogLevel::INFO) << " - computing quality energy gradient..." << std::endl;
    quality_g.resize(fe_problem.x().size());
    quality_g.setZero();
    fe_barrier.accumulate_gradient(fe_problem.x().data(), quality_g.data());
    log(LogLevel::INFO) << " -> quality energy gradient norm = "<<quality_g.norm() << std::endl;
    energy_record.final_quality_gradient_norm_ = quality_g.norm();
    energy_record.final_quality_gradient_max_coeff_ = quality_g.maxCoeff();

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

    log(LogLevel::IMPORTANT)<<" ================================================================================================================================ "<<std::endl;

    return 0;
}

