
#include <iostream>


#include "NonConformingHelpers.hh"
#include "RemeshingOptimizerV2.hh"

#include "VoxelGridMeshGen.hh"
#include "ParametricMeshGen.hh"
#include "SimpleMeshGen.hh"

//#include "CrossFieldTestingT.hh"



using namespace OpenVolumeMesh;


int parse_args(int argc, char** argv, argparse::ArgumentParser& parser);

void short_tests(argparse::ArgumentParser& parser);

void extract_orientation_stats(const TetrahedralMesh& mesh,
                           const std::string& mesh_name,
                           const std::string& stats_filepath);

//here because using the furch's knot generator
void setup_torsion_cube_mesh(TetrahedralMesh &ref_mesh,
                              TetrahedralMesh &mesh,
                              int cube_len,
                              double torsion);



/*
int setup_args(const argparse::ArgumentParser& parser,
               std::string&   domain_path,
               std::string& codomain_path,
               std::string& boundary_path);*/



bool quick_and_dirty_OBJ_load(const std::string &path, TetrahedralMesh& mesh);

std::string encode_vertex_trajectories_as_OBJ(const std::vector<std::vector<Vec3d>>& trajectories);

std::vector<std::string> split_string(const std::string& s) {
    std::istringstream iss(s);
    std::vector<std::string> result;
    std::string word;
    while (iss >> word) {
        result.push_back(word);
    }
    return result;
}


#warning TODO: move this to helpers
void copy_boundary_conditions_to_mesh(TetrahedralMesh& mesh,
                                      const std::string& boundary_conditions_filepath) {

    std::ifstream file(boundary_conditions_filepath);
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
// Read idx
            int idx(-1);
            iss >> idx;

            Vec3d pos;
            for (int i = 0;i < 3;++i) {
                iss >> pos[i];
            }
            mesh.set_vertex(VertexHandle(idx), pos);
        }
        file.close();
    } else {
        std::cerr << "Couldn't open boundary file " << boundary_conditions_filepath << std::endl;
        exit(EXIT_FAILURE);
    }
}


void export_boundary_conditions_from_mesh(TetrahedralMesh &mesh,
                                          const std::string& boundary_conditions_filepath) {

    std::ofstream file(boundary_conditions_filepath);
    if (file.is_open()) {
        std::string line;

        for (auto vh: mesh.vertices()) {
            if (mesh.is_boundary(vh)) {
                auto& pos = mesh.vertex(vh);

                file <<vh;
                for (int i(0); i<3; i++) {
                    file << " "<<std::setprecision(30)<<pos[i];
                }
                file << std::endl;
            }
        }

        /*while (std::getline(file, line)) {
            std::istringstream iss(line);
            // Read idx
            int idx(-1);
            iss >> idx;

            Vec3d pos;
            for (int i = 0;i < 3;++i) {
                iss >> pos[i];
            }
            mesh.set_vertex(VertexHandle(idx), pos);
        }*/
        file.close();
    } else {
        std::cerr << "Couldn't open boundary file " << boundary_conditions_filepath<<" for writing" << std::endl;
        exit(EXIT_FAILURE);
    }
}

template<typename MeshT>
void setup_minimal_mesh(MeshT& ref_mesh, MeshT& mesh);

template<>
void setup_minimal_mesh<TetrahedralMesh>(TetrahedralMesh& ref_mesh, TetrahedralMesh& mesh){

    //tri-tet for testing
    auto v0 = mesh.add_vertex({-1,-1,0});
    auto v1 = mesh.add_vertex({1,-1,0});
    auto v2 = mesh.add_vertex({0,1,0});
    auto v3 = mesh.add_vertex({0,0, -1});
    auto v4 = mesh.add_vertex({0,0,1});

    mesh.add_cell({v0,v1,v3,v4});
    mesh.add_cell({v1,v2,v3,v4});
    mesh.add_cell({v2,v0,v3,v4});

    mesh.split_edge(mesh.edge_handle(mesh.find_halfedge(v3, v4)));

    mesh.collect_garbage();

    ref_mesh = mesh;
}

template<>
void setup_minimal_mesh<NonConformingTetrahedralMesh>(NonConformingTetrahedralMesh& ref_mesh, NonConformingTetrahedralMesh& mesh){

    //tri-tet for testing
    auto v0 = mesh.add_vertex({-1,-1,0});
    auto v1 = mesh.add_vertex({1,-1,0});
    auto v2 = mesh.add_vertex({0,1,0});
    auto v3 = mesh.add_vertex({0,0, -1});
    auto v4 = mesh.add_vertex({0,0,1});

    mesh.add_cell({v0,v1,v3,v4});
    mesh.add_cell({v1,v2,v3,v4});
    mesh.add_cell({v2,v0,v3,v4});

    mesh.collect_garbage();

    ref_mesh = mesh;
}


void write_validation_file(std::string filename){


    std::ofstream validation_file;
    validation_file.open(filename);

    if(!validation_file.is_open()){
        log(LogLevel::WARNING)<<" WARNING - COULD NOT WRITE FILE VALIDATION FILE AT "<<filename<<std::endl;
        return;
    }
    log(LogLevel::INFO)<<" wrote validation file at "<<filename<<std::endl;
    validation_file.close();

}


void remove_validation_file(std::string filename){

    if(!remove(filename.c_str())){
        log(LogLevel::INFO)<<" --> removed validation file at "<<filename<<std::endl;
    }else{
        log(LogLevel::WARNING)<<" --> WARNING - FAILED TO REMOVE VALIDATION FILE AT "<<filename<<std::endl;
    }

}




int run_3d_boundary_deformation_map(TetrahedralMesh& domain_mesh,
                          TetrahedralMesh& codomain_mesh,
                          const std::string& mesh_name,
                          bool skip_pre_processing,
                          const std::string& json_stats_path,
                          argparse::ArgumentParser& parser,
                          nlohmann::json& json_data){

    if (domain_mesh.n_logical_vertices() != codomain_mesh.n_logical_vertices() ||
    domain_mesh.n_logical_edges() != codomain_mesh.n_logical_edges() ||
    domain_mesh.n_logical_faces() != codomain_mesh.n_logical_faces() ||
    domain_mesh.n_logical_cells() != codomain_mesh.n_logical_cells()) {
        log(LogLevel::ERROR)<<" ERROR - domain and codomain meshes don't have the same size"<<std::endl;
        exit(EXIT_FAILURE);
    }

    log(LogLevel::IMPORTANT)<<" - exporting remeshing optimization stats to "<<json_stats_path<<std::endl;

    LightWeightStopWatch sw;

    std::string log_level = json_data["log_level"];
    LogConfig::set_min_level(log_level);

    OVM::IO::FileManager fm;

    NonConformingTetrahedralMesh ref_mesh, mesh;
    ref_mesh = domain_mesh;
    mesh = codomain_mesh;


    int initial_v_count = mesh.n_vertices();
    int initial_e_count = mesh.n_edges();
    int initial_f_count = mesh.n_faces();
    int initial_c_count = mesh.n_cells();

    log(LogLevel::INFO)<<" --- initial domain mesh boundary stats"<<std::endl;
    print_boundary_quality_stats(ref_mesh);

    double min_boundary_dist, avg_boundary_dist, max_boundary_dist;
    check_boundary_SDE(ref_mesh, mesh, min_boundary_dist, avg_boundary_dist, max_boundary_dist);
    log(LogLevel::INFO)<<" ---------------------- ∂ distortion: "<<std::endl;
    log(LogLevel::INFO)<<" -        min: "<<min_boundary_dist<<std::endl;
    log(LogLevel::INFO)<<" -        avg: "<<avg_boundary_dist<<std::endl;
    log(LogLevel::INFO)<<" -        max: "<<max_boundary_dist<<std::endl;
    log(LogLevel::INFO)<<" ------------------------------------"<<std::endl;

    double ref_mesh_tot_vol = ref_mesh.total_volume();

    if(json_data["enable_initial_elements_analysis"]) {
        log(LogLevel::INFO) << " ---------------------- analyzing initial ref mesh elements " << std::endl;
        analyze_all_elements(ref_mesh, ref_mesh, 1.0, false,
            ref_mesh_tot_vol/ ref_mesh.n_logical_elements(), CONFORMING,
                                            mesh_name + "_domain_initial_elements.csv");
        log(LogLevel::INFO) << " ------------------------------------" << std::endl;
    }


    double input_min_quality = std::numeric_limits<double>::max();
    double input_min_rel_vol = std::numeric_limits<double>::max();
    for (auto ch: ref_mesh.cells()) {
        auto qual = ref_mesh.compute_quality_vector(ch, VOLUME_LENGTH, CONFORMING)[0];
        double vol = compute_volume(ref_mesh, ch);
        input_min_quality = std::min(input_min_quality, qual);
        input_min_rel_vol = std::min(input_min_rel_vol, vol);
    }
    input_min_rel_vol *= (double)initial_c_count / ref_mesh_tot_vol;

    double min_codomain_boundary_edge_length(0), max_codomain_boundary_edge_length(0), avg_codomain_boundary_edge_length(0);
    compute_boundary_edge_lengths(mesh, min_codomain_boundary_edge_length, max_codomain_boundary_edge_length, avg_codomain_boundary_edge_length);
    double input_min_rel_edge_len = compute_min_edge_length(ref_mesh) / min_codomain_boundary_edge_length;


    fm.writeFile("input_ref_tet_mesh.ovm", ref_mesh);
    export_boundary_as_triangle_mesh(mesh, "input_tet_mesh_boundary.obj");
    export_boundary_as_triangle_mesh(ref_mesh, "input_ref_tet_mesh_boundary.obj");
    log(LogLevel::IMPORTANT) << " - input meshes written" << std::endl;


    if (!skip_pre_processing) {

        std::vector<std::string> pre_remeshing_steps = json_data["pre_remeshing_steps"];

        for (auto step: pre_remeshing_steps) {

            if (step == "decimation") {
                if(json_data["enable_initial_mesh_decimation"]) {

                    if (json_data["decimation_mode"] == "interleaving") {
                        decimate_mesh_with_interleaving(ref_mesh, mesh, json_data);
                    }else {
                        decimate_mesh_with_phases(ref_mesh, mesh, json_data);
                    }
                    ref_mesh.collect_garbage();
                    mesh.collect_garbage();

                    if(json_data["enable_post_ref_mesh_decimation_elements_analysis"]) {
                        log(LogLevel::INFO) << " ---------------------- analyzing elements after decimation " << std::endl;
                        analyze_all_elements(ref_mesh, ref_mesh, 1.0, false,
                                             ref_mesh_tot_vol / ref_mesh.n_logical_elements(), CONFORMING,
                                             mesh_name + "_domain_elements_post_decimation.csv");
                        log(LogLevel::INFO) << " ------------------------------------" << std::endl;
                    }

                    fm.writeFile("input_ref_tet_mesh_post_decimation.ovm", ref_mesh);

                }
            }else if (step == "hill_climbing") {
                if(json_data["enable_ref_mesh_pre_remeshing"]){

                    if (json_data["pre_remeshing_mode"] == "interleaving") {
                        hill_climbing_ref_mesh_remeshing_with_interleaving(ref_mesh, mesh, json_data);
                    } else {
                        hill_climbing_ref_mesh_remeshing_with_phases(ref_mesh, mesh, json_data);
                    }

                    ref_mesh.collect_garbage();
                    mesh.collect_garbage();

                    if(json_data["enable_post_ref_mesh_pre_remeshing_elements_analysis"]) {
                        log(LogLevel::INFO) << " ---------------------- analyzing elements after pre-remeshing " << std::endl;
                        analyze_all_elements(ref_mesh, ref_mesh, 1.0, false,
                            ref_mesh_tot_vol / ref_mesh.n_logical_elements(), CONFORMING,
                                                            mesh_name + "_domain_elements_post_pre_remeshing.csv");
                        log(LogLevel::INFO) << " ------------------------------------" << std::endl;
                    }
                    fm.writeFile("input_ref_tet_mesh_post_remeshing.ovm", ref_mesh);
                }
            }else {
                log(LogLevel::ERROR)<<" ERROR - unknown pre-remeshing step "<<step<<std::endl;
            }
        }
    }else {
        log(LogLevel::WARNING)<<" WARNING - SKIPPING ALL PRE-PROCESSING, BYPASSING JSON PARAMETERS (--skip-pre-processing flag was set)"<<std::endl;
    }


    if(ref_mesh.contains_degenerate_or_flipped_elements(CONFORMING)){
        log(LogLevel::ERROR)<<" ERROR - ref mesh contains degenerate or flipped elements after setup"<<std::endl;

        for(auto elh: ref_mesh.elements()){
            if(ref_mesh.is_flipped_or_degenerate(elh, CONFORMING)){
                log(LogLevel::ERROR)<<" --> element "<<ref_mesh.element_to_str(elh)<<" with volume "<<CGAL_wrapper::exact_volume(ref_mesh, elh).to_double()<<std::endl;
            }
        }
        return -1;
    }

    if (check_boundary_holonomy(ref_mesh)) {
        log(LogLevel::ERROR) << " ERROR - ref mesh holonomy issue after mesh setup" << std::endl;
        return -1;
    }

    log(LogLevel::INFO) << " --- domain mesh boundary stats before optimization" << std::endl;
    print_boundary_quality_stats(ref_mesh);

    std::string remeshed_domain_mesh_path = parser.get<std::string>("--remeshed-domain-mesh-filepath");
    if(remeshed_domain_mesh_path != "") {

        std::string remeshed_boundary_conditions_path = parser.get<std::string>("--remeshed-boundary-conditions-filepath");
        if (remeshed_boundary_conditions_path == "") {
            log(LogLevel::ERROR)<<" ERROR - if you want to export the remeshed mesh you have to provide a path for the updated boundary conditions (with --remeshed-boundary-conditions-filepath)"<<std::endl;
            log(LogLevel::ERROR)<<" NOTE: even if the remeshing does NOT change the boundary conditions, it reorganizes the mesh's indexing, so the boundary conditions' vertex indices have to be updated accoridngly"<<std::endl;
            exit(EXIT_FAILURE);
        }
        log(LogLevel::IMPORTANT) << " - exporting pre-remeshed domain mesh to " << remeshed_domain_mesh_path <<std::endl<<
                     " and boundary conditions to "<< remeshed_boundary_conditions_path<< std::endl;

        export_boundary_conditions_from_mesh(mesh, remeshed_boundary_conditions_path);
        fm.writeFile(remeshed_domain_mesh_path, ref_mesh);
    }


    if(json_data["stop_after_mesh_setup"]){
        log(LogLevel::IMPORTANT)<<" -> stopping after mesh setup"<<std::endl;
        remove_validation_file(mesh_name + ".val");
        exit(EXIT_FAILURE);
    }


    int opti_result(-1);
    std::vector<std::pair<int, int>> mesh_size_history;
    std::vector<double> constraints_violation_history;

    if(json_data["solver_name"] != "truncated_newton_PCG") {

        log(LogLevel::ERROR) << " ERROR - solver " << json_data["solver_name"]
                  << " provided by opti-parameters.json is not compatible with mappings from files" << std::endl;
        exit(EXIT_FAILURE);
    }

    double setup_time_s = sw.total_duration();
    auto setup_time_h = setup_time_s / 3600.0;
    log(LogLevel::IMPORTANT)<<" - done with mesh setup in "<<setup_time_h<<" hours"<<std::endl;

    double max_allocated_optimization_time_s = std::max(0.0, 3600.0 * ((double)json_data["timeout_time_h"] - setup_time_h));

    OptimizerStats opti_stats;
    std::vector<std::vector<Vec3d>> vertex_trajectories;
    std::vector<std::vector<Vec3d>> vertex_target_paths;



    if(parser.get<bool>("--run-positional-opti")){
        //NOTE: the choice of virtual elements doesn't matter since we're not gonna use them anyway
        BoundaryGlueingRemeshingOptimizerT<IVL_B_3D_AD, AUG_SD_3D_AD> optimizer(ref_mesh, mesh, json_data, mesh_name);
        opti_result = optimizer.run_positional(max_allocated_optimization_time_s);
        opti_stats = optimizer.opti_stats();
        vertex_trajectories = optimizer.boundary_vertices_trajectories();
        vertex_target_paths = optimizer.boundary_vertices_target_paths();

    }else {
        if (json_data["global_wrapping_virtual_element"] == "AUG_SDE") {
            BoundaryGlueingRemeshingOptimizerT<IVL_B_3D_AD, AUG_SD_3D_AD> optimizer(ref_mesh, mesh, json_data, mesh_name);
            opti_result = optimizer.run(max_allocated_optimization_time_s);
            opti_stats = optimizer.opti_stats();
            vertex_trajectories = optimizer.boundary_vertices_trajectories();

        } else if (json_data["global_wrapping_virtual_element"] == "SDE") {
            BoundaryGlueingRemeshingOptimizerT<IVL_B_3D_AD, SD_3D_AD> optimizer(ref_mesh, mesh, json_data, mesh_name);
            opti_result = optimizer.run(max_allocated_optimization_time_s);
            opti_stats = optimizer.opti_stats();
            vertex_trajectories = optimizer.boundary_vertices_trajectories();

        } else if (json_data["global_wrapping_virtual_element"] == "AUG_FOF") {
            BoundaryGlueingRemeshingOptimizerT<IVL_B_3D_AD, AUG_FOF_NO_EPS_3D_AD> optimizer(ref_mesh, mesh, json_data, mesh_name);
            opti_result = optimizer.run(max_allocated_optimization_time_s);
            opti_stats = optimizer.opti_stats();
            vertex_trajectories = optimizer.boundary_vertices_trajectories();

        } else if (json_data["global_wrapping_virtual_element"] == "FOF") {
            BoundaryGlueingRemeshingOptimizerT<IVL_B_3D_AD, FOF_NO_EPS_3D_AD> optimizer(ref_mesh, mesh, json_data, mesh_name);
            opti_result = optimizer.run(max_allocated_optimization_time_s);
            opti_stats = optimizer.opti_stats();
            vertex_trajectories = optimizer.boundary_vertices_trajectories();

        } else if (json_data["global_wrapping_virtual_element"] == "AUG_DA") {
            BoundaryGlueingRemeshingOptimizerT<IVL_B_3D_AD, AUG_ANGLE_3D_PH_AD> optimizer(ref_mesh, mesh, json_data, mesh_name);
            opti_result = optimizer.run(max_allocated_optimization_time_s);
            opti_stats = optimizer.opti_stats();
            vertex_trajectories = optimizer.boundary_vertices_trajectories();

        } else {
            log(LogLevel::ERROR) << " ERROR - unhandled distortion element " << json_data["wrapping_distortion_element"] <<
                    std::endl;
            exit(EXIT_FAILURE);
        }
    }

    std::vector<int> n_verts_history, n_elems_history;
    for(auto p: mesh_size_history){
        n_verts_history.push_back(p.first);
        n_elems_history.push_back(p.second);
    }

    double runtime_s = sw.total_duration();

    ref_mesh.collect_garbage();
    mesh.collect_garbage();

    int final_v_count = mesh.n_logical_vertices();
    int final_e_count = mesh.n_logical_edges();
    int final_f_count = mesh.n_logical_faces();
    int final_c_count = mesh.n_logical_cells();

    double min_SDE, avg_SDE, max_SDE;
    check_full_SDE(mesh, ref_mesh, min_SDE, avg_SDE, max_SDE);

    double min_conf_dist, avg_conf_dist, max_conf_dist;
    check_full_conformal_distortion(mesh, ref_mesh, min_conf_dist, avg_conf_dist, max_conf_dist);

    double post_SDE_opti_min_SDE, post_SDE_opti_avg_SDE, post_SDE_opti_max_SDE;
    double post_conf_opti_min_conf_dist, post_conf_opti_avg_conf_dist, post_conf_opti_max_conf_dist;

    if(true) {
        int max_inner_iters(1000);
        auto target_vol = compute_total_volume(ref_mesh) / ref_mesh.n_logical_cells();
        // SDE opti
        optimize_mesh<NonConformingTetrahedralMesh, SD_3D_AD>(ref_mesh, target_vol, max_inner_iters, true, nullptr, false, json_data, LogConfig::get_min_level() > LogLevel::INFO);
        optimize_mesh<NonConformingTetrahedralMesh, SD_3D_AD>(mesh, ref_mesh, max_inner_iters, LogConfig::get_min_level() > LogLevel::INFO);
        check_full_SDE(mesh, ref_mesh, post_SDE_opti_min_SDE, post_SDE_opti_avg_SDE, post_SDE_opti_max_SDE);

        // conf dist opti
        optimize_mesh<NonConformingTetrahedralMesh, FOF_NO_EPS_3D_AD>(ref_mesh, target_vol, max_inner_iters, true, nullptr, false, json_data,
                                                                      LogConfig::get_min_level() > LogLevel::INFO);
        optimize_mesh<NonConformingTetrahedralMesh, FOF_NO_EPS_3D_AD>(mesh, ref_mesh, max_inner_iters, LogConfig::get_min_level() > LogLevel::INFO);
        check_full_conformal_distortion(mesh, ref_mesh, post_conf_opti_min_conf_dist, post_conf_opti_avg_conf_dist, post_conf_opti_max_conf_dist);
    }


    std::string domain_mesh_path = parser.get<std::string>("--output-domain-mesh-filepath");
    std::string codomain_mesh_path = parser.get<std::string>("--output-codomain-mesh-filepath");
    log(LogLevel::IMPORTANT)<<" - exporting resulting domain mesh to "<<domain_mesh_path<<std::endl;
    log(LogLevel::IMPORTANT)<<" -               and codomain mesh to "<<codomain_mesh_path<<std::endl;
    fm.writeFile(domain_mesh_path, ref_mesh);
    fm.writeFile(codomain_mesh_path, mesh);

    nlohmann::json output_json;
    output_json["mesh_name"] = mesh_name;
    output_json["opti_result"]     = opti_result;
    output_json["wrap_iterations"] = opti_stats.outer_iters_;
    output_json["first_snap_it"]   = opti_stats.first_snap_iter_;
    output_json["backup_mode"]   = opti_stats.backup_mode_;

    output_json["total_split_count"]    = opti_stats.total_split_count_;
    output_json["total_swap_count"]     = opti_stats.total_swap_count_;
    output_json["total_collapse_count"] = opti_stats.total_collapse_count_;

    output_json["runtime_s"] = runtime_s;
    output_json["total_opti_time_s"]            = opti_stats.total_opti_time_s_;
    output_json["total_upkeep_time_s"]          = opti_stats.total_upkeep_time_s_;
    output_json["total_remeshing_time_s"]       = opti_stats.total_remeshing_time_s_;
    output_json["total_unaccounted_for_time_s"] = opti_stats.total_unaccounted_for_time_s_;
    output_json["pre_processing_time_s"] = setup_time_s;

    output_json["n_init_vertices"] = initial_v_count;
    output_json["n_init_edges"] = initial_e_count;
    output_json["n_init_faces"] = initial_f_count;
    output_json["n_init_cells"] = initial_c_count;

    output_json["n_final_vertices"] = final_v_count;
    output_json["n_final_edges"] = final_e_count;
    output_json["n_final_faces"] = final_f_count;
    output_json["n_final_cells"] = final_c_count;

    output_json["min_boundary_SDE"] = min_boundary_dist;
    output_json["avg_boundary_SDE"] = avg_boundary_dist;
    output_json["max_boundary_SDE"] = max_boundary_dist;

    output_json["min_final_SDE"] = min_SDE;
    output_json["avg_final_SDE"] = avg_SDE;
    output_json["max_final_SDE"] = max_SDE;

    output_json["min_post_opti_SDE"] = post_SDE_opti_min_SDE;
    output_json["avg_post_opti_SDE"] = post_SDE_opti_avg_SDE;
    output_json["max_post_opti_SDE"] = post_SDE_opti_max_SDE;

    output_json["min_final_conf_dist"] = min_conf_dist;
    output_json["avg_final_conf_dist"] = avg_conf_dist;
    output_json["max_final_conf_dist"] = max_conf_dist;

    output_json["min_post_opti_conf_dist"] = post_conf_opti_min_conf_dist;
    output_json["avg_post_opti_conf_dist"] = post_conf_opti_avg_conf_dist;
    output_json["max_post_opti_conf_dist"] = post_conf_opti_max_conf_dist;

    output_json["deep_remeshing_iterations"] = opti_stats.deep_remeshing_iterations_;
    output_json["max_dihedral_angle_diff_deg_history"] = opti_stats.max_angle_diff_deg_history_;
    output_json["avg_dihedral_angle_diff_deg_history"] = opti_stats.avg_angle_diff_deg_history_;
    output_json["virtual_weight_scaling_history"] = opti_stats.virtual_weight_scaling_history_;
    output_json["opti_time_s_history"] = opti_stats.opti_time_s_history_;
    output_json["wrapped_ratio_history"] = opti_stats.wrapped_ratio_history_;
    output_json["newton_iters_history"] = opti_stats.newton_iters_history_;
    output_json["worst_quality_history"] = opti_stats.worst_quality_history_;
    output_json["min_relative_volume_history"] = opti_stats.min_relative_volume_history_;
    output_json["min_relative_edge_length_history"] = opti_stats.min_relative_edge_length_history_;


    output_json["input_min_worst_quality"] = input_min_quality;
    output_json["input_min_relative_volume"] = input_min_rel_vol;
    output_json["input_min_relative_edge_length"] = input_min_rel_edge_len;

    output_json["init_min_worst_quality"] = opti_stats.init_worst_quality_;
    output_json["init_min_relative_volume"] = opti_stats.init_min_relative_volume_;
    output_json["init_min_relative_edge_length"] = opti_stats.init_min_relative_edge_length_;

    output_json["final_min_worst_quality"] = opti_stats.final_worst_quality_;
    output_json["final_min_relative_volume"] = opti_stats.final_min_relative_volume_;
    output_json["final_min_relative_edge_length"] = opti_stats.final_min_relative_edge_length_;

    output_json["min_quality_before_remeshing_history"] = opti_stats.min_quality_before_remeshing_history_;
    output_json["avg_quality_before_remeshing_history"] = opti_stats.avg_quality_before_remeshing_history_;
    output_json["min_quality_after_remeshing_history"]  = opti_stats.min_quality_after_remeshing_history_;
    output_json["avg_quality_after_remeshing_history"]  = opti_stats.avg_quality_after_remeshing_history_;

    output_json["unweighted_initial_quality_energy_history"] = opti_stats.unweighted_initial_quality_energy_history_;
    output_json["unweighted_initial_primary_virtual_energy_history"] = opti_stats.unweighted_initial_primary_virtual_energy_history_;
    output_json["unweighted_initial_secondary_virtual_energy_history"] = opti_stats.unweighted_initial_secondary_virtual_energy_history_;

    output_json["unweighted_final_quality_energy_history"] = opti_stats.unweighted_final_quality_energy_history_;
    output_json["unweighted_final_primary_virtual_energy_history"] = opti_stats.unweighted_final_primary_virtual_energy_history_;
    output_json["unweighted_final_secondary_virtual_energy_history"] = opti_stats.unweighted_final_secondary_virtual_energy_history_;

    output_json["unweighted_final_quality_gradient_norm_history"] = opti_stats.unweighted_final_quality_gradient_norm_history_;
    output_json["unweighted_final_primary_virtual_gradient_norm_history"] = opti_stats.unweighted_final_primary_virtual_gradient_norm_history_;
    output_json["unweighted_final_secondary_virtual_gradient_norm_history"] = opti_stats.unweighted_final_secondary_virtual_gradient_norm_history_;

    output_json["final_total_gradient_max_coeff_history"] = opti_stats.final_total_gradient_max_coeff_history_;
    output_json["final_total_Hessian_max_coeff_history"] = opti_stats.final_total_Hessian_max_coeff_history_;

    output_json["full_quality_history"] = opti_stats.quality_history_;
    output_json["full_final_virtual_energies_history"] = opti_stats.final_virtual_energies_history_;

    output_json["n_verts_history"] = n_verts_history;
    output_json["n_elems_history"] = n_elems_history;


    export_boundary_as_triangle_mesh(mesh, "output_tet_mesh_boundary.obj");
    export_boundary_as_triangle_mesh(ref_mesh, "output_ref_tet_mesh_boundary.obj");


    log(LogLevel::IMPORTANT)<<" - exporting remeshing optimization stats to "<<json_stats_path<<std::endl;
    std::ofstream of(json_stats_path);
    if(!of.is_open()){
        log(LogLevel::ERROR)<<" --> couldn't open JSON file..."<<std::endl;
        exit(EXIT_FAILURE);
    }
    of<<output_json.dump()<<std::endl;
    of.close();

    remove_validation_file(mesh_name + ".val");

    if(json_data["enable_vertex_trajectory_mesh"]){
        std::string trajectories_filename("vertex_trajectories.obj");
        std::string targets_filename("vertex_targets.obj");

        log(LogLevel::IMPORTANT)<<" - exporting vertex trajectories to "<<trajectories_filename<<std::endl;
        log(LogLevel::IMPORTANT)<<" - exporting vertex targets to "<<targets_filename<<std::endl;

        std::ofstream(trajectories_filename) << encode_vertex_trajectories_as_OBJ(vertex_trajectories);
        std::ofstream(targets_filename) << encode_vertex_trajectories_as_OBJ(vertex_target_paths);
    }

    log(LogLevel::IMPORTANT)<<" ****************************************************** "<<std::endl;

    domain_mesh = ref_mesh;
    codomain_mesh = mesh;

    return opti_result;

}


void run_map_on_input_from_files(argparse::ArgumentParser& parser, nlohmann::json& json_data){
    log(LogLevel::IMPORTANT)<<" ****************************************************** "<<std::endl;
    log(LogLevel::IMPORTANT)<<" Running mapping method on file from mesh..."<<std::endl;


    OVM::IO::FileManager fm;
    fm.setVerbosityLevel(1);

    TetrahedralMesh domain_mesh, codomain_mesh;
    auto domain_mesh_filepath = parser.get<std::string>("--domain-mesh");

    size_t dot_loc = domain_mesh_filepath.find_last_of(".");
    std::string domain_mesh_extension = domain_mesh_filepath.substr(dot_loc + 1);

    if (domain_mesh_extension == "obj") {

        log(LogLevel::IMPORTANT)<<" --> loading OBJ mesh at "<<domain_mesh_filepath<<std::endl;
        if(!quick_and_dirty_OBJ_load(domain_mesh_filepath, domain_mesh)){
            log(LogLevel::ERROR)<<" ERROR - couldn't read domain mesh at "<<domain_mesh_filepath<<std::endl;
            exit(EXIT_FAILURE);
        }

    }else if (domain_mesh_extension == "ovm") {
        if(!fm.readFile(domain_mesh_filepath, domain_mesh)){
            log(LogLevel::ERROR)<<" ERROR - couldn't read domain mesh at "<<domain_mesh_filepath<<std::endl;
            exit(EXIT_FAILURE);
        }
    }else {
        log(LogLevel::ERROR)<<" ERROR - unhandled file extension "<<domain_mesh_extension<<std::endl;
        exit(EXIT_FAILURE);
    }

    std::string mesh_name = parser.get<std::string>("--mesh-name");
    log(LogLevel::IMPORTANT)<<" - loaded mesh "<<mesh_name<<std::endl;

    bool skip_pre_processing = parser.get<bool>("--skip-pre-processing");

    int bo_count = count_fully_boundary_elements(domain_mesh);
    log(LogLevel::IMPORTANT)<<" - boundary-only elements count before: "<<bo_count<<std::endl;
    //split_fully_boundary_elements(domain_mesh);
    if (!skip_pre_processing) {
        split_cord_edges_and_faces(domain_mesh);
        bo_count = count_fully_boundary_elements(domain_mesh);
        log(LogLevel::IMPORTANT)<<" - boundary-only elements count after cleanup: "<<bo_count<<std::endl;
        if(bo_count){
            log(LogLevel::ERROR)<<" ERROR - mesh still contains boundary-only elements"<<std::endl;

            std::string json_stats_path = parser.get<std::string>("--json-stats-filepath");
            log(LogLevel::ERROR)<<" - writing emergency error file instead of json stats"<<std::endl;
            std::ofstream of(json_stats_path);
            if(!of.is_open()){
                log(LogLevel::ERROR)<<" --> couldn't open JSON file..."<<std::endl;
                exit(EXIT_FAILURE);
            }
            of<<"boundary sanitizing failure"<<std::endl;
            of.close();
            exit(EXIT_FAILURE);
        }
    }

    codomain_mesh = domain_mesh;

    auto boundary_filepath = parser.get<std::string>("--boundary-txt");
    if (boundary_filepath != "") {
        std::cerr << " copying boundary conditions..." << std::endl;
        copy_boundary_conditions_to_mesh(codomain_mesh, boundary_filepath);
    } else {

        TetrahedralMesh temp_codomain_mesh;
        log(LogLevel::IMPORTANT)<<" - boundary conditions not specified with option '--boundary-txt', checking option '--codomain-mesh'"<<std::endl;
        auto codomain_mesh_filepath = parser.get<std::string>("--codomain-mesh");

        size_t dot_loc = codomain_mesh_filepath.find_last_of(".");
        std::string codomain_mesh_extension = codomain_mesh_filepath.substr(dot_loc + 1);

        if (codomain_mesh_extension == "obj") {

            log(LogLevel::IMPORTANT) << " --> loading OBJ codomain mesh at " << codomain_mesh_filepath << std::endl;
            if (!quick_and_dirty_OBJ_load(codomain_mesh_filepath, temp_codomain_mesh)) {
                log(LogLevel::ERROR) << " ERROR - couldn't read codomain mesh at " << codomain_mesh_filepath << std::endl;
                exit(EXIT_FAILURE);
            }

        }else if(!fm.readFile(codomain_mesh_filepath, temp_codomain_mesh)){
            log(LogLevel::ERROR)<<" ERROR - couldn't read codomain mesh at "<<codomain_mesh_filepath<<std::endl;
            exit(EXIT_FAILURE);
        }

        for (auto vh: temp_codomain_mesh.vertices()) {
            if (temp_codomain_mesh.is_boundary(vh)) {
                codomain_mesh.set_vertex(vh, temp_codomain_mesh.vertex(vh));
            }
        }
    }

    log(LogLevel::IMPORTANT) << " - read codomain file" << std::endl;


    if(parser.get<bool>("--invert-input")){
        invert_mesh(domain_mesh, false);
        invert_mesh(codomain_mesh, true);
    }



    domain_mesh.collect_garbage();
    codomain_mesh.collect_garbage();

    if(parser.get<bool>("--invert-problem")){
        if(invert_problem(domain_mesh, codomain_mesh)){
            log(LogLevel::ERROR)<<" ERROR - couldn't invert problem"<<std::endl;
            exit(EXIT_FAILURE);
        }
    }

    write_validation_file(mesh_name + ".val");

    run_3d_boundary_deformation_map(domain_mesh, codomain_mesh,
                         parser.get<std::string>("--mesh-name"),
                         parser.get<bool>("--skip-pre-processing"),
                         parser.get<std::string>("--json-stats-filepath"),
                         parser, json_data);
}


void run_map_on_test_mesh(argparse::ArgumentParser& parser, nlohmann::json& json_data){

    log(LogLevel::IMPORTANT)<<" ****************************************************** "<<std::endl;
    log(LogLevel::IMPORTANT)<<" testing non-conforming tetrahedral mesh stuff..."<<std::endl;


    NonConformingTetrahedralMesh ref_mesh, mesh;

    int length(10);
    length = parser.get<int>("--mesh-len");

    auto test_mesh(parser.get<std::string>("--test-mesh"));
    double torsion = parser.get<double>("--torsion");
    log(LogLevel::IMPORTANT)<<" testing opti with a "<<test_mesh<<" mesh of length "<<length<<" and torsion "<<torsion<<std::endl;

    bool compress_ref_mesh_to_cube = parser.get<bool>("--compress-ref-mesh-to-cube");

    ref_mesh = tet_sweep::SimpleMeshGen::generate_rod_mesh(length, 1.0, torsion);

    if(test_mesh == "minimal"){
        setup_minimal_mesh(ref_mesh, mesh);
    }else if(test_mesh == "minimal-non-star-shaped"){
        mesh = tet_sweep::SimpleMeshGen::generate_minimal_non_star_shaped_mesh();
        ref_mesh = tet_sweep::SimpleMeshGen::generate_minimal_non_star_shaped_domain_mesh();
    }else if(test_mesh == "straight"){
        mesh = tet_sweep::SimpleMeshGen::generate_rod_mesh(length, 1.0, torsion);

    }else if(test_mesh == "sine"){

        double period(1.0), x_scale(1.0), y_scale(1.0), z_scale(1.0);
        try {
            period = parser.get<double>("--sine-period");
            x_scale = parser.get<double>("--sine-x-scale");
            y_scale = parser.get<double>("--sine-y-scale");
            z_scale = parser.get<double>("--sine-z-scale");
        }catch(std::exception e){
            std::cout<<" parser exception: "<<e.what()<<std::endl;
            return;
        }

        mesh = tet_sweep::SimpleMeshGen::generate_sine_mesh(length, x_scale, torsion, period, z_scale);

    }else if(test_mesh == "spiral") {
        mesh = tet_sweep::SimpleMeshGen::generate_spiral_mesh(length, torsion,
                                                    parser.get<double>("--spiral-x-scale"), parser.get<double>("--spiral-turns"));

    }else if(test_mesh == "trefoil-knot") {

        mesh = tet_sweep::SimpleMeshGen::generate_trefoil_knot_mesh(length, torsion, parser.get<double>("--parametric-range"));

    } else if (test_mesh == "") {

        log(LogLevel::ERROR)<<" ERROR - no test mesh provided"<<std::endl;
        exit(EXIT_FAILURE);
    }else{
        log(LogLevel::ERROR)<<" ERROR - unhandled test mesh type "<<test_mesh<<std::endl;
        exit(EXIT_FAILURE);
        //setup_straight_mesh(ref_mesh, mesh, length);
    }

    log(LogLevel::IMPORTANT)<<" - mesh setup"<<std::endl;

    run_3d_boundary_deformation_map(ref_mesh, mesh,
                         parser.get<std::string>("--mesh-name"),
                         parser.get<bool>("--skip-pre-processing"),
                         parser.get<std::string>("--json-stats-filepath"),
                         parser, json_data);

}


int main(int argc, char** argv) {

    argparse::ArgumentParser parser;
    if(parse_args(argc, argv, parser)){
        std::cout<<" - arguments parsing failed"<<std::endl;
        return 1;
    }

    //only print error messages by default
    LogConfig::set_min_level(LogLevel::ERROR);


    try {

        std::ifstream f("opti-parameters.json");
        if(!f.is_open()){
            std::cout<<" --> couldn't open JSON file..."<<std::endl;
            exit(EXIT_FAILURE);
        }
        nlohmann::json json_data = nlohmann::json::parse(f);
        //std::cout<<" - parsed json file"<<std::endl;

        std::string log_level = json_data["log_level"];
        LogConfig::set_min_level(log_level);

        std::string success_stats_filepath = parser.get<std::string>("--orientation-stats-filepath");
        if(success_stats_filepath != ""){
            std::string codomain_mesh_filepath = parser.get<std::string>("--codomain-mesh");
            if(codomain_mesh_filepath == "") {
                std::cout<<" --> ERROR - success stats will be computed on the mesh given by --codomain-mesh"<<std::endl;
                exit(EXIT_FAILURE);
            }

            TetrahedralMesh codomain_mesh;
            IO::FileManager fm;
            fm.setVerbosityLevel(1);
    
            if(!fm.readFile(codomain_mesh_filepath, codomain_mesh)){
                std::cout<<" ERROR - couldn't read codomain mesh at "<<codomain_mesh_filepath<<std::endl;
                exit(EXIT_FAILURE);
            }

            std::string mesh_name = parser.get<std::string>("--mesh-name");
            extract_orientation_stats(codomain_mesh, mesh_name,success_stats_filepath);

        }else{

            if (parser.get<std::string>("--test-mesh") != "") {
                run_map_on_test_mesh(parser, json_data);
            
            }else{
                run_map_on_input_from_files(parser, json_data);
            }
        }
    }catch ( nlohmann::json_abi_v3_12_0::detail::type_error e) {
        std::cerr<<" parser error: "<<e.what()<<std::endl;
    }
    //system("leaks boundary_deformation_maps");

    /*std::cout<<"done. Type in anything to exit"<<std::endl;
    std::string ans;
    std::cin>>ans;*/
    return 0;

}




int parse_args(int argc, char** argv, argparse::ArgumentParser& parser){


    parser.add_argument("-h", "--help")
            .help("prints this help message")
            .flag();

    parser.add_argument("--skip-pre-processing")
            .help("skips all pre-processing stuff like cord-splitting, decimation, remeshing, etc.")
            .default_value(false)
            .flag();

    parser.add_argument("--domain-mesh")
            .help("specify the mesh to use as domain. If this is set, it will superseed the `--test-mesh` mode (see below).")
            .default_value("");

    parser.add_argument("--boundary-txt")
            .help("specify the boundary conditions to use for the codomain. Requires `--domain-mesh` to be set as well.")
            .default_value("");

    parser.add_argument("--codomain-mesh")
            .help("specify the mesh to use as codomain. Will only be checked if --boundary-txt is not set.")
            .default_value("");


    parser.add_argument("--invert-input")
            .help("geometrically inverts the input mesh (useful if input is left-hand-oriented)")
            .flag();


    parser.add_argument("--invert-problem")
            .help("inverts the codomain and domain mesh, and uses a single interior vertex for the codomain (that is then the domain of the optimization problem)")
            .flag();

    parser.add_argument("--mesh-name")
            .help("can be used to specify the name written in the json output.")
            .default_value("no_name");


    parser.add_argument("--run-positional-opti")
            .help("instead of the virtual-elements based deformation, this runs the boundary positions, equality-constrained version")
            .default_value(false)
            .flag();


    parser.add_argument("--test-mesh")
            .help("specify the test mesh boundary conditions. All options have a single interior vertex. Options are:\n"
                  " -------------------------------------------------------------- \n"
                  " - 'straight' a rectangular mesh of variable length \n"
                  " - 'minimal-non-star-shaped' five triangles forming a non-star-shaped space \n"
                  " - 'sine' the straight mesh shaped into a sine wave \n"
                  " - 'spiral' the straight mesh shaped into a spiral \n"
                  " -------------------------------------------------------------- \n"
            )
            .default_value("");


    parser.add_argument("--mesh-len")
            .help("specify the test mesh integer length")
            .default_value(10)
            .nargs(1)
            .scan<'i', int>();

    parser.add_argument("--compress-ref-mesh-to-cube")
            .help("makes the reference mesh a geometric cube, instead of a straight parallelepiped of length --mesh-len. WARNING: only works for tetrahedral spirals for now")
            .flag();

    parser.add_argument("--torsion")
            .help("specify the test mesh z-axis torsion (in degrees)")
            .default_value(0.0)
            .nargs(1)
            .scan<'g', double>();

    //------------------------------------------------------------------------------------ test mesh-specific parameters

    parser.add_argument("--sine-x-scale")
            .help("specify the sine test mesh horizontal scaling")
            .default_value(10.0)
            .nargs(1)
            .scan<'g', double>();

    parser.add_argument("--sine-y-scale")
            .help("specify the sine test mesh vertical scaling")
            .default_value(1.0)
            .nargs(1)
            .scan<'g', double>();

    parser.add_argument("--sine-z-scale")
            .help("specify the sine test mesh depth scaling")
            .default_value(1.0)
            .nargs(1)
            .scan<'g', double>();

    parser.add_argument("--sine-period")
            .help("specify the sine test mesh period")
            .default_value(1.0)
            .nargs(1)
            .scan<'g', double>();


    parser.add_argument("--spiral-turns")
            .help("specify the spiral test mesh number of turns")
            .default_value(1.0)
            .nargs(1)
            .scan<'g', double>();

    parser.add_argument("--spiral-x-scale")
            .help("specify the spiral test mesh x scaling (only for tet meshes)")
            .default_value(0.0)
            .nargs(1)
            .scan<'g', double>();


    parser.add_argument("--clover-center-radius")
            .help("specify the clover test mesh center radius")
            .default_value(1.0)
            .nargs(1)
            .scan<'g', double>();

    parser.add_argument("--clover-shift")
            .help("specify the clover test mesh shift between the bottom and top (mesh will be star-shaped if set to 0)")
            .default_value(0.0)
            .nargs(1)
            .scan<'g', double>();

    parser.add_argument("--clover-leaves")
            .help("specify the clover test mesh number of 'leaves'")
            .default_value(2)
            .nargs(1)
            .scan<'i', int>();

    parser.add_argument("--mouth-teeth-count")
            .help("specify the (rough) number of 'teeth' in the mouth test mesh")
            .default_value(5.0)
            .nargs(1)
            .scan<'g', double>();

    parser.add_argument("--mouth-opening")
            .help("specify how much the mouth test mesh should be 'open' (capped at 1)")
            .default_value(0.1)
            .nargs(1)
            .scan<'g', double>();

    parser.add_argument("--paraspiral-center-radius")
            .help("specify the radius of the center of the parametric spiral test mesh")
            .default_value(1.0)
            .nargs(1)
            .scan<'g', double>();

    parser.add_argument("--paraspiral-branches")
            .help("specify the number of branches of the parametric spiral test mesh")
            .default_value(3)
            .nargs(1)
            .scan<'i', int>();

    parser.add_argument("--paraspiral-branch-len")
            .help("specify the (conceptual) length of the parametric spiral test mesh")
            .default_value(1.0)
            .nargs(1)
            .scan<'g', double>();

    parser.add_argument("--compression-mode")
            .help("specify the number of branches of the parametric spiral test mesh")
            .default_value(0)
            .nargs(1)
            .scan<'i', int>();

    parser.add_argument("--i-to-b-v-count-ratio")
            .help("specify the ratio between interior to boundary vertices (for all compression tests)")
            .default_value(2.0)
            .nargs(1)
            .scan<'g', double>();

    parser.add_argument("--compression-disc-ratio")
            .help("specify the horizontal scaling for the compression disc test mesh")
            .default_value(1.0)
            .nargs(1)
            .scan<'g', double>();


    parser.add_argument("--parametric-range")
            .help("specify the the range of the parametric curve mesh should be generated, in [0, 1]")
            .default_value(0.9)
            .nargs(1)
            .scan<'g', double>();

    parser.add_argument("--parametric-width")
            .help("specify the the width of the parametric curve mesh. 0 will give a fully degenerate mesh, while values around 1 will probably generate some self-intersections. The best value here probably depends on the curve.")
            .default_value(0.5)
            .nargs(1)
            .scan<'g', double>();

    parser.add_argument("--parametric-z-scale")
            .help("specify the depth scaling for the parametric tet meshes. Setting it to 0 will leave it in the XY plane.")
            .default_value(1.0)
            .nargs(1)
            .scan<'g', double>();

    parser.add_argument("--complex-spiral-a")
            .help("specify the a parameter for the complex spiral (see curve function for details).")
            .default_value(5.0)
            .nargs(1)
            .scan<'g', double>();

    parser.add_argument("--complex-spiral-b")
            .help("specify the b parameter for the complex spiral (see curve function for details).")
            .default_value(60.0)
            .nargs(1)
            .scan<'g', double>();

    parser.add_argument("--complex-spiral-c")
            .help("specify the c parameter for the complex spiral (see curve function for details).")
            .default_value(20.0)
            .nargs(1)
            .scan<'g', double>();


    //------------------------------------------------------------------------------------------------- other parameters

    parser.add_argument("--json-stats-filepath")
            .help("specifiy the path of a json file that will contain information on the remeshing optimization process")
            .default_value("default_output.json")
            .nargs(1);


    parser.add_argument("--output-domain-mesh-filepath")
            .help("specifiy the path for the output domain mesh")
            .default_value("output_domain_mesh.ovm")
            .nargs(1);

    parser.add_argument("--output-codomain-mesh-filepath")
            .help("specify the path for the output codomain mesh")
            .default_value("output_codomain_mesh.ovm")
            .nargs(1);

    parser.add_argument("--remeshed-domain-mesh-filepath")
            .help("specify the path for the remeshed domain mesh (can be used to 'save' the higher quality domain meshes)")
            .default_value("")
            .nargs(1);

    parser.add_argument("--remeshed-boundary-conditions-filepath")
            .help("specify the path for the remeshed codomain boundary conditions (used with --remeshed-domain-mesh-filepath)")
            .default_value("")
            .nargs(1);

    parser.add_argument("--orientation-stats-filepath")
            .help("specify the path for a file containing very basic element orientation information about a mesh given by --codomain-mesh-filepath.")
            .default_value("")
            .nargs(1);

    parser.add_argument("-tet", "--volume")
            .help("[DEPRECATED] runs the remeshing optimization on a volumetric (tet) mesh")
            .flag();

    parser.add_argument("-v2")
            .help("[DEPRECATED] runs the v2 of the method (now the only available version)")
            .flag();



    try {
        parser.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << parser;
        std::exit(1);
    }

    if (parser["--help"] == true) {
        std::cerr<<parser;
        exit(1);
    }

    return 0;

}


void extract_orientation_stats(const TetrahedralMesh& mesh,
                           const std::string& mesh_name,
                           const std::string& stats_filepath) {


    std::ofstream output_file;
    output_file.open(stats_filepath);

    if(!output_file.is_open()){
        std::cout<<" ERROR - couldn't open stats file "<<stats_filepath<<std::endl;
        return;
    }

    auto bad_tets = CGAL_wrapper::find_degenerate_or_flipped_tets(mesh);

    std::cout<<" --> #degenerate cells: "<<bad_tets.degenerate.size()<<std::endl;
    std::cout<<" -->    #flipped cells: "<<bad_tets.flipped.size()<<std::endl;

    output_file << mesh_name << std::endl;
    output_file << mesh.n_logical_vertices() << std::endl;
    output_file << mesh.n_logical_edges() << std::endl;
    output_file << mesh.n_logical_faces() << std::endl;
    output_file << mesh.n_logical_cells() << std::endl;
    output_file << bad_tets.degenerate.size() << std::endl;
    output_file << bad_tets.flipped.size() << std::endl;
    output_file << 0 << std::endl;
    output_file.close();
}


/** Reads the OBJ and creates a tet mesh with a single interior vertex connected to all surface triangles */
bool quick_and_dirty_OBJ_load(const std::string& path, TetrahedralMesh& mesh) {

    std::cout<<" - loading .obj mesh from "<<path<<std::endl;
    std::vector<Vec3d> positions;
    std::vector<std::array<VertexHandle, 3>> faces;

    std::ifstream file(path);
    if (!file.is_open()) {
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::istringstream ss(line);
        std::string token;
        ss >> token;

        if (token == "v") {
            Vec3d p;
            ss >> p[0] >> p[1] >> p[2];
            positions.push_back(p);

        } else if (token == "f") {
            std::array<VertexHandle, 3> face;
            for (int i = 0; i < 3; ++i) {
                std::string entry;
                ss >> entry;
                int idx = std::stoi(entry.substr(0, entry.find('/')));
                face[i] = VertexHandle((idx < 0) ? (int)positions.size() + idx + 1: idx);
            }
            faces.push_back(face);
        }
    }

    std::cout<<" - found surface mesh with "<<positions.size()<<" vertices and "<<faces.size()<<" faces"<<std::endl;

    mesh.clear();

    auto interior_v = mesh.add_vertex({0,0,0});
    Vec3d centroid = {0,0,0};
    for(const auto& p: positions){
        mesh.add_vertex(p);
        centroid += p;
    }
    centroid /= positions.size();
    mesh.set_vertex(interior_v, centroid);

    for(const auto& f: faces){
        mesh.add_cell({interior_v, f[0], f[1], f[2]});
    }
    std::cout<<" -> tet mesh contains "<<mesh.n_vertices()<<" vertices and "<<mesh.n_cells()<<" cells"<<std::endl;

    return true;
}

/** Generated with Claude AI */
std::string encode_vertex_trajectories_as_OBJ(const std::vector<std::vector<Vec3d>>& trajectories){

    std::ostringstream oss;

    oss << "# Vertex trajectories exported as Wavefront OBJ\n";
    oss << "# " << trajectories.size() << " trajectory(ies)\n\n";

    // --- Pass 1: emit all vertex positions ---
    // OBJ indices are 1-based and global across the whole file.
    for (const auto& traj : trajectories)
        for (const auto& p : traj)
            oss << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";

    // --- Pass 2: emit polyline ('l') elements per trajectory ---
    oss << "\n";
    int vertex_offset = 1; // running 1-based index
    for (std::size_t t = 0; t < trajectories.size(); ++t)
    {
        const auto& traj = trajectories[t];
        if (traj.size() < 2)
        {
            // A single point can't form an edge; skip the 'l' record but
            // still advance the offset so indices stay consistent.
            vertex_offset += static_cast<int>(traj.size());
            continue;
        }

        oss << "g trajectory_" << t << "\n";
        oss << "l";
        for (std::size_t i = 0; i < traj.size(); ++i)
            oss << " " << (vertex_offset + static_cast<int>(i));
        oss << "\n";

        vertex_offset += static_cast<int>(traj.size());
    }

    return oss.str();
}


void extract_boundary_conditions_from_mesh(argparse::ArgumentParser& parser) {

    TetrahedralMesh domain_mesh;
    auto domain_mesh_filepath = parser.get<std::string>("--domain-mesh");

    IO::FileManager fm;
    fm.setVerbosityLevel(1);
    if(!fm.readFile(domain_mesh_filepath, domain_mesh)){
        std::cout<<" ERROR - couldn't read domain mesh at "<<domain_mesh_filepath<<std::endl;
        exit(EXIT_FAILURE);
    }

    auto boundary_filepath = parser.get<std::string>("--boundary-txt");
    export_boundary_conditions_from_mesh(domain_mesh, boundary_filepath);
}
