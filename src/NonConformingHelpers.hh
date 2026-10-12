#pragma once

#include "BadTetAnalysis.hh"

struct RemeshingStats {
    //-1: error, 0: no remeshing, 1: remeshed something
    int remeshed_something = -1;
    int split_count = 0;
    int swap_count = 0;
    int collapse_count = 0;
};


//This should probably be made a template (on the dimension) at some point
struct CustomArray{ //inheritance from Vec3d for operators

    double primary   = 0.0;
    double secondary = 0.0;
    double tertiary  = 0.0;

    CustomArray();
    CustomArray(const double& x);
    CustomArray(const double& p, const double& s, const double& t);
    CustomArray(const std::vector<double>& vec);

    // component-wise operators
    CustomArray operator*(const CustomArray& other);
    void operator*=(const CustomArray& other);

    CustomArray operator/(const CustomArray& other);
    void operator/=(const CustomArray& other);

    CustomArray operator+(const CustomArray& other);

    CustomArray operator-(const CustomArray& other);

    //replaces NaN and inf with 0 (based on std::isfinite)
    void make_finite();

    /*WeightArray& operator=(const WeightArray& other);
    WeightArray& operator=(const std::array<double,3>& other);*/
};

using EnergyArray = CustomArray;
using WeightArray = CustomArray;


std::ostream& operator<< (std::ostream& os, const CustomArray& w);


void collapse_relatively_short_boundary_edges(NonConformingTetrahedralMesh& domain_mesh,
                                              NonConformingTetrahedralMesh& codomain_mesh,
                                              double length_ratio_thresh);



template <class MeshT>
void
setup_regular_ref_element(double target_vol,
                          std::vector<typename MeshT::EigenVec>& P_ref){

    const int topo_dim(MeshT::topo_dim), geo_dim(MeshT::geo_dim);

    double ref_volume(1.0);

    P_ref.resize(MeshT::topo_dim);

    if(MeshT::topo_dim == 3){
        P_ref[0] << 0,0;
        P_ref[1] << 1,0;
        P_ref[2] << 0.5, 0.5 * std::sqrt(3.0);

        //NOTE: this is using the position matrix directly to compute the volume because the first vertex is at the origin
        //So basically the two other points' positions are the same as the edge vectors
        ref_volume = P_ref[1][0] * P_ref[2][1] - P_ref[2][0] * P_ref[1][1];

    }else{

        P_ref[0] <<  1, 1, 1;
        P_ref[1] << -1, 1,-1;
        P_ref[2] <<  1,-1,-1;
        P_ref[3] << -1,-1, 1;

        Eigen::Matrix<double, MeshT::geo_dim, MeshT::geo_dim> E_ref;
        for(int i(0); i<3; i++){
            E_ref.col(i) = P_ref[i+1] - P_ref[0];
        }
        ref_volume = E_ref.determinant() / 6.0;
    }
    //log(LogLevel::INFO)<<" - initial ref element volume = "<<ref_volume<<std::endl;
    //log(LogLevel::INFO)<<" -  target ref element volume = "<<target_vol<<std::endl;

    double scale_factor = std::pow(target_vol / ref_volume, 1./(double)(geo_dim));
    for(int i(0); i<P_ref.size(); i++){
        P_ref[i] *= scale_factor;
    }
}


template<typename VELEM_T>
int setup_parameters(const nlohmann::json& json_data, std::vector<double>& parameters) {

    if(VELEM_T::element_code() == "FOF_NO_EPS" || VELEM_T::element_code() == "AUG_FOF_NO_EPS") {
        parameters = {json_data["global_wrapping_virtual_FOF_conf_w"], json_data["global_wrapping_virtual_FOF_vol_w"] };

    }else if(VELEM_T::element_code() == "SD"){
        parameters = {json_data["global_wrapping_virtual_SDE_w"], 1.0, 0.0};

    }else if(VELEM_T::element_code() == "AUG_SD"){
        parameters = {json_data["global_wrapping_virtual_SDE_w"], 1.0, json_data["global_wrapping_virtual_SDE_sqrt_reg_eps"]};

    }else if(VELEM_T::element_code() == "AUG_DA"){
        parameters = {
            json_data["global_wrapping_virtual_DA_w"],
            0.0, //will be set later, a different value for each element
            json_data["global_wrapping_virtual_DA_sqrt_w"],
            json_data["global_wrapping_virtual_DA_pow"],
            json_data["global_wrapping_virtual_DA_bound_w"],
            0.0, //will be set later, a different value for each element
        };
        log(LogLevel::INFO)<<" - dihedral angle element regularization alpha: "<<parameters[2]<<std::endl;

    }else{
        log(LogLevel::ERROR)<<" ERROR - unhandled virtual element type "<<VELEM_T::element_name()<<std::endl;
        return -1;
    }

    return 0;
}


#warning TODO: move this to RemeshingOptimizer
/** return values
     * -2 virtual vertex position computation failure
     * -1 other errors
 * NOTE: mesh is not const because we temporarily shift some boundary vertices to compute the target virtual elements shape
 * NOTE: the mesh pointers are used for debugging / visualization */
template<typename VELEM_T_A, typename VELEM_T_B, typename VELEM_T_C>
int setup_holonomy_virtual_elements(TetrahedralMesh &mesh,
                                    const VertexPropertyT<Vec3d>& target_pos_prop,
                                    const OVM::VertexPropertyT<int> &global_vid_prop,
                                    const WeightArray& weight_scaling,
                                    const nlohmann::json &json_data,
                                    OVM::EdgePropertyT<Vec3d>& virtual_vertex_pos_prop,
                                    OVM::EdgePropertyT<double>& target_dihedral_angle_prop,
                                    COMISO::FiniteElementProblem &fe_problem,
                                    COMISO::FiniteElementSet<VELEM_T_A>& fe_virtual_elements_A,
                                    COMISO::FiniteElementSet<VELEM_T_B>& fe_virtual_elements_B,
                                    COMISO::FiniteElementSet<VELEM_T_C>& fe_virtual_elements_C,
                                    double dihedral_angle_energy_diff_bound_deg,
                                    TetrahedralMesh* p_focus_virtual_elements_mesh = nullptr,
                                    TetrahedralMesh* p_ref_focus_virtual_elements_mesh = nullptr) {

    log(LogLevel::INFO)<<" --------------------------------------------------------------------------------- "<<std::endl;
    log(LogLevel::INFO)<<" setting up holonomy virtual elements"<<std::endl;
    log(LogLevel::INFO)<<" - A: "<<VELEM_T_A::element_name()<<std::endl;
    log(LogLevel::INFO)<<" - B: "<<VELEM_T_B::element_name()<<std::endl;
    log(LogLevel::INFO)<<" - C: "<<VELEM_T_C::element_name()<<std::endl;
    log(LogLevel::INFO)<<" - weight scaling: "<<weight_scaling<<std::endl;


    int DA_param_idx = VELEM_T_A::element_code() == "AUG_DA" ? 0 : VELEM_T_B::element_code() == "AUG_DA" ? 1 : VELEM_T_C::element_code() == "AUG_DA" ? 2 : -1;

    log(LogLevel::INFO)<<" - dihedral angle energy diff bound: "<<dihedral_angle_energy_diff_bound_deg<<"°"<<std::endl;
    log(LogLevel::INFO)<<" - dihedral angle energy set idx: "<<DA_param_idx<<std::endl;


    auto log_level = LogLevel::DEBUG;

    const int virtual_vertices_index_start = fe_problem.n_unknowns();

    int boundary_edges_count(0);
    for (auto eh: mesh.edges()) {
        if (mesh.is_boundary(eh)) {
            boundary_edges_count++;
        }
    }

    double target_dihedral_angle_sum(0);
    for (auto eh: mesh.edges()) {
        if (mesh.is_boundary(eh)) {
            target_dihedral_angle_sum += std::abs(180 - target_dihedral_angle_prop[eh]);
        }
    }
    log(LogLevel::INFO)<<" - total target dihedral angle sum: "<<target_dihedral_angle_sum<<std::endl;

    int updated_n_vars = fe_problem.n_unknowns() + 3 * boundary_edges_count;
    fe_problem = COMISO::FiniteElementProblem(updated_n_vars);
    log(log_level)<<" --> updated n_vars as "<<virtual_vertices_index_start<<" + 3 * "<<boundary_edges_count<<" = "<<fe_problem.n_unknowns()<<std::endl;

    TetrahedralMesh virtual_elements_mesh, ref_virtual_elements_mesh;
    /*auto hinge_edge_prop = focus_virtual_elements_mesh.request_edge_property<int>("hinge_edge", -1);
    auto ref_hinge_edge_prop = ref_focus_virtual_elements_mesh.request_edge_property<int>("hinge_edge", -1);
    focus_virtual_elements_mesh.set_persistent(hinge_edge_prop);
    ref_focus_virtual_elements_mesh.set_persistent(ref_hinge_edge_prop);*/


    std::vector<std::vector<double>> parameters(3);
    if (setup_parameters<VELEM_T_A>(json_data, parameters[0]) ||
        setup_parameters<VELEM_T_B>(json_data, parameters[1]) ||
        setup_parameters<VELEM_T_C>(json_data, parameters[2])) {

        log(LogLevel::ERROR)<<" - error while setting up parameters"<<std::endl;
        return -1;
    }

    //global weight scaling
    double w = 0.5 / (double) boundary_edges_count;
    std::vector<double> normalized_w_sum = {0,0,0};

    log(LogLevel::INFO)<<" - global weight scaling: "<<w<<std::endl;

    //removed for now because not really checked anymore and no time to adapt to multi-FE sets
    /*double min_unw_e(std::numeric_limits<double>::max()), avg_unw_e(0.0), max_unw_e(0.0), max_edge_unw_e(0.0);
    double min_e(std::numeric_limits<double>::max()), avg_e(0.0), max_e(0.0), max_edge_e(0.0);
    EdgeHandle     max_unw_e_eh(-1),  max_e_eh(-1);
    HalfFaceHandle max_unw_e_hfh(-1), max_e_hfh(-1);*/

    int backup_virtual_pos_count(0);

    int virtual_vertex_idx(0);
    for (auto eh: mesh.edges()) {
        if (mesh.is_boundary(eh)) {
            auto from_v = mesh.edge(eh).from_vertex();
            auto to_v = mesh.edge(eh).to_vertex();

            /*if (!focus_edge && fixed_vertices[from_v] && fixed_vertices[to_v]) {
                continue;
            }*/

            log(log_level) <<
                           " ----------------------------------------------------------------------------------------------------------------------------"
                           << std::endl;
            log(log_level) << " ------------------------------- handling virtual vertex "<<virtual_vertex_idx<<" for boundary edge (" << mesh.edge(eh).from_vertex()
                           << ", " << mesh.edge(eh).to_vertex() << ")" << std::endl;

            std::vector<OVM::HalfFaceHandle> incident_boundary_halffaces;
            for (auto ehf_it = mesh.ehf_iter(eh); ehf_it.valid(); ehf_it++) {
                if (mesh.is_boundary(*ehf_it)) {
                    incident_boundary_halffaces.push_back(*ehf_it);
                }
            }
            auto &left_hf = incident_boundary_halffaces[0];
            auto &right_hf = incident_boundary_halffaces[1];

            log(log_level) << " -  left hf: " << halfface_to_str(mesh, left_hf) << std::endl;
            log(log_level) << " - right hf: " << halfface_to_str(mesh, right_hf) << std::endl;

            if(incident_boundary_halffaces.size() != 2){
                log(LogLevel::ERROR)<<" ERROR - boundary surface is not 2-manifold at boundary edge ("<<mesh.edge(eh).from_vertex()
                           << ", " << mesh.edge(eh).to_vertex() << "): incident to "<<incident_boundary_halffaces.size()<<" boundary faces"<<std::endl;
                    return -2;
            }

            std::unordered_map<VertexHandle, Vec3d> original_pos;
            //temporarily move vertices to their target positions
            for (auto hf: incident_boundary_halffaces) {
                for (auto vh: mesh.get_halfface_vertices(hf)) {
                    if (!original_pos.contains(vh)) {
                        original_pos.insert({vh, mesh.vertex(vh)});
                        mesh.set_vertex(vh, target_pos_prop[vh]);
                    }
                }
            }

#warning TODO: pre compute those and pass them as another prop
            //compute the target pos for the virtual elements
            auto ref_x = compute_optimal_virtual_vertex_pos(mesh, eh, {left_hf, right_hf});
            log(log_level) << " -  rotational ref mid-point: "<<ref_x << std::endl;

            //then move the vertices back to their current position
            for (auto hf: incident_boundary_halffaces) {
                for (auto vh: mesh.get_halfface_vertices(hf)) {
                    mesh.set_vertex(vh, original_pos[vh]);
                }
            }

            double edge_unw_e(0), edge_e(0);

            //----------------------------------------------------------------------------------------------------
            auto x = virtual_vertex_pos_prop[eh];
            log(log_level) << " -  virtual vertex pos: "<<x << std::endl;


            for (auto hf: incident_boundary_halffaces) {
                //auto verts = mesh.get_halfface_vertices(hf);
                log(log_level) << " --------------------------------------  handling halfface: "<<halfface_to_str(mesh, hf) << std::endl;


                //order vertices so the two end-vertices of the edge are at the beginning
                std::vector<VertexHandle> verts;
                for (auto hfhe_it = mesh.hfhe_iter(hf); hfhe_it.valid(); hfhe_it++) {
                    if (mesh.edge_handle(*hfhe_it) == eh) {
                        verts = {
                            mesh.from_vertex_handle(*hfhe_it),
                            mesh.to_vertex_handle(*hfhe_it),
                            mesh.to_vertex_handle(mesh.next_halfedge_in_halfface(*hfhe_it, hf))};
                    }
                }


                Vec3d c,n;
                compute_centroid_and_normal(mesh, hf, c, n);


                //NOTE: this is assuming that all elements use the same number of variables
                typename VELEM_T_A::VecV vx_virtual;
                typename VELEM_T_A::VecI vi_virtual;

                typename VELEM_T_A::VecC vc_A;
                typename VELEM_T_B::VecC vc_B;
                typename VELEM_T_C::VecC vc_C;

                std::vector<Eigen::Vector3d> P_ref;

                std::vector<VertexHandle> virtual_vertices;
                log(log_level)<<" - positions: "<<std::endl;
                std::vector<OVM::VertexHandle> local_verts;
                for (int i(0); i < verts.size(); i++) {
                    auto v = verts[i];
                    auto pos = mesh.vertex(v);
                    auto ref_pos = target_pos_prop[v];

                    log(log_level)<<" -     "<<pos<<std::endl;
                    log(log_level)<<" - ref "<<ref_pos<<std::endl;

                    for (int j(0); j < 3; j++) {
                        int local_idx = 3 * i + j;
                        int global_idx = 3 * global_vid_prop[v] + j;

                        vi_virtual[local_idx] = global_idx;
                        fe_problem.x()[global_idx] = pos[j];
                        vx_virtual[local_idx] = fe_problem.x()[global_idx];
                    }
                    P_ref.push_back(vec2vec(ref_pos));

                    virtual_vertices.push_back(ref_virtual_elements_mesh.add_vertex(ref_pos));
                    virtual_elements_mesh.add_vertex(pos);
                }


                for (int j(0); j < 3; j++) {
                    int local_idx = 3 * 3 + j;
                    int global_idx = virtual_vertices_index_start + 3 * virtual_vertex_idx + j;

                    vi_virtual[local_idx] = global_idx;
                    fe_problem.x()[global_idx] = x[j];
                    vx_virtual[local_idx] = fe_problem.x()[global_idx];
                }
                P_ref.push_back(vec2vec(ref_x));

                virtual_vertices.push_back(ref_virtual_elements_mesh.add_vertex(ref_x));
                virtual_elements_mesh.add_vertex(x);
                auto virtual_ch = ref_virtual_elements_mesh.add_cell(virtual_vertices);
                virtual_elements_mesh.add_cell(virtual_vertices);


                if(CGAL_wrapper::is_flipped(ref_virtual_elements_mesh, virtual_ch) || CGAL_wrapper::is_degenerate(ref_virtual_elements_mesh, virtual_ch) ){
                    log(LogLevel::ERROR)<<" ERROR - ref virtual element is either flipped or degenerate (detected by CGAL). Volume: "<<CGAL_wrapper::exact_volume(ref_virtual_elements_mesh, virtual_ch).to_double()<<std::endl;
                    return -2;
                }
                if(CGAL_wrapper::is_flipped(virtual_elements_mesh, virtual_ch) || CGAL_wrapper::is_degenerate(virtual_elements_mesh, virtual_ch) ){
                    log(LogLevel::ERROR)<<" ERROR - virtual element is either flipped or degenerate (detected by CGAL). Volume: "<<CGAL_wrapper::exact_volume(virtual_elements_mesh, virtual_ch).to_double()<<std::endl;
                    return -2;
                }
                double ref_virtual_cell_vol = compute_volume(ref_virtual_elements_mesh, virtual_ch);
                double virtual_cell_vol = compute_volume(virtual_elements_mesh, virtual_ch);
                if(ref_virtual_cell_vol <= 0 || virtual_cell_vol <= 0){
                    log(LogLevel::ERROR)<<" ERROR - either mesh or ref virtual element is either flipped or degenerate."<<std::endl;
                    log(LogLevel::ERROR)<<" - ref volume = "<<ref_virtual_cell_vol<<std::endl;
                    log(LogLevel::ERROR)<<" -     volume = "<<virtual_cell_vol<<std::endl;
                    return -2;
                }

                log(log_level)<<" -     x "<<x<<std::endl;
                log(log_level)<<" - ref x "<<ref_x<<std::endl;
                log(log_level)<<" -->    edge target dihedral angle: "<<target_dihedral_angle_prop[eh]<<std::endl;
                log(log_level)<<" --> element target dihedral angle: "<<(0.5*target_dihedral_angle_prop[eh])<<std::endl;


                if(DA_param_idx >= 0) {
                    parameters[DA_param_idx][1] = deg_to_rad(0.5 * target_dihedral_angle_prop[eh]);
                    log(log_level)<<" --> updated parameter set "<<DA_param_idx<<", index 1 (target theta) as "<<parameters[DA_param_idx][1]<<" rad = "<<rad_to_deg(parameters[DA_param_idx][1])<<"°"<<std::endl;

                    //without bound for base energy computation
                    parameters[DA_param_idx][4] = 0.0;
                    log(log_level)<<" --> updated parameter set "<<DA_param_idx<<", index 4 (bound weight) as "<<parameters[DA_param_idx][4]<<std::endl;
                    //log(LogLevel::IMPORTANT)<<" w update to zero"<<std::endl;

                }

                VELEM_T_A::compute_constants(parameters[0], P_ref, vc_A);
                VELEM_T_B::compute_constants(parameters[1], P_ref, vc_B);
                VELEM_T_C::compute_constants(parameters[2], P_ref, vc_C);

                VELEM_T_A::scale_weights(w * weight_scaling.primary,   vc_A);
                VELEM_T_B::scale_weights(w * weight_scaling.secondary, vc_B);
                VELEM_T_C::scale_weights(w * weight_scaling.tertiary,  vc_C);

                //VELEM_T::scale_weights(std::abs(180 - target_dihedral_angle_prop[eh]) / target_dihedral_angle_sum, vc_virtual);
                log(log_level)<<" --> final vc_A: "<<vc_A.transpose()<<std::endl;
                log(log_level)<<" --> final vc_B: "<<vc_B.transpose()<<std::endl;
                log(log_level)<<" --> final vc_A: "<<vc_C.transpose()<<std::endl;

                normalized_w_sum[0] += vc_A[0];
                normalized_w_sum[1] += vc_B[0];
                normalized_w_sum[2] += vc_C[0];

                VELEM_T_A elem_A;
                VELEM_T_B elem_B;
                VELEM_T_C elem_C;

                std::vector<double> energies = {
                    vc_A[0] > 0.0 ? elem_A.eval_f(vx_virtual, vc_A) : 0.0,
                    vc_B[0] > 0.0 ? elem_B.eval_f(vx_virtual, vc_B) : 0.0,
                    vc_C[0] > 0.0 ? elem_C.eval_f(vx_virtual, vc_C) : 0.0
                };

                log(log_level)<<" -->            virtual energy A = "<<energies[0]<<std::endl;
                log(log_level)<<" -->            virtual energy B = "<<energies[1]<<std::endl;
                log(log_level)<<" -->            virtual energy C = "<<energies[2]<<std::endl;
                log(log_level)<<" --> unweighted virtual energy A = "<<(energies[0] / vc_A[0])<<std::endl;
                log(log_level)<<" --> unweighted virtual energy B = "<<(energies[1] / vc_B[0])<<std::endl;
                log(log_level)<<" --> unweighted virtual energy C = "<<(energies[2] / vc_C[0])<<std::endl;

                if (dihedral_angle_energy_diff_bound_deg > 0) {

                    if (DA_param_idx != 1) {
                        log(LogLevel::ERROR)<<" ERROR - enabling the dihedral energy bound assumes it's the secondary virtual element set"<<std::endl;
                        return -1;
                    }

                    log(log_level)<<" ------------------------------------------------"<<std::endl;
                    double unw_DA_e = energies[DA_param_idx] / vc_B[0];
                    log(log_level)<<" -->   weighted energy = "<<std::setprecision(8)<<energies[DA_param_idx]<<std::endl;
                    log(log_level)<<" --> unweighted energy = "<<unw_DA_e<<std::endl;

                    double dihedral_sum_deg = compute_dihedral_angle_sum_deg(mesh, eh);
                    double current_edge_dihedral_angle_deg = 360.0 - dihedral_sum_deg;
                    log(log_level)<<" --> current   edge dihedral angle = "<<current_edge_dihedral_angle_deg<<std::endl;

                    Vec3d a = {vx_virtual[0], vx_virtual[1], vx_virtual[2]};
                    Vec3d b = {vx_virtual[3], vx_virtual[4], vx_virtual[5]};
                    Vec3d virtual_n = (x-a).cross(b-a).normalized();
                    log(log_level)<<" --> boundary face normal: "<<n<<std::endl;
                    log(log_level)<<" -->  virtual face normal: "<<virtual_n<<std::endl;
                    double current_element_dihedral_angle_deg = rad_to_deg(std::acos(std::max(-1.0, std::min(1.0, n.dot(-virtual_n)))));
                    log(log_level)<<" --> current element dihedral angle = "<<current_element_dihedral_angle_deg<<std::endl;


                    vc_B[4] = vc_B[0] * (double)json_data["global_wrapping_virtual_DA_bound_w"];
                    vc_B[5] =  std::abs(current_element_dihedral_angle_deg - 0.5 * target_dihedral_angle_prop[eh]) + dihedral_angle_energy_diff_bound_deg;

                    log(log_level)<<" --> diff to target = "<<std::abs(current_element_dihedral_angle_deg - 0.5 * target_dihedral_angle_prop[eh])<<std::endl;

                    log(log_level)<<" --> bound set to = "<<vc_B[5]<<std::endl;

                    energies[1] = elem_B.eval_f(vx_virtual, vc_B);
                    log(log_level)<<" --> updated (weighted) dihedral angle energy (with bound) = "<<std::setprecision(8)<<energies[1]<<std::endl;

                    //auto t = elem_B.max_feasible_step(vx_virtual, -vx_virtual, vc_B);
                    //log(LogLevel::IMPORTANT)<<" --> max feasible step = "<<t<<std::endl;

                    log(log_level)<<" ------------------------------------------------"<<std::endl;
                }

                log(log_level)<<" ----------------------------------------------"<<std::endl;


                /*auto t = elem.max_feasible_step(vx_virtual, -vx_virtual, vc_virtual);
                log(log_level)<<" --> max feasible step = "<<t<<std::endl;*/

                if (!std::isfinite(energies[0]) || !std::isfinite(energies[1]) || !std::isfinite(energies[2])) {
                    log(LogLevel::ERROR) << " ERROR - infinite energy for virtual element " << virtual_vertex_idx <<
                            " for halfface " <<
                            halfface_to_str(mesh, hf) << " incident to edge " << edge_to_str(mesh, eh) << std::endl;
                    //log(LogLevel::ERROR) << " ERROR - infinite energy for virtual element "<<virtual_vertex_idx<<" for halfface " <<
                    //                     halfface_to_str(mesh, hf) <<" incident to edge "<<edge_to_str(mesh, eh)<<", even after recomputing virtual vertex pos" <<std::endl;

                    auto optimal_x = compute_optimal_virtual_vertex_pos(mesh, eh, {left_hf, right_hf});
                    auto old_x = compute_rotational_mid_point(mesh, eh, {left_hf, right_hf});
                    std::cerr << " -    current x: " << virtual_vertex_pos_prop[eh] << std::endl;
                    std::cerr << " -    optimal x: " << optimal_x << std::endl;
                    std::cerr << " - rotational x: " << old_x << std::endl;

                    log(LogLevel::ERROR) << " -->            virtual energy A = " << energies[0] << std::endl;
                    log(LogLevel::ERROR) << " -->            virtual energy B = " << energies[1] << std::endl;
                    log(LogLevel::ERROR) << " -->            virtual energy C = " << energies[2] << std::endl;
                    log(LogLevel::ERROR) << " --> unweighted virtual energy A = " << (energies[0] / vc_A[0]) << std::endl;
                    log(LogLevel::ERROR) << " --> unweighted virtual energy B = " << (energies[1] / vc_B[0]) << std::endl;
                    log(LogLevel::ERROR) << " --> unweighted virtual energy C = " << (energies[2] / vc_C[0]) << std::endl;
                    log(LogLevel::ERROR) << " ----------------------------------------------" << std::endl;

                    return -1;


                    /*log(LogLevel::WARNING) << " WARNING - infinite energy for virtual element "<<virtual_vertex_idx<<" incident to halfface " <<
                                           halfface_to_str(mesh, hf) <<", recomputing virtual vertex pos" <<std::endl;*/

                    auto backup_x = compute_optimal_virtual_vertex_pos(mesh, eh, {left_hf, right_hf});
                    virtual_vertex_pos_prop[eh] = backup_x;
                    backup_virtual_pos_count++;
                    break;

                }

                if (vc_A[0] > 0.0) {
                    fe_virtual_elements_A.instances().add_element(vi_virtual, vc_A);
                }
                if (vc_B[0] > 0.0) {
                    fe_virtual_elements_B.instances().add_element(vi_virtual, vc_B);
                }
                if (vc_C[0] > 0.0) {
                    fe_virtual_elements_C.instances().add_element(vi_virtual, vc_C);
                }
                log(log_level) << " -  added elements to their respective sets " << std::endl;


                /*instances_to_add.push_back({vi_virtual, vc_virtual});

                fe_virtual_reg_elements.instances().add_element(vi_reg, vc_reg);

                //NOTE: that's an approximation, assuming that the conformal part is dominant if the virtual element is FOF
                auto unw_e = e / vc_virtual[0];

                log(log_level) << " ->  virtual element energy = " << e << ", unweighted = " << unw_e << std::endl;

                min_unw_e = std::min(min_unw_e, unw_e);
                min_e = std::min(min_e, e);

                avg_unw_e += unw_e;
                avg_e += e;

                if (unw_e >= max_unw_e) {
                    max_unw_e = unw_e;
                    max_unw_e_hfh = hf;
                }
                if (e >= max_e) {
                    max_e = e;
                    max_e_hfh = hf;
                }

                edge_unw_e += unw_e;
                edge_e += e;*/

            }
            //----------------------------------------------------------------------------------------------------



            /*if (instances_to_add.size() != 2) {
                log(LogLevel::ERROR)<<" ERROR - there should be two instances to add but there are "<<instances_to_add.size()<<std::endl;
                return -1;
            }
            for (auto& pair : instances_to_add) {
                fe_virtual_elements.instances().add_element(pair.first, pair.second);
            }

            if (edge_unw_e >= max_edge_unw_e) {
                max_edge_unw_e = edge_unw_e;
                max_unw_e_eh = eh;
            }
            if (edge_e >= max_edge_e) {
                max_edge_e = edge_e;
                max_e_eh = eh;
            }*/

            virtual_vertex_idx++;
        }
    }


    //avg_unw_e /= fe_virtual_elements.n_instances();
    //avg_e     /= fe_virtual_elements.n_instances();

    log(LogLevel::INFO)<<"---------------------------------------------"<<std::endl;
    log(LogLevel::INFO)<<" -> holonomy virtual elements setup done"<<std::endl;
    log(LogLevel::INFO)<<" ->            normalized sum weight for set A = "<<normalized_w_sum[0]<<std::endl;
    log(LogLevel::INFO)<<" ->            normalized sum weight for set B = "<<normalized_w_sum[1]<<std::endl;
    log(LogLevel::INFO)<<" ->            normalized sum weight for set C = "<<normalized_w_sum[2]<<std::endl;
    log(LogLevel::INFO)<<" -> recomputed virtual vertex position "<<backup_virtual_pos_count<<" times"<<std::endl;
    log(LogLevel::INFO)<<" ---------------------------------------------"<<std::endl;
    /*log(LogLevel::INFO)<<" ------------- edge elements stats -----------"<<std::endl;
    log(LogLevel::INFO)<<" ---------------------------------------------"<<std::endl;
    log(LogLevel::INFO)<<" NOTE: unweighted energies are not accurate if the virtual element is FOF"<<std::endl;
    log(LogLevel::INFO)<<"       min / max / avg unweighted energy: "<<min_unw_e<<" / "<<max_unw_e<<" / "<<avg_unw_e<<std::endl;
    log(LogLevel::INFO)<<" halfface with max unweighted energy was: "<<halfface_to_str(mesh, max_unw_e_hfh)<<std::endl;
    log(LogLevel::INFO)<<"     edge with max unweighted energy was: "<<edge_to_str(mesh, max_unw_e_eh)<<" with sum = "<<max_edge_unw_e<<std::endl;
    log(LogLevel::INFO)<<" ---------------------------------------------"<<std::endl;
    log(LogLevel::INFO)<<"       min / max / avg energy: "<<min_e<<" / "<<max_e<<" / "<<avg_e<<std::endl;
    log(LogLevel::INFO)<<" halfface with max energy was: "<<halfface_to_str(mesh, max_e_hfh)<<std::endl;
    log(LogLevel::INFO)<<"     edge with max energy was: "<<edge_to_str(mesh, max_e_eh)<<" with sum = "<<max_edge_e<<std::endl;
    log(LogLevel::INFO)<<" ---------------------------------------------"<<std::endl;*/

    virtual_elements_mesh.collect_garbage();
    ref_virtual_elements_mesh.collect_garbage();

    if (p_focus_virtual_elements_mesh) {
        *p_focus_virtual_elements_mesh = virtual_elements_mesh;
    }
    if (p_ref_focus_virtual_elements_mesh) {
        *p_ref_focus_virtual_elements_mesh = ref_virtual_elements_mesh;
    }

    if (fe_virtual_elements_A.n_instances()) {
        fe_problem.add_set(&fe_virtual_elements_A);
    }
    if (fe_virtual_elements_B.n_instances()) {
        fe_problem.add_set(&fe_virtual_elements_B);
    }
    if (fe_virtual_elements_C.n_instances()) {
        fe_problem.add_set(&fe_virtual_elements_C);
    }

    //manual_break("breaking after holo setup");
    log(LogLevel::INFO)<<" --------------------------------------------------------------------------------- "<<std::endl;

    //manual_stop("for debugging");
    return 0;

}





using BOUNDARY_EDGE_ELEM = EV3D_AD;
/** NOTE: constraints are NOT set for now */
int setup_boundary_edge_elements(TetrahedralMesh &mesh,
                                 const VertexPropertyT<Vec3d> &target_pos_prop,
                                 const OVM::VertexPropertyT<int> &global_vid_prop,
                                 const double &edge_w,
                                 COMISO::FiniteElementProblem &fe_problem,
                                 COMISO::FiniteElementSet<BOUNDARY_EDGE_ELEM> &fe_edge_elements,
                                 std::vector<COMISO::LinearConstraint> &edge_eq_constraints);




template<typename MeshT, typename ElementT>
int optimize_mesh(MeshT &mesh,
                  MeshT &ref_mesh,
                  int max_inner_iters,
                  bool run_silently = true){

    using EigenVec = typename MeshT::EigenVec;
    const int geo_dim = MeshT::geo_dim;

    std::vector<EigenVec> P_ref;

    const double Dirichlet_w(1.0), Dirichlet_scaling(1.0), Dirichlet_eps(0.0);
    const double conf_w(1.0), vol_w(0.0);
    int n_vars(mesh.n_logical_vertices() * geo_dim);

    std::string ref_opti_element_type = ElementT::element_code();

    log(LogLevel::INFO) << " - using element type: " << ElementT::element_name() << std::endl;

    double ref_to_mesh_scaling_factor = std::pow(compute_total_volume(mesh) / compute_total_volume(ref_mesh), 1.0 / 3.0);
    //std::cout<<" - ref to mesh scaling factor = "<<ref_to_mesh_scaling_factor<<std::endl;

    auto global_vid_prop = mesh. template request_vertex_property<int>();
    int global_vid(0);
    for(auto v: mesh.vertices()){
        global_vid_prop[v] = global_vid;
        global_vid++;
    }

    COMISO::FiniteElementProblem fe_problem(n_vars);


    std::vector<COMISO::LinearConstraint> eq_constraints;
    eq_constraints.reserve(n_vars);
    //COMISO::FiniteElementSet <BaseElemT> fe_barrier("barrier elements");
    COMISO::FiniteElementSet <ElementT> fe_barrier(ElementT::element_name());


    double max_energy(0), max_w(0);
    double w_sum(0);

    double w = 1; //todo here
    //w = mesh.element_volume(elh) / target_vol_per_element;
    w /= mesh.n_logical_elements();

    std::vector<double> parameters;
    if(ref_opti_element_type == "FOF_NO_EPS"){
        parameters = {conf_w, vol_w};
    }else if (ref_opti_element_type == "SD") {
        parameters = {Dirichlet_w, Dirichlet_scaling, Dirichlet_eps};
    }else {
        log(LogLevel::ERROR)<<" ERROR - unhandled relaxation element "<<ElementT::element_name()<<std::endl;
        return -1;
    }

    //add the barrier elements
    for (auto elh: mesh.elements()) {

        //log(LogLevel::INFO)<<" --------------------------------------------------- handling element "<<mesh.element_to_str(elh)<<std::endl;
        std::vector<EigenVec> P_ref(4);

        w_sum += w;
        max_w = std::max(w, max_w);

        typename ElementT::VecI vi_barrier;
        typename ElementT::VecC vc_barrier;
        typename ElementT::VecV vx_barrier;

        auto verts = mesh.get_element_vertices(elh);

        for (int local_vid(0); local_vid < verts.size(); local_vid++) {
            const auto &vh = verts[local_vid];
            P_ref[local_vid] = ref_to_mesh_scaling_factor * vec2vec(ref_mesh.vertex(vh));

            int global_id = geo_dim * global_vid_prop[vh];
            for (int i(0); i < geo_dim; i++) {
                vi_barrier[geo_dim * local_vid + i] = global_id + i;
            }

            for (int i(0); i < geo_dim; i++) {
                fe_problem.x()[global_id + i] = mesh.vertex(vh)[i];
                vx_barrier[geo_dim * local_vid + i] = fe_problem.x()[global_id + i];
            }
        }

        auto setup_res = ElementT::compute_constants(parameters, P_ref, vc_barrier);

        if (setup_res) {
            log(LogLevel::ERROR) << " ERROR: constants couldn't be computed for face: " << elh << std::endl;
            return -1;
        }

        ElementT elem;
        double barrier_e = elem.eval_f(vx_barrier, vc_barrier);

        auto vol = compute_volume(mesh, elh);

        if (!std::isfinite(barrier_e)) {
            log(LogLevel::ERROR) << " ERROR: infinite energy for element "<<cell_to_str(mesh, elh)<<std::endl;
            return -1;
        }

        max_energy = std::max(max_energy, barrier_e);

        fe_barrier.instances().add_element(vi_barrier,
                                           vc_barrier);
    }

    fe_problem.add_set(&fe_barrier);

    int constrained_count(0);
    //constrain the fixed vertices
    for(auto v: mesh.vertices()) {
        if(mesh.is_boundary(v)) {

            auto pos = mesh.vertex(v);
            int var_idx = geo_dim * global_vid_prop[v];

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

    std::vector<COMISO::NConstraintInterface*> constraints_pointers;
    constraints_pointers.reserve(constraints_pointers.size());
    for( auto& c : eq_constraints) {
        constraints_pointers.push_back(&c);
    }

    if(!run_silently) {
        log(LogLevel::INFO) << " energies after optimization: " << std::endl;
        fe_problem.print_objectives();
    }

    double initial_energy = fe_problem.eval_f(fe_problem.x().data());
    if(initial_energy == std::numeric_limits<double>::infinity()){
        log(LogLevel::INFO)<<" -- infinite initial energy, stopping there."<<std::endl;
        return -1;
    }

    bool use_truncated_newton(true);
    bool converged(false);
    bool feasible_solution_found(false);
    if(use_truncated_newton){
        COMISO::TruncatedNewtonPCG tn(1e-4);
        tn.always_update_preconditioner() = true;
        tn.max_iters() = max_inner_iters;
        tn.max_pcg_iters() = 100;
        tn.allow_warmstart() = false;
        tn.tolerance_gdx() = 1e-4;
        tn.set_silent(run_silently);
        tn.set_eps_constraints_violation(1e-9);//1e-9

        tn.solve_projected_normal_equation(&fe_problem, constraints_pointers);
        converged = tn.converged();
        feasible_solution_found = tn.feasible_solution_found();
    }else {
        COMISO::NewtonSolver ns(1e-1, 1e-9, max_inner_iters);
        ns.solve_infeasible_start(&fe_problem, constraints_pointers);
        converged = ns.converged();
    }

    if(!run_silently) {
        log(LogLevel::INFO) << " energies after optimization: " << std::endl;
        fe_problem.print_objectives();
    }

    for(auto v: mesh.vertices()){
        if(!mesh.is_boundary(v)) {
            int global_id = geo_dim * global_vid_prop[v];
            typename MeshT::PointT new_pos;
            for (int i(0); i < MeshT::geo_dim; i++) {
                new_pos[i] = fe_problem.x()[global_id + i];
            }
            mesh.set_vertex(v, new_pos);
        }
    }

    return 0;
}




template<typename MeshT, typename ElementT>
int optimize_mesh(MeshT &mesh,
                  double target_vol_per_element,
                  int max_inner_iters,
                  const VertexPropertyT<bool> &constrained_prop,
                  const VertexPropertyT<Vec3d>* p_target_pos_prop,
                  bool enable_boundary_edge_elements,
                  const nlohmann::json &json_data,
                  bool run_silently = true){

    using EigenVec = typename MeshT::EigenVec;
    const int geo_dim = MeshT::geo_dim;
    const int topo_dim = MeshT::topo_dim;

    /*if(MeshT::topo_dim != 3){
        log(LogLevel::INFO)<<" ERROR - ref mesh opti is not implemented for tet meshes yet"<<std::endl;
        exit(EXIT_FAILURE);
    }*/

    std::vector<EigenVec> P_ref;
    setup_regular_ref_element<MeshT>(target_vol_per_element, P_ref);

    const double Dirichlet_w(1), Dirichlet_barrier_w(0), Dirichlet_barrier_min(1e1);
    int n_vars(mesh.n_logical_vertices() * geo_dim);

    double fof_epsilon(0);
    double conf_w(json_data["fof_conf_w"]);
    double vol_w(json_data["fof_vol_w"]);

    std::string ref_opti_element_type = ElementT::element_code();

    if(!run_silently) {
        log(LogLevel::INFO) << " - using element type: " << ElementT::element_name() << std::endl;
        if(p_target_pos_prop){
            log(LogLevel::INFO) << " - using target pos prop as constraints (if any) " << std::endl;
        }
    }


    auto global_vid_prop = mesh. template request_vertex_property<int>();
    int global_vid(0);
    for(auto v: mesh.vertices()){
        global_vid_prop[v] = global_vid;
        global_vid++;
    }

    COMISO::FiniteElementProblem fe_problem(n_vars);


    std::vector<COMISO::LinearConstraint> eq_constraints;
    eq_constraints.reserve(n_vars);
    //COMISO::FiniteElementSet <BaseElemT> fe_barrier("barrier elements");
    COMISO::FiniteElementSet <ElementT> fe_barrier(ElementT::element_name());


    double max_energy(0), max_w(0);
    double w_sum(0);

    double w = 1; //todo here
    //w = mesh.element_volume(elh) / target_vol_per_element;
    w /= mesh.n_logical_elements();

    bool enable_smoothing_IVL_bounds = json_data["enable_smoothing_IVL_bounds"];

    std::vector<double> parameters = {w * Dirichlet_w, Dirichlet_barrier_w, Dirichlet_barrier_min};

    if(ref_opti_element_type == "SD"){
        parameters = {w * 1.0, 1.0, 0.0};
    }else if(ref_opti_element_type == "FOF_NO_EPS"){
        parameters = {w * 1.0, 0.0};
    }else if (ref_opti_element_type == "IVL_B") {

        auto vol_bound = (double)json_data["global_wrapping_IVL_rel_vol_bound"] * target_vol_per_element;

        if(enable_smoothing_IVL_bounds){

            if (!run_silently) {
                log(LogLevel::INFO) << " --> using IVL bounds" << std::endl;
                log(LogLevel::INFO) << "        IVL bound = " << json_data["global_wrapping_IVL_bound"] << std::endl;
                log(LogLevel::INFO) << " rel volume bound = " << json_data["global_wrapping_IVL_rel_vol_bound"] << std::endl;
                log(LogLevel::INFO) << "       avg volume = " << target_vol_per_element << std::endl;
                log(LogLevel::INFO) << "     volume bound = " << vol_bound << std::endl;
            }

            parameters = {
                    json_data["wrapping_relaxation_IVL_w"],
                    json_data["global_wrapping_IVL_bound_w"],
                    json_data["global_wrapping_IVL_bound"],
                    json_data["global_wrapping_IVL_vol_w"],
                    vol_bound,
                    target_vol_per_element
            };

        }else {
            parameters = {
                    json_data["wrapping_relaxation_IVL_w"],
                    0, //json_data["global_wrapping_IVL_bound_w"],
                    1, //json_data["global_wrapping_IVL_bound"],
                    0, //json_data["global_wrapping_IVL_vol_w"],
                    1, //vol_bound
                    1 //target_vol_per_element
            };
        }
        //log(LogLevel::INFO)<<" - IVL weight = "<<parameters[0]<<std::endl;
    }else {
        log(LogLevel::ERROR)<<" ERROR - unhandled relaxation element "<<ElementT::element_name()<<std::endl;
        return -1;
    }

    int backup_bound_use_count(0), backup_vol_bound_use_count(0);
    double max_bound(0);

    //add the barrier elements
    for (auto elh: mesh.elements()) {

        //log(LogLevel::INFO)<<" --------------------------------------------------- handling element "<<mesh.element_to_str(elh)<<std::endl;

        w_sum += w;
        max_w = std::max(w, max_w);

        typename ElementT::VecI vi_barrier;
        typename ElementT::VecC vc_barrier;
        typename ElementT::VecV vx_barrier;

        auto verts = mesh.get_element_vertices(elh);

        for (int local_vid(0); local_vid < verts.size(); local_vid++) {
            const auto &vh = verts[local_vid];

            int global_id = geo_dim * global_vid_prop[vh];
            for (int i(0); i < geo_dim; i++) {
                vi_barrier[geo_dim * local_vid + i] = global_id + i;
            }

            for (int i(0); i < geo_dim; i++) {
                fe_problem.x()[global_id + i] = mesh.vertex(vh)[i];
                vx_barrier[geo_dim * local_vid + i] = fe_problem.x()[global_id + i];
            }
        }

        auto setup_res = ElementT::compute_constants(parameters, P_ref, vc_barrier);
        /*auto setup_res = BaseElemT::compute_constants(w * Dirichlet_w, Dirichlet_barrier_w, Dirichlet_barrier_min,
                                                      P_ref,
                                                      vc_barrier);*/

        //log(LogLevel::INFO)<<" - constants vector: "<<vc_barrier.transpose()<<std::endl;

        if (setup_res) {
            log(LogLevel::ERROR) << " ERROR: constants couldn't be computed for face: " << elh << std::endl;
            return -1;
        }

        ElementT elem;
        double barrier_e = elem.eval_f(vx_barrier, vc_barrier);

        auto vol = compute_volume(mesh, elh);
        /*log(LogLevel::INFO)<<" -     volume = "<<vol<<std::endl;
        log(LogLevel::INFO)<<" - rel volume = "<<(vol / target_vol_per_element) <<std::endl;
        log(LogLevel::INFO)<<" -     energy = "<<barrier_e<<std::endl;*/

        if (!std::isfinite(barrier_e)) {


            double VL_qual = mesh.compute_quality_vector(elh, VOLUME_LENGTH, CONFORMING)[0];
            double IVL_qual = 1.0 / VL_qual;


            double bound_w_save = vc_barrier[1];
            double vol_w_save = vc_barrier[3];
            //disable the bounds to compute the base value
            vc_barrier[1] = 0;
            vc_barrier[3] = 0;

            double unw_barrier_e = elem.eval_f(vx_barrier, vc_barrier) / vc_barrier[0];

            if(unw_barrier_e >= vc_barrier[2]) {

                double bound_scale_factor = json_data["global_wrapping_IVL_bound_scaling"];

                vc_barrier[2] = unw_barrier_e * bound_scale_factor;

                backup_bound_use_count++;
                max_bound = std::max(max_bound, vc_barrier[2]);
            }

            vc_barrier[1] = bound_w_save;
            vc_barrier[3] = vol_w_save;

            barrier_e = elem.eval_f(vx_barrier, vc_barrier);

            if (!run_silently) {
                log(LogLevel::DEBUG) << " ----------------------------------------------------------------------------" <<
                                    std::endl;
                log(LogLevel::DEBUG) << " distortion energy is infinite for element " << cell_to_str(mesh, elh) <<std::endl;
                log(LogLevel::DEBUG) << " ->    base energy = " << unw_barrier_e <<std::endl;
                log(LogLevel::DEBUG) << " -> updated bound to " << vc_barrier[2] << std::endl;
                log(LogLevel::DEBUG) << " -> updated energy = " << barrier_e << std::endl;
            }

            if (!std::isfinite(barrier_e)) {
                auto current_vol = compute_volume(mesh, elh);

                double vol_bound_scale_factor = json_data["global_wrapping_IVL_vol_scaling"];

                vc_barrier[4] = current_vol / vol_bound_scale_factor;
                barrier_e = elem.eval_f(vx_barrier, vc_barrier);
                backup_vol_bound_use_count++;


                if (!run_silently) {
                    log(LogLevel::DEBUG) << " ------------------------------------------------------------" << std::endl;
                    log(LogLevel::DEBUG) << " barrier energy is STILL infinite, using a lower volume bound" << std::endl;
                    log(LogLevel::DEBUG) << " -> current volume = " << current_vol << " (relative to average = " << (
                            current_vol / target_vol_per_element) << ")" << std::endl;
                    log(LogLevel::DEBUG) << " -> updated energy = " << barrier_e << std::endl;
                    log(LogLevel::DEBUG) << " -> updated vol bound to = " << vc_barrier[4] << std::endl;
                }

                if (!std::isfinite(barrier_e)) {

                    vc_barrier[1] = 0; //disable the bound to compute the actual value
                    vc_barrier[3] = 0;

                    //unw_barrier_e = elem.eval_f(vx_barrier, vc_barrier) / vc_barrier[0];

                    log(LogLevel::ERROR) << " ------------------------------------------------------------" << std::endl;
                    log(LogLevel::ERROR) << " barrier energy is STILL infinite after using lower volume bound" << std::endl;
                    log(LogLevel::ERROR) << " -> current volume = " << current_vol << " (relative to average = " << (
                            current_vol / target_vol_per_element) << ")" << std::endl;
                    log(LogLevel::ERROR) << " -> unweighted energy = " << unw_barrier_e << std::endl;

                }
            }
        }

        max_energy = std::max(max_energy, barrier_e);

        fe_barrier.instances().add_element(vi_barrier,
                                           vc_barrier);
    }

    if(!run_silently) {
        log(LogLevel::INFO) << " ---> done with setup " << std::endl;
        log(LogLevel::INFO) << " -> max energy = " << max_energy << std::endl;
        log(LogLevel::INFO) << " -> max weight = " << max_w << std::endl;
        log(LogLevel::INFO) << " -> weight sum = " << w_sum << std::endl;
        log(LogLevel::INFO) << " -> #times a backup     bound was used: " << backup_bound_use_count << std::endl;
        log(LogLevel::INFO) << " -> #times a backup vol bound was used: " << backup_vol_bound_use_count << std::endl;
    }

    fe_problem.add_set(&fe_barrier);

    int constrained_count(0);
    //constrain the fixed vertices
    for(auto v: mesh.vertices()) {
        if(constrained_prop[v]) {

            auto pos = mesh.vertex(v);

            if(p_target_pos_prop){
                pos = (*p_target_pos_prop)[v];
            }

            int var_idx = geo_dim * global_vid_prop[v];

            for(int i(0); i<geo_dim; i++) {
                COMISO::LinearConstraint::SVectorNC coeffs;
                coeffs.resize(n_vars);

                coeffs.coeffRef(var_idx + i) = 1;
                eq_constraints.push_back(COMISO::LinearConstraint(coeffs, -pos[i]));
            }
            //log(LogLevel::INFO) << " - constrained vertex " << v << " to " << pos<<" (current position: "<<mesh.vertex(v)<<")" << std::endl;
            /*if(!constrain_boundary){
                log(LogLevel::INFO)<<" --> only constraining first boundary vertex, the rest is free"<<std::endl;
                break;
            }*/
            constrained_count++;
        }
    }
    if (!constrained_count) {
        log(LogLevel::ERROR)<<" ERROR - No vertex constrained at all. You should constrain at least one "<<std::endl;
        return -1;
    }

    std::vector<COMISO::NConstraintInterface*> constraints_pointers;
    constraints_pointers.reserve(constraints_pointers.size());
    for( auto& c : eq_constraints) {
        constraints_pointers.push_back(&c);
    }

    if(!run_silently) {
        log(LogLevel::INFO)<<" energies before optimization: "<<std::endl;
        fe_problem.print_objectives();
    }

    double initial_energy = fe_problem.eval_f(fe_problem.x().data());
    if(initial_energy == std::numeric_limits<double>::infinity()){
        log(LogLevel::INFO)<<" -- infinite initial energy, stopping there."<<std::endl;
        return -1;
    }

    bool use_truncated_newton(true);
    bool converged(false);
    bool feasible_solution_found(false);
    if(use_truncated_newton){
        COMISO::TruncatedNewtonPCG tn(json_data["smoothing_convergence_eps"]);
        tn.always_update_preconditioner() = true;
        tn.max_iters() = max_inner_iters;
        tn.max_pcg_iters() = 100;
        tn.allow_warmstart() = false;
        tn.tolerance_gdx() = json_data["smoothing_tolerance_gdx"];
        tn.set_silent(run_silently);
        tn.set_eps_constraints_violation(json_data["smoothing_constr_eps"]);//1e-9
        //tn.max_preconditioner_range() = json_data["preconditioner_range"];
        //tn.stop_if_progress_is_lower_than_thresh() = stop_if_progress_is_lower_than_thresh && !relaxation_pass;
        //tn.constraint_progress_thresh() = constraint_progress_thresh;

        tn.solve_projected_normal_equation(&fe_problem, constraints_pointers);
        converged = tn.converged();
        feasible_solution_found = tn.feasible_solution_found();
    }else {
        COMISO::NewtonSolver ns(1e-1, 1e-9, max_inner_iters);

        ns.solve_infeasible_start(&fe_problem, constraints_pointers);

        converged = ns.converged();
    }

    //solver_iters += ns.iters();


    if(!run_silently) {
        log(LogLevel::INFO) << " energies after optimization: " << std::endl;
        fe_problem.print_objectives();
    }
    //copy problem data back to mesh
    /*if(!feasible_solution_found){
        log(LogLevel::ERROR)<<" ERROR - couldn't optimize reference mesh"<<std::endl;
        exit(EXIT_FAILURE);
    }*/

    for(auto v: mesh.vertices()){
        if(!constrained_prop[v] || constrained_prop) {
            int global_id = geo_dim * global_vid_prop[v];
            typename MeshT::PointT new_pos;
            for (int i(0); i < MeshT::geo_dim; i++) {
                new_pos[i] = fe_problem.x()[global_id + i];
            }
            mesh.set_vertex(v, new_pos);
        }
    }

    return 0;
}



template<typename MeshT, typename ElementT>
int optimize_mesh(MeshT& mesh,
                  double target_vol_per_element,
                  int max_inner_iters,
                  bool constrain_boundary,
                  const OVM::VertexPropertyT<Vec3d>* p_target_pos_prop,
                  bool enable_boundary_edge_elements,
                  const nlohmann::json& json_data,
                  bool run_silently = true) {

    /*if (constrain_boundary && p_target_pos_prop) {
        log(LogLevel::ERROR) <<
                " ERROR - boundary is constrained, but virtual elements are enabled, which doesn't make sense" <<
                std::endl;
        return -1;
    }*/

    auto constrained_prop = mesh.template request_vertex_property<bool>();
    if (constrain_boundary) {
        for (auto vh: mesh.vertices()) {
            constrained_prop[vh] = mesh.is_boundary(vh);
        }
    } else {
        //if boundary is not constrained, we just fix the first boundary vertex
        for (auto vh: mesh.vertices()) {
            if (!mesh.is_deleted(vh) && mesh.is_boundary(vh)) {
                log(LogLevel::INFO)<<" --> mesh is not constrained globally, only constraining first boundary vertex "<<vh<<std::endl;
                constrained_prop[vh] = true;
                break;
            }
        }
    }


    return optimize_mesh<MeshT, ElementT>(mesh, target_vol_per_element, max_inner_iters, constrained_prop, p_target_pos_prop, enable_boundary_edge_elements, json_data, run_silently);
}


/** Note: only copies the CONFORMING positions for now. You'll have to implement the nonconforming positions copy later */
int setup_mesh_as_one_ring_submesh(const NonConformingTetrahedralMesh& mesh,
                                    const OpenVolumeMesh::EdgeHandle& he,
                                    VERTEX_POSITION_SOURCE pos_source,
                                    NonConformingTetrahedralMesh& submesh,
                                    OpenVolumeMesh::EdgeHandle& submesh_he);

int setup_mesh_as_one_ring_submesh(const NonConformingTetrahedralMesh& mesh,
                                    const OpenVolumeMesh::VertexHandle& vh,
                                    NonConformingTetrahedralMesh& submesh,
                                    OpenVolumeMesh::VertexHandle& submesh_vh);

void setup_priority_queue(NonConformingTetrahedralMesh &mesh,
                          const nlohmann::json &json_data,
                          int valence_filter,
                          EntityPriorityQueue<OVM::EdgeHandle, double>& queue);

void decimate_mesh_with_phases(NonConformingTetrahedralMesh& ref_mesh,
                   NonConformingTetrahedralMesh& mesh,
                   const nlohmann::json& json_data);


void decimate_mesh_with_interleaving(NonConformingTetrahedralMesh &ref_mesh,
                                     NonConformingTetrahedralMesh &mesh,
                                     const nlohmann::json &json_data);


double compute_decimation_priority(NonConformingTetrahedralMesh &ref_mesh,
                                   const OVM::EdgeHandle &eh,
                                   const nlohmann::json &json_data);


/** NOTE: as a hack, this can be called with (mesh, ref_mesh, json_data, true) to remesh the codomain mesh instead of the domain one
 * The first of the two meshes is the only one that is checked for quality improvement. */
RemeshingStats hill_climbing_ref_mesh_remeshing_with_phases(NonConformingTetrahedralMesh &ref_mesh,
                                                            NonConformingTetrahedralMesh &mesh,
                                                            const nlohmann::json &json_data,
                                                            bool also_check_mesh_for_validity = false);

/* NOTE: same comment as above. */
RemeshingStats focused_remeshing(NonConformingTetrahedralMesh &ref_mesh,
                                 NonConformingTetrahedralMesh &mesh,
                                 bool enable_focus_splitting,
                                 const nlohmann::json &json_data,
                                 bool enable_refill);

RemeshingStats hill_climbing_ref_mesh_remeshing_with_interleaving(NonConformingTetrahedralMesh &ref_mesh,
                                                                  NonConformingTetrahedralMesh &mesh,
                                                                  const nlohmann::json &json_data,
                                                                  bool also_check_mesh_for_validity = false);



/** return value:
 * -1: error
 *  0: success
 *  1: worse metric
 *  2: cell flip
 *  3: all neighbors are boundary
 *  4: now deleted
 *  5: boundary
 *  6: good enough quality */
int swap_if_improves_metric(NonConformingTetrahedralMesh &ref_mesh,
                            NonConformingTetrahedralMesh &mesh,
                            const nlohmann::json &json_data,
                            const OVM::EdgeHandle &to_swap,
                            double quality_threshold,
                            double quality_vector_rel_eps,
                            double volume_threshold,
                            double swap_priority,
                            bool also_check_mesh_for_validity,
                            QUALITY_METRIC metric_code,
                            std::vector<EdgeHandle>& impacted_edges);



/** return value:
 * -1: error
 *  0: success
 *  1: worse metric
 *  2: resulting min volume too small
 *  3: good enough quality */
int split_if_improves_metric(NonConformingTetrahedralMesh &ref_mesh,
                             NonConformingTetrahedralMesh &mesh,
                             const nlohmann::json &json_data,
                             const OVM::EdgeHandle &to_split,
                             double quality_threshold,
                             bool enable_boundary_split,
                             bool enable_pre_split_opti,
                             double quality_vector_rel_eps,
                             double volume_threshold,
                             double split_priority, //just used for debugging/logging
                             bool also_check_mesh_for_validity,
                             QUALITY_METRIC metric_code,
                             std::vector<EdgeHandle>& impacted_edges);


/** return value:
 * -1: error
 *  0: success
 *  1: worse metric
 *  2: boundary incident
 *  3: invalid in domain
 *  4: invalid in codomain
 *  5: worse volume
 *  6: good enough quality*/
int collapse_if_improves_metric(NonConformingTetrahedralMesh &ref_mesh,
                                NonConformingTetrahedralMesh &mesh,
                                const nlohmann::json &json_data,
                                const OVM::EdgeHandle &to_collapse,
                                double quality_threshold,
                                double quality_vector_rel_eps,
                                double volume_threshold,
                                double collapse_priority,
                                bool also_check_mesh_for_validity,
                                QUALITY_METRIC metric_code,
                                std::vector<EdgeHandle> &impacted_edges);
