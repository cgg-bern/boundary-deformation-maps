#include "NonConformingHelpers.hh"

#include "LightWeightStopWatch.hh"


CustomArray::CustomArray() {
    primary   = 0.0;
    secondary = 0.0;
    tertiary  = 0.0;
}

CustomArray::CustomArray(const double& x) {
    primary   = x;
    secondary = x;
    tertiary  = x;
}

CustomArray::CustomArray(const std::vector<double>& vec) {
    primary   = vec[0];
    secondary = vec[1];
    tertiary  = vec[2];
}

CustomArray::CustomArray(const double& p, const double& s, const double& t) {
    primary   = p;
    secondary = s;
    tertiary  = t;
}

CustomArray CustomArray::operator*(const CustomArray& other) {
    CustomArray new_t(
        this -> primary * other.primary,
        this -> secondary * other.secondary,
        this -> tertiary * other.tertiary
        );
    return new_t;
}

void CustomArray::operator*=(const CustomArray& other) {
    *this = (*this) * other;
}

CustomArray CustomArray::operator/(const CustomArray& other) {
    CustomArray new_t(
        this -> primary / other.primary,
        this -> secondary / other.secondary,
        this -> tertiary / other.tertiary
        );
    return new_t;
}

void CustomArray::operator/=(const CustomArray& other) {
    *this = (*this) / other;
}

CustomArray CustomArray::operator+(const CustomArray& other) {
    CustomArray new_t(
        this -> primary + other.primary,
        this -> secondary + other.secondary,
        this -> tertiary + other.tertiary
        );
    return new_t;
}

CustomArray CustomArray::operator-(const CustomArray& other) {
    CustomArray new_t(
        this -> primary - other.primary,
        this -> secondary - other.secondary,
        this -> tertiary - other.tertiary
        );
    return new_t;
}


void CustomArray::make_finite(){
    primary = std::isfinite(primary) ? primary : 0;
    secondary = std::isfinite(secondary) ? secondary : 0;
    tertiary = std::isfinite(tertiary) ? tertiary : 0;
}

std::ostream& operator<< (std::ostream& os, const CustomArray& t) {
    os<<"["<<t.primary<<"; "<<t.secondary<<"; "<<t.tertiary<<"]";
    return os;
}

void collapse_relatively_short_boundary_edges(NonConformingTetrahedralMesh &domain_mesh,
                                              NonConformingTetrahedralMesh &codomain_mesh,
                                              double length_ratio_thresh) {
    log(LogLevel::INFO) <<
            " ------------------------------------------ collapsing edges whose domain length (relative to the average) is below "
            << length_ratio_thresh << "... " << std::endl;

    double K_thresh = std::numeric_limits<double>::max();

    domain_mesh.initialize_props(K_thresh);
    codomain_mesh.initialize_props(K_thresh);

    for (auto v: domain_mesh.vertices()) {
        domain_mesh.set_as_original_vertex(v, false);
        codomain_mesh.set_as_original_vertex(v, false);
    }

    auto length_prop = domain_mesh.request_edge_property<double>();

    double total_length(0);
    int boundary_edges_count(0);
    for (auto e: domain_mesh.edges()) {
        if (domain_mesh.is_boundary(e)) {
            auto from_v = domain_mesh.edge(e).from_vertex();
            auto to_v = domain_mesh.edge(e).to_vertex();

            auto from_v_pos = domain_mesh.vertex(from_v);
            auto to_v_pos = domain_mesh.vertex(to_v);

            auto length = (from_v_pos - to_v_pos).norm();
            length_prop[e] = length;

            total_length += length;
            boundary_edges_count++;
        }
    }
    double avg_length = total_length / boundary_edges_count;
    log(LogLevel::INFO) << " -> average length = " << avg_length << std::endl;

    std::vector<OVM::EdgeHandle> edges_to_collapse;
    for (auto e: domain_mesh.edges()) {
        if (domain_mesh.is_boundary(e)) {
            if (length_prop[e] < length_ratio_thresh * avg_length) {
                log(LogLevel::INFO) << " -- edge (" << domain_mesh.edge(e).from_vertex() << ", "
                        << domain_mesh.edge(e).to_vertex() << ") has length " << length_prop[e]
                        << " which is less than " << length_ratio_thresh << " of the average" << std::endl;
                edges_to_collapse.push_back(e);
            }
        }
    }

    log(LogLevel::INFO) << " - found " << edges_to_collapse.size() << " edges to collapse" << std::endl;

    int collapse_count(0);

    for (auto to_collapse: edges_to_collapse) {
        Eigen::Vector3d surviving_v_pos;
        auto heh = domain_mesh.is_collapse_ok(to_collapse, surviving_v_pos);
        if (heh.is_valid()) {
            auto surviving_v = domain_mesh.collapse_edge_and_update_props(heh);
            codomain_mesh.collapse_edge_and_update_props(heh);
            domain_mesh.set_vertex(surviving_v, vec2vec(surviving_v_pos));
            collapse_count++;
        }
    }

    log(LogLevel::INFO) << " -> done, collapsed " << collapse_count << " / " << edges_to_collapse.size() <<
            " boundary edges" << std::endl;

    log(LogLevel::INFO) << " ------------------------------------------ " << std::endl;
}






int setup_boundary_edge_elements(TetrahedralMesh &mesh,
                                 const VertexPropertyT<Vec3d> &target_pos_prop,
                                 const OVM::VertexPropertyT<int> &global_vid_prop,
                                 const double &edge_w,
                                 COMISO::FiniteElementProblem &fe_problem,
                                 COMISO::FiniteElementSet<BOUNDARY_EDGE_ELEM> &fe_edge_elements,
                                 std::vector<COMISO::LinearConstraint> &edge_eq_constraints){


    log(LogLevel::INFO)<<" --------------------------------------------------------------------------------- "<<std::endl;
    log(LogLevel::INFO)<<" setting up boundary edge elements..."<<std::endl;
    log(LogLevel::INFO)<<" -  weight: "<<edge_w<<std::endl;
    auto log_level = LogLevel::DEBUG;

    if(edge_w <= 0.0){
        log(LogLevel::WARNING)<<" edge weight <= 0, skipping setup"<<std::endl;
        return 1;
    }


    double min_unw_e(std::numeric_limits<double>::max()), avg_unw_e(0.0), max_unw_e(0.0);
    double min_e(std::numeric_limits<double>::max()), avg_e(0.0), max_e(0.0);

    EdgeHandle     max_unw_e_eh(-1),  max_e_eh(-1);


    for (auto eh: mesh.edges()) {
        if (mesh.is_boundary(eh)) {

            auto from_v = mesh.edge(eh).from_vertex();
            auto to_v = mesh.edge(eh).to_vertex();

            log(log_level) <<
                           " ----------------------------------------------------------------------------------------------------------------------------"
                           << std::endl;
            log(log_level) << " ------------------------------- handling boundary edge (" << mesh.edge(eh).from_vertex()
                           << ", " << mesh.edge(eh).to_vertex() << ")" << std::endl;


            BOUNDARY_EDGE_ELEM::VecV vx_edge;
            BOUNDARY_EDGE_ELEM::VecI vi_edge;
            BOUNDARY_EDGE_ELEM::VecC vc_edge;

            auto from_v_idx = global_vid_prop[from_v];
            auto   to_v_idx = global_vid_prop[  to_v];

            auto from_pos = mesh.vertex(from_v);
            auto   to_pos = mesh.vertex(  to_v);

            auto target_from_pos = target_pos_prop[from_v];
            auto   target_to_pos = target_pos_prop[  to_v];

            auto target_edge_vec = target_to_pos - target_from_pos;

            vx_edge << from_pos[0], from_pos[1], from_pos[2],
                         to_pos[0],   to_pos[1],   to_pos[2];


            vi_edge << 3 * from_v_idx + 0, 3 * from_v_idx + 1, 3 * from_v_idx + 2,
                       3 *   to_v_idx + 0, 3 *   to_v_idx + 1, 3 *   to_v_idx + 2;


            double target_length = target_edge_vec.norm();

            vc_edge << edge_w / target_length, target_edge_vec[0], target_edge_vec[1], target_edge_vec[2];

            fe_edge_elements.instances().add_element(vi_edge, vc_edge);

            BOUNDARY_EDGE_ELEM elem;
            auto e = elem.eval_f(vx_edge, vc_edge);
            auto unw_e = e / edge_w;

            min_e = std::min(min_e, e);
            min_unw_e = std::min(min_unw_e, unw_e);

            if(e >= max_e){
                max_e_eh = eh;
                max_e = e;
            }
            if(unw_e >= max_unw_e){
                max_unw_e_eh = eh;
                max_unw_e = unw_e;
            }

            avg_e += e;
            avg_unw_e += unw_e;

        }
    }

    fe_problem.add_set(&fe_edge_elements);

    log(LogLevel::INFO)<<"---------------------------------------------"<<std::endl;
    log(LogLevel::INFO)<<" -> boundary edge elements setup done"<<std::endl;
    log(LogLevel::INFO)<<" ---------------------------------------------"<<std::endl;
    log(LogLevel::INFO)<<"       min / max / avg unweighted energy: "<<min_unw_e<<" / "<<max_unw_e<<" / "<<avg_unw_e<<std::endl;
    log(LogLevel::INFO)<<"     edge with max unweighted energy was: "<<edge_to_str(mesh, max_unw_e_eh)<<std::endl;
    log(LogLevel::INFO)<<" ---------------------------------------------"<<std::endl;
    log(LogLevel::INFO)<<"       min / max / avg energy: "<<min_e<<" / "<<max_e<<" / "<<avg_e<<std::endl;
    log(LogLevel::INFO)<<"     edge with max energy was: "<<edge_to_str(mesh, max_e_eh)<<std::endl;
    log(LogLevel::INFO)<<" ---------------------------------------------"<<std::endl;
    log(LogLevel::INFO)<<" --------------------------------------------------------------------------------- "<<std::endl;

    return 0;


}


int setup_mesh_as_one_ring_submesh(const NonConformingTetrahedralMesh &mesh,
                                   const OpenVolumeMesh::EdgeHandle &eh,
                                   VERTEX_POSITION_SOURCE pos_source,
                                   NonConformingTetrahedralMesh &submesh,
                                   OpenVolumeMesh::EdgeHandle &submesh_he) {
    std::vector<OVM::VertexHandle> local_v_prop(mesh.n_vertices(), OVM::VertexHandle(-1));

    auto global_v_prop = submesh.request_vertex_property<VertexHandle>("global_vid", OVM::VertexHandle(-1));
    submesh.set_persistent(global_v_prop);


    //submesh.initialize_props();

    for (auto ec_it = mesh.ec_iter(eh); ec_it.valid(); ec_it++) {
        auto elem_verts = mesh.get_cell_vertices(*ec_it);
        std::vector<OVM::VertexHandle> local_verts;
        for (int i(0); i < elem_verts.size(); i++) {
            auto v = elem_verts[i];
            if (!local_v_prop[v.idx()].is_valid()) {
                NonConformingTetrahedralMesh::EigenVec pos;
                switch (pos_source) {
                    case CONFORMING: {
                        // "double jeopardy" here but that's ok
                        pos = vec2vec(mesh.vertex(v));
                        break;
                    }
                    case NON_CONFORMING: {
                        pos = mesh.per_vertex_pos_prop(*ec_it)[i];
                        break;
                    }
                    case REFERENCE: {
                        pos = mesh.per_vertex_ref_pos_prop(*ec_it)[i];
                        break;
                    }
                }

                auto local_v = submesh.add_vertex(vec2vec(pos));
                local_v_prop[v.idx()] = local_v;
                global_v_prop[local_v] = v;
            }
            local_verts.push_back(local_v_prop[v.idx()]);
        }
        submesh.add_cell(local_verts);
    }

    if (submesh.n_cells() != (mesh.valence(eh) - mesh.is_boundary(eh))) {
        log(LogLevel::ERROR) << " ERROR - 1-ring mesh has " << submesh.n_cells() <<
                " cells but (interior) edge to split is valence " << mesh.valence(eh) << std::endl;
        return -1;
    }

    auto from_v = mesh.edge(eh).from_vertex();
    auto to_v = mesh.edge(eh).to_vertex();

    auto local_from_v = local_v_prop[from_v.idx()];
    auto local_to_v = local_v_prop[to_v.idx()];
    if (!local_from_v.is_valid() || !local_to_v.is_valid()) {
        log(LogLevel::ERROR) << " ERROR - local end vertices are (" << local_from_v << ", " << local_to_v << ")" <<
                std::endl;
        return -1;
    }

    auto local_to_split_heh = submesh.find_halfedge(local_from_v, local_to_v);
    if (!local_to_split_heh.is_valid()) {
        log(LogLevel::ERROR) << " ERROR - couldn't find edge to split in 1-ring mesh" << std::endl;
        return -1;
    }

    submesh_he = submesh.edge_handle(local_to_split_heh);

    return 0;
}


int setup_mesh_as_one_ring_submesh(const NonConformingTetrahedralMesh &mesh,
                                   const OpenVolumeMesh::VertexHandle &vh,
                                   NonConformingTetrahedralMesh &submesh,
                                   OpenVolumeMesh::VertexHandle &submesh_vh) {

    std::vector<OVM::VertexHandle> local_v_prop(mesh.n_vertices(), OVM::VertexHandle(-1));
    auto global_v_prop = submesh.request_vertex_property<VertexHandle>("global_vid", OVM::VertexHandle(-1));
    submesh.set_persistent(global_v_prop);

    for (auto vc_it = mesh.vc_iter(vh); vc_it.valid(); vc_it++) {
        auto elem_verts = mesh.get_cell_vertices(*vc_it);
        std::vector<OVM::VertexHandle> local_verts;
        for (int i(0); i < elem_verts.size(); i++) {
            auto cv = elem_verts[i];
            if (!local_v_prop[cv.idx()].is_valid()) {

                auto pos = mesh.vertex(cv);

                auto local_v = submesh.add_vertex(pos);
                local_v_prop[cv.idx()] = local_v;
                global_v_prop[local_v] = cv;

                if (cv == vh) {
                    submesh_vh = local_v;
                }
            }
            local_verts.push_back(local_v_prop[cv.idx()]);
        }
        submesh.add_cell(local_verts);
    }

    return 0;
}



void setup_priority_queue(NonConformingTetrahedralMesh &mesh,
                          const nlohmann::json &json_data,
                          int valence_filter,
                          EntityPriorityQueue<OVM::EdgeHandle, double>& queue) {

    for (auto eh: mesh.edges()) {

        if (valence_filter > 0 && mesh.valence(eh) != valence_filter) {
            continue;
        }
        auto from_v = mesh.edge(eh).from_vertex();
        auto to_v = mesh.edge(eh).to_vertex();

        if (!mesh.is_boundary(from_v) && !mesh.is_boundary(to_v)) {
            auto from_v_pos = mesh.vertex(from_v);
            auto to_v_pos = mesh.vertex(to_v);

            auto len = (from_v_pos - to_v_pos).norm();


            auto to_skip_prop = mesh.request_cell_property<bool>("", true);
            for (auto ec_it = mesh.ec_iter(eh); ec_it.valid(); ec_it++) {
                to_skip_prop[*ec_it] = false;
            }

            queue.push(eh, compute_decimation_priority(mesh, eh, json_data));
        }
    }
}


void decimate_mesh_with_phases(NonConformingTetrahedralMesh &ref_mesh,
                               NonConformingTetrahedralMesh &mesh,
                               const nlohmann::json &json_data) {

    if (&ref_mesh == &mesh) {
        log(LogLevel::ERROR) << " ERROR - both ref mesh and mesh are the same"<<std::endl;
        return;
    }
    LightWeightStopWatch sw;

    const int max_iterations = json_data["decimation_max_iterations"];
    const bool enable_smoothing = json_data["decimation_enable_smoothing"];
    const bool enable_swaps = json_data["decimation_enable_swaps"];

    double collapse_metric_thresh = json_data["decimation_collapse_metric_thresh"];
    QUALITY_METRIC metric_code = metric_name_to_enum(json_data["decimation_quality_metric"]);

    double rel_volume_thresh = json_data["decimation_rel_volume_thresh"];

    const int init_n_vertices(mesh.n_logical_vertices()),
        init_n_edges(mesh.n_logical_edges()),
        init_n_faces(mesh.n_logical_faces()),
        init_n_cells(mesh.n_logical_cells());

    log(LogLevel::INFO) <<
            " =================================================================================================================="
            << std::endl;
    log(LogLevel::INFO) << " DECIMATING MESH FOR AT MOST " << max_iterations <<
            " ITERATIONS, USING A METRIC MIN THRESHOLD OF "
            << collapse_metric_thresh << " AND WITH BETWEEN-STEP OPTIMIZATION "
            << (enable_smoothing ? "ENABLED" : "DISABLED") << std::endl;
    log(LogLevel::INFO) << " EDGE SWAPS ARE " << (enable_swaps ? "ENABLED" : "DISABLED") << std::endl;

    double ref_mesh_tot_vol = ref_mesh.total_volume();
    double mesh_tot_volume = ref_mesh.total_volume();
    double ref_to_mesh_scale_factor = std::pow(mesh_tot_volume / ref_mesh_tot_vol, 1. / 3.);

    double volume_thresh = rel_volume_thresh * ref_mesh_tot_vol / ref_mesh.n_logical_elements();

    log(LogLevel::INFO) << " -->    quality threshold: " <<collapse_metric_thresh<< std::endl;
    ref_mesh.print_mesh_size();
    mesh.print_mesh_size();

    log(LogLevel::INFO) << " -->    quality threshold: " <<collapse_metric_thresh<< std::endl;
    log(LogLevel::INFO) << " --> rel volume threshold: " <<rel_volume_thresh<< std::endl;
    log(LogLevel::INFO) << " -->     volume threshold: " <<volume_thresh<< std::endl;

    int total_collapse_count(0);
    int total_swap_count(0);

    for (int i(0); i < max_iterations; i++) {
        log(LogLevel::INFO) << " -------------------------------------------------------- RUNNING DECIMATION ITERATION "
                << i << std::endl;

        LightWeightStopWatch it_sw;
        int swap_count(0);

        if (enable_swaps) {
            int no_longer_valence_3_count(0), worse_angles_count(0), cell_flip_count(0),
                    all_boundary_neighbors_count(0);

            auto to_skip_prop = ref_mesh.request_cell_property<bool>("", false);

            /*std::vector<EdgeHandle> to_swap;
            for (auto eh: ref_mesh.edges()) {
                if (!ref_mesh.is_boundary(eh) && ref_mesh.valence(eh) == 3) {
                    to_swap.push_back(eh);
                }
            }*/

            EntityPriorityQueue<OVM::EdgeHandle, double> to_swap;
            setup_priority_queue(ref_mesh, json_data, 3, to_swap);
            log(LogLevel::INFO) << " --> found " << to_swap.size() << " valence 3 edges to swap" << std::endl;

            while (!to_swap.empty()) {

                auto eh = to_swap.front();
                to_swap.pop();

                if (!ref_mesh.is_deleted(eh) && ref_mesh.valence(eh) == 3) {
                    auto from_v = ref_mesh.edge(eh).from_vertex();
                    auto to_v = ref_mesh.edge(eh).to_vertex();
                    log(LogLevel::DEBUG)<<" ---------------- swapping edge "<<edge_to_str(ref_mesh, eh)<<std::endl;

                    //check that the edge can be swapped
                    //log(LogLevel::INFO)<<" - opposite vertices: "<<std::endl;
                    std::vector<VertexHandle> opposite_vertices;
                    HalfEdgeHandle heh = ref_mesh.halfedge_handle(eh, 0);
                    int boundary_count(0);
                    for (auto hehf_it = ref_mesh.hehf_iter(heh); hehf_it.valid(); hehf_it++) {
                        opposite_vertices.push_back(
                            ref_mesh.to_vertex_handle(ref_mesh.next_halfedge_in_halfface(heh, *hehf_it)));
                        //log(LogLevel::INFO)<<" - "<<opposite_vertices.back()<<std::endl;
                        boundary_count += ref_mesh.is_boundary(opposite_vertices.back());
                    }

                    if (boundary_count != 3) {
                        NonConformingTetrahedralMesh temp_mesh;
                        auto v0 = temp_mesh.add_vertex(ref_mesh.vertex(opposite_vertices[0]));
                        auto v1 = temp_mesh.add_vertex(ref_mesh.vertex(opposite_vertices[1]));
                        auto v2 = temp_mesh.add_vertex(ref_mesh.vertex(opposite_vertices[2]));
                        auto vb = temp_mesh.add_vertex(ref_mesh.vertex(from_v));
                        auto vt = temp_mesh.add_vertex(ref_mesh.vertex(to_v));
                        auto cb = temp_mesh.add_cell(v0, v2, v1, vb);
                        auto ct = temp_mesh.add_cell(v0, v1, v2, vt);


                        if (!CGAL_wrapper::is_flipped(temp_mesh, cb) &&
                            !CGAL_wrapper::is_degenerate(temp_mesh, cb) &&
                            !CGAL_wrapper::is_flipped(temp_mesh, ct) &&
                            !CGAL_wrapper::is_degenerate(temp_mesh, ct)) {
                            auto temp_interior_hfh = temp_mesh.face_handle(temp_mesh.halfface({v0, v1, v2}));


                            /*auto angles_before = ref_mesh.compute_angles_in_neighborhood(eh, to_skip_prop, CONFORMING);
                            auto angles_after = temp_mesh.compute_angles_in_neighborhood(temp_interior_hfh, CONFORMING);


                            double min_diff_before = std::min(angles_before.min_dihedral_angle_deg,
                                                              180.0 - angles_before.max_dihedral_angle_deg);

                            double min_diff_after = std::min(angles_after.min_dihedral_angle_deg,
                                                             180.0 - angles_after.max_dihedral_angle_deg);*/
                            auto min_volume_before = ref_mesh.compute_min_volume_in_neighborhood(eh);
                            auto min_volume_after = temp_mesh.compute_min_volume_in_neighborhood(temp_interior_hfh);

                            auto quality_before = ref_mesh.compute_quality_vector_in_neighborhood(eh, metric_code, to_skip_prop, CONFORMING);
                            auto quality_after = temp_mesh.compute_quality_vector_in_neighborhood(temp_interior_hfh, metric_code, to_skip_prop, CONFORMING);


                            log(LogLevel::DEBUG) << " -> quality before: " << quality_before << std::endl;
                            log(LogLevel::DEBUG) << " -> quality after: " << quality_after << std::endl;


                            //log(LogLevel::INFO)<<" -> min angle diff before: "<<min_diff_before<<", min angle diff after: "<<min_diff_after<<std::endl;

                            /*if (min_diff_after > min_diff_before ||
                                min_diff_after > swap_angle_thresh_deg) {*/
                            if((quality_before.compare(quality_after) < 0 ||
                               quality_after[0] > collapse_metric_thresh) &&
                                    (min_volume_after > volume_thresh || min_volume_after > min_volume_before)){

                                log(LogLevel::INFO) << " -> quality vector changed from " << std::endl
                                                    << quality_before << std::endl
                                                    << " to " << std::endl << quality_after << std::endl
                                                    << " and min volume changed from "<<min_volume_before<<" to "<<min_volume_after
                                                    <<" --> performing swap of edge "<<edge_to_str(ref_mesh, eh)<< std::endl;

                                auto mid_v = ref_mesh.split_edge(eh);
                                mesh.split_edge(eh);
                                HalfEdgeHandle to_collapse(-1);
                                for (auto out_he: ref_mesh.outgoing_halfedges(mid_v)) {
                                    auto neighbor = ref_mesh.to_vertex_handle(out_he);
                                    log(LogLevel::DEBUG)<<" -> trying to collapse to neighbor "<<neighbor<<std::endl;
                                    if (neighbor != from_v && neighbor != to_v) {
                                        to_collapse = out_he;
                                        log(LogLevel::DEBUG)<<" --> can collapse on neighbor "<<neighbor<<std::endl;
                                        break;
                                    }
                                }

                                if (to_collapse.is_valid()) {
                                    log(LogLevel::DEBUG)<<" --> trying to collapse out halfedge "<<to_collapse<<": "<<mid_v<<" -> "<<ref_mesh.to_vertex_handle(to_collapse)<<std::endl;

                                    ref_mesh.collapse_edge(to_collapse);
                                    log(LogLevel::DEBUG)<<" --> ref collapse ok "<<std::endl;

                                    mesh.collapse_edge(to_collapse);
                                    log(LogLevel::DEBUG)<<" --> mesh collapse ok "<<std::endl;
                                }else{
                                    log(LogLevel::ERROR) <<
                                            " ERROR - couldn't find edge to collapse after splitting edge (" << from_v
                                            << ", " << to_v << ")" << std::endl;
                                    exit(EXIT_FAILURE);
                                }
                                log(LogLevel::DEBUG)<<" -> swapped edge!"<<std::endl;
                                swap_count++;
                            } else {
                                //log(LogLevel::INFO)<<" -> worse angles, skipping swap"<<std::endl;
                                worse_angles_count++;
                            }
                        } else {
                            //log(LogLevel::INFO) << " -> either cell is flipped, cannot swap edge" << std::endl;
                            cell_flip_count++;
                        }
                    } else {
                        //log(LogLevel::INFO)<<" -> opposite vertices are all boundary, skipping swap"<<std::endl;
                        all_boundary_neighbors_count++;
                    }
                } else {
                    //log(LogLevel::INFO)<<" -> edge is either deleted or no longer valence 3"<<std::endl;
                    no_longer_valence_3_count++;
                }
            }

            log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;
            log(LogLevel::INFO) << " --> edges swapped: " << swap_count << std::endl;
            log(LogLevel::INFO) << " --------------------" << std::endl;
            log(LogLevel::INFO) << " -- deleted or no-longer valence 3: " << no_longer_valence_3_count << std::endl;
            log(LogLevel::INFO) << " -- all-boundary opposite vertices: " << all_boundary_neighbors_count << std::endl;
            log(LogLevel::INFO) << " --    flipped or degenerate cells: " << cell_flip_count << std::endl;
            log(LogLevel::INFO) << " --         worse resulting angles: " << worse_angles_count << std::endl;
            log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;
            total_swap_count += swap_count;
        }


        //1. gather all interior edges, sorted by length
        EntityPriorityQueue<OVM::EdgeHandle, double> edge_queue;
        setup_priority_queue(ref_mesh, json_data, -1, edge_queue);

        log(LogLevel::INFO) << " -> found " << edge_queue.size() << " interior edges to collapse" << std::endl;

        //2. try to collapse as many of them as possible. Collapsible means
        //a. it satisfies the link condition
        //b. it doesn't invert any element
        //c. the dihedral angles after collapse are above the threshold (or better)
        int collapse_count(0);
        while (!edge_queue.empty()) {
            auto to_collapse = edge_queue.front();
            double to_collapse_priority = edge_queue.front_priority();
            edge_queue.pop();

            auto from_v = ref_mesh.edge(to_collapse).from_vertex();
            auto to_v = ref_mesh.edge(to_collapse).to_vertex();

            if (!ref_mesh.is_deleted(from_v) && !ref_mesh.is_deleted(to_v)) {
                /*log(LogLevel::INFO) << " --------- trying to collapse edge (" << mesh.edge(to_collapse).from_vertex() << ", "
                          << mesh.edge(to_collapse).to_vertex() << ")" << std::endl;*/
                bool lc = link_condition(ref_mesh, to_collapse);
                //log(LogLevel::INFO)<<" -> lc: "<<lc<<std::endl;
                if (lc) {
                    //log(LogLevel::INFO) << " - interior edge, checking Chebyshev centroid" << std::endl;
                    double cheb_r(0);
                    auto surviving_v_pos = vec2vec(ref_mesh.solve_edge_Chebyshev(to_collapse, cheb_r));

                    if (cheb_r > 1e-3) {
                        //log(LogLevel::INFO) << " - Chebyshev radius = " << cheb_r << ", edge is geometrically collapsible" << std::endl;


                        auto init_from_v_pos = ref_mesh.vertex(from_v);
                        auto init_to_v_pos = ref_mesh.vertex(to_v);

                        auto edge_incident_prop = mesh.request_cell_property<bool>("", false);
                        auto min_volume_before = ref_mesh.compute_min_volume_in_neighborhood(to_collapse, edge_incident_prop);
                        auto quality_before = ref_mesh.compute_quality_vector_in_neighborhood(to_collapse, metric_code, edge_incident_prop, CONFORMING);


                        ref_mesh.set_vertex(from_v, surviving_v_pos);
                        ref_mesh.set_vertex(to_v, surviving_v_pos);

                        for (auto ec_it = mesh.ec_iter(to_collapse); ec_it.valid(); ec_it++) {
                            edge_incident_prop[*ec_it] = true;
                            //log(LogLevel::INFO) << " - skipping edge-incident element " << ref_mesh.element_to_str(*ec_it) << std::endl;
                        }

                        auto min_volume_after = ref_mesh.compute_min_volume_in_neighborhood(to_collapse, edge_incident_prop);
                        auto quality_after = ref_mesh.compute_quality_vector_in_neighborhood(to_collapse, metric_code, edge_incident_prop, CONFORMING);



                        log(LogLevel::DEBUG) << " -> quality before: " << quality_before << std::endl;
                        log(LogLevel::DEBUG) << " -> quality after: " << quality_after << std::endl;

                    if((quality_after[0] > collapse_metric_thresh || quality_before.compare(quality_after) < 0) &&
                            (min_volume_after >= volume_thresh || min_volume_after >= min_volume_before)){

                        log(LogLevel::INFO) << " -> quality vector changed from " << std::endl
                                            << quality_before << std::endl
                                            << " to " << std::endl << quality_after << std::endl
                                            << " and min volume changed from "<<min_volume_before<<" to "<<min_volume_after
                                            << " --> performing collapse of edge (" << from_v << ", "
                                            << to_v
                                            << ") with priority " << to_collapse_priority
                                            << std::endl;


                        auto to_collapse_he = ref_mesh.halfedge_handle(to_collapse, 0);
                            auto surviving_v = ref_mesh.collapse_edge(to_collapse_he);

                            ref_mesh.set_vertex(surviving_v, surviving_v_pos);

                            mesh.collapse_edge(to_collapse_he);

                            collapse_count++;
                        } else {
                            ref_mesh.set_vertex(from_v, init_from_v_pos);
                            ref_mesh.set_vertex(to_v, init_to_v_pos);
                        }
                    }
                }
            }
        }

        log(LogLevel::INFO) << " --> swapped " << swap_count << " edges and collapsed " << collapse_count <<
                " edges during iteration " << i << " in "<<it_sw.total_duration()<<" seconds, moving on to optimization" << std::endl;
        log(LogLevel::INFO) << " ---------------------------------------- " << std::endl;

        total_collapse_count += collapse_count;

        if (!collapse_count && !swap_count) {
            log(LogLevel::INFO) << " -> nothing was swapped or collapsed during this iteration, stopping here" <<
                    std::endl;
            break;
        }

        //3. optimize the mesh to try to collapse new edges
        int opti_result(-1);
        if (json_data["ref_opti_element_type"] == "fof") {
            opti_result = optimize_mesh<NonConformingTetrahedralMesh, FOF_B_3D_AD>(
                ref_mesh, ref_mesh_tot_vol / ref_mesh.n_logical_elements(),
                100, true, nullptr, false, json_data);
        } else if (json_data["ref_opti_element_type"] == "SDE") {
            opti_result = optimize_mesh<NonConformingTetrahedralMesh, SD_B_3D_AD>(
                ref_mesh, ref_mesh_tot_vol / ref_mesh.n_logical_elements(),
                100, true, nullptr, false,json_data);
        } else if (json_data["ref_opti_element_type"] == "IVL") {
            opti_result = optimize_mesh<NonConformingTetrahedralMesh, IVL_B_3D_AD>(
                ref_mesh, ref_mesh_tot_vol / ref_mesh.n_logical_elements(),
                100, true, nullptr, false, json_data);
        } else {
            log(LogLevel::ERROR) << " ERROR - unsupported ref mesh opti element name " << json_data[
                "ref_opti_element_type"] << std::endl;
            exit(EXIT_FAILURE);
        }

        if (opti_result) {
            log(LogLevel::ERROR) << " ERROR - couldn't optimize mesh at iteration " << i << std::endl;
            exit(EXIT_FAILURE);
        }


        log(LogLevel::INFO) << " -> done with iteration " << i << ", updated mesh size: " << std::endl;
        mesh.print_mesh_size();
        int interior_count(0);
        for (auto v: ref_mesh.vertices()) {
            interior_count += !ref_mesh.is_boundary(v);
        }
        log(LogLevel::INFO) << " --> down to " << interior_count << " interior vertices" << std::endl;
        log(LogLevel::INFO) <<
                " ----------------------------------------------------------------------------------------------------------------"
                << std::endl;
    }

    log(LogLevel::INFO) << " done! Final mesh size after " << total_swap_count << " swaps and " << total_collapse_count
            << " collapses in "<<sw.total_duration()<<" seconds: " << std::endl;
    mesh.print_mesh_size();


    log(LogLevel::INFO) << " - decimation rate (V, E, F, C): ("
              << ((double) mesh.n_logical_vertices() / (double) init_n_vertices) << ", "
              << ((double) mesh.n_logical_edges() / (double) init_n_edges) << ", "
              << ((double) mesh.n_logical_faces() / (double) init_n_faces) << ", "
              << ((double) mesh.n_logical_cells() / (double) init_n_cells) << ") " << std::endl;

    /*log(LogLevel::INFO) << " ---------------------- analyzing elements after decimation " << std::endl;
    analyze_all_elements(ref_mesh, ref_mesh,
                         1.0, false,
                         ref_mesh_tot_vol / ref_mesh.n_logical_elements(), CONFORMING,
                         "elements_after_decimation.csv");
    log(LogLevel::INFO) << " ------------------------------------" << std::endl;*/

    log(LogLevel::INFO) <<
            " =================================================================================================================="
            << std::endl;
}

double compute_decimation_priority(NonConformingTetrahedralMesh &ref_mesh,
                                   const OVM::EdgeHandle &eh,
                                   const nlohmann::json &json_data) {

    auto from_v = ref_mesh.edge(eh).from_vertex();
    auto to_v = ref_mesh.edge(eh).to_vertex();

    auto from_v_pos = ref_mesh.vertex(from_v);
    auto to_v_pos = ref_mesh.vertex(to_v);


    auto to_skip_prop = ref_mesh.request_cell_property<bool>("", true);
    for (auto ec_it = ref_mesh.ec_iter(eh); ec_it.valid(); ec_it++) {
        to_skip_prop[*ec_it] = false;
    }

    std::string priority_criterion = json_data["decimation_priority_criterion"];

    if (priority_criterion == "volume_length") {
        return ref_mesh.compute_quality_vector_in_neighborhood(eh, VOLUME_LENGTH, to_skip_prop, CONFORMING)[0];

    }else {
        auto len = (from_v_pos - to_v_pos).norm();

        auto angles = ref_mesh.compute_angles_in_neighborhood(eh, to_skip_prop, CONFORMING);

        if (priority_criterion == "min_angle_diff") {
            //diff to either zero or 180
            return std::min(angles.min_dihedral_angle_deg, 180.0 - angles.max_dihedral_angle_deg);
        } else if (priority_criterion == "longest_edge") {
            return -len;
        }else if (priority_criterion == "shortest_edge") {
            return len;
        }else if (priority_criterion == "max_dihedral_angle") {
            return -angles.max_dihedral_angle_deg;
        }else if (priority_criterion == "min_dihedral_angle") {
            return -angles.min_dihedral_angle_deg;
        }else{
            log(LogLevel::ERROR) << " ERROR - unhandled decimation priority criterion " << priority_criterion << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}

void decimate_mesh_with_interleaving(NonConformingTetrahedralMesh &ref_mesh,
                                     NonConformingTetrahedralMesh &mesh,
                                     const nlohmann::json &json_data) {

    const int max_iterations = json_data["decimation_max_iterations"];
    const bool enable_smoothing = json_data["decimation_enable_smoothing"];
    const bool enable_swaps = json_data["decimation_enable_swaps"];
    const double collapse_angle_thresh_deg = json_data["decimation_collapse_angle_thresh_deg"];

    std::cerr<<" ERROR - this isn't up-to-date anymore. use the 'phases' version instead"<<std::endl;
    exit(EXIT_FAILURE);

    const int init_n_vertices(mesh.n_logical_vertices()),
        init_n_edges(mesh.n_logical_edges()),
        init_n_faces(mesh.n_logical_faces()),
        init_n_cells(mesh.n_logical_cells());

    log(LogLevel::INFO) <<
            " =================================================================================================================="
            << std::endl;
    log(LogLevel::INFO) << " DECIMATING MESH WITH INTERLEAVING, USING A COLLAPSE DIHEDRAL ANGLE THRESHOLD OF "
            << collapse_angle_thresh_deg << "° AND WITH BETWEEN-STEP OPTIMIZATION "
            << (enable_smoothing ? "ENABLED" : "DISABLED") << std::endl;
    log(LogLevel::INFO) << " EDGE SWAPS ARE " << (enable_swaps ? "ENABLED" : "DISABLED") <<std::endl;

    double ref_mesh_tot_vol = ref_mesh.total_volume();
    double mesh_tot_volume = mesh.total_volume();
    double ref_to_mesh_scale_factor = std::pow(mesh_tot_volume / ref_mesh_tot_vol, 1. / 3.);


    double tot_vol = ref_mesh.total_volume();

    int collapse_count(0);
    int swap_count(0);

    EntityPriorityQueue<OVM::EdgeHandle, double> edge_queue;
    auto in_queue_prop = ref_mesh.request_edge_property<bool>("", false);
    for (auto eh: ref_mesh.edges()) {
        auto from_v = ref_mesh.edge(eh).from_vertex();
        auto to_v = ref_mesh.edge(eh).to_vertex();

        if (!ref_mesh.is_boundary(from_v) && !mesh.is_boundary(to_v)) {
            edge_queue.push(eh, compute_decimation_priority(ref_mesh, eh, json_data));
            in_queue_prop[eh] = true;
        }
    }

    log(LogLevel::INFO) << " --> found " << edge_queue.size() << " interior edges to decimate " << std::endl;
    int not_valence_3_count(0), deleted_count(0), boundary_count(0), worse_angles_count(0), cell_flip_count(0),
        all_boundary_neighbors_count(0);

    int last_op_tot_count(-1);
    int remeshing_period(100);
    while (!edge_queue.empty()) {

        int op_tot_count = swap_count + collapse_count;
        if (!(op_tot_count % remeshing_period) && op_tot_count != last_op_tot_count) {
            last_op_tot_count = op_tot_count;

            log(LogLevel::INFO) << " ------------------------------------------ EDGES LEFT IN QUEUE: "<< edge_queue.size() << std::endl;

            //3. optimize the mesh to try to collapse new edges
            int opti_result(-1);
            if (json_data["ref_opti_element_type"] == "fof") {
                opti_result = optimize_mesh<NonConformingTetrahedralMesh, FOF_B_3D_AD>(
                    ref_mesh, tot_vol / ref_mesh.n_logical_elements(),
                    100, true, nullptr,false, json_data);
            } else if (json_data["ref_opti_element_type"] == "SDE") {
                opti_result = optimize_mesh<NonConformingTetrahedralMesh, SD_B_3D_AD>(
                    ref_mesh, tot_vol / ref_mesh.n_logical_elements(),
                    100, true, nullptr, false, json_data);
            } else if (json_data["ref_opti_element_type"] == "IVL") {
                opti_result = optimize_mesh<NonConformingTetrahedralMesh, IVL_B_3D_AD>(
                    ref_mesh, tot_vol / ref_mesh.n_logical_elements(),
                    100, true, nullptr, false, json_data);
            } else {
                log(LogLevel::ERROR) << " ERROR - unsupported ref mesh opti element name " << json_data[
                    "ref_opti_element_type"] << std::endl;
                exit(EXIT_FAILURE);
            }

            if (opti_result) {
                log(LogLevel::ERROR) << " ERROR - couldn't optimize mesh " << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        auto to_remesh_eh = edge_queue.front();
        edge_queue.pop();
        in_queue_prop[to_remesh_eh] = false;


        std::vector<EdgeHandle> impacted_edges;

        auto from_v = ref_mesh.edge(to_remesh_eh).from_vertex();
        auto to_v = ref_mesh.edge(to_remesh_eh).to_vertex();

        bool swapped_edge(false);
        if (enable_swaps) {
            auto to_skip_prop = ref_mesh.request_cell_property<bool>("", false);

            if (!ref_mesh.is_deleted(to_remesh_eh)) {
                if (ref_mesh.valence(to_remesh_eh) == 3) {
                    if (!ref_mesh.is_boundary(to_remesh_eh)) {
                        //log(LogLevel::INFO)<<" ---------------- swapping edge ("<<from_v<<", "<<to_v<<")"<<std::endl;

                        //check that the edge can be swapped
                        //log(LogLevel::INFO)<<" - opposite vertices: "<<std::endl;
                        std::vector<VertexHandle> opposite_vertices;
                        HalfEdgeHandle heh = ref_mesh.halfedge_handle(to_remesh_eh, 0);
                        int boundary_count(0);
                        for (auto hehf_it = ref_mesh.hehf_iter(heh); hehf_it.valid(); hehf_it++) {
                            opposite_vertices.push_back(
                                ref_mesh.to_vertex_handle(ref_mesh.next_halfedge_in_halfface(heh, *hehf_it)));
                            //log(LogLevel::INFO)<<" - "<<opposite_vertices.back()<<std::endl;
                            boundary_count += ref_mesh.is_boundary(opposite_vertices.back());
                        }

                        if (boundary_count != 3) {
                            NonConformingTetrahedralMesh temp_mesh;
                            auto v0 = temp_mesh.add_vertex(ref_mesh.vertex(opposite_vertices[0]));
                            auto v1 = temp_mesh.add_vertex(ref_mesh.vertex(opposite_vertices[1]));
                            auto v2 = temp_mesh.add_vertex(ref_mesh.vertex(opposite_vertices[2]));
                            auto vb = temp_mesh.add_vertex(ref_mesh.vertex(from_v));
                            auto vt = temp_mesh.add_vertex(ref_mesh.vertex(to_v));
                            auto cb = temp_mesh.add_cell(v0, v2, v1, vb);
                            auto ct = temp_mesh.add_cell(v0, v1, v2, vt);


                            if (!CGAL_wrapper::is_flipped(temp_mesh, cb) &&
                                !CGAL_wrapper::is_degenerate(temp_mesh, cb) &&
                                !CGAL_wrapper::is_flipped(temp_mesh, ct) &&
                                !CGAL_wrapper::is_degenerate(temp_mesh, ct)) {
                                auto temp_interior_hfh = temp_mesh.face_handle(temp_mesh.halfface({v0, v1, v2}));


                                auto angles_before = ref_mesh.compute_angles_in_neighborhood(
                                    to_remesh_eh, to_skip_prop, CONFORMING);
                                auto angles_after = temp_mesh.compute_angles_in_neighborhood(
                                    temp_interior_hfh, CONFORMING);


                                double min_diff_before = std::min(angles_before.min_dihedral_angle_deg,
                                                                  180.0 - angles_before.max_dihedral_angle_deg);

                                double min_diff_after = std::min(angles_after.min_dihedral_angle_deg,
                                                                 180.0 - angles_after.max_dihedral_angle_deg);

                                //log(LogLevel::INFO)<<" -> min angle diff before: "<<min_diff_before<<", min angle diff after: "<<min_diff_after<<std::endl;


                                if (min_diff_after > min_diff_before /*||
                                    min_diff_after > swap_angle_thresh_deg*/) {
                                    auto mid_v = ref_mesh.split_edge(to_remesh_eh);
                                    mesh.split_edge(to_remesh_eh);

                                    /*log(LogLevel::INFO) << " - mid vertex " << mid_v << " neighbors: " << std::endl;
                                    for (auto out_he: ref_mesh.outgoing_halfedges(mid_v)) {
                                        auto neighbor = ref_mesh.to_vertex_handle(out_he);
                                        log(LogLevel::INFO) << " -- " << neighbor << std::endl;
                                    }*/

                                    bool collapsed_new_edge(false);
                                    for (auto out_he: ref_mesh.outgoing_halfedges(mid_v)) {
                                        auto neighbor = ref_mesh.to_vertex_handle(out_he);
                                        //log(LogLevel::INFO)<<" -> trying to collapse to neighbor "<<neighbor<<std::endl;
                                        if (neighbor != from_v && neighbor != to_v) {
                                            ref_mesh.collapse_edge(out_he);
                                            mesh.collapse_edge(out_he);
                                            collapsed_new_edge = true;
                                            //log(LogLevel::INFO)<<" --> collapsed on neighbor "<<neighbor<<std::endl;
                                            break;
                                        }
                                    }
                                    if (!collapsed_new_edge) {
                                        log(LogLevel::ERROR) <<
                                                " ERROR - couldn't find edge to collapse after splitting edge (" <<
                                                from_v
                                                << ", " << to_v << ")" << std::endl;
                                        exit(EXIT_FAILURE);
                                    }
                                    //log(LogLevel::INFO)<<" -> swapped edge!"<<std::endl;
                                    swap_count++;
                                    swapped_edge = true;

                                    std::vector<VertexHandle> end_vertices = {vt, vb};
                                    std::vector<VertexHandle> eq_vertices = {v0, v1, v2};

                                    //pole edges
                                    for (auto end_v: end_vertices) {
                                        for (auto eq_v: eq_vertices) {
                                            auto new_eh = ref_mesh.edge_handle(ref_mesh.find_halfedge(end_v, eq_v));

                                            if (!in_queue_prop[new_eh]  && !ref_mesh.is_boundary(new_eh)) {
                                                impacted_edges.push_back(new_eh);
                                            }
                                        }
                                    }

                                    //equatorial edges
                                    for (int i(0); i<eq_vertices.size() - 1; i++) {
                                        for (int j(i+1); j<eq_vertices.size(); j++) {
                                            auto new_eh = ref_mesh.edge_handle(ref_mesh.find_halfedge(eq_vertices[i], eq_vertices[j]));
                                            if (!in_queue_prop[new_eh] && !ref_mesh.is_boundary(new_eh)) {
                                                impacted_edges.push_back(new_eh);
                                            }
                                        }
                                    }

                                } else {
                                    worse_angles_count++;
                                }
                            } else {
                                cell_flip_count++;
                            }
                        } else {
                            all_boundary_neighbors_count++;
                        }
                    } else {
                        boundary_count++;
                    }
                } else {
                    not_valence_3_count++;
                }
            } else {
                deleted_count++;
            }
        }


        if (!swapped_edge) {

            //2. try to collapse as many of them as possible. Collapsible means
            //a. it satisfies the link condition
            //b. it doesn't invert any element
            //c. the dihedral angles after collapse are above the threshold (or better)

            if (!ref_mesh.is_deleted(from_v) && !ref_mesh.is_deleted(to_v)) {
                /*log(LogLevel::INFO) << " --------- trying to collapse edge (" << mesh.edge(to_collapse).from_vertex() << ", "
                          << mesh.edge(to_collapse).to_vertex() << ")" << std::endl;*/
                bool lc = link_condition(ref_mesh, to_remesh_eh);
                //log(LogLevel::INFO)<<" -> lc: "<<lc<<std::endl;
                if (lc) {
                    //log(LogLevel::INFO) << " - interior edge, checking Chebyshev centroid" << std::endl;
                    double cheb_r(0);
                    auto surviving_v_pos = vec2vec(ref_mesh.solve_edge_Chebyshev(to_remesh_eh, cheb_r));

                    if (cheb_r > 1e-3) {
                        //log(LogLevel::INFO) << " - Chebyshev radius = " << cheb_r << ", edge is geometrically collapsible" << std::endl;


                        auto init_from_v_pos = ref_mesh.vertex(from_v);
                        auto init_to_v_pos = ref_mesh.vertex(to_v);

                        auto edge_incident_prop = mesh.request_cell_property<bool>("", false);

                        auto angles_before = ref_mesh.compute_angles_in_neighborhood(to_remesh_eh,
                            edge_incident_prop,
                            CONFORMING);
                        //log(LogLevel::INFO) << " - min dihedral angle before: " << min_angle_deg_before << std::endl;

                        ref_mesh.set_vertex(from_v, surviving_v_pos);
                        ref_mesh.set_vertex(to_v, surviving_v_pos);

                        for (auto ec_it = mesh.ec_iter(to_remesh_eh); ec_it.valid(); ec_it++) {
                            edge_incident_prop[*ec_it] = true;
                            //log(LogLevel::INFO) << " - skipping edge-incident element " << ref_mesh.element_to_str(*ec_it) << std::endl;
                        }
                        auto angles_after = ref_mesh.compute_angles_in_neighborhood(to_remesh_eh,
                            edge_incident_prop,
                            CONFORMING);
                        /*log(LogLevel::INFO) << " - dihedral angle changed from " << min_angle_deg_before << " to "
                                  << min_angle_deg_after << "°" << std::endl;*/

                        double min_diff_before = std::min(angles_before.min_dihedral_angle_deg,
                                                          180.0 - angles_before.max_dihedral_angle_deg);

                        double min_diff_after = std::min(angles_after.min_dihedral_angle_deg,
                                                         180.0 - angles_after.max_dihedral_angle_deg);


                        /*if (angles_after.min_dihedral_angle_deg > min_dihedral_angle_deg_thresh ||
                            angles_after.min_dihedral_angle_deg > angles_before.min_dihedral_angle_deg) {*/
                        if (min_diff_after > min_diff_before ||
                            min_diff_after > collapse_angle_thresh_deg) {
                            //log(LogLevel::INFO) << " -> performing collapse of edge with priority "<<to_collapse_priority << std::endl;

                            auto to_collapse_he = ref_mesh.halfedge_handle(to_remesh_eh, 0);
                            auto surviving_v = ref_mesh.collapse_edge(to_collapse_he);

                            ref_mesh.set_vertex(surviving_v, surviving_v_pos);

                            mesh.collapse_edge(to_collapse_he);

                            collapse_count++;

                            // add impacted edges to list
                            for (auto out_he_it: ref_mesh.outgoing_halfedges(surviving_v)) {
                                impacted_edges.push_back(ref_mesh.edge_handle(out_he_it));
                            }

                            for (auto vc_it = ref_mesh.vc_iter(surviving_v); vc_it.valid(); vc_it++) {
                                auto op_hf = ref_mesh.vertex_opposite_halfface(*vc_it, surviving_v);

                                for (auto hfe_it = ref_mesh.hfe_iter(op_hf); hfe_it.valid(); hfe_it++) {
                                    if (!in_queue_prop[*hfe_it] && !ref_mesh.is_boundary(*hfe_it)) {
                                        impacted_edges.push_back(*hfe_it);
                                    }
                                }
                            }
                        } else {
                            ref_mesh.set_vertex(from_v, init_from_v_pos);
                            ref_mesh.set_vertex(to_v, init_to_v_pos);
                        }
                    }
                }
            }
        }

        if (!impacted_edges.empty()) {
            int added_count(0);
            for (auto eh: impacted_edges) {
                if (!in_queue_prop[eh] && !ref_mesh.is_deleted(eh)) {
                    auto new_from_v = ref_mesh.edge(eh).from_vertex();
                    auto new_to_v = ref_mesh.edge(eh).to_vertex();

                    if (!ref_mesh.is_boundary(new_from_v) && !ref_mesh.is_boundary(new_to_v)) {
                        in_queue_prop[eh] = true;
                        edge_queue.push(eh, compute_decimation_priority(ref_mesh, eh, json_data));
                        added_count++;
                    }
                }
            }
            log(LogLevel::INFO) << " - added " << added_count << " / " << impacted_edges.size() <<
                    " impacted edges to queue " << std::endl;
        }
    }


    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;
    log(LogLevel::INFO) << " --> edges swapped: " << swap_count << std::endl;
    log(LogLevel::INFO) << " --------------------" << std::endl;
    log(LogLevel::INFO) << " --                          deleted: " << deleted_count << std::endl;
    log(LogLevel::INFO) << " --                         boundary: " << boundary_count << std::endl;
    log(LogLevel::INFO) << " --                    not valence 3: " << not_valence_3_count << std::endl;
    log(LogLevel::INFO) << " -- all-boundary equatorial vertices: " << all_boundary_neighbors_count << std::endl;
    log(LogLevel::INFO) << " --      flipped or degenerate cells: " << cell_flip_count << std::endl;
    log(LogLevel::INFO) << " --           worse resulting angles: " << worse_angles_count << std::endl;
    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;

    log(LogLevel::INFO) << " done! Final mesh size after " << swap_count << " swaps and " << collapse_count
            << " collapses: " << std::endl;
    mesh.print_mesh_size();


    log(LogLevel::INFO) << " - decimation rate (V, E, F, C): ("
              << ((double) mesh.n_logical_vertices() / (double) init_n_vertices) << ", "
              << ((double) mesh.n_logical_edges() / (double) init_n_edges) << ", "
              << ((double) mesh.n_logical_faces() / (double) init_n_faces) << ", "
              << ((double) mesh.n_logical_cells() / (double) init_n_cells) << ") " << std::endl;

    log(LogLevel::INFO) << " ---------------------- analyzing elements after decimation " << std::endl;
    analyze_all_elements(ref_mesh, ref_mesh,
                         1.0, false,
                         ref_mesh_tot_vol / ref_mesh.n_logical_elements(), CONFORMING,
                         "elements_after_decimation.csv");
    log(LogLevel::INFO) << " ------------------------------------" << std::endl;

    log(LogLevel::INFO) <<
            " =================================================================================================================="
            << std::endl;
}



RemeshingStats hill_climbing_ref_mesh_remeshing_with_phases(NonConformingTetrahedralMesh &ref_mesh,
                                                            NonConformingTetrahedralMesh &mesh,
                                                            const nlohmann::json &json_data,
                                                            bool also_check_mesh_for_validity) {

    LightWeightStopWatch sw;
    log(LogLevel::INFO) << " ================================================================== " << std::endl;
    log(LogLevel::INFO) << " - improving ref mesh metric with hill climbing remeshing in successive phases " << std::endl;

    const bool enable_smoothing = json_data["ref_mesh_pre_remeshing_enable_smoothing"];
    const bool enable_pre_split_opti = json_data["ref_mesh_pre_remeshing_enable_pre_split_local_opti"];
    const double target_quality_threshold = json_data["ref_mesh_pre_remeshing_target_quality_threshold"];
    const double candidate_quality_threshold = json_data["ref_mesh_pre_remeshing_candidate_quality_threshold"];
    const double quality_vector_rel_eps = json_data["ref_mesh_pre_remeshing_quality_relative_epsilon"];
    const bool enable_boundary_splits = json_data["ref_mesh_pre_remeshing_enable_boundary_splits"];
    const double rel_volume_threshold = json_data["ref_mesh_pre_remeshing_relative_volume_threshold"];
    const bool enable_early_stop = json_data["ref_mesh_pre_remeshing_enable_early_stop"];

    const double early_stop_progress_thresh = 1.0;


    double ref_mesh_tot_vol = ref_mesh.total_volume();
    double     mesh_tot_vol = mesh.total_volume();

    double volume_threshold = rel_volume_threshold * ref_mesh_tot_vol / ref_mesh.n_logical_elements();
    log(LogLevel::INFO) << " -->       target quality threshold: " <<target_quality_threshold<< std::endl;
    log(LogLevel::INFO) << " -->    candidate quality threshold: " <<candidate_quality_threshold<< std::endl;
    log(LogLevel::INFO) << " -->                    quality eps: " <<quality_vector_rel_eps<< std::endl;
    log(LogLevel::INFO) << " -->           rel volume threshold: " <<rel_volume_threshold<< std::endl;
    log(LogLevel::INFO) << " -->               volume threshold: " <<volume_threshold<< std::endl;

    log(LogLevel::DEBUG) << " - ref mesh before remeshing: " << std::endl;
    if (LogConfig::get_min_level() <= LogLevel::DEBUG) { ref_mesh.print_mesh_size(); }

    int total_swap_count(0), total_collapse_count(0), total_split_count(0);
    int boundary_split_count(0);
    int interior_to_boundary_collapse_count(0);

    const int max_iterations(json_data["ref_mesh_pre_remeshing_max_iterations"]);

    log(LogLevel::INFO) << " --> running at most " << max_iterations << " iterations" << std::endl;

    QUALITY_METRIC metric_code = metric_name_to_enum(json_data["ref_mesh_pre_remeshing_quality_metric"]);

    log(LogLevel::INFO) << " --> using metric: " << json_data["ref_mesh_pre_remeshing_quality_metric"] << ", with code "
            << metric_code << std::endl;

    log(LogLevel::INFO) << " --> also checking mesh for validity: "<< also_check_mesh_for_validity << std::endl;



    enum PHASE { COLLAPSE, SWAP, SPLIT, PHASE_COUNT};
    std::vector<std::string> phase_names = {
        "collapses", "swaps", "splits"
    };

    std::vector<int> phases = json_data["ref_mesh_pre_remeshing_phases"];
    log(LogLevel::INFO) << " --> remeshing will consist of the following phases:";
    for (auto phase: phases) {
        if (phase > PHASE_COUNT) {
            log(LogLevel::ERROR) << " ERROR - unknown phase code " << phase << std::endl;
            return {};
        }
        log(LogLevel::INFO) << " [" << phase_names[phase] << "]";
    }
    log(LogLevel::INFO) << std::endl;

    double last_min_qual = 0.0;
    double last_min_vol = 0.0;

    int last_low_qual_count = 0;
    int last_low_vol_count = 0;


    double IVL_bound_scaling = json_data["global_wrapping_IVL_bound_scaling"];
    double min_quality_before(0.0);
    if (json_data["enable_all_safety_checks"]) {
        min_quality_before = std::numeric_limits<double>::max();
        for (auto ch: ref_mesh.cells()) {
            auto qual = ref_mesh.compute_quality_vector(ch, VOLUME_LENGTH, CONFORMING)[0];
            min_quality_before = std::min(min_quality_before, qual);
        }
        log(LogLevel::INFO) << " --> min quality before remeshing: "<<min_quality_before << std::endl;
    }


    for (int it(0); it < max_iterations; it++) {
        LightWeightStopWatch it_sw;

        log(LogLevel::INFO) << " -------------------------------------------------------- RUNNING REMESHING ITERATION "<< it << std::endl;

        double min_qual(std::numeric_limits<double>::max());
        double min_vol = min_qual;
        int low_qual_tet_count(0), low_vol_tet_count(0);
        bool found_low_qual_tet(false);
        for (auto ch: ref_mesh.cells()) {
            auto qual = ref_mesh.compute_quality_vector(ch, VOLUME_LENGTH, CONFORMING)[0];
            auto vol = compute_volume(ref_mesh, ch);

            if (qual < target_quality_threshold) {
                low_qual_tet_count++;
                //found_low_qual_tet = true;
                //break;
            }
            min_qual = std::min(min_qual, qual);

            if (vol < volume_threshold) {
                low_vol_tet_count++;
            }
            min_vol = std::min(min_vol, vol);
        }
        double min_vol_improvement  = (min_vol  - last_min_vol) / last_min_vol;
        double min_qual_improvement = (min_qual - last_min_qual) / last_min_qual;

        log(LogLevel::INFO) << " --------------------------------------------------------------- " << std::endl;
        log(LogLevel::INFO)<<" --> found "<<std::setw(5)<<low_qual_tet_count<<" elements with quality below target threshold "<<target_quality_threshold<<std::endl;
        log(LogLevel::INFO)<<" -->        last min quality: "<<last_min_qual<<std::endl;
        log(LogLevel::INFO)<<" -->     updated min quality: "<<min_qual<<std::endl;
        log(LogLevel::INFO)<<" --> min quality improvement: "<<min_qual_improvement<<std::endl;
        log(LogLevel::INFO) << " --------------------------------------------- " << std::endl;
        log(LogLevel::INFO)<<" --> found "<<std::setw(5)<<low_vol_tet_count<<" elements with volume below min threshold "<<volume_threshold<<std::endl;
        log(LogLevel::INFO)<<" -->        last min  volume: "<<last_min_vol<<std::endl;
        log(LogLevel::INFO)<<" -->     udpated min  volume: "<<min_vol<<std::endl;
        log(LogLevel::INFO)<<" --> min  volume improvement: "<<min_vol_improvement<<std::endl;

        if (!low_qual_tet_count) {
            log(LogLevel::INFO)<<" --> no element with quality below target threshold "<<target_quality_threshold<<", stopping here "<<std::endl;
            break;
        }

        if (enable_early_stop &&
            low_qual_tet_count >= last_low_qual_count &&
            low_vol_tet_count >= last_low_vol_count &&
            min_qual_improvement < early_stop_progress_thresh &&
            min_vol_improvement <= early_stop_progress_thresh) {
            log(LogLevel::INFO) << " --> min quality and min volume haven't improved enough, and no fewer low-quality or low-volume elements, stopping here " <<std::endl;
            break;
        }
        last_min_qual = min_qual;
        last_low_qual_count = low_qual_tet_count;

        last_min_vol = min_vol;
        last_low_vol_count = low_vol_tet_count;

        log(LogLevel::INFO) << " --------------------------------------------------------------- " << std::endl;


        bool only_split_hinge_equator_edges(false);
        bool only_swap_hinge_equator_edges(false);

        EntityPriorityQueue<OVM::EdgeHandle, double> edge_queue;

        int swap_count(0), collapse_count(0), split_count(0);

        for (auto phase: phases) {

            switch (phase) {
                /*case HINGE_EQUATOR_SWAP:{
                    only_swap_hinge_equator_edges = true;
                }*/
                case SWAP: {
                    sw.lap_duration();

                    log(LogLevel::INFO) <<
                                        " ------------------------------------------------------------------------ entering " << (
                                                only_swap_hinge_equator_edges ? "hinge equator" : "") << " swap phase..." << std::endl;

                    for (auto eh: ref_mesh.edges()) {
                        if (!ref_mesh.is_boundary(eh)) {
                            edge_queue.push(eh, 0.0);
                        }
                    }

                    //gather_fully_glued_edges(false, 1, 3, true, edge_queue);

                    log(LogLevel::INFO) << " --> found " << edge_queue.size() << " interior edges to swap" << std::endl;

                    std::array<int, 7> swap_result_histogram;
                    for (auto &x: swap_result_histogram) { x = 0; }

                    while (!edge_queue.empty()) {
                        auto to_swap = edge_queue.front();
                        double to_swap_priority = edge_queue.front_priority();
                        edge_queue.pop();

                        std::vector<EdgeHandle> impacted_edges;
                        int swap_result = swap_if_improves_metric(ref_mesh, mesh, json_data, to_swap,
                                                                  candidate_quality_threshold, quality_vector_rel_eps,
                                                                  volume_threshold,
                                                                  to_swap_priority, also_check_mesh_for_validity,
                                                                  metric_code, impacted_edges);
                        if (swap_result < 0) {
                            log(LogLevel::ERROR) << " error while trying to swap" << std::endl;
                            return {};
                        }
                        swap_result_histogram[swap_result]++;
                    }

                    swap_count = swap_result_histogram[0];

                    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;
                    log(LogLevel::INFO) << " --> edges swapped: " << swap_count << std::endl;
                    log(LogLevel::INFO) << " ----------------------------------" << std::endl;
                    log(LogLevel::INFO) << " --                       boundary: " << swap_result_histogram[5] <<
                                        std::endl;
                    log(LogLevel::INFO) << " --                        deleted: " << swap_result_histogram[4] <<
                                        std::endl;
                    log(LogLevel::INFO) << " -- all-boundary opposite vertices: " << swap_result_histogram[3] <<
                                        std::endl;
                    log(LogLevel::INFO) << " --    flipped or degenerate cells: " << swap_result_histogram[2] <<
                                        std::endl;
                    log(LogLevel::INFO) << " --         worse resulting metric: " << swap_result_histogram[1] <<
                                        std::endl;
                    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;
                    log(LogLevel::INFO) << " swap phase duration: " << sw.lap_duration() << "s" << std::endl;
                    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;

                    break;
                }

                case COLLAPSE: {
                    bool collapse_debug(true);

                    sw.lap_duration();

                    log(LogLevel::INFO)
                            << " ------------------------------------------------------------------------ entering collapse phase..."
                            << std::endl;

                    for (auto eh: ref_mesh.edges()) {
                        auto from_v = ref_mesh.edge(eh).from_vertex();
                        auto to_v = ref_mesh.edge(eh).to_vertex();
                        if (!ref_mesh.is_boundary(from_v) && !ref_mesh.is_boundary(to_v)) {
                            //if (!ref_mesh.is_boundary(eh)){
                            edge_queue.push(eh, 0.0);
                        }
                    }

                    log(LogLevel::INFO) << " --> found " << edge_queue.size() << " interior edges to collapse " << std::endl;

                    std::array<int, 7> collapse_result_histogram;
                    for (auto &x: collapse_result_histogram) { x = 0; }
                    while (!edge_queue.empty()) {
                        auto to_collapse = edge_queue.front();
                        double to_collapse_priority = edge_queue.front_priority();
                        edge_queue.pop();

                        auto from_v = ref_mesh.edge(to_collapse).from_vertex();
                        auto to_v = ref_mesh.edge(to_collapse).to_vertex();
                        bool interior_to_boundary_collapse = ref_mesh.is_boundary(from_v) != ref_mesh.is_boundary(to_v);

                        std::vector<EdgeHandle> impacted_edges;
                        auto collapse_result = collapse_if_improves_metric(
                                ref_mesh, mesh, json_data, to_collapse, candidate_quality_threshold, quality_vector_rel_eps, volume_threshold, to_collapse_priority,
                                also_check_mesh_for_validity, metric_code, impacted_edges);
                        if (collapse_result < 0) {
                            log(LogLevel::ERROR) << " error while trying to collapse" << std::endl;
                            return {};
                        }
                        interior_to_boundary_collapse_count += !collapse_result;
                        collapse_result_histogram[collapse_result]++;
                    }
                    collapse_count = collapse_result_histogram[0];

                    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;
                    log(LogLevel::INFO) << " --> edges collapsed: " << collapse_count << std::endl;
                    log(LogLevel::INFO) << " ----------------------------------" << std::endl;
                    log(LogLevel::INFO) << " --      boundary-incident: " << collapse_result_histogram[2] << std::endl;
                    log(LogLevel::INFO) << " --    invalid in   domain: " << collapse_result_histogram[3] << std::endl;
                    log(LogLevel::INFO) << " --    invalid in codomain: " << collapse_result_histogram[4] << std::endl;
                    log(LogLevel::INFO) << " -- worse resulting volume: " << collapse_result_histogram[5] << std::endl;
                    log(LogLevel::INFO) << " -- worse resulting metric: " << collapse_result_histogram[1] << std::endl;
                    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;
                    log(LogLevel::INFO) << " collapse phase duration: " << sw.lap_duration() << "s" << std::endl;
                    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;

                    break;
                }

                case SPLIT: {

                    sw.lap_duration();

                    log(LogLevel::INFO) <<
                                        " ------------------------------------------------------------------------ entering " << (
                                                only_split_hinge_equator_edges ? "hinge equator" : "") << " split phase..." <<
                                        std::endl;
                    log(LogLevel::INFO) << " --> with pre-split optimization "
                                        << (enable_pre_split_opti ? "enabled" : "disabled") << std::endl;

                    log(LogLevel::INFO) << " --> using quality vector comparison relative epsilon: " << quality_vector_rel_eps <<
                                        std::endl;

                    for (auto eh: ref_mesh.edges()) {
                        if (enable_boundary_splits || !ref_mesh.is_boundary(eh)) {
                            edge_queue.push(eh, 0.0);
                        }
                    }

                    std::array<int, 4> split_result_histogram;
                    for (auto &x: split_result_histogram) { x = 0; }

                    log(LogLevel::INFO) << " -> found " << edge_queue.size() << " interior edges" <<
                                        std::endl;
                    while (!edge_queue.empty()) {

                        auto to_split = edge_queue.front();
                        double to_split_priority = edge_queue.front_priority();
                        edge_queue.pop();

                        bool boundary_edge = ref_mesh.is_boundary(to_split);
                        std::vector<EdgeHandle> impacted_edges;
                        int split_result = split_if_improves_metric(ref_mesh, mesh, json_data,
                                                                    to_split,
                                                                    candidate_quality_threshold,
                                                                    enable_boundary_splits,
                                                                    enable_pre_split_opti,
                                                                    quality_vector_rel_eps,
                                                                    volume_threshold,
                                                                    to_split_priority,
                                                                    also_check_mesh_for_validity,
                                                                    metric_code,
                                                                    impacted_edges);
                        if (split_result < 0) {
                            log(LogLevel::ERROR) << " error while trying to split" << std::endl;
                            return {};
                        }
                        split_result_histogram[split_result]++;
                        if (!split_result) {
                            split_count++;
                            boundary_split_count += boundary_edge;
                        }
                    }

                    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;
                    log(LogLevel::INFO) << " --> edges split: " << split_count << std::endl;
                    log(LogLevel::INFO) << " ----------------------------------" << std::endl;
                    log(LogLevel::INFO) << " --         worse resulting volume: " << split_result_histogram[2] << std::endl;
                    log(LogLevel::INFO) << " --         worse resulting metric: " << split_result_histogram[1] << std::endl;
                    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;
                    log(LogLevel::INFO) << " split phase duration: " << sw.lap_duration() << "s" << std::endl;
                    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;

                    break;
                }
            }

            if (json_data["enable_all_safety_checks"]) {
                CellHandle min_qual_ch(-1);
                double min_quality_after = std::numeric_limits<double>::max();
                for (auto ch: ref_mesh.cells()) {
                    auto qual = ref_mesh.compute_quality_vector(ch, VOLUME_LENGTH, CONFORMING)[0];
                    if(qual < min_quality_after){
                        min_qual_ch = ch;
                    }
                    min_quality_after = std::min(min_quality_after, qual);
                }
                std::cout<<" --> min quality updated to "<<min_quality_after<<" (was "<<min_quality_before<<" before) after phase "<<phase<<std::endl;
                std::cout<<" --> min quality cell is now "<<cell_to_str(ref_mesh, min_qual_ch)<<std::endl;

                if(min_quality_after < (min_quality_before / IVL_bound_scaling) * 0.9){
                    log(LogLevel::ERROR)<<" ERROR - quality got worse after phase "<<phase<<std::endl;
                    exit(EXIT_FAILURE);
                }
            }

            if (enable_smoothing /*&& phase == SPLIT*/) {

                LightWeightStopWatch smoothing_sw;
                log(LogLevel::INFO) << " ----------------------------------------------------------------------------------------------------------------" << std::endl;
                log(LogLevel::INFO) << " optimizing reference mesh after remeshing..." << std::endl;
                if (optimize_mesh<NonConformingTetrahedralMesh, IVL_B_3D_AD>(ref_mesh, ref_mesh_tot_vol / ref_mesh.n_logical_elements(), 100, true, nullptr, false, json_data,
                                                                             LogConfig::get_min_level() > LogLevel::INFO)) {
                    log(LogLevel::ERROR) << " error while optimizing reference mesh" << std::endl;
                    return {};
                }

                if (also_check_mesh_for_validity) {
                    log(LogLevel::INFO) << " ----------------------------------------------------------------------------------------------------------------" << std::endl;
                    log(LogLevel::INFO) << " optimizing mesh after remeshing..." << std::endl;
                    if (optimize_mesh<NonConformingTetrahedralMesh, IVL_B_3D_AD>(mesh, mesh_tot_vol / mesh.n_logical_elements(), 100, true, nullptr, false, json_data,
                                                                                 LogConfig::get_min_level() > LogLevel::INFO)) {
                        log(LogLevel::WARNING) << " error while optimizing mesh" << std::endl;
                        //return {};
                    }
                    log(LogLevel::INFO) << " ----------------------------------------------------------------------------------------------------------------" << std::endl;
                }
                log(LogLevel::INFO) << " -> smoothed mesh in " << smoothing_sw.total_duration() << " seconds " << std::endl;

                if (json_data["enable_all_safety_checks"]) {
                    double min_quality_after = std::numeric_limits<double>::max();
                    CellHandle min_qual_ch(-1);
                    for (auto ch: ref_mesh.cells()) {
                        auto qual = ref_mesh.compute_quality_vector(ch, VOLUME_LENGTH, CONFORMING)[0];
                        if(qual < min_quality_after){
                            min_qual_ch = ch;
                        }
                        min_quality_after = std::min(min_quality_after, qual);
                    }
                    std::cout<<" --> min quality updated to "<<min_quality_after<<" (was "<<min_quality_before<<" before) after smoothing after phase "<<phase<<std::endl;
                    std::cout<<" --> min quality cell is now "<<cell_to_str(ref_mesh, min_qual_ch)<<std::endl;

                    if(min_quality_after < (min_quality_before / IVL_bound_scaling) * 0.9){
                        log(LogLevel::ERROR)<<" ERROR - quality got worse after smoothing after phase "<<phase<<std::endl;
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }



        log(LogLevel::INFO) <<" ----------------------------------------------------------------------------------------------------------------"<< std::endl;
        log(LogLevel::INFO) << " --> swapped " << swap_count << ", collapsed " << collapse_count << " and split " <<
                   split_count << " edges during iteration " << it << std::endl;

        total_swap_count += swap_count;
        total_collapse_count += collapse_count;
        total_split_count += split_count;

        if (!collapse_count && !swap_count && !split_count) {
            log(LogLevel::INFO) << " -> nothing was remeshed during this iteration, stopping here" << std::endl;
            break;
        }

        log(LogLevel::INFO) << " -> done with iteration " << it << " in "<<it_sw.total_duration()<<"seconds. Updated mesh size: " << std::endl;
        if (LogConfig::get_min_level() <= LogLevel::INFO) { ref_mesh.print_mesh_size(); };

        int interior_count(0);
        for (auto v: ref_mesh.vertices()) {
            interior_count += !ref_mesh.is_boundary(v);
        }
        log(LogLevel::INFO) << " --> down to " << interior_count << " interior vertices" << std::endl;

        log(LogLevel::INFO) <<" ----------------------------------------------------------------------------------------------------------------"<< std::endl;
    }

    log(LogLevel::INFO) <<" ----------------------------------------------------------------------------------------------------------------"<< std::endl;
    log(LogLevel::INFO) << " done with remeshing in "<<sw.total_duration()<<" seconds!"<<std::endl;
    log(LogLevel::INFO)<<" Final mesh size after " << total_swap_count << " swaps, " << total_collapse_count <<
            " collapses and " << total_split_count << " splits: "<<std::endl;
    log(LogLevel::INFO) << " - boundary splits: "<<boundary_split_count<<std::endl;
    log(LogLevel::INFO) << " - interior to boundary collapses: "<<interior_to_boundary_collapse_count<<std::endl;
    if (LogConfig::get_min_level() <= LogLevel::INFO) { ref_mesh.print_mesh_size(); }
    log(LogLevel::INFO) <<" ----------------------------------------------------------------------------------------------------------------"<< std::endl;

    log(LogLevel::INFO) << " ================================================================== " << std::endl;

    return {total_swap_count || total_collapse_count || total_split_count, total_split_count, total_swap_count, total_collapse_count};
}


RemeshingStats focused_remeshing(NonConformingTetrahedralMesh &ref_mesh,
                                 NonConformingTetrahedralMesh &mesh,
                                 bool enable_focus_splitting,
                                 const nlohmann::json &json_data,
                                 bool enable_refill) {

    LightWeightStopWatch sw;
    log(LogLevel::INFO) << " ================================================================== " << std::endl;
    log(LogLevel::INFO) << " - improving ref mesh metric with hill climbing remeshing by focusing on elements with quality below threshold " << std::endl;

    /* NOTE: these parameters can be high-jacked by updating them inside the json. That's a bit hacky but it's just leftovers from a previous approach*/
    const bool enable_smoothing = json_data["ref_mesh_pre_remeshing_enable_smoothing"];
    const bool enable_pre_split_opti = json_data["ref_mesh_pre_remeshing_enable_pre_split_local_opti"];
    const bool enable_boundary_splits = json_data["ref_mesh_pre_remeshing_enable_boundary_splits"];
    //const double candidate_quality_threshold = json_data["ref_mesh_pre_remeshing_candidate_quality_threshold"];
    const double target_quality_threshold = json_data["ref_mesh_pre_remeshing_target_quality_threshold"];
    const double quality_vector_rel_eps = json_data["ref_mesh_pre_remeshing_quality_relative_epsilon"];
    const double rel_volume_threshold = json_data["ref_mesh_pre_remeshing_relative_volume_threshold"];

    const double ref_mesh_tot_vol = ref_mesh.total_volume();
    const double     mesh_tot_vol = ref_mesh.total_volume();

    const double volume_threshold = rel_volume_threshold * ref_mesh_tot_vol / ref_mesh.n_logical_elements();

    bool also_check_mesh_for_validity = true;

    //log(LogLevel::INFO) << " --> candidate quality threshold: " <<candidate_quality_threshold<< std::endl;
    log(LogLevel::INFO) << " -->    target quality threshold: " <<target_quality_threshold<< std::endl;
    log(LogLevel::INFO) << " -->                 quality eps: " <<quality_vector_rel_eps<< std::endl;
    log(LogLevel::INFO) << " -->        rel volume threshold: " <<rel_volume_threshold<< std::endl;
    log(LogLevel::INFO) << " -->            volume threshold: " <<volume_threshold<< std::endl;

    log(LogLevel::INFO) << " - ref mesh before remeshing: " << std::endl;
    if (LogConfig::get_min_level() <= LogLevel::INFO) { ref_mesh.print_mesh_size(); }

    int total_swap_count(0), total_collapse_count(0), total_split_count(0);
    int boundary_split_count(0);
    int interior_to_boundary_collapse_count(0);

    QUALITY_METRIC metric_code = metric_name_to_enum(json_data["ref_mesh_pre_remeshing_quality_metric"]);

    log(LogLevel::INFO) << " --> using metric: " << json_data["ref_mesh_pre_remeshing_quality_metric"] << ", with code "
            << metric_code << std::endl;


    enum PHASE { COLLAPSE, SWAP, SPLIT, PHASE_COUNT};
    std::vector<std::string> phase_names = {
        "collapses", "swaps", "splits"
    };

    std::vector<int> phases = json_data["ref_mesh_pre_remeshing_phases"];
    log(LogLevel::INFO) << " --> remeshing will consist of the following phases:";
    for (auto phase: phases) {
        if (phase > PHASE_COUNT) {
            log(LogLevel::ERROR) << " ERROR - unknown phase code " << phase << std::endl;
            return {};
        }
        log(LogLevel::INFO) << " [" << phase_names[phase] << "]";
    }
    log(LogLevel::INFO) << std::endl;

    int below_target_count(0);
    //gather edges on bad quality elements
    EntityPriorityQueue<CellHandle, double> cells_to_remesh;
    for (auto ch: mesh.cells()) {
        auto qual = ref_mesh.compute_quality_vector(ch, VOLUME_LENGTH, CONFORMING)[0];
        if (qual <= target_quality_threshold) {
            cells_to_remesh.push(ch, qual);
        }
    }
    int initial_cells_to_remesh_count = cells_to_remesh.size();
    log(LogLevel::INFO) << " --> initially found "<<initial_cells_to_remesh_count<<" cells to remesh"<<std::endl;

    std::array<int, 7> swap_result_histogram;
    for (auto &x: swap_result_histogram) { x = 0; }

    std::array<int, 6> collapse_result_histogram;
    for (auto &x: collapse_result_histogram) { x = 0; }

    std::array<int, 4> split_result_histogram;
    for (auto &x: split_result_histogram) { x = 0; }

    int focused_splitting_added_cells(0);
    int added_later_cells_count(0);
    int attempted_cells_count(0);
    int remeshed_cell_count(0);

    while (!cells_to_remesh.empty()) {
        auto ch = cells_to_remesh.front();
        auto priority = cells_to_remesh.front_priority();
        cells_to_remesh.pop();

        if (ref_mesh.is_deleted(ch)) {
            continue;
        }

        log(LogLevel::DEBUG) << " ------------ trying to remesh cell " << cell_to_str(ref_mesh, ch) << " with priority "
                << priority << std::endl;

        attempted_cells_count++;

        int new_cells_start_idx = ref_mesh.n_cells();

        bool remeshed_edge(false);
        for (auto ce_it = ref_mesh.ce_iter(ch); ce_it.valid(); ce_it++) {
            auto eh = *ce_it;
            log(LogLevel::DEBUG) << " -------- trying to remesh edge " << edge_to_str(ref_mesh, *ce_it) << std::endl;

            if (!mesh.is_boundary(eh)) {
                auto from_v = mesh.edge(eh).from_vertex();
                auto   to_v = mesh.edge(eh).to_vertex();

                for (auto phase: phases) {
                    if (phase == SWAP) {
                        std::vector<EdgeHandle> impacted_edges;
                        int swap_result = swap_if_improves_metric(ref_mesh, mesh, json_data, eh,
                                                                  target_quality_threshold, quality_vector_rel_eps,
                                                                  volume_threshold,
                                                                  priority, also_check_mesh_for_validity,
                                                                  metric_code, impacted_edges);
                        if (swap_result < 0) {
                            log(LogLevel::ERROR) << " error while trying to swap" << std::endl;
                            return {};
                        }
                        swap_result_histogram[swap_result]++;

                        if (!swap_result) {
                            //log(LogLevel::INFO) << " -> edge was swapped, moving on to next one"<<std::endl;
                            total_swap_count++;
                            break;
                        }
                    } else if (phase == COLLAPSE) {
                        auto from_v = ref_mesh.edge(eh).from_vertex();
                        auto to_v = ref_mesh.edge(eh).to_vertex();
                        bool interior_to_boundary_collapse = ref_mesh.is_boundary(from_v) != ref_mesh.is_boundary(to_v);

                        std::vector<EdgeHandle> impacted_edges;
                        auto collapse_result = collapse_if_improves_metric(
                            ref_mesh, mesh, json_data, eh, target_quality_threshold, quality_vector_rel_eps,
                            volume_threshold, priority, also_check_mesh_for_validity, metric_code, impacted_edges);
                        if (collapse_result < 0) {
                            log(LogLevel::ERROR) << " error while trying to collapse" << std::endl;
                            return {};
                        }
                        interior_to_boundary_collapse_count += !collapse_result;
                        collapse_result_histogram[collapse_result]++;


                        if (!collapse_result) {
                            //log(LogLevel::INFO) << " -> edge was collapsed, moving on to next one"<<std::endl;
                            total_collapse_count++;
                            break;
                        }
                    } else if (phase == SPLIT) {
                        bool boundary_edge = ref_mesh.is_boundary(eh);
                        std::vector<EdgeHandle> impacted_edges;
                        int split_result = split_if_improves_metric(ref_mesh, mesh, json_data,
                                                                    eh,
                                                                    target_quality_threshold,
                                                                    enable_boundary_splits,
                                                                    enable_pre_split_opti,
                                                                    quality_vector_rel_eps,
                                                                    volume_threshold,
                                                                    priority,
                                                                    also_check_mesh_for_validity,
                                                                    metric_code,
                                                                    impacted_edges);
                        if (split_result < 0) {
                            log(LogLevel::ERROR) << " error while trying to split" << std::endl;
                            return {};
                        }

                        split_result_histogram[split_result]++;

                        if (!split_result) {
                            total_split_count++;
                            boundary_split_count += boundary_edge;
                        }
                    }
                }

                if (json_data["enable_all_safety_checks"]) {
                    if (ref_mesh.contains_degenerate_or_flipped_elements(CONFORMING) ||
                        (also_check_mesh_for_validity && mesh.contains_degenerate_or_flipped_elements(CONFORMING))) {
                        log(LogLevel::ERROR) <<
                                             " ERROR - either ref mesh or mesh contains degenerate or flipped elements after dealing with edge ("<<from_v<<", "<<to_v<<")" <<
                                             std::endl;

                        for (auto elh: ref_mesh.elements()) {
                            if (ref_mesh.is_flipped_or_degenerate(elh, CONFORMING)) {
                                log(LogLevel::ERROR) << " --> ref element " << ref_mesh.element_to_str(elh) << " with volume " <<
                                                     ref_mesh.element_volume(elh) << std::endl;
                            }
                            if (mesh.is_flipped_or_degenerate(elh, CONFORMING)) {
                                log(LogLevel::ERROR) << " --> element " << mesh.element_to_str(elh) << " with volume " <<
                                                     mesh.element_volume(elh) << std::endl;
                            }
                        }
                        return {};
                    }
                }

                if (remeshed_edge) {
                    log(LogLevel::DEBUG) << " -> remeshed edge " << eh << " on low quality cell " << ch << std::endl;
                    remeshed_cell_count++;
                    break;
                }
            }
        }

        if (enable_refill) {
            below_target_count = 0;
            for (int idx(new_cells_start_idx); idx < ref_mesh.n_cells(); idx++) {
                CellHandle new_ch(idx);
                if (!ref_mesh.is_deleted(new_ch)) {
                    auto qual = ref_mesh.compute_quality_vector(new_ch, VOLUME_LENGTH, CONFORMING)[0];
                    if (qual < target_quality_threshold) {
                        cells_to_remesh.push(new_ch, qual);
                        added_later_cells_count++;
                        log(LogLevel::DEBUG) << " - added new cell " << cell_to_str(ref_mesh, new_ch) <<
                                "  with quality " << qual << " to queue" << std::endl;
                        below_target_count += qual < target_quality_threshold;
                    }
                }
            }
        }
    }



    if (enable_smoothing) {
        log(LogLevel::INFO) <<
                " ----------------------------------------------------------------------------------------------------------------"
                << std::endl;
        log(LogLevel::INFO) << " optimizing reference mesh after remeshing..." << std::endl;
        if (optimize_mesh<NonConformingTetrahedralMesh, IVL_B_3D_AD>(
            ref_mesh, ref_mesh_tot_vol, 10, true, nullptr, false, json_data,
            LogConfig::get_min_level() > LogLevel::INFO)) {
            log(LogLevel::ERROR) << " error while optimizing reference mesh" << std::endl;
            return {};
        }

        log(LogLevel::INFO) <<
                " ----------------------------------------------------------------------------------------------------------------"
                << std::endl;
        log(LogLevel::INFO) << " optimizing mesh after remeshing..." << std::endl;
        if (optimize_mesh<NonConformingTetrahedralMesh, IVL_B_3D_AD>(
            mesh, mesh_tot_vol, 10, true, nullptr, false, json_data,
            LogConfig::get_min_level() > LogLevel::INFO)) {
            log(LogLevel::ERROR) << " error while optimizing mesh" << std::endl;
            return {};
        }
        log(LogLevel::INFO) <<
                " ----------------------------------------------------------------------------------------------------------------"
                << std::endl;
    }

    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;
    log(LogLevel::INFO) << " --> edges swapped: " << total_swap_count << std::endl;
    log(LogLevel::INFO) << " ----------------------------------" << std::endl;
    log(LogLevel::INFO) << " --                       boundary: " << swap_result_histogram[5] << std::endl;
    log(LogLevel::INFO) << " --                        deleted: " << swap_result_histogram[4] << std::endl;
    log(LogLevel::INFO) << " -- all-boundary opposite vertices: " << swap_result_histogram[3] << std::endl;
    log(LogLevel::INFO) << " --    flipped or degenerate cells: " << swap_result_histogram[2] << std::endl;
    log(LogLevel::INFO) << " --         worse resulting metric: " << swap_result_histogram[1] << std::endl;
    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;

    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;
    log(LogLevel::INFO) << " --> edges collapsed: " << total_collapse_count << std::endl;
    log(LogLevel::INFO) << " ----------------------------------" << std::endl;
    log(LogLevel::INFO) << " --      boundary-incident: " << collapse_result_histogram[2] << std::endl;
    log(LogLevel::INFO) << " --    invalid in   domain: " << collapse_result_histogram[3] << std::endl;
    log(LogLevel::INFO) << " --    invalid in codomain: " << collapse_result_histogram[4] << std::endl;
    log(LogLevel::INFO) << " -- worse resulting volume: " << collapse_result_histogram[5] << std::endl;
    log(LogLevel::INFO) << " -- worse resulting metric: " << collapse_result_histogram[1] << std::endl;
    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;

    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;
    log(LogLevel::INFO) << " --> edges split: " << total_split_count << std::endl;
    log(LogLevel::INFO) << " ----------------------------------" << std::endl;
    log(LogLevel::INFO) << " -- worse resulting volume: " << split_result_histogram[2] << std::endl;
    log(LogLevel::INFO) << " -- worse resulting metric: " << split_result_histogram[1] << std::endl;
    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;

    log(LogLevel::INFO) <<" ----------------------------------------------------------------------------------------------------------------"<< std::endl;
    log(LogLevel::INFO) << " done with remeshing in " << sw.total_duration() << " seconds!" << std::endl;
    log(LogLevel::INFO) << " Final mesh size after " << total_swap_count << " swaps, " << total_collapse_count <<
            " collapses and " << total_split_count << " splits: " << std::endl;
    if (LogConfig::get_min_level() <= LogLevel::INFO) { ref_mesh.print_mesh_size(); }
    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;
    log(LogLevel::INFO) << " - #cells remeshed: " << remeshed_cell_count << std::endl;
    log(LogLevel::INFO) << " - #cells attempted: " << attempted_cells_count << std::endl;
    log(LogLevel::INFO) << " - #cells added to queue later: " << added_later_cells_count << std::endl;
    log(LogLevel::INFO) << " - #cells added after focused splitting: " << focused_splitting_added_cells << std::endl;
    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;
    log(LogLevel::INFO) << " - boundary splits: " << boundary_split_count << std::endl;
    log(LogLevel::INFO) << " - interior to boundary collapses: " << interior_to_boundary_collapse_count << std::endl;
    log(LogLevel::INFO) <<" ----------------------------------------------------------------------------------------------------------------"<< std::endl;

    log(LogLevel::INFO) << " ================================================================== " << std::endl;

    /*if (focused_splitting_added_cells) {
        manual_stop("stopping after focused splitting for debugging");
    }*/


    return {total_swap_count || total_collapse_count || total_split_count, total_split_count, total_swap_count, total_collapse_count};
}


RemeshingStats hill_climbing_ref_mesh_remeshing_with_interleaving(NonConformingTetrahedralMesh &ref_mesh,
                                                       NonConformingTetrahedralMesh &mesh,
                                                       const nlohmann::json &json_data,
                                                       bool also_check_mesh_for_validity) {

    log(LogLevel::INFO) << " ================================================================== " << std::endl;
    log(LogLevel::INFO) << " - improving ref mesh metric with hill climbing with operation interleaving " << std::endl;

    const bool enable_smoothing = json_data["ref_mesh_pre_remeshing_enable_smoothing"];
    const bool enable_pre_split_opti = json_data["ref_mesh_pre_remeshing_enable_pre_split_local_opti"];

    const double quality_threshold = json_data["ref_mesh_pre_remeshing_candidate_quality_threshold"];
    const double quality_vector_rel_eps = json_data["ref_mesh_pre_remeshing_quality_relative_epsilon"];
    const double relative_volume_threshold = json_data["ref_mesh_pre_remeshing_relative_volume_threshold"];

    const bool enable_boundary_splits = json_data["ref_mesh_pre_remeshing_enable_boundary_splits"];

    double ref_mesh_tot_vol = ref_mesh.total_volume();

    double volume_threshold = relative_volume_threshold * ref_mesh_tot_vol / ref_mesh.n_logical_elements();

    log(LogLevel::DEBUG) << " - ref mesh before remeshing: " << std::endl;
    if (LogConfig::get_min_level() <= LogLevel::DEBUG) { ref_mesh.print_mesh_size(); }

    int total_swap_count(0), total_collapse_count(0), total_split_count(0);

    QUALITY_METRIC metric_code = metric_name_to_enum(json_data["ref_mesh_pre_remeshing_quality_metric"]);

    log(LogLevel::INFO) << " --> using metric: " << json_data["ref_mesh_pre_remeshing_quality_metric"] << ", with code "
            << metric_code << std::endl;


    enum PHASE { COLLAPSE, SWAP, SPLIT, HINGE_EQUATOR_SWAP, HINGE_EQUATOR_SPLIT, PHASE_COUNT };
    std::vector<std::string> phase_names = {
        "collapses", "swaps", "splits", "hinge_equator_swap", "hinge_equator_split"
    };

    std::vector<int> phases = json_data["ref_mesh_pre_remeshing_phases"];
    log(LogLevel::INFO) << " --> remeshing will consist of the following phases:";
    for (auto phase: phases) {
        if (phase > PHASE_COUNT) {
            log(LogLevel::ERROR) << " ERROR - unknown phase code " << phase << std::endl;
            return {};
        }
        log(LogLevel::INFO) << " [" << phase_names[phase] << "]";
    }
    log(LogLevel::INFO) << std::endl;

    // results histograms setup
    std::array<int, 7> swap_result_histogram;
    for (auto &x: swap_result_histogram) { x = 0; }

    std::array<int, 6> collapse_result_histogram;
    for (auto &x: collapse_result_histogram) { x = 0; }


    auto in_queue_prop = ref_mesh.request_edge_property<bool>("", false);
    EntityPriorityQueue<OVM::EdgeHandle, double> edge_queue;
    for (auto eh: ref_mesh.edges()) {
        if (enable_boundary_splits || !ref_mesh.is_boundary(eh)) {
            edge_queue.push(eh, 0.0);
            in_queue_prop[eh] = true;
        }
    }
    log(LogLevel::INFO) << " --> found " << edge_queue.size() << " interior edges to remesh initially" << std::endl;

    int remeshing_period(100);
    int last_op_tot_count(-1);

    while (!edge_queue.empty()) {
        int op_tot_count = total_swap_count + total_split_count + total_collapse_count;

        if (!(op_tot_count % remeshing_period) && op_tot_count != last_op_tot_count) {
            log(LogLevel::INFO) << " -------------------------------------------------------- EDGES LEFT IN QUEUE: " <<
                    edge_queue.size() << std::endl;

            last_op_tot_count = op_tot_count;

            if (optimize_mesh<NonConformingTetrahedralMesh, IVL_B_3D_AD>(ref_mesh, ref_mesh_tot_vol, 100, true, nullptr, false,json_data,
                                                                      LogConfig::get_min_level() > LogLevel::INFO)) {
                log(LogLevel::ERROR) << " error while optimizing reference mesh" << std::endl;
                return {};
            }
        }

        auto to_remesh_priority = edge_queue.front_priority();
        auto to_remesh = edge_queue.front();
        edge_queue.pop();
        in_queue_prop[to_remesh] = false;

        if (ref_mesh.is_deleted(to_remesh) || (!enable_boundary_splits && ref_mesh.is_boundary(to_remesh))) {
            continue;
        }

        std::vector<EdgeHandle> impacted_edges;

        for (auto phase: phases) {
            switch (phase) {
                case SWAP: {
                    int swap_result = swap_if_improves_metric(ref_mesh, mesh, json_data, to_remesh, quality_threshold, quality_vector_rel_eps, volume_threshold,
                                                              to_remesh_priority, also_check_mesh_for_validity,
                                                              metric_code,
                                                              impacted_edges);
                    if (swap_result < 0) {
                        log(LogLevel::ERROR) << " error while trying to swap" << std::endl;
                        return {};
                    }
                    if (swap_result >= swap_result_histogram.size()) {
                        log(LogLevel::ERROR) << " ERROR - unhandled swap result "<<swap_result << std::endl;
                        return {};
                    }
                    swap_result_histogram[swap_result]++;
                    total_swap_count += !swap_result;

                    break;
                }

                case COLLAPSE: {
                    if (!ref_mesh.is_deleted(to_remesh)) {
                        auto collapse_result = collapse_if_improves_metric(
                            ref_mesh, mesh, json_data, to_remesh, quality_threshold, quality_vector_rel_eps, volume_threshold,
                            to_remesh_priority, also_check_mesh_for_validity, metric_code, impacted_edges);

                        if (collapse_result < 0) {
                            log(LogLevel::ERROR) << " error while trying to collapse" << std::endl;
                            return {};
                        }
                        if (collapse_result >= collapse_result_histogram.size()) {
                            log(LogLevel::ERROR) << " ERROR - unhandled collapse result "<<collapse_result << std::endl;
                            return {};
                        }
                        collapse_result_histogram[collapse_result]++;
                        total_collapse_count += !collapse_result;
                    }
                    break;
                }

                case SPLIT: {
                    if (!ref_mesh.is_deleted(to_remesh)) {
                        int split_result = split_if_improves_metric(ref_mesh, mesh, json_data,
                                                                    to_remesh,
                                                                    quality_threshold,
                                                                    enable_boundary_splits,
                                                                    enable_pre_split_opti,
                                                                    quality_vector_rel_eps,
                                                                    volume_threshold,
                                                                    to_remesh_priority,
                                                                    also_check_mesh_for_validity,
                                                                    metric_code,
                                                                    impacted_edges);
                        if (split_result < 0) {
                            log(LogLevel::ERROR) << " error while trying to split" << std::endl;
                            return {};
                        } else if (!split_result) {
                            total_split_count++;
                        }
                    }
                    break;
                }
            }
        }


        if (!impacted_edges.empty()) {
            int added_count(0);
            for (auto eh: impacted_edges) {
                if (!in_queue_prop[eh] && !ref_mesh.is_boundary(eh) && !ref_mesh.is_boundary(eh)) {
                    edge_queue.push(eh, 0.0);
                    in_queue_prop[eh] = true;
                    added_count++;
                }
            }
            log(LogLevel::INFO) << " - added "<<added_count<<" / "<<impacted_edges.size()<<" impacted edges to queue "<< std::endl;
        }
    }


        /*if (enable_smoothing) {
            if(optimize_mesh<NonConformingTetrahedralMesh, IVL3D_AD>(ref_mesh, ref_mesh_tot_vol , 500, true, json_data, LogConfig::get_min_level()> LogLevel::INFO)){
                log(LogLevel::ERROR)<<" error while optimizing reference mesh"<<std::endl;
                return -1;
            }
        }*/


    log(LogLevel::INFO) <<" ----------------------------------------------------------------------------------------------------------------"<< std::endl;
    log(LogLevel::INFO) << " done! Final mesh size after " << total_swap_count << " swaps, " << total_collapse_count <<
            " collapses and " << total_split_count << " splits: " << std::endl;
    if (LogConfig::get_min_level() <= LogLevel::INFO) { ref_mesh.print_mesh_size(); }
    log(LogLevel::INFO) <<" ----------------------------------------------------------------------------------------------------------------"<< std::endl;


    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;
    log(LogLevel::INFO) << " --> edges swapped: " << total_swap_count << std::endl;
    log(LogLevel::INFO) << " ----------------------------------" << std::endl;
    log(LogLevel::INFO) << " --                       boundary: " << swap_result_histogram[5] <<
            std::endl;
    log(LogLevel::INFO) << " --                        deleted: " << swap_result_histogram[4] <<
            std::endl;
    log(LogLevel::INFO) << " -- all-boundary opposite vertices: " << swap_result_histogram[3] <<
            std::endl;
    log(LogLevel::INFO) << " --    flipped or degenerate cells: " << swap_result_histogram[2] <<
            std::endl;
    log(LogLevel::INFO) << " --         worse resulting metric: " << swap_result_histogram[1] <<
            std::endl;
    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;


    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;
    log(LogLevel::INFO) << " --> edges collapsed: " << total_collapse_count << std::endl;
    log(LogLevel::INFO) << " ----------------------------------" << std::endl;
    log(LogLevel::INFO) << " -- worse resulting metric: " << collapse_result_histogram[1] << std::endl;
    log(LogLevel::INFO) << " --      boundary-incident: " << collapse_result_histogram[2] << std::endl;
    log(LogLevel::INFO) << " --    invalid in   domain: " << collapse_result_histogram[3] << std::endl;
    log(LogLevel::INFO) << " --    invalid in codomain: " << collapse_result_histogram[4] << std::endl;
    log(LogLevel::INFO) << " ---------------------------------------------" << std::endl;


    log(LogLevel::INFO) << " -     ref mesh after remeshing: " << std::endl;
    if (LogConfig::get_min_level() <= LogLevel::INFO) { ref_mesh.print_mesh_size(); }

    log(LogLevel::INFO) << " ================================================================== " << std::endl;

    return {total_swap_count || total_collapse_count || total_split_count, total_split_count, total_swap_count, total_collapse_count};
}


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
                                std::vector<EdgeHandle>& impacted_edges) {


    auto from_v = ref_mesh.edge(to_collapse).from_vertex();
    auto to_v = ref_mesh.edge(to_collapse).to_vertex();

    if (/*!ref_mesh.is_boundary(to_collapse) && !collapse_would_create_cord(ref_mesh, to_collapse)*/
        !ref_mesh.is_boundary(from_v) && !ref_mesh.is_boundary(to_v) &&
        !ref_mesh.is_deleted(from_v) && !ref_mesh.is_deleted(to_v)) {
            log(LogLevel::DEBUG)
                    << " ----------------------------------------------------------------------"
                    << std::endl;
            log(LogLevel::DEBUG) << " --------------------- trying to collapse edge (" << from_v << ", "
                    << to_v << ")" << std::endl;
        log(LogLevel::DEBUG) << " - quality threshold: " << quality_threshold << std::endl;


        if (ref_mesh.is_boundary(from_v) || ref_mesh.is_boundary(to_v)) {
            log(LogLevel::WARNING) << " WARNING - collapsing to or from boundary" << std::endl;
            //exit(EXIT_FAILURE);
        }

        NonConformingTetrahedralMesh::EigenVec ref_mesh_surviving_v_pos, mesh_surviving_v_pos;
        OVM::HalfEdgeHandle collapse_heh = ref_mesh.is_collapse_ok(to_collapse, ref_mesh_surviving_v_pos);
        if (collapse_heh.is_valid()) {
            log(LogLevel::DEBUG) << " -> collapse valid in ref mesh, checking in mesh..."
                    << std::endl;

            if (also_check_mesh_for_validity) {
                collapse_heh = mesh.is_collapse_ok(to_collapse, mesh_surviving_v_pos);
            }

            if (collapse_heh.is_valid()) {
                log(LogLevel::DEBUG) << " -> collapse also valid in mesh, computing quality"<<std::endl;


                if (ref_mesh.is_boundary(ref_mesh.from_vertex_handle(collapse_heh))) {
                    log(LogLevel::ERROR) << " ERROR - collapsing from the boundary" << std::endl;
                    return -1;
                }

                auto init_from_v_pos = ref_mesh.vertex(from_v);
                auto init_to_v_pos = ref_mesh.vertex(to_v);
                //log(LogLevel::INFO)<<" - conforming positions are currently "<<init_from_v_pos<<" and "<<init_to_v_pos<<std::endl;

                auto edge_incident_prop = ref_mesh.template request_cell_property<bool>("", false);
                //std::vector<typename MeshT::ElemHandleT> edge_neighbors;
                //ref_mesh.gather_neighbor_elements(to_collapse, edge_neighbors);

                /*auto angles_before = ref_mesh.compute_angles_in_neighborhood(to_collapse,
                                                                             edge_incident_prop,
                                                                             CONFORMING);*/

                auto quality_before = ref_mesh.compute_quality_vector_in_neighborhood(
                    to_collapse,
                    metric_code,
                    edge_incident_prop,
                    CONFORMING);
                //log(LogLevel::INFO) << " - min dihedral angle before: " << min_angle_deg_before << std::endl;

                /*log(LogLevel::DEBUG) << " - min / max dihedral angles before collapse: "
                        << std::setprecision(5)
                        << angles_before.min_dihedral_angle_deg << " / "
                        << angles_before.max_dihedral_angle_deg << std::endl;*/

                if (quality_before[0] > quality_threshold) {
                    log(LogLevel::DEBUG) << " -> worst element has quality above threshold, skipping" << std::endl;
                    return 5;
                }


                double ref_min_volume_before = ref_mesh.compute_min_volume_in_neighborhood(to_collapse);
                double     min_volume_before = mesh.compute_min_volume_in_neighborhood(to_collapse);

                ref_mesh.set_vertex(from_v, vec2vec(ref_mesh_surviving_v_pos));
                ref_mesh.set_vertex(to_v, vec2vec(ref_mesh_surviving_v_pos));

                //mesh.set_vertex(from_v, vec2vec(mesh_surviving_v_pos));
                //mesh.set_vertex(to_v, vec2vec(mesh_surviving_v_pos));

                //skip the "virtually" deleted cells
                for (auto ec_it = ref_mesh.ec_iter(to_collapse); ec_it.valid(); ec_it++) {
                    edge_incident_prop[*ec_it] = true;
                    //log(LogLevel::INFO) << " - skipping edge-incident element " << ref_ref_mesh.element_to_str(*ec_it) << std::endl;
                }

                //log(LogLevel::INFO)<<" - set end vertices duplicates to "<<mesh_surviving_v_pos.transpose()<<" to simulate collapse and compute new angles"<<std::endl;
                /*auto angles_after = ref_mesh.compute_angles_in_neighborhood(to_collapse,
                                                                            edge_incident_prop,
                                                                            CONFORMING);*/

                auto quality_after = ref_mesh.compute_quality_vector_in_neighborhood(
                    to_collapse,
                    metric_code,
                    edge_incident_prop,
                    CONFORMING);

                auto ref_min_volume_after = ref_mesh.compute_min_volume_in_neighborhood(to_collapse, edge_incident_prop);
                auto     min_volume_after = mesh.compute_min_volume_in_neighborhood(to_collapse, edge_incident_prop);

                /*log(LogLevel::DEBUG) << " - min / max dihedral angles after collapse: "
                        << std::setprecision(5)
                        << angles_after.min_dihedral_angle_deg << " / "
                        << angles_after.max_dihedral_angle_deg << std::endl;*/


                if ((ref_min_volume_after > volume_threshold || ref_min_volume_after > ref_min_volume_before) &&
                        (!also_check_mesh_for_validity || min_volume_after > volume_threshold || min_volume_after > min_volume_before)) {

                    if (quality_before.compare(quality_after, quality_vector_rel_eps) < 0) {

                        log(LogLevel::INFO) << " -> quality vector changed from " << std::endl
                                << quality_before << std::endl
                                << " to " << std::endl << quality_after << std::endl
                                << " and min volume changed from "<<min_volume_before<<" to "<<min_volume_after
                                << " --> performing collapse of edge (" << from_v << ", " << to_v
                                << ") with priority " << collapse_priority
                                << std::endl;


                        // auto to_collapse_he = ref_mesh.halfedge_handle(to_collapse, 0);
                        auto surviving_v = ref_mesh.collapse_edge(collapse_heh);
                        mesh.collapse_edge(collapse_heh);

                        ref_mesh.set_vertex(surviving_v, vec2vec(ref_mesh_surviving_v_pos));
                        mesh.set_vertex(surviving_v, vec2vec(mesh_surviving_v_pos));

                        // add impacted edges to list
                        for (auto out_he_it: ref_mesh.outgoing_halfedges(surviving_v)) {
                            impacted_edges.push_back(ref_mesh.edge_handle(out_he_it));
                        }

                        for (auto vc_it = ref_mesh.vc_iter(surviving_v); vc_it.valid(); vc_it++) {
                            auto op_hf = ref_mesh.vertex_opposite_halfface(*vc_it, surviving_v);

                            for (auto hfe_it = ref_mesh.hfe_iter(op_hf); hfe_it.valid(); hfe_it++) {
                                impacted_edges.push_back(*hfe_it);
                            }
                        }

                        return 0;
                    } else {
                        ref_mesh.set_vertex(from_v, init_from_v_pos);
                        ref_mesh.set_vertex(to_v,   init_to_v_pos);
                        return 1;
                    }
                }else {
                    ref_mesh.set_vertex(from_v, init_from_v_pos);
                    ref_mesh.set_vertex(to_v,   init_to_v_pos);
                    return 5;
                }
            } else {
                return 4;
            }
        } else {
            return 3;
        }
    } else {
        return 2;
    }
}


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
                            std::vector<EdgeHandle>& impacted_edges) {

    if (!ref_mesh.is_boundary(to_swap)) {
        if (!ref_mesh.is_deleted(to_swap)) {
            auto from_v = ref_mesh.edge(to_swap).from_vertex();
            auto to_v = ref_mesh.edge(to_swap).to_vertex();

            log(LogLevel::DEBUG) << " ---------------------------------------------- trying to swap edge (" << from_v << ", " << to_v << ")" <<
                    std::endl;
            log(LogLevel::DEBUG) << " - quality threshold: " << quality_threshold << std::endl;
            if (also_check_mesh_for_validity) {
                log(LogLevel::DEBUG) << " -> also checking mesh for validity of swap" << std::endl;
            }



            //check that the swap won't create cords
            log(LogLevel::DEBUG) << " - equatorial vertices: "<<std::endl;
            HalfEdgeHandle to_swap_heh = ref_mesh.halfedge_handle(to_swap, 0);
            std::vector<VertexHandle> equatorial_vertices;
            int boundary_count(0);
            for (auto hehf_it = ref_mesh.hehf_iter(to_swap_heh); hehf_it.valid(); hehf_it++) {
                equatorial_vertices.push_back(
                    ref_mesh.to_vertex_handle(ref_mesh.next_halfedge_in_halfface(to_swap_heh, *hehf_it)));
                log(LogLevel::DEBUG) << " - "<<equatorial_vertices.back()<<" at "<<ref_mesh.vertex(equatorial_vertices.back())<<std::endl;
                boundary_count += ref_mesh.is_boundary(equatorial_vertices.back());
            }
            log(LogLevel::DEBUG) << " -> that's " << equatorial_vertices.size() << " equatorial vertices" << std::endl;

            //auto collapse_target_vh = equatorial_vertices.back();
            if (boundary_count != equatorial_vertices.size()) {

                const int n_equ(equatorial_vertices.size());

                NonConformingTetrahedralMesh one_ring_mesh, one_ring_ref_mesh;
                OVM::EdgeHandle local_to_swap;
                setup_mesh_as_one_ring_submesh(ref_mesh, to_swap, CONFORMING,
                one_ring_ref_mesh, local_to_swap);

                if (also_check_mesh_for_validity) {
                    setup_mesh_as_one_ring_submesh(mesh, to_swap, CONFORMING,
                                                   one_ring_mesh, local_to_swap);
                }

                log(LogLevel::DEBUG) << " -> local edge to swap: " << local_to_swap << std::endl;

                auto global_v_prop = one_ring_ref_mesh.template request_vertex_property<VertexHandle>("global_vid");

                auto to_skip_prop = one_ring_ref_mesh.template request_cell_property<bool>("", false);

                auto quality_before = one_ring_ref_mesh.compute_quality_vector_in_neighborhood(local_to_swap,
                    metric_code,
                    to_skip_prop,
                    CONFORMING);
                log(LogLevel::DEBUG) << " -> quality before: " << quality_before << std::endl;

                if (quality_before[0] > quality_threshold) {
                    log(LogLevel::DEBUG) << " -> worst element has quality above threshold, skipping" << std::endl;
                    return 6;
                }

                auto local_mid_v = one_ring_ref_mesh.split_edge(local_to_swap);
                if (also_check_mesh_for_validity) {
                    one_ring_mesh.split_edge(local_to_swap);
                }


                auto local_from_v = one_ring_ref_mesh.edge(local_to_swap).from_vertex();
                auto local_to_v = one_ring_ref_mesh.edge(local_to_swap).to_vertex();

                log(LogLevel::DEBUG) << " - 1-ring mesh after split: " << std::endl;
                if(LogConfig::get_min_level() <= LogLevel::DEBUG){one_ring_ref_mesh.print_mesh_size();}

                log(LogLevel::DEBUG) << " - local from, to and mid-vertices: " << local_from_v << ", " << local_to_v <<
                        ", " << local_mid_v << std::endl;


                VertexHandle collapse_target_vh(-1);
                QualityVector best_quality = quality_before;
                bool found_valid_target(false);

                for (auto local_v: one_ring_ref_mesh.vertices()) {
                    if (local_v != local_from_v &&
                        local_v != local_to_v &&
                        local_v != local_mid_v) {

                        auto global_v = global_v_prop[local_v];
                        log(LogLevel::DEBUG) << " -- trying to collapse on equatorial vertex " << local_v << " -> " <<
                                global_v
                                << std::endl;
                        if (!ref_mesh.is_boundary(global_v)) {
                            auto local_to_collapse = one_ring_ref_mesh.find_halfedge(local_mid_v, local_v);
                            log(LogLevel::DEBUG) << " -- that's local eddge " << local_to_collapse << std::endl;

                            NonConformingTetrahedralMesh one_ring_ref_mesh_copy(one_ring_ref_mesh);
                            NonConformingTetrahedralMesh     one_ring_mesh_copy(one_ring_mesh);

                            log(LogLevel::DEBUG) << " -- 1-ring mesh copy: " << std::endl;
                            if(LogConfig::get_min_level() <= LogLevel::DEBUG){one_ring_ref_mesh.print_mesh_size();}

                            double ref_min_vol_before = one_ring_ref_mesh.min_element_volume();
                            double     min_vol_before = one_ring_mesh.min_element_volume();

                            auto local_surviving_v = one_ring_ref_mesh_copy.collapse_edge(local_to_collapse);
                            if (also_check_mesh_for_validity) {
                                one_ring_mesh_copy.collapse_edge(local_to_collapse);
                            }

                            double ref_min_vol_after = one_ring_ref_mesh_copy.min_element_volume();
                            double     min_vol_after = one_ring_mesh_copy.min_element_volume();

                            if ((ref_min_vol_after > volume_threshold || ref_min_vol_after > ref_min_vol_before)  &&
                                !one_ring_ref_mesh_copy.contains_degenerate_or_flipped_elements(CONFORMING) &&
                                (!also_check_mesh_for_validity ||
                                    ((min_vol_after > volume_threshold || min_vol_after > min_vol_before) &&
                                    !one_ring_mesh_copy.contains_degenerate_or_flipped_elements(CONFORMING)))) {

                                found_valid_target = true;
                                log(LogLevel::DEBUG) << " ---> valid collapse" << std::endl;


                                auto to_skip_prop_copy = one_ring_ref_mesh_copy.template request_cell_property<bool>("", false);
                                auto quality_after = one_ring_ref_mesh_copy.compute_quality_vector_in_neighborhood(
                                    local_surviving_v,
                                    metric_code,
                                    to_skip_prop_copy,
                                    CONFORMING);

                                if (best_quality.compare(quality_after, quality_vector_rel_eps) < 0) {
                                    best_quality = quality_after;
                                    collapse_target_vh = global_v;
                                    log(LogLevel::DEBUG) << " -> updated best quality vector to " << best_quality
                                            << " after collapsing on equatorial vertex " << global_v << std::endl;

                                    if (n_equ == 3) {
                                        log(LogLevel::DEBUG)
                                                << " -> swapping valence-3 edge, so no need to test other equatorial vertices"
                                                << std::endl;
                                        break;
                                    }
                                }
                            }
                        } else {
                            log(LogLevel::DEBUG) << " -- equatorial vertex is boundary, skipping" << std::endl;
                        }
                    }
                }


                if ((metric_code == CONSTANT_METRIC && collapse_target_vh.is_valid()) ||
                    quality_before.compare(best_quality) < 0) {
                    //do the swap in the actual mesh

                    log(LogLevel::DEBUG) << " -> swap was successful in 1-ring mesh, reproducing it in actual mesh" <<
                            std::endl;
                    log(LogLevel::INFO) << " -> quality vector changed from " << std::endl << quality_before
                            << std::endl
                            << " to " << std::endl << best_quality << std::endl
                            << " --> performing swap of edge "<<edge_to_str(ref_mesh, to_swap)<<" -> " << collapse_target_vh << " with priority " << swap_priority<<" and "<<equatorial_vertices.size()<<" equatorial vertices "
                            << std::endl;

                    auto mid_v = ref_mesh.split_edge(to_swap);
                    mesh.split_edge(to_swap);
                    log(LogLevel::DEBUG) << " - mid-vertex: " << mid_v << std::endl;

                    /*if (json_data["enable_all_safety_checks"]) {
                        if (ref_mesh.contains_degenerate_or_flipped_elements(CONFORMING) ||
                            (also_check_mesh_for_validity && mesh.contains_degenerate_or_flipped_elements(CONFORMING))) {
                            log(LogLevel::ERROR) <<
                                    " ERROR - either ref mesh or mesh contains degenerate or flipped elements right after splitting (for swap)" <<
                                    std::endl;

                            for (auto elh: ref_mesh.elements()) {
                                if (ref_mesh.is_flipped_or_degenerate(elh, CONFORMING)) {
                                    log(LogLevel::ERROR) << " --> ref element " << ref_mesh.element_to_str(elh) << " with volume " <<
                                            ref_mesh.element_volume(elh) << std::endl;
                                }
                                if (mesh.is_flipped_or_degenerate(elh, CONFORMING)) {
                                    log(LogLevel::ERROR) << " --> element " << mesh.element_to_str(elh) << " with volume " <<
                                            mesh.element_volume(elh) << std::endl;
                                }
                            }
                            return -1;
                        }
                    }*/

                    auto to_collapse_heh = ref_mesh.find_halfedge(mid_v, collapse_target_vh);
                    log(LogLevel::DEBUG) << " - edge to collapse: " << to_collapse_heh << " (" << mid_v << ", " <<
                            collapse_target_vh << ")" << std::endl;

                    auto surviving_v = ref_mesh.collapse_edge(to_collapse_heh);
                    if (surviving_v != collapse_target_vh) {
                        log(LogLevel::ERROR) << " ERROR - collapse target vertex is "
                                << collapse_target_vh << " but surviving vertex is "
                                << surviving_v << std::endl;
                        return -1;
                    }
                    mesh.collapse_edge(to_collapse_heh);

                    if (json_data["enable_all_safety_checks"]) {
                        if (ref_mesh.contains_degenerate_or_flipped_elements(CONFORMING) ||
                            (also_check_mesh_for_validity && mesh.contains_degenerate_or_flipped_elements(CONFORMING))) {
                            log(LogLevel::ERROR) <<
                                    " ERROR - either ref mesh or mesh contains degenerate or flipped elements right after collapsing (for swap)" <<
                                    std::endl;

                            for (auto elh: ref_mesh.elements()) {
                                if (ref_mesh.is_flipped_or_degenerate(elh, CONFORMING)) {
                                    log(LogLevel::ERROR) << " --> element " << ref_mesh.element_to_str(elh) << " with volume " <<
                                            ref_mesh.element_volume(elh) << std::endl;
                                }
                                if (mesh.is_flipped_or_degenerate(elh, CONFORMING)) {
                                    log(LogLevel::ERROR) << " --> element " << mesh.element_to_str(elh) << " with volume " <<
                                            mesh.element_volume(elh) << std::endl;
                                }
                            }
                            return -1;
                        }
                    }

                    // add impacted edges to list
                    std::vector<VertexHandle> end_vertices = {from_v, to_v};
                    for (auto end_v: end_vertices) {
                        for (auto eq_v: equatorial_vertices) {
                            auto heh = ref_mesh.find_halfedge(end_v, eq_v);
                            if (!heh.is_valid()) {
                                log(LogLevel::ERROR) << " couldn't find pole edge (" << end_v << ", " <<eq_v<<")"<< std::endl;
                                return -1;
                            }
                            impacted_edges.push_back(ref_mesh.edge_handle(heh));
                        }
                    }

                    for (int i(0); i < equatorial_vertices.size() - 1; i++) {
                        for (int j(i+1); j < equatorial_vertices.size(); j++) {

                            auto eq_from_v = equatorial_vertices[i];
                            auto eq_to_v = equatorial_vertices[j];
                            auto heh = ref_mesh.find_halfedge(eq_from_v, eq_to_v);

                            if (heh.is_valid()) {
                                log(LogLevel::DEBUG) << " - added new edge (" << eq_from_v << ", " << eq_to_v <<
                                        ") to list" << std::endl;
                                impacted_edges.push_back(ref_mesh.edge_handle(heh));
                            }
                        }
                    }

                    return 0;
                }

                if (found_valid_target) {
                    //log(LogLevel::INFO)<<" -> worse angles, skipping swap"<<std::endl;
                    return 1;
                } else {
                    //log(LogLevel::INFO) << " -> couldn't collapse to any equatorial vertex, cannot swap edge" << std::endl;
                    return 2;
                }
            } else {
                //log(LogLevel::INFO)<<" -> opposite vertices are all boundary, skipping swap"<<std::endl;
                return 3;
            }
        } else {
            //log(LogLevel::INFO)<<" -> edge is either deleted or no longer valence 3"<<std::endl;
            return 4;
        }
    } else {
        return 5;
    }
}


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
                            std::vector<EdgeHandle>& impacted_edges) {

    std::vector<VertexHandle> equatorial_vertices;
    auto to_split_heh = ref_mesh.halfedge_handle(to_split, 0);
    for (auto hehf_it = ref_mesh.hehf_iter(to_split_heh); hehf_it.valid(); hehf_it++) {
        equatorial_vertices.push_back(
            ref_mesh.to_vertex_handle(ref_mesh.next_halfedge_in_halfface(to_split_heh, *hehf_it)));
    }

    auto from_v = ref_mesh.edge(to_split).from_vertex();
    auto to_v = ref_mesh.edge(to_split).to_vertex();

    if ((enable_boundary_split || (!ref_mesh.is_boundary(from_v) || !ref_mesh.is_boundary(to_v))) &&
            !ref_mesh.is_deleted(from_v) && !ref_mesh.is_deleted(to_v)) {
        log(LogLevel::DEBUG)
                << " ----------------------------------------------------------------------"
                << std::endl;
        log(LogLevel::DEBUG) << " --------------------- trying to split edge (" << from_v << ", "
                             << to_v
                             << ") with priority " << split_priority << std::endl;
        log(LogLevel::DEBUG) << " - quality threshold: " << quality_threshold << std::endl;
        log(LogLevel::DEBUG) << " - volume threshold: " << volume_threshold << std::endl;
        log(LogLevel::DEBUG) << " - also checking mesh for validity: " << also_check_mesh_for_validity << std::endl;


        //gather equatorial vertices, for the glueing set upkeeping
        auto first_elh = ref_mesh.ec_iter(to_split);
        std::map<VertexHandle, int> equatorial_vertex_to_glueing_id_map;

        NonConformingTetrahedralMesh one_ring_mesh, one_ring_ref_mesh;
        OVM::EdgeHandle local_to_split;
        setup_mesh_as_one_ring_submesh(ref_mesh, to_split, CONFORMING, one_ring_ref_mesh,
                                       local_to_split);
        setup_mesh_as_one_ring_submesh(mesh, to_split, CONFORMING, one_ring_mesh,
                                       local_to_split);


        if (json_data["enable_all_safety_checks"]) {
            if (one_ring_ref_mesh.contains_degenerate_or_flipped_elements(CONFORMING)) {
                log(LogLevel::ERROR) <<
                        " ERROR - 1-ring mesh contains degenerate or flipped elements righ after setup" <<
                        std::endl;

                for (auto elh: one_ring_ref_mesh.elements()) {
                    if (one_ring_ref_mesh.is_flipped_or_degenerate(elh, CONFORMING)) {
                        log(LogLevel::ERROR) << " --> element " << one_ring_ref_mesh.element_to_str(elh) <<
                                             " with volume " << one_ring_ref_mesh.element_volume(elh) << std::endl;
                    }
                }
                return -1;
            }
        }

        //auto to_skip_prop = one_ring_ref_mesh.template request_cell_property<bool>("", false);
        auto to_skip_prop = one_ring_ref_mesh.template request_cell_property<bool>("", false);

        auto angles_before = one_ring_ref_mesh.compute_angles_in_neighborhood(local_to_split,
                                                                              to_skip_prop,
                                                                              CONFORMING);

        auto quality_before = one_ring_ref_mesh.compute_quality_vector_in_neighborhood(local_to_split,
                                                                                       metric_code,
                                                                                       to_skip_prop,
                                                                                       CONFORMING);

        auto ref_min_volume_before = one_ring_ref_mesh.compute_min_volume_in_neighborhood(local_to_split);
        auto min_volume_before = one_ring_mesh.compute_min_volume_in_neighborhood(local_to_split);

            log(LogLevel::DEBUG) << " - min / max dihedral angles before split: "
                    << std::setprecision(4)
                    << angles_before.min_dihedral_angle_deg << " / "
                    << angles_before.max_dihedral_angle_deg << std::endl;

        if (quality_before[0] > quality_threshold) {
            log(LogLevel::DEBUG) << " -> worst element has quality above threshold, skipping" << std::endl;
            return 1;
        }

        auto local_mid_v = one_ring_ref_mesh.split_edge(local_to_split);
        one_ring_mesh.split_edge(local_to_split);

        //note: this works but makes split happen "too often".
        // It should probably be used after the split is done
        // That or only if the max angle is above a certain threshold
        auto angles_after_split = one_ring_ref_mesh.compute_angles_in_neighborhood(local_to_split,
                                                                                   to_skip_prop,
                                                                                   CONFORMING);

        auto quality_after_split = one_ring_ref_mesh.compute_quality_vector_in_neighborhood(local_to_split,
                                                                                            metric_code,
                                                                                            to_skip_prop,
                                                                                            CONFORMING);

            log(LogLevel::DEBUG) << " - min / max dihedral angles after split: "
                    << std::setprecision(5)
                    << angles_after_split.min_dihedral_angle_deg << " / "
                    << angles_after_split.max_dihedral_angle_deg << std::endl;



        auto initial_pos = one_ring_ref_mesh.vertex(local_mid_v);
        //weird bug but I don't have time to deal with it now.
        // Hopefully it will manifest in an easier to figure out way later.
        Vec3d mid_v_opti_pos;
        if (enable_pre_split_opti) {
            if (optimize_mesh<NonConformingTetrahedralMesh, IVL_B_3D_AD>(one_ring_ref_mesh,
                                                                         one_ring_ref_mesh.total_volume() /
                                                                         one_ring_ref_mesh.n_cells(),
                                                                         10,
                                                                         true,
                                                                         nullptr,
                                                                         false,
                                                                         json_data,
                                                                         true)) {
                log(LogLevel::ERROR) << " error while optimizing 1-ring mesh" << std::endl;
                return -1;
            }
            mid_v_opti_pos = one_ring_ref_mesh.vertex(local_mid_v);
        }

        auto angles_after = one_ring_ref_mesh.compute_angles_in_neighborhood(local_mid_v,
                                                                             to_skip_prop,
                                                                             CONFORMING);

        auto quality_after = one_ring_ref_mesh.compute_quality_vector_in_neighborhood(local_mid_v,
                                                                                      metric_code,
                                                                                      to_skip_prop,
                                                                                      CONFORMING);

        auto ref_min_volume_after = one_ring_ref_mesh.compute_min_volume_in_neighborhood(local_mid_v);
        auto     min_volume_after = one_ring_mesh.compute_min_volume_in_neighborhood(local_mid_v);


        log(LogLevel::DEBUG) << " - min / max dihedral angles after mid-vertex opti: "
                    << std::setprecision(6)
                    << angles_after.min_dihedral_angle_deg << " / "
                    << angles_after.max_dihedral_angle_deg << std::endl;

            log(LogLevel::DEBUG) << " -> quality vector changed from " << std::endl
                    << quality_before << std::endl
                    << " to " << std::endl << quality_after << std::endl;

            log(LogLevel::DEBUG) << " -> quality rel epsilon rel:  " << quality_vector_rel_eps << std::endl;

        log(LogLevel::DEBUG) << " -> ref mesh volume change: " << ref_min_volume_before <<" -> "<<ref_min_volume_after << std::endl;
        log(LogLevel::DEBUG) << " ->     mesh volume change: " << min_volume_before <<" -> "<<min_volume_after << std::endl;
        log(LogLevel::DEBUG) << " ->     mesh volume after >          threshold: " << (min_volume_after > volume_threshold) << std::endl;
        log(LogLevel::DEBUG) << " ->     mesh volume after > mesh volume before: " << (min_volume_after > min_volume_before) << std::endl;
        log(LogLevel::DEBUG) << " ->                            first condition: " << (ref_min_volume_after > volume_threshold || ref_min_volume_after > ref_min_volume_before) << std::endl;
        log(LogLevel::DEBUG) << " ->                           second condition: " << (!also_check_mesh_for_validity || min_volume_after > volume_threshold || min_volume_after > min_volume_before) << std::endl;


        if((ref_min_volume_after > volume_threshold || ref_min_volume_after > ref_min_volume_before) &&
           (!also_check_mesh_for_validity || min_volume_after > volume_threshold || min_volume_after > min_volume_before)) {
            if ((metric_code == CONSTANT_METRIC || quality_before.compare(quality_after, quality_vector_rel_eps) < 0)) {

                log(LogLevel::INFO) << " ----------------------------------------------------------------------" << std::endl;
                log(LogLevel::INFO) << " -> quality vector changed from " << std::endl << quality_before << std::endl
                        << " to " << std::endl << quality_after << std::endl
                        << " and min volume changed from "<<ref_min_volume_before<<" to "<<ref_min_volume_after<<" --> performing split of edge "<<to_split<<": (" << from_v << ", " << to_v
                        << ") with priority " << split_priority
                        << std::endl;


                auto mid_v = ref_mesh.split_edge(to_split);
                mesh.split_edge(to_split);

                log(LogLevel::DEBUG) << " new mid-vertex is " << mid_v << std::endl;

                if (json_data["enable_all_safety_checks"]) {
                    if (ref_mesh.contains_degenerate_or_flipped_elements(CONFORMING)) {
                        log(LogLevel::ERROR) <<
                                " ERROR - ref mesh contains degenerate or flipped elements right after splitting"
                                << std::endl;

                        for (auto elh: ref_mesh.elements()) {
                            if (ref_mesh.is_flipped_or_degenerate(elh, CONFORMING)) {
                                log(LogLevel::ERROR) << " --> element " << ref_mesh.element_to_str(elh) << " with volume "
                                        << ref_mesh.element_volume(elh) << std::endl;
                            }
                        }
                        return -1;
                    }
                }

                if (enable_pre_split_opti) {
                    ref_mesh.set_vertex(mid_v, mid_v_opti_pos);

                    if (json_data["enable_all_safety_checks"]) {
                        std::vector<CellHandle> neighbors;
                        ref_mesh.gather_neighbor_elements(mid_v, neighbors);
                        for (auto elh: neighbors) {
                            if (ref_mesh.is_flipped_or_degenerate(elh, CONFORMING)) {
                                log(LogLevel::ERROR) << " ERROR - moving mid-vertex duplicates to "
                                        << mid_v_opti_pos
                                        << " flipped or degenerated element "
                                        << ref_mesh.element_to_str(elh) << std::endl;
                                return -1;
                            }
                        }
                    }
                }

                if (json_data["enable_all_safety_checks"]) {
                    if (ref_mesh.contains_degenerate_or_flipped_elements(CONFORMING)) {
                        log(LogLevel::ERROR) <<
                                " ERROR - ref mesh contains degenerate or flipped elements after updating glueing props and stuff"
                                << std::endl;

                        for (auto elh: ref_mesh.elements()) {
                            if (ref_mesh.is_flipped_or_degenerate(elh, CONFORMING)) {
                                log(LogLevel::ERROR) << " --> element " << ref_mesh.element_to_str(elh) << " with volume "
                                        << ref_mesh.element_volume(elh) << std::endl;
                            }
                        }
                        return -1;
                    }
                }


                // add impacted edges to list
                for (auto out_he_it: ref_mesh.outgoing_halfedges(mid_v)) {
                    impacted_edges.push_back(ref_mesh.edge_handle(out_he_it));
                }

                for (int i(0); i < equatorial_vertices.size() - 1; i++) {
                    for (int j(i+1); j < equatorial_vertices.size(); j++) {

                        auto eq_from_v = equatorial_vertices[i];
                        auto eq_to_v = equatorial_vertices[j];
                        auto heh = ref_mesh.find_halfedge(eq_from_v, eq_to_v);

                        if (heh.is_valid()) {
                            log(LogLevel::DEBUG) << " - added new edge (" << eq_from_v << ", " << eq_to_v <<
                                    ") to list" << std::endl;
                            impacted_edges.push_back(ref_mesh.edge_handle(heh));
                        }
                    }
                }

                return 0;
            }
            return 1;
        }
        return 2;
    }
    return 3;
}
