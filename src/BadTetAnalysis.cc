#include "BadTetAnalysis.hh"
#include "NonConformingHelpers.hh"

using namespace OpenVolumeMesh;

void set_mesh_as_regular_tet(TetrahedralMesh& mesh){
    mesh.add_vertex({0, 0, 0});
    mesh.add_vertex({1, 0, 0});
    mesh.add_vertex({0.5, 0.5 * std::sqrt(3.0), 0});
    mesh.add_vertex({0.5, std::sqrt(3.0) / 6.0, std::sqrt(6.0) / 3.0});

    mesh.add_cell({VertexHandle(0),
                   VertexHandle(1),
                   VertexHandle(2),
                   VertexHandle(3)});
}

void run_cap_tet_experiments(bool change_domain_or_codomain,
                             std::vector<PerElementStats>& stats){

    TetrahedralMesh ref_mesh;
    set_mesh_as_regular_tet(ref_mesh);

    TetrahedralMesh  mesh = ref_mesh;

    VertexHandle moving_vertex(3);
    auto initial_pos = mesh.vertex(moving_vertex);

    double h_min = 1e-30;
    double h_max = 1e0;
    int step_count = 50;

    double pow_range = std::log10(h_max / h_min);
    log(LogLevel::INFO)<<" pow range = "<<pow_range<<std::endl;

    IO::FileManager fm;


    for(int i(0);i < step_count;i++){

        double n = std::log10(h_min) + i * pow_range / step_count;
        double h = std::pow(10, n);

        log(LogLevel::INFO)<<" ------- step "<<i<<", pow = "<<n<<", h = "<<h<<std::endl;

        auto pos = initial_pos;
        pos[2] *= h;
        mesh.set_vertex(moving_vertex, pos);

        fm.writeFile("cap_tet_h" + std::to_string(h) + ".ovm", mesh);

        optimize_tet_and_fill_stats(ref_mesh, mesh, CellHandle(0), change_domain_or_codomain, false, 1.0, stats);
        stats.back().param1 = h;
    }

}


void run_sliver_tet_experiments(bool change_domain_or_codomain,
                                std::vector<PerElementStats>& stats){

    TetrahedralMesh ref_mesh;
    set_mesh_as_regular_tet(ref_mesh);
    TetrahedralMesh  mesh = ref_mesh;

    std::vector<VertexHandle> moving_vertices = {VertexHandle(2), VertexHandle(3)};

    // π/3 rotation around the x axis to get the rotating edge on the z axis
    Eigen::Matrix3d init_rot = Eigen::AngleAxis<double>(M_PI / 3.0, Eigen::Vector3d::UnitX()).toRotationMatrix();
    for(auto v: moving_vertices){
        Eigen::Vector3d pos = vec2vec(mesh.vertex(v));
        pos = (init_rot * pos).eval();
        mesh.set_vertex(v, vec2vec(pos));
    }

    std::vector<Vec3d> initial_pos = {mesh.vertex(VertexHandle(2)), mesh.vertex(VertexHandle(3))};

    double h_min = 1e-30;
    double h_max = 1e0;
    int step_count = 50;

    double pow_range = std::log10(h_max / h_min);

    IO::FileManager fm;

    for(int i(0);i < step_count;i++){

        double n = std::log10(h_min) + i * pow_range / step_count;
        double h = std::pow(10, n);

        log(LogLevel::INFO)<<" ------- step "<<i<<", pow = "<<n<<", h = "<<h<<std::endl;

        for(int i(0); i<2; i++){
            auto pos = initial_pos[i];
            pos[2] *= h;
            mesh.set_vertex(moving_vertices[i], pos);
        }

        fm.writeFile("sliver_tet_h" + std::to_string(h) + ".ovm", mesh);

        optimize_tet_and_fill_stats(ref_mesh, mesh, CellHandle(0), change_domain_or_codomain, false, 1.0, stats);
        stats.back().param1 = h;
    }

}


void run_spade_tet_experiments(bool change_domain_or_codomain,
                               std::vector<PerElementStats>& stats){

    std::vector<VertexHandle> moving_vertices = {VertexHandle(2), VertexHandle(3)};

    double spade_epsilon = 1e-10;

    IO::FileManager fm;

    TetrahedralMesh ref_mesh;
    set_mesh_as_regular_tet(ref_mesh);
    TetrahedralMesh  mesh = ref_mesh;

    // π/3 rotation around the x axis to get the rotating edge on the z axis
    Eigen::Matrix3d init_rot = Eigen::AngleAxis<double>(0.5 * M_PI - std::acos(1.0 / 3.0) - spade_epsilon,
                                                        Eigen::Vector3d::UnitX()).toRotationMatrix();
    for (auto v: moving_vertices) {
        Eigen::Vector3d pos = vec2vec(mesh.vertex(v));
        pos = (init_rot * pos).eval();
        mesh.set_vertex(v, vec2vec(pos));
    }
    fm.writeFile("init_spade_eps_" + std::to_string(spade_epsilon) + ".ovm", mesh);

    std::vector <Vec3d> initial_pos = {mesh.vertex(VertexHandle(2)), mesh.vertex(VertexHandle(3))};

    double h_min = 1e-30;
    double h_max = 1e0;
    int step_count = 50;

    double pow_range = std::log10(h_max / h_min);


    for (int i(0); i < step_count; i++) {

        double n = std::log10(h_min) + i * pow_range / step_count;
        double h = std::pow(10, n);

        log(LogLevel::INFO) << " ------- step " << i << ", pow = " << n << ", h = " << h << std::endl;

        for (int i(0); i < 2; i++) {
            auto pos = initial_pos[i];
            pos[2] *= h;
            mesh.set_vertex(moving_vertices[i], pos);
        }

        fm.writeFile("spade_tet_h" + std::to_string(h) + ".ovm", mesh);

        optimize_tet_and_fill_stats(ref_mesh, mesh, CellHandle(0), change_domain_or_codomain, false, 1.0, stats);
        stats.back().param1 = h;
    }
}

//not sure about the name
void run_triangle_tet_experiments(bool change_domain_or_codomain,
                                  std::vector<PerElementStats>& stats){

    TetrahedralMesh ref_mesh;
    set_mesh_as_regular_tet(ref_mesh);

    TetrahedralMesh  mesh = ref_mesh;

    VertexHandle moving_vertex(3);
    double triangle_epsilon = 1e-10;

    auto init_pos = mesh.vertex(moving_vertex);
    mesh.set_vertex(moving_vertex, mesh.vertex(VertexHandle(2)) + Vec3d(0, -triangle_epsilon, init_pos[2]));

    auto initial_pos = mesh.vertex(moving_vertex);


    double h_min = 1e-30;
    double h_max = 1e0;
    int step_count = 50;

    double pow_range = std::log10(h_max / h_min);
    log(LogLevel::INFO)<<" pow range = "<<pow_range<<std::endl;

    IO::FileManager fm;
    fm.writeFile("init_triangle_tet.ovm", mesh);

    for(int i(0);i < step_count;i++){

        double n = std::log10(h_min) + i * pow_range / step_count;
        double h = std::pow(10, n);

        log(LogLevel::INFO)<<" ------- step "<<i<<", pow = "<<n<<", h = "<<h<<std::endl;

        auto pos = initial_pos;
        pos[2] *= h;
        mesh.set_vertex(moving_vertex, pos);

        fm.writeFile("triangle_tet_h" + std::to_string(h) + ".ovm", mesh);

        optimize_tet_and_fill_stats(ref_mesh, mesh, CellHandle(0), change_domain_or_codomain, false, 1.0, stats);
        stats.back().param1 = h;
    }
}


void run_sliver_tet_rotation_experiments(bool change_domain_or_codomain,
                                         std::vector<PerElementStats>& stats){

    TetrahedralMesh ref_mesh;
    set_mesh_as_regular_tet(ref_mesh);
    TetrahedralMesh  mesh = ref_mesh;

    std::vector<VertexHandle> moving_vertices = {VertexHandle(2), VertexHandle(3)};

    // π/3 rotation around the x axis to get the rotating edge on the z axis
    Eigen::Matrix3d init_rot = Eigen::AngleAxis<double>(M_PI / 3.0, Eigen::Vector3d::UnitX()).toRotationMatrix();
    for(auto v: moving_vertices){
        Eigen::Vector3d pos = vec2vec(mesh.vertex(v));
        pos = (init_rot * pos).eval();
        mesh.set_vertex(v, vec2vec(pos));
    }

    std::vector<Vec3d> initial_pos = {mesh.vertex(VertexHandle(2)), mesh.vertex(VertexHandle(3))};
    double h = 1e-5;
    int step_count = 50;


    for(int i(0); i<2; i++){
        auto pos = initial_pos[i];
        pos[2] *= h;
        mesh.set_vertex(moving_vertices[i], pos);
    }


    double alpha_incr = 2 * M_PI / step_count;
    Eigen::Matrix3d incr_rot = Eigen::AngleAxis<double>(alpha_incr, Eigen::Vector3d::UnitX()).toRotationMatrix();


    IO::FileManager fm;

    double alpha(0);
    for(int i(0);i < step_count;i++){

        for(auto v: moving_vertices){
            Eigen::Vector3d pos = vec2vec(mesh.vertex(v));
            pos = (incr_rot * pos).eval();
            mesh.set_vertex(v, vec2vec(pos));
        }

        fm.writeFile("sliver_tet_angle_" + std::to_string(alpha) + ".ovm", mesh);

        optimize_tet_and_fill_stats(ref_mesh, mesh, CellHandle(0), change_domain_or_codomain, false, 1.0, stats);
        stats.back().param1 = alpha;
        alpha += alpha_incr;
    }

}

void run_flat_tet_experiments(bool change_domain_or_codomain,
                              std::vector<PerElementStats>& stats){

    TetrahedralMesh ref_mesh;
    set_mesh_as_regular_tet(ref_mesh);

    TetrahedralMesh  mesh = ref_mesh;

    VertexHandle moving_vertex(3);

    std::vector<double> h_values = {1e-1, 1e-5, 1e-10, 1e-15, 1e-20};


    Vec2d bl_corner(-1, -1);
    Vec2d ur_corner(1, 1);

    auto centroid = mesh.vertex(moving_vertex);
    centroid[2] = 0;
    log(LogLevel::INFO)<<" - centroid at "<<centroid<<std::endl;
    log(LogLevel::INFO)<<" - vertices moved to "<<std::endl;
    for(auto v: mesh.vertices()){
        mesh.set_vertex(v, mesh.vertex(v) - centroid);
        log(LogLevel::INFO)<<" -- v"<<v<<" at "<<mesh.vertex(v)<<std::endl;
    }

    int grid_dim = 32;

    bool export_meshes = true;

    IO::FileManager fm;
    for(auto h: h_values) {
        for (int i(0); i < grid_dim; i++) {
            for (int j(0); j < grid_dim; j++) {

                double x = bl_corner[0] + i * (ur_corner - bl_corner)[0] / (grid_dim-1);
                double y = bl_corner[1] + j * (ur_corner - bl_corner)[1] / (grid_dim-1);
                Vec3d pos = {x, y, h};
                log(LogLevel::INFO)<<" - grid point ("<<i<<", "<<j<<"), pos = "<<pos<<std::endl;

                mesh.set_vertex(moving_vertex, pos);

                if(export_meshes) {
                    fm.writeFile("flat_tet_" + std::to_string(i) + "_" + std::to_string(j) + ".ovm", mesh);
                }

                optimize_tet_and_fill_stats(ref_mesh, mesh, CellHandle(0), change_domain_or_codomain, false, 1.0, stats);
                stats.back().param1 = x;
                stats.back().param2 = y;
                stats.back().param3 = h;
            }
        }
    }
}

void run_flat_tet_exp_stretch_experiments(bool change_domain_or_codomain,
                                          std::vector<PerElementStats>& stats){

    TetrahedralMesh ref_mesh;
    set_mesh_as_regular_tet(ref_mesh);

    TetrahedralMesh  mesh = ref_mesh;

    VertexHandle moving_vertex(3);

    std::vector<double> h_values = {1e0, 1e-1, 1e-5, 1e-10, 1e-15, 1e-20};
    h_values = {1e0, 1e-1};

    double neg_pow = -4;
    double pos_pow = 4;
    double pow_range = pos_pow - neg_pow; //
    log(LogLevel::INFO)<<" pow range = "<<pow_range<<std::endl;

    auto centroid = mesh.vertex(moving_vertex);
    centroid[2] = 0;
    log(LogLevel::INFO)<<" - centroid at "<<centroid<<std::endl;
    log(LogLevel::INFO)<<" - vertices moved to "<<std::endl;
    for(auto v: mesh.vertices()){
        mesh.set_vertex(v, mesh.vertex(v) - centroid);
        log(LogLevel::INFO)<<" -- v"<<v<<" at "<<mesh.vertex(v)<<std::endl;
    }

    int grid_dim = 32;

    bool export_meshes = true;

    IO::FileManager fm;
    for(auto h: h_values) {
        for (int i(0); i < grid_dim; i++) {

            double n = neg_pow + i * pow_range / (grid_dim - 1);
            log(LogLevel::INFO)<<" - n = "<<n<<std::endl;
            double x = (n < 0 ? -1 : 1) * std::pow(10, std::abs(n));
            log(LogLevel::INFO)<<" - x = "<<x<<std::endl;

            for (int j(0); j < grid_dim; j++) {

                double m = neg_pow + j * pow_range / (grid_dim - 1);
                log(LogLevel::INFO)<<" - m = "<<m<<std::endl;
                double y = (m < 0 ? -1 : 1) * std::pow(10, std::abs(m));
                log(LogLevel::INFO)<<" - y = "<<y<<std::endl;

                Vec3d pos = {x, y, h};
                log(LogLevel::INFO)<<" - grid point ("<<i<<", "<<j<<") pos = "<<pos<<std::endl;

                mesh.set_vertex(moving_vertex, pos);

                if(export_meshes) {
                    fm.writeFile("flat_tet_exp_stretch_" + std::to_string(i) + "_" + std::to_string(j) + ".ovm", mesh);
                }

                optimize_tet_and_fill_stats(ref_mesh, mesh, CellHandle(0), change_domain_or_codomain, false, 1.0, stats);
                stats.back().param1 = x;
                stats.back().param2 = y;
                stats.back().param3 = h;
            }
        }
    }
}




void swap_meshes(TetrahedralMesh& domain_mesh,
                 TetrahedralMesh& codomain_mesh){
    TetrahedralMesh temp = codomain_mesh;
    codomain_mesh = domain_mesh;
    domain_mesh = temp;
}

void optimize_tet_and_fill_stats(const TetrahedralMesh& ref_mesh,
                                 const TetrahedralMesh& mesh,
                                 OpenVolumeMesh::CellHandle ch,
                                 bool change_domain_or_codomain,
                                 bool constrain_boundary,
                                 double scale_factor,
                                 std::vector<PerElementStats>& stats){

    PerElementStats stat;
    stat.domain_or_codomain_deformation = change_domain_or_codomain;
    stat.min_dihedral_angle_deg = std::numeric_limits<double>::max();
    for(auto ce_it = mesh.ce_iter(ch); ce_it.valid(); ce_it++){
        double angle = compute_dihedral_angle_deg(mesh, ch, *ce_it);
        stat.min_dihedral_angle_deg = std::min(stat.min_dihedral_angle_deg, angle);
        stat.avg_dihedral_angle_deg += angle;
        stat.max_dihedral_angle_deg = std::max(stat.max_dihedral_angle_deg, angle);
    }
    stat.avg_dihedral_angle_deg /= 6.0;
    stat.volume = compute_volume(mesh, ch);

    if(change_domain_or_codomain){
        //log(LogLevel::INFO)<<" -> domain mesh was deformed "<<std::endl;
        optimize_tet(mesh, ref_mesh, ch, constrain_boundary, scale_factor, stat);
    } else {
        //log(LogLevel::INFO)<<" -> codomain mesh was deformed "<<std::endl;
        optimize_tet(ref_mesh, mesh, ch, constrain_boundary, scale_factor,  stat);
    }

    stats.push_back(stat);
    //log(LogLevel::INFO)<<" --------------------------------------------------------------------------- "<<std::endl;
}


void optimize_tet(const TetrahedralMesh& ref_mesh,
                  const TetrahedralMesh& mesh,
                  OpenVolumeMesh::CellHandle ch,
                  bool constrain_boundary,
                  double scale_factor,
                  PerElementStats& stat){

    std::array<Vec3d, 4> domain_positions;
    std::array<Vec3d, 4> codomain_positions;
    std::array<bool, 4> constraints_mask = {false, false, false, false};

    if(constrain_boundary){
        //log(LogLevel::INFO)<<" -> constraining boundary "<<std::endl;
    }

    auto verts = mesh.get_cell_vertices(ch);

    for (int local_vid(0);local_vid < verts.size();local_vid++) {
        const auto &vh = verts[local_vid];

        //log(LogLevel::INFO)<<" - v"<<vh<<" at "<<mesh.vertex(vh)<<std::endl;

        domain_positions[local_vid] = ref_mesh.vertex(vh);
        domain_positions[local_vid] *= scale_factor;

        codomain_positions[local_vid] = mesh.vertex(vh);

        if(constrain_boundary && mesh.is_boundary(vh)) {
            constraints_mask[local_vid] = true;
        }

        stat.n_boundary_vertices += mesh.is_boundary(vh);
    }


    optimize_tet(domain_positions, codomain_positions, constraints_mask, stat);

}


void optimize_tet(const std::array<Vec3d, 4>& domain_positions,
                  const std::array<Vec3d, 4>& codomain_positions,
                  const std::array<bool, 4>& constraints_mask,
                  PerElementStats& stat){

    // ---- opti
    const int n_vars(3 * 4);

    using ElemT = SD_B_3D_AD;
    COMISO::FiniteElementProblem fe_problem(n_vars);
    std::vector<COMISO::LinearConstraint> eq_constraints;
    eq_constraints.reserve(n_vars);
    COMISO::FiniteElementSet <ElemT> fe_barrier("barrier elements");

    int global_vid(0);

    double w = 1;

    typename ElemT::VecI vi_barrier;
    typename ElemT::VecC vc_barrier;
    typename ElemT::VecV vx_barrier;


    //log(LogLevel::INFO)<<" -- current vertices positions: "<<std::endl;
    std::vector<Eigen::Vector3d> P_ref;//only works for tets
    for (int local_vid(0); local_vid < 4; local_vid++) {

        int vid = 3 * local_vid;
        //auto v_pos = mesh.vertex(vh);
        auto v_pos = codomain_positions[local_vid];

        for (int i(0);i < 3;i++) {
            vi_barrier[vid + i] = vid + i;

            fe_problem.x()[vid + i] = v_pos[i];
            vx_barrier[vid + i] = fe_problem.x()[vid + i];
        }

        if(constraints_mask[local_vid] ||
           (local_vid == 3 && eq_constraints.empty())) { //constrain at least one vertex to work with the Newton solver

            //log(LogLevel::INFO)<<" - adding constraint "<<std::endl;
            for (int i(0);i < 3;i++) {
                COMISO::LinearConstraint::SVectorNC coeffs;
                coeffs.resize(n_vars);
                coeffs.coeffRef(vid + i) = 1.0;
                eq_constraints.emplace_back(COMISO::LinearConstraint(coeffs, -v_pos[i]));
            }
        }

        P_ref.push_back(vec2vec(domain_positions[local_vid]));
    }

    if(eq_constraints.empty()){
        log(LogLevel::ERROR)<<" ERROR - EMPTY CONSTRAINTS"<<std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<double> parameters = {1.0, 0.0, std::numeric_limits<double>::max()};
    /*if(json_data_["base_element_type"] == "fof"){
        parameters = {fof_epsilon, conf_w, vol_w, Dirichlet_barrier_w, Dirichlet_barrier_thresh};
    }*/
    auto setup_res = ElemT::compute_constants(parameters, P_ref, vc_barrier);


    ElemT elem;
    double barrier_e = elem.eval_f(vx_barrier, vc_barrier);

    fe_barrier.instances().add_element(vi_barrier,
                                       vc_barrier);

    fe_problem.add_set(&fe_barrier);

    std::vector<COMISO::NConstraintInterface*> constraints_pointers;
    constraints_pointers.reserve(eq_constraints.size());
    for( auto& c : eq_constraints) {
        constraints_pointers.emplace_back(&c);
    }

    double initial_energy = fe_problem.eval_f(fe_problem.x().data());
    if(initial_energy == std::numeric_limits<double>::infinity()){
        //log(LogLevel::INFO)<<" -- infinite initial energy, stopping there."<<std::endl;
        stat.final_SDE = std::numeric_limits<double>::infinity();;
        stat.converged = false;
        stat.opti_iters = -1;
    }else {
        //log(LogLevel::INFO) << " - initial SDE = " << initial_energy << std::endl;
        stat.initial_SDE = initial_energy;

        if (false) {
            bool use_truncated_newton(true);
            bool converged(false);
            int max_iters(1000);
            if (use_truncated_newton) {
                COMISO::TruncatedNewtonPCG tn(1e-8);
                tn.always_update_preconditioner() = true;
                tn.max_iters() = max_iters;
                tn.max_pcg_iters() = 100;
                tn.allow_warmstart() = false;
                tn.tolerance_gdx() = 1e-8;//previously 1e-1
                tn.set_silent(LogConfig::get_min_level() >= LogLevel::INFO);
                tn.set_eps_constraints_violation(1e-9);//1e-9
                //tn.max_preconditioner_range() = json_data_["preconditioner_range"];
                //tn.stop_if_progress_is_lower_than_thresh() = stop_if_progress_is_lower_than_thresh && !relaxation_pass;
                //tn.constraint_progress_thresh() = constraint_progress_thresh;
                tn.solve_projected_normal_equation(&fe_problem, constraints_pointers);
                stat.converged = tn.converged();
                stat.opti_iters = tn.n_iterations_used();
            } else {
                COMISO::NewtonSolver ns(1e-1, 1e-9, max_iters);
                ns.solve_infeasible_start(&fe_problem, constraints_pointers);
                converged = ns.converged();
            }

            if(!stat.converged){
                // log(LogLevel::INFO)<<" --> solver didn't converge"<<std::endl;
            }

            //solver_iters += ns.iters();


            //log(LogLevel::INFO) << " energies after optimization: " << std::endl;
            double final_energy = fe_problem.eval_f(fe_problem.x().data());

            stat.final_SDE = final_energy;
        }else {
            stat.final_SDE = -1.0;
        }
    }

}


void export_stats_to_csv(const std::vector<PerElementStats>& stats,
                         const std::string& csv_path){

        std::ofstream out_file(csv_path);
        if(!out_file.is_open()){
            log(LogLevel::INFO)<<" -> couldn't open file "<<csv_path<<std::endl;
            exit(EXIT_FAILURE);
        }

        out_file<<"parameter1;parameter2;parameter3;"
                  "n_boundary_vertices;deformation;opti_iters;converged;initial_SDE;final_SDE;ref_SDE;"
                  "min_dihedral_angle_deg;avg_dihedral_angle_deg;max_dihedral_angle_deg;"
                  "min_interior_angle_deg;avg_interior_angle_deg;max_interior_angle_deg;"
                  "volume;volume_variation;"
                  "bbox_diagonal_norm;"
                  "boundary_size_variation;height_variation"<<std::endl;
        for(auto stat: stats){
            out_file<<std::setprecision(20)<<stat.param1<<";"<<stat.param2<<";"<<stat.param3<<";"
                    <<stat.n_boundary_vertices<<";"<<stat.domain_or_codomain_deformation<<";"
                    <<stat.opti_iters<<";"<<stat.converged<<";"<<stat.initial_SDE<<";"<<stat.final_SDE<<";"<<stat.ref_SDE<<";"
                    <<stat.min_dihedral_angle_deg<<";"<<stat.avg_dihedral_angle_deg<<";"<<stat.max_dihedral_angle_deg<<";"
                    <<stat.min_interior_angle_deg<<";"<<stat.avg_interior_angle_deg<<";"<<stat.max_interior_angle_deg<<";"
                    <<stat.volume<<";"<<stat.relative_volume<<";"
                    <<stat.bounding_box_diagonal_norm<<";"
                    <<stat.boundary_size_variation<<";"<<stat.relative_height<<std::endl;
        }
        out_file.close();
};


void write_experiment_results(TetOptiExperimentFunction test_function, std::string test_name){

    std::vector<bool> mesh_to_change = {true, false};

    std::vector<PerElementStats> stats;

    for(auto option: mesh_to_change){
        test_function(option, stats);
    }

    std::string file_name = test_name + ".csv";
    export_stats_to_csv(stats, file_name);


}


void analyze_all_elements(const NonConformingTetrahedralMesh &ref_mesh,
                          const NonConformingTetrahedralMesh &mesh,
                          double ref_to_mesh_scale_factor,
                          bool use_ref_mesh_as_reference,
                          double reg_tet_target_vol,
                          VERTEX_POSITION_SOURCE pos_source,
                          std::string csv_file_path){

    log(LogLevel::INFO)<<std::endl;
    log(LogLevel::INFO)<<" +--------------------------------------------------------------------+ "<<std::endl;
    if (use_ref_mesh_as_reference) {
        log(LogLevel::INFO)<<" |--------------- analyzing all elements map distortion --------------| "<<std::endl;
    }else {
        log(LogLevel::INFO)<<" |------------------ analyzing all elements quality ------------------| "<<std::endl;
    }


    std::vector<PerElementStats> stats;

    bool constrain_boundary = true;

    double mesh_tot_vol = mesh.total_volume();
    double mesh_avg_vol = mesh.total_volume()/ref_mesh.n_logical_elements();

    double ref_mesh_tot_vol = ref_mesh.total_volume();
    double ref_mesh_avg_vol = ref_mesh.total_volume()/ref_mesh.n_logical_elements();

    double unit_to_ref_scale_factor = std::pow(ref_mesh_avg_vol, 1.0 / 3.0);

    double global_codomain_bbox_diagonal_norm = mesh.global_bounding_box_diagonal().norm();
    //log(LogLevel::INFO)<<" - global codomain mesh (boundary) bounding box diagonal: "<<global_codomain_bbox_diagonal_norm<<std::endl;

    std::vector<NonConformingTetrahedralMesh::EigenVec> reg_elem_ref;
    setup_regular_ref_element<NonConformingTetrahedralMesh>(ref_mesh_avg_vol, reg_elem_ref);

    double max_bbox_diagonal_norm(0);
    OVM::CellHandle max_bbox_diagonal_norm_elem(-1);

    IO::FileManager fm;

    MeshAngles angles;

    double max_initial_SDE(0);
    OVM::CellHandle max_initial_SDE_element(-1);
    PerElementStats max_initial_SDE_element_stats;

    double max_final_SDE(0);
    OVM::CellHandle max_final_SDE_element(-1);
    PerElementStats max_final_SDE_element_stats;

    double min_volume_length(std::numeric_limits<double>::max());
    OVM::CellHandle min_volume_length_element(-1);
    PerElementStats min_volume_length_element_stats;

    //{min dihedral angle, max dihedral angle diff, relative volume, SDE, volume-length}
    std::vector<int> min_powers = {-2, -2, -8, -1, -8};
    std::vector<int> max_powers = {1, 1, 2, 5, 0};
    std::vector<int> histo_sizes;
    for (int i(0); i<min_powers.size(); i++) {
        histo_sizes.push_back(max_powers[i] - min_powers[i] + 2);
    }

    //int histo_size(max_power - min_power + 2);
    // axis 0: either edges connected to the boundary or cells with at least one boundary edge
    // axis 1: either fully interior edges, or cells with at most one vertex connected to the boundary
    std::vector<std::vector<int>>
            min_dihedral_angle_power_histo(      histo_sizes[0], {0, 0}),
            max_dihedral_angle_diff_power_histo( histo_sizes[1], {0, 0}),
            relative_volume_power_histo(         histo_sizes[2], {0, 0}),
            SDE_power_histo(                     histo_sizes[3], {0, 0}),
            volume_length_power_histo(           histo_sizes[4], {0, 0});

    for (auto ch: mesh.cells()){

        //log(LogLevel::INFO)<<" ------------------------------------------------------------------ analyzing element "<<ch<<std::endl;

        auto verts = mesh.get_cell_vertices(ch);


        PerElementStats stat;
        stat.domain_or_codomain_deformation = false;
        stat.param1 = ch.idx();

        std::array<Vec3d, 4> codomain_positions;
        std::array<Vec3d, 4> domain_positions;
        std::array<bool, 4> constraints_mask = {false, false, false, false};

        for (int local_vid(0);local_vid < verts.size();local_vid++) {
            const auto &vh = verts[local_vid];

            if (use_ref_mesh_as_reference) {
                domain_positions[local_vid] = ref_mesh.vertex(vh);
                domain_positions[local_vid] *= ref_to_mesh_scale_factor;

                switch (pos_source) {
                    case REFERENCE: {
                codomain_positions[local_vid] = vec2vec(mesh.per_vertex_ref_pos_prop(ch)[local_vid]);
                        break;
                    }
                    case CONFORMING: {
                codomain_positions[local_vid] = mesh.vertex(vh);

                        break;
                    }
                        case NON_CONFORMING: {
                codomain_positions[local_vid] = vec2vec(mesh.per_vertex_pos_prop(ch)[local_vid]);
                        break;
                    }
                }
            }else {
                domain_positions[local_vid] = vec2vec(reg_elem_ref[local_vid]);
                codomain_positions[local_vid] = mesh.vertex(vh);
            }

            if(ref_mesh.is_boundary(vh)){
                if((codomain_positions[local_vid] - mesh.vertex(vh)).norm() > 1e-5){
                    log(LogLevel::ERROR)<<" ERROR - boundary vertex "<<vh<<" is at "<<mesh.vertex(vh)<<" in the codomain but duplicate is at "<<codomain_positions[local_vid]<<std::endl;
                    exit(EXIT_FAILURE);
                }
            }


            if(constrain_boundary && mesh.is_boundary(vh)) {
                constraints_mask[local_vid] = true;
            }

            stat.n_boundary_vertices += mesh.is_boundary(vh);
        }

        if (use_ref_mesh_as_reference) {
            std::vector<Eigen::Vector3d> P_ref;
            setup_boundary_adapted_ref_positions(ref_mesh, mesh, ch, ref_to_mesh_scale_factor, P_ref);
            for (int i(0); i < 4; i++) {
                domain_positions[i] = vec2vec(P_ref[i]);
                //ref_cell_mesh.set_vertex(VertexHandle(i), domain_positions[i]);
            }
        }

        // just one boundary vertex counts as interior in this context
        bool full_interior_element = stat.n_boundary_vertices <= 1;


        // ----------------------------------------------------------------------------------------- SDE stuff

        optimize_tet(domain_positions, codomain_positions, constraints_mask, stat);

        if(stat.initial_SDE > 1e7){
            //fm.writeFile("non_conforming_tet_" + std::to_string(ch.idx()) + ".ovm", cell_mesh);
        }
        if(stat.opti_iters < 0){
            //fm.writeFile("non_conforming_failed_tet_" + std::to_string(ch.idx()) + ".ovm", cell_mesh);
        }

        if(stat.initial_SDE > max_final_SDE){
            max_initial_SDE = stat.initial_SDE;
            max_initial_SDE_element = ch;
            max_initial_SDE_element_stats = stat;
        }

        if(stat.final_SDE > max_final_SDE){
            max_final_SDE = stat.final_SDE;
            max_final_SDE_element = ch;
            max_final_SDE_element_stats = stat;
        }

        stat.ref_SDE = cell_SDE(ref_mesh, ch, unit_to_ref_scale_factor);

        int SDE_power = std::min(max_powers[3], std::max(min_powers[3] - 1, (int)std::floor(std::log10(stat.initial_SDE))));
        /*std::cout<<" ----------------------- "<<std::endl;
        std::cout<<" -         SDE power: "<<int(std::floor(std::log10(stat.initial_SDE)))<<std::endl;
        std::cout<<" - clipped SDE power: "<<SDE_power<<std::endl;
        std::cout<<" -         min power: "<<min_powers[3]<<std::endl;
        std::cout<<" -         histo idx: "<<(SDE_power - min_powers[3] + 1)<<std::endl;*/
        SDE_power_histo[SDE_power - min_powers[3] + 1][full_interior_element]++;

        // ----------------------------------------------------------------------------------------- volume-length stuff

        auto vol_length_qv =  mesh.compute_quality_vector(ch, VOLUME_LENGTH, pos_source);
        if (vol_length_qv.size() != 1) {
            log(LogLevel::ERROR)<<" ERROR - volume length quality vector has size "<<vol_length_qv.size()<<std::endl;
            return;
        }
        stat.volume_length = vol_length_qv[0];

        int volume_length_power = std::min(max_powers[4], std::max(min_powers[4] - 1, (int)std::floor(std::log10(stat.volume_length))));
        volume_length_power_histo[volume_length_power - min_powers[4] + 1][full_interior_element]++;

        if (stat.volume_length < min_volume_length) {
            min_volume_length = stat.volume_length;
            min_volume_length_element = ch;
            min_volume_length_element_stats = stat;
        }
        /*if(full_interior_element) {
            volume_length_power_histo[volume_length_power - min_powers[4] + 1].second++;
        }else {
            volume_length_power_histo[volume_length_power - min_powers[4] + 1].first++;
        }*/

        // ----------------------------------------------------------------------------------------- angles stuff

        stat.min_dihedral_angle_deg = std::numeric_limits<double>::max();
        for(auto ce_it = mesh.ce_iter(ch); ce_it.valid(); ce_it++){
            double angle_deg = mesh.compute_dihedral_angle_deg(ch, *ce_it, pos_source);
            stat.min_dihedral_angle_deg = std::min(stat.min_dihedral_angle_deg, angle_deg);
            stat.avg_dihedral_angle_deg += angle_deg;
            stat.max_dihedral_angle_deg = std::max(stat.max_dihedral_angle_deg, angle_deg);

            auto from_v = mesh.edge(*ce_it).from_vertex();
            auto   to_v = mesh.edge(*ce_it).to_vertex();
            bool full_interior = !mesh.is_boundary(from_v) && !mesh.is_boundary(to_v);
            if(angle_deg > 90.0){
                double angle_diff = 180.0 - angle_deg;
                int angle_power = std::min(max_powers[0], std::max(min_powers[0] - 1, (int)std::floor(std::log10(angle_diff))));

                max_dihedral_angle_diff_power_histo[angle_power - min_powers[0] + 1][full_interior]++;

                /*if(full_interior) {
                    max_dihedral_angle_diff_power_histo[angle_power - min_powers[0] + 1].second++;
                }else{
                    max_dihedral_angle_diff_power_histo[angle_power - min_powers[0] + 1].first++;
                }*/
            }else{
                int angle_power = std::min(max_powers[0], std::max(min_powers[0] - 1, (int)std::floor(std::log10(angle_deg))));
                min_dihedral_angle_power_histo[angle_power - min_powers[0] + 1][full_interior]++;

                /*if(full_interior) {
                    min_dihedral_angle_power_histo[angle_power - min_powers[0] + 1].second++;
                }else{
                    min_dihedral_angle_power_histo[angle_power - min_powers[0] + 1].first++;
                }*/
            }
        }
        stat.avg_dihedral_angle_deg /= 6.0;

        for(auto cf_it = mesh.cf_iter(ch); cf_it.valid(); cf_it++){
            for(auto fv_it = mesh.fv_iter(*cf_it); fv_it.valid(); fv_it++){
                double angle_deg = mesh.compute_interior_angle_deg(*cf_it, *fv_it, pos_source, &ch);

                stat.min_interior_angle_deg = std::min(stat.min_interior_angle_deg, angle_deg);
                stat.avg_interior_angle_deg += angle_deg;
                stat.max_interior_angle_deg = std::max(stat.max_interior_angle_deg, angle_deg);
            }
        }
        stat.avg_interior_angle_deg /= 6.0;

        angles.min_dihedral_angle_deg = std::min(angles.min_dihedral_angle_deg, stat.min_dihedral_angle_deg);
        angles.max_dihedral_angle_deg = std::max(angles.max_dihedral_angle_deg, stat.max_dihedral_angle_deg);

        angles.min_interior_angle_deg = std::min(angles.min_interior_angle_deg, stat.min_interior_angle_deg);
        angles.max_interior_angle_deg = std::max(angles.max_interior_angle_deg, stat.max_interior_angle_deg);


        // ----------------------------------------------------------------------------------------- bounding box stuff

        if (pos_source == CONFORMING) {
            stat.volume = mesh.element_volume(ch);
        stat.bounding_box_diagonal_norm = mesh.element_bounding_box_diagonal(ch).norm();
        }else if (pos_source == NON_CONFORMING){
            stat.volume = mesh.non_conforming_element_volume(ch);
        stat.bounding_box_diagonal_norm = mesh.nonconforming_element_bounding_box_diagonal(ch).norm();
        }

        //log(LogLevel::INFO)<<" - codomain non-conforming bbox diagonal norm = "<<stat.bounding_box_diagonal_norm<<std::endl;
        stat.bounding_box_diagonal_norm /= global_codomain_bbox_diagonal_norm;
        //log(LogLevel::INFO)<<" - ratio to global diagonal norm = "<<stat.bounding_box_diagonal_norm<<std::endl;
        if(stat.bounding_box_diagonal_norm > 1.1){
            log(LogLevel::INFO)<<" --> found element "<<ch<<" with large bounding box diagonal: "<<stat.bounding_box_diagonal_norm<<std::endl;
        }
        if(stat.bounding_box_diagonal_norm > max_bbox_diagonal_norm){
            max_bbox_diagonal_norm = stat.bounding_box_diagonal_norm;
            max_bbox_diagonal_norm_elem = ch;
        }


        // ----------------------------------------------------------------------------------------- volume stuff

        if(stat.volume <= 0){
            log(LogLevel::INFO)<<" --> found cell with volume "<<stat.volume<<std::endl;
            //fm.writeFile("zero_vol_tet_" + std::to_string(ch.idx()) + ".ovm", cell_mesh);
            //exit(EXIT_FAILURE);
        }
        double ref_volume = (compute_volume(ref_mesh, ch) * std::pow(ref_to_mesh_scale_factor, 3.0));
        /*if (use_ref_mesh_as_reference) {
            stat.relative_volume = stat.volume / ref_volume;
        }else {
            stat.relative_volume = stat.volume / mesh_avg_vol;
        }*/
        stat.relative_volume = stat.volume / mesh_avg_vol;


        int relative_volume_power = std::min(max_powers[2], std::max(min_powers[2] - 1, (int)std::floor(std::log10(stat.relative_volume))));
        relative_volume_power_histo[relative_volume_power - min_powers[2] + 1][full_interior_element]++;

        /*if(full_interior_element) {
            relative_volume_power_histo[relative_volume_power - min_powers[2] + 1].second++;
        }else {
            relative_volume_power_histo[relative_volume_power - min_powers[2] + 1].first++;
        }*/

        // ----------------------------------------------------------------------------------------- boundary area stuff

        double ref_boundary_face_area, boundary_face_area;
        for(auto chf_it = mesh.chf_iter(ch); chf_it.valid(); chf_it++){
            if(mesh.is_boundary(mesh.opposite_halfface_handle(*chf_it))){
                ref_boundary_face_area = compute_area(ref_mesh, *chf_it) * ref_to_mesh_scale_factor * ref_to_mesh_scale_factor;
                    boundary_face_area = compute_area(    mesh, *chf_it);
            }
        }


        // ----------------------------------------------------------------------------------------- height stuff

        double ref_height =  ref_volume / ref_boundary_face_area;
        double     height = stat.volume / boundary_face_area;
        stat.relative_height = height / ref_height;


        stat.boundary_size_variation = boundary_size_variation(ref_mesh, mesh, ch, ref_to_mesh_scale_factor);
        //log(LogLevel::INFO)<<" - boundary size variation = "<<stat.boundary_size_variation<<std::endl;


        stats.push_back(stat);
    }

    log(LogLevel::INFO) << " +--------------------------------------------------------------------+ " << std::endl;
    log(LogLevel::INFO) << " | -> done with analysis. Angles stats:                               |" << std::endl;
    log(LogLevel::INFO) << " | - min dihedral angle = " << std::setprecision(5) << std::setw(10) << angles.min_dihedral_angle_deg << "                                  |" << std::endl;
    log(LogLevel::INFO) << " | - max dihedral angle = " << std::setprecision(5) << std::setw(10) << angles.max_dihedral_angle_deg << "                                  |" << std::endl;
    log(LogLevel::INFO) << " | - min interior angle = " << std::setprecision(5) << std::setw(10) << angles.min_interior_angle_deg << "                                  |" << std::endl;
    log(LogLevel::INFO) << " | - max interior angle = " << std::setprecision(5) << std::setw(10) << angles.max_interior_angle_deg << "                                  |" << std::endl;

    std::vector<std::string> labels = {
        " |     dihedral angle < 90° histogram (boundary / fully interior)     |",
        " | (180 - dihedral angle > 90°) histogram (boundary / fully interior) |",
        " |                  relative volume (to average vol)                  |",
        " |        (unweighted) SDE (either to ref mesh or regular tet)        |",
        " |                            volume-length                           |"
    };
    std::vector<std::vector<std::vector<int> > *> histograms = {
        &min_dihedral_angle_power_histo,
        &max_dihedral_angle_diff_power_histo,
        &relative_volume_power_histo,
        &SDE_power_histo,
        &volume_length_power_histo
    };

    for(int i(0); i<histograms.size(); i++) {
        log(LogLevel::INFO) << " +--------------------------------------------------------------------+ " << std::endl;
        log(LogLevel::INFO)<<labels[i]<<std::endl;
        auto& histo = *histograms[i];
        for (int j(0); j < histo.size(); j++) {
            if (!j) {
                log(LogLevel::INFO) << " |        < 1e" << (min_powers[i]) << ":"
                          << std::setw(6) << histo[j][0]
                          << " / "
                          << std::setw(6) << histo[j][1]
                          << "                                      |" << std::endl;
            } else if (j + 1 == histo.size()) {
                log(LogLevel::INFO) << " |         > 1e" << (max_powers[i]) << ":"
                          << std::setw(6) << histo[j][0]
                          << " / "
                          << std::setw(6) << histo[j][1]
                          << "                                      |" << std::endl;
            } else {
                int lower_power = min_powers[i] + j - 1;
                int upper_power = lower_power + 1;
                std::string lower_shift(lower_power >= 0 ? " " : "");
                std::string upper_shift(upper_power >= 0 ? " " : "");
                log(LogLevel::INFO) << " |  [1e" << lower_power << "," << lower_shift << (min_powers[i] >= -1 ? "" : " ") << "1e"
                          << (min_powers[i] + j) << upper_shift << "]:"
                          << std::setw(6) << histo[j][0]
                          << " / "
                          << std::setw(6) << histo[j][1]
                          << "                                      |" << std::endl;
            }
        }
    }


    // other metrics
    log(LogLevel::INFO) << " +--------------------------------------------------------------------+ " << std::endl;
    log(LogLevel::INFO) << " | -> element with largest bounding box was " << std::setw(7) << max_bbox_diagonal_norm_elem
              << "                   |" << std::endl;
    log(LogLevel::INFO) << " |    with diagonal norm = " << std::setw(6) << std::setprecision(2) << max_bbox_diagonal_norm
              << " times the mesh's bounding box       | " << std::endl;

    log(LogLevel::INFO) << " | -> element with max init. SDE was " << std::setw(6) << max_initial_SDE_element
                        << " with init. SDE = " << std::setw(8) << std::setprecision(3) << max_initial_SDE << " | " << std::endl;
    log(LogLevel::INFO) << " |    it has " << max_initial_SDE_element_stats.n_boundary_vertices
                        << " boundary vertices and its final SDE was " << std::setw(6)
                        << max_initial_SDE_element_stats.final_SDE << "        | " << std::endl;

    log(LogLevel::INFO) << " | -> element with max final SDE was " << std::setw(6) << max_final_SDE_element
              << " with final SDE = " << std::setw(8) << std::setprecision(3) << max_final_SDE << " | " << std::endl;
    log(LogLevel::INFO) << " |    it has " << max_final_SDE_element_stats.n_boundary_vertices
              << " boundary vertices and its initial SDE was " << std::setw(7)
    << max_final_SDE_element_stats.initial_SDE << "      | " << std::endl;

    log(LogLevel::INFO) << " | -> element with min volume-length was " << std::setw(6) << min_volume_length_element
              << " with VL = " << std::setw(11) << std::setprecision(3) << min_volume_length << " | " << std::endl;
    log(LogLevel::INFO) << " |    it has " << min_volume_length_element_stats.n_boundary_vertices
              << " boundary vertices                                      | " << std::endl;
    log(LogLevel::INFO) << " +--------------------------------------------------------------------+ " << std::endl;
    if(csv_file_path != "") {
        export_stats_to_csv(stats, csv_file_path);
        log(LogLevel::INFO) << " | --> exported stats to "<<std::setw(44)<<csv_file_path<<" | "<<std::endl;
    }
    log(LogLevel::INFO) << " +--------------------------------------------------------------------+ " << std::endl;
    log(LogLevel::INFO)<<std::endl;

}
