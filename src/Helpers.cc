//
// Created by Valentin Nigolian on 14.12.2024.
//

#include "Helpers.hh"

#include "NonConformingHelpers.hh"

namespace OVM = OpenVolumeMesh;


QUALITY_METRIC metric_name_to_enum(const std::string& metric_name) {
    if (metric_name == "dihedral_angle") {
        return DIHEDRAL_ANGLE;
    } else if (metric_name == "negative_dihedral_angle"){
        return NEGATIVE_DIHEDRAL_ANGLE;
    } else if(metric_name == "dihedral_angle_diff") {
        return DIHEDRAL_ANGLE_DIFF;
    }else if (metric_name == "volume_length"){
        return VOLUME_LENGTH;
    }else if (metric_name == "wicke"){
        return WICKE;
    }else {
        return UNKNOWN_METRIC;
    }
}



std::ostream& operator<<(std::ostream& os, const std::vector<OpenVolumeMesh::VertexHandle>& verts){
    for(auto v: verts){
        os<<" "<<v;
    }
    return os;
}


void split_boundary_edges_opposite_to_high_angles(TetrahedralMesh& domain_mesh,
                                                  TetrahedralMesh& codomain_mesh,
                                                  double angle_threshold_deg,
                                                  int max_iterations,
                                                  bool only_allow_safe_splits){

    log(LogLevel::INFO)<<" ------------------------------------------------------------------------------------- "<<std::endl;
    log(LogLevel::INFO)<<" splitting boundary edges to minimize maximum angle... "<<std::endl;
    if(only_allow_safe_splits){
        log(LogLevel::INFO)<<" --> splits will strictly decrease the max angle"<<std::endl;
    }else{
        log(LogLevel::INFO)<<" --> WARNING - splits might increase the max angle"<<std::endl;
    }

    //not very efficient but just to test some stuff for now
    double max_steps = 1000;


    EntityPriorityQueue<OVM::EdgeHandle, double> edges_to_split;
    int iteration_count(0);
    int total_split_count(0);
    int split_count(0);

    do{
        log(LogLevel::INFO)<<" ----------- running round "<<iteration_count<<std::endl;
        edges_to_split.clear();
        split_count = 0;

        for(auto fh: domain_mesh.faces()){
            if(domain_mesh.is_boundary(fh)){
                double max_angle;
                OVM::EdgeHandle op_eh;
                max_face_angle(domain_mesh, fh, max_angle, op_eh);

                if(max_angle > angle_threshold_deg){
                    auto from_v = domain_mesh.edge(op_eh).from_vertex();
                    auto to_v   = domain_mesh.edge(op_eh).to_vertex();

                    auto length = (domain_mesh.vertex(from_v) - domain_mesh.vertex(to_v)).norm();

                    log(LogLevel::INFO)<<" - face "<<fh<<" has max angle "<<max_angle<<" > "<<angle_threshold_deg<<", adding opposite edge ("<<from_v<<", "<<to_v<<") with length "<<length<<" to list"<<std::endl;
                    edges_to_split.push(op_eh, -length);
                }
            }
        }

        int to_split_count = edges_to_split.size();
        log(LogLevel::INFO)<<" - found "<<edges_to_split.size()<<" edges to split"<<std::endl;

        while(!edges_to_split.empty()){

            auto eh = edges_to_split.front();
            auto length = -edges_to_split.front_priority();
            edges_to_split.pop();

            if(domain_mesh.is_deleted(eh)){
                continue;
            }

            auto from_v = domain_mesh.edge(eh).from_vertex();
            auto to_v   = domain_mesh.edge(eh).to_vertex();

            log(LogLevel::INFO)<<" ----- splitting edge ("<<from_v<<", "<<to_v<<") with length "<<length<<std::endl;

            double initial_max_max_angle(0);
            if(only_allow_safe_splits) {
                for (auto ef_it = domain_mesh.ef_iter(eh); ef_it.valid(); ef_it++) {
                    if (domain_mesh.is_boundary(*ef_it)) {
                        double max_angle(0);
                        OVM::EdgeHandle op_eh;
                        max_face_angle(domain_mesh, *ef_it, max_angle, op_eh);
                        //log(LogLevel::INFO)<<" - max angle is "<<max_angle<<" opposite to edge "<<domain_mesh.edge(op_eh)<<std::endl;

                        initial_max_max_angle = std::max(initial_max_max_angle, max_angle);
                    }
                }
            }

            auto from_v_pos = domain_mesh.vertex(from_v);
            auto   to_v_pos = domain_mesh.vertex(to_v);

            auto d = to_v_pos - from_v_pos;

            auto mid_v = domain_mesh.split_edge(eh);
            codomain_mesh.split_edge(eh);

            double best_t(0);
            double min_max_angle(std::numeric_limits<double>::max());
            double prev_max_max(180.0);
            for(int i(1); i < (int)max_steps; i++){

                double t = (double) i / max_steps;

                auto mid_v_pos = from_v_pos + t * d;
                domain_mesh.set_vertex(mid_v, mid_v_pos);

                double max_max_angle(0);
                OVM::EdgeHandle max_max_angle_op_eh(-1);
                for(auto vf_it = domain_mesh.vf_iter(mid_v); vf_it.valid(); vf_it++){
                    if(domain_mesh.is_boundary(*vf_it)){
                        double max_angle;
                        OVM::EdgeHandle op_eh;
                        max_face_angle(domain_mesh, *vf_it, max_angle, op_eh);

                        if(max_angle > max_max_angle){
                            max_max_angle = max_angle;
                            max_max_angle_op_eh = op_eh;
                        }
                    }
                }

                if(max_max_angle < min_max_angle){

                    auto op_from_v = domain_mesh.edge(max_max_angle_op_eh).from_vertex();
                    auto op_to_v   = domain_mesh.edge(max_max_angle_op_eh).to_vertex();

                    min_max_angle = max_max_angle;
                    best_t = t;
                    //log(LogLevel::INFO) << " -- updated best t to " << best_t << " with min max angle = " << min_max_angle
                    //          << " (opposite to edge (" << op_from_v << ", " << op_to_v << ")" << std::endl;
                }

                if(max_max_angle > prev_max_max){
                    log(LogLevel::INFO)<<" -> angle increased, stopping there"<<std::endl;
                    break;
                }
                prev_max_max = max_max_angle;
            }

            log(LogLevel::INFO)<<" --> final best t is "<<best_t<<" giving max angle "<<min_max_angle<<std::endl;

            if(!only_allow_safe_splits || min_max_angle < initial_max_max_angle){
                auto mid_v_pos = from_v_pos + best_t * d;
                domain_mesh.set_vertex(mid_v, mid_v_pos);
                split_count++;
            }else{
                log(LogLevel::INFO)<<" new max angle is "<<min_max_angle<<" > "<<initial_max_max_angle<<", cancelling split"<<std::endl;
                auto to_collapse = domain_mesh.find_halfedge(mid_v, to_v);
                  domain_mesh.collapse_edge(to_collapse);
                codomain_mesh.collapse_edge(to_collapse);
            }
        }

        total_split_count += split_count;
        log(LogLevel::INFO)<<" -> split "<<split_count<<" / "<<to_split_count<<" edges at iteration "<<iteration_count<<", updating total to "<<total_split_count<<std::endl;


        iteration_count++;
    }while(iteration_count < max_iterations && split_count > 0);


    log(LogLevel::INFO)<<" ------------------------------------------------------------------------------------- "<<std::endl;
}


void split_edges_to_minimize_max_boundary_angle_with_terminal_edge_propagation(TetrahedralMesh& domain_mesh,
                                                                               TetrahedralMesh& codomain_mesh,
                                                                               double angle_threshold_deg,
                                                                               int max_iterations){

    log(LogLevel::INFO)<<" ------------------------------------------------------------------------------------- "<<std::endl;
    log(LogLevel::INFO)<<" splitting boundary edges to maximize minimum angle... "<<std::endl;

    //not very efficient but just to test some stuff for now
    double max_steps = 1000;


    EntityPriorityQueue<OVM::EdgeHandle, double> edges_to_split;
    int iteration_count(0);
    int total_split_count(0);

    bool split_something(true);
    //v1: find an edge, backtrack to terminal edge, split in opposite order, then repeat.
    //not very efficient but only to test hypothesis for now

    while(iteration_count < max_iterations && split_something){
        log(LogLevel::INFO)<<" =================================================================="<<std::endl;
        log(LogLevel::INFO)<<" - looking for "<<iteration_count<<"-th edge to split..."<<std::endl;

        split_something = false;
        for (auto fh: domain_mesh.faces()) {
            if (domain_mesh.is_boundary(fh)) {

                double max_angle(0);

                OVM::EdgeHandle op_eh;
                max_face_angle(domain_mesh, fh, max_angle, op_eh);

                if (max_angle > angle_threshold_deg) {
                    int split_count(0);

                    auto from_v = domain_mesh.edge(op_eh).from_vertex();
                    auto to_v = domain_mesh.edge(op_eh).to_vertex();

                    auto op_length = (domain_mesh.vertex(from_v) - domain_mesh.vertex(to_v)).norm();

                    log(LogLevel::INFO) << " - face " << fh << " has max angle " << max_angle << " > " << angle_threshold_deg
                              << ", splitting opposite edge (" << from_v << ", " << to_v << ") with length "
                              << op_length
                              << std::endl;


                    OVM::EdgeHandle terminal_edge(-1);

                    do {

                        terminal_edge = find_terminal_boundary_edge(domain_mesh, op_eh);

                        auto from_v = domain_mesh.edge(terminal_edge).from_vertex();
                        auto to_v = domain_mesh.edge(terminal_edge).to_vertex();

                        double length = edge_length(domain_mesh, terminal_edge);
                        log(LogLevel::INFO) << " ----- splitting edge (" << from_v << ", " << to_v << ") with length "
                                  << length
                                  << std::endl;

                        auto from_v_pos = domain_mesh.vertex(from_v);
                        auto to_v_pos = domain_mesh.vertex(to_v);

                        auto d = to_v_pos - from_v_pos;

                        auto mid_v = domain_mesh.split_edge(terminal_edge);
                        codomain_mesh.split_edge(terminal_edge);

                        double best_t(0);
                        double min_max_angle(std::numeric_limits<double>::max());
                        double prev_max_max(180.0);
                        for (int i(1); i <= (int) max_steps; i++) {

                            double t = (double) i / (max_steps + 1);

                            auto mid_v_pos = from_v_pos + t * d;
                            domain_mesh.set_vertex(mid_v, mid_v_pos);

                            double max_max_angle(0);
                            OVM::EdgeHandle max_max_angle_op_eh(-1);
                            for (auto vf_it = domain_mesh.vf_iter(mid_v); vf_it.valid(); vf_it++) {
                                if (domain_mesh.is_boundary(*vf_it)) {
                                    double max_angle;
                                    OVM::EdgeHandle op_eh;
                                    max_face_angle(domain_mesh, *vf_it, max_angle, op_eh);

                                    if (max_angle > max_max_angle) {
                                        max_max_angle = max_angle;
                                        max_max_angle_op_eh = op_eh;
                                    }
                                }
                            }

                            if (max_max_angle < min_max_angle) {

                                auto op_from_v = domain_mesh.edge(max_max_angle_op_eh).from_vertex();
                                auto op_to_v = domain_mesh.edge(max_max_angle_op_eh).to_vertex();

                                min_max_angle = max_max_angle;
                                best_t = t;
                                //log(LogLevel::INFO) << " -- updated best t to " << best_t << " with min max angle = " << min_max_angle
                                //          << " (opposite to edge (" << op_from_v << ", " << op_to_v << ")" << std::endl;
                            }

                            if (max_max_angle > prev_max_max) {
                                log(LogLevel::INFO) << " -> angle increased, stopping there" << std::endl;
                                break;
                            }
                            prev_max_max = max_max_angle;
                        }

                        log(LogLevel::INFO) << " --> final best t is " << best_t << " giving max angle " << min_max_angle
                                  << std::endl;
                        auto mid_v_pos = from_v_pos + best_t * d;
                        domain_mesh.set_vertex(mid_v, mid_v_pos);

                        log(LogLevel::INFO) << " --> moved mid-vertex " << mid_v << " to " << mid_v_pos << std::endl;

                        split_count++;
                        split_something = true;

                    } while (terminal_edge != op_eh);

                    total_split_count += split_count;
                    log(LogLevel::INFO) << " -> split " << split_count << " edges at iteration " << iteration_count
                              << ", updating total to " << total_split_count << std::endl;

                    //log(LogLevel::INFO)<<" -> edge "<<domain_mesh.edge(op_eh)<<" is now split"<<std::endl;


                    iteration_count++;
                    break;
                }
            }
        }
    }

    log(LogLevel::INFO)<<" ------------------------------------------------------------------------------------- "<<std::endl;
}


OpenVolumeMesh::EdgeHandle find_terminal_boundary_edge(TetrahedralMesh& mesh,
                                                       const OpenVolumeMesh::EdgeHandle& eh){

    OVM::EdgeHandle terminal_edge(-1);

    if(!mesh.is_boundary(eh)){
        log(LogLevel::INFO)<<" WARNING - asking for terminal boundary edge on a non-boundary edge"<<std::endl;
        return terminal_edge;
    }

    //log(LogLevel::INFO)<<" --------------------- looking for terminal edges starting from "<<mesh.edge(eh)<<std::endl;

    auto candidate_prop = mesh.request_edge_property<bool>("", false);
    std::queue<OVM::EdgeHandle> candidates;
    candidate_prop[eh] = true;
    candidates.push(eh);


    while(!candidates.empty()){
        auto candidate = candidates.front();
        candidates.pop();

        double current_len = edge_length(mesh, candidate);
        //log(LogLevel::INFO)<<" ---------- checking candidate "<<mesh.edge(candidate)<<" with length "<<current_len<<std::endl;

        bool is_terminal(true);

        for(auto ef_it = mesh.ef_iter(candidate); ef_it.valid(); ef_it++){
            if(mesh.is_boundary(*ef_it)){
                double max_len(0);
                OVM::EdgeHandle max_len_eh(-1);

                for(auto fe_it = mesh.fe_iter(*ef_it); fe_it.valid(); fe_it++){
                    if(!candidate_prop[*fe_it]){
                        auto len = edge_length(mesh, *fe_it);
                        //log(LogLevel::INFO)<<" - testing neighbor edge "<<mesh.edge(*fe_it)<<" with len "<<len<<std::endl;
                        if(len > current_len && len > max_len){
                            max_len = len;
                            max_len_eh = *fe_it;
                        }
                    }
                }

                if(max_len_eh.is_valid()) {
                    //log(LogLevel::INFO)<<" --> found neighbor edge "<<mesh.edge(max_len_eh)<<" with max len "<<max_len<<" --> candidate "<<mesh.edge(candidate)<<" is not terminal"<<std::endl;
                    candidate_prop[max_len_eh] = true;
                    candidates.push(max_len_eh);
                    is_terminal = false;
                }
            }
        }

        if(is_terminal){
            //log(LogLevel::INFO)<<" --> edge "<<mesh.edge(candidate)<<" is terminal"<<std::endl;
            terminal_edge = candidate;
        }
    }

    return terminal_edge;
}

double edge_length(const TetrahedralMesh& mesh, const OVM::EdgeHandle& eh){
    return (mesh.vertex(mesh.edge(eh).from_vertex()) - mesh.vertex(mesh.edge(eh).to_vertex())).norm();
}


double compute_average_incident_edge_length(const TetrahedralMesh &mesh,
                                            const OpenVolumeMesh::VertexHandle &vh,
                                            bool only_boundary) {
    double avg(0.0);
    int val(0);
    for (auto out_he: mesh.outgoing_halfedges(vh)) {
        avg += ((!only_boundary || mesh.is_boundary(out_he)) ? edge_length(mesh, mesh.edge_handle(out_he)) : 0.0);
        val++;
    }

    return avg / val;
}

void max_face_angle(TetrahedralMesh& mesh,
                    OpenVolumeMesh::FaceHandle fh,
                    double& max_angle,
                    OpenVolumeMesh::EdgeHandle& op_eh){

    auto verts = mesh.get_halfface_vertices(mesh.halfface_handle(fh, 0));

    max_angle = 0;
    for(int i(0); i<3; i++){

        int idx_m1 = (3 + i-1) % 3;
        int idx = i;
        int idx_p1 = (i+1) % 3;

        auto vm1 = verts[idx_m1];
        auto v   = verts[idx];
        auto vp1 = verts[idx_p1];

        auto pm1 = mesh.vertex(vm1);
        auto p   = mesh.vertex(v);
        auto pp1 = mesh.vertex(vp1);

        auto em1 = (pm1 - p).normalized();
        auto ep1 = (pp1 - p).normalized();

        double angle = (180.0 / M_PI) * std::acos(em1.dot(ep1));

        if(angle > max_angle){
            max_angle = angle;
            op_eh = mesh.edge_handle(mesh.find_halfedge(vm1, vp1));
        }
    }

}


void invert_mesh(TetrahedralMesh& mesh,
                 bool move_interior_to_centroid){
    log(LogLevel::INFO)<<" ------------------------------------------ inverting mesh... "<<std::endl;

    Vec3d centroid = {0,0,0};

    int b_count(0);
    for(auto v: mesh.vertices()){
        if(mesh.is_boundary(v)) {
            centroid += mesh.vertex(v);
            b_count++;
        }
    }
    centroid /= b_count;
    log(LogLevel::INFO)<<" centroid at "<<centroid<<std::endl;

    for(auto v: mesh.vertices()){
        if(move_interior_to_centroid && !mesh.is_boundary(v)){
            mesh.set_vertex(v, centroid);
        }
        mesh.set_vertex(v, 2 * centroid - mesh.vertex(v));
    }

    log(LogLevel::INFO)<<" ------------------------------------------ "<<std::endl;
}

int invert_problem(TetrahedralMesh& domain_mesh,
                   TetrahedralMesh& codomain_mesh){

    log(LogLevel::INFO)<<" ------------------------------------------ inverting problem... "<<std::endl;

    TetrahedralMesh new_domain_mesh, new_codomain_mesh;

    auto old_to_new_prop = codomain_mesh.request_vertex_property<VertexHandle>();
    auto new_to_old_prop = new_domain_mesh.request_vertex_property<VertexHandle>();

    std::vector<Eigen::Vector3d> centroids;
    std::vector<Eigen::Vector3d> normals;

    for(auto hfh: codomain_mesh.halffaces()){
        if(codomain_mesh.is_boundary(hfh)){
            Vec3d c,n;
            compute_centroid_and_normal(codomain_mesh, hfh, c, n);

            centroids.push_back(vec2vec(c));
            normals.push_back(vec2vec(n));
        }
    }
    log(LogLevel::INFO)<<" - setup Chebyshev input for "<<centroids.size()<<" boundary halffaces"<<std::endl;

    double radius(0);
    auto cheb_center = solve_Chebyshev(normals, centroids, radius, false);

    if(radius <= 0.0){
        std::cerr<<" ERROR - codomain is not star-shaped"<<std::endl;
        return 1;
    }
    log(LogLevel::INFO)<<" - Chebyshev centroid at "<<cheb_center.transpose()<<std::endl;


    auto int_vh = new_domain_mesh.add_vertex(vec2vec(cheb_center));

    for(auto vh: codomain_mesh.vertices()){
        if(codomain_mesh.is_boundary(vh)) {
            auto new_vh = new_domain_mesh.add_vertex(codomain_mesh.vertex(vh));

            new_to_old_prop[new_vh] = vh;
            old_to_new_prop[vh] = new_vh;
        }
    }
    log(LogLevel::INFO)<<" - added "<<new_domain_mesh.n_vertices()<<" boundary vertices to the new domain mesh"<<std::endl;


    for(auto hfh: codomain_mesh.halffaces()){
        if(codomain_mesh.is_boundary(hfh)){

            auto op_hfh = codomain_mesh.opposite_halfface_handle(hfh);
            auto verts = codomain_mesh.halfface_vertices(op_hfh);

            std::vector<VertexHandle> new_verts;
            for(auto vh: verts){
                new_verts.push_back(old_to_new_prop[vh]);
            }

            new_verts.push_back(int_vh);

            new_domain_mesh.add_cell(new_verts);
        }
    }
    log(LogLevel::INFO)<<" - added "<<new_domain_mesh.n_cells()<<" cells to the new domain mesh"<<std::endl;

    if(new_domain_mesh.n_cells() != centroids.size()){
        log(LogLevel::ERROR)<<" ERROR - not as many cellls in new domain mesh as there are boundary faces in the input codomain mesh "<<std::endl;
        return -1;
    }

    //copy the new domain mesh
    new_codomain_mesh = new_domain_mesh;

    //and copy the domain's boundary positions to the new codomain
    for(auto vh: domain_mesh.vertices()){
        if(domain_mesh.is_boundary(vh)){
            auto new_vh = old_to_new_prop[vh];

            new_codomain_mesh.set_vertex(new_vh, domain_mesh.vertex(vh));
        }
    }

    //and finally swap the domain and codomain
    domain_mesh   = new_domain_mesh;
    codomain_mesh = new_codomain_mesh;

    log(LogLevel::INFO)<<" ------------------------------------------ "<<std::endl;

    return 0;
}


int translate_and_rotate_to_target(TetrahedralMesh& mesh,
                                   const OpenVolumeMesh::HalfFaceHandle& hfh,
                                   const Vec3d& target_p0,
                                   const Vec3d& target_n,
                                   const Vec3d& target_e1) {

    using EigenVec = Eigen::Vector3d;

    Vec3d tmp_source_c, source_n;
    compute_centroid_and_normal(mesh, hfh, tmp_source_c, source_n);


    auto verts = mesh.get_halfface_vertices(hfh);

    Vec3d source_p0 = mesh.vertex(verts[0]);
    Vec3d source_e1 = (mesh.vertex(verts[1]) - mesh.vertex(verts[0])).normalized();

    log(LogLevel::INFO)<<" - source anchor: "<<source_p0<<", normal: "<<source_n<<" and edge vector: "<<source_e1<<std::endl;
    log(LogLevel::INFO)<<" - target anchor: "<<target_p0<<", normal: "<<target_n<<" and edge vector: "<<target_e1<<std::endl;


    auto n_diff = (target_n - source_n).norm();
    log(LogLevel::INFO)<<" - normals diff norm = "<<n_diff<<std::endl;
    if (n_diff > DBL_EPSILON) {


        Vec3d rot_axis;
        double theta;
        if (source_n.cross(target_n).norm() > DBL_EPSILON) {
            rot_axis = source_n.cross(target_n).normalized();
            theta = std::acos(std::min(std::max(source_n.dot(target_n),-1.0), 1.0));

        }else {
            log(LogLevel::INFO)<<" - domain and codomain normals are nearly identical, but opposite. Rotating around the edge vector instead"<<std::endl;

            rot_axis = source_e1;
            theta = M_PI;
        }



        log(LogLevel::INFO)<<" - rotation axis: "<<rot_axis<<std::endl;
        log(LogLevel::INFO)<<" -         theta: "<<theta<<std::endl;

        Eigen::Matrix3d rot_mat = Eigen::AngleAxisd(theta, vec2vec(rot_axis)).toRotationMatrix();

        for (auto vh: mesh.vertices()) {
            auto source_p = mesh.vertex(vh);
            auto rel_p = source_p - source_p0;
            EigenVec p_rot = rot_mat * vec2vec(rel_p);
            mesh.set_vertex(vh, vec2vec(p_rot) + target_p0);
        }
    }else {
        log(LogLevel::INFO)<<" - source and target normals are nearly identical, only applying translation (before rotation around normal)"<<std::endl;

        for (auto vh: mesh.vertices()) {
            mesh.set_vertex(vh, mesh.vertex(vh) + target_p0 - source_p0);
        }
    }


    Vec3d tmp_updated_source_c, updated_source_n;
    compute_centroid_and_normal(mesh, hfh, tmp_updated_source_c, updated_source_n);


    if ((updated_source_n - target_n).norm() > 1e-7) {
        log(LogLevel::ERROR)<<" - source and target normals are not aligned properly. Diff between normals: "<<(updated_source_n - target_n).norm()<<std::endl;
        return -1;
    }

    Vec3d current_e1 = (mesh.vertex(verts[1]) - mesh.vertex(verts[0])).normalized();

    auto e1_diff = (target_e1 - current_e1).norm();
    log(LogLevel::INFO)<<" - edge vec diff norm = "<<n_diff<<std::endl;
    if (e1_diff > DBL_EPSILON) {
        double dot_prod = current_e1.dot(target_e1);
        double theta2 = std::acos(std::min(1.0, std::max(-1.0, dot_prod)));
        Vec3d rot_axis2 = current_e1.cross(target_e1).normalized();
        //log(LogLevel::INFO)<<" - dot prod: "<<std::setprecision(20)<<dot_prod<<std::endl;
        log(LogLevel::INFO)<<" - rotation axis 2: "<<rot_axis2<<std::endl;
        log(LogLevel::INFO)<<" -     face normal: "<<target_n<<std::endl;
        log(LogLevel::INFO)<<" -         theta 2: "<<theta2<<std::endl;

        Eigen::Matrix3d rot_mat2 = Eigen::AngleAxisd(theta2, vec2vec(rot_axis2)).toRotationMatrix();

        for (auto vh: mesh.vertices()) {
            auto current_p = mesh.vertex(vh);
            auto rel_p = current_p - target_p0;
            EigenVec p_rot = rot_mat2 * vec2vec(rel_p);
            mesh.set_vertex(vh, vec2vec(p_rot) + target_p0);
        }
    }else {
        log(LogLevel::INFO)<<" - current and target edge vectors are nearly identical, done."<<std::endl;
    }

    Vec3d updated_e1 = (mesh.vertex(verts[1]) - mesh.vertex(verts[0])).normalized();
    if ((updated_e1 - target_e1).norm() > 1e-7) {
        log(LogLevel::ERROR)<<" - target edge vector is not aligned with updated one. Diff = "<<(updated_e1 - target_e1).norm()<<std::endl;
        return -1;
    }
    return 0;
}


void translate_to_pos(TetrahedralMesh &mesh,
                      const OpenVolumeMesh::VertexHandle &anchor_vh,
                      const Vec3d& target_pos){

    auto anchor_p = mesh.vertex(anchor_vh);
    auto t = target_pos - anchor_p;

    for (auto vh: mesh.vertices()) {
        mesh.set_vertex(vh, mesh.vertex(vh) + t);
    }
}

void scale_to_volume(TetrahedralMesh& mesh,
                    const OpenVolumeMesh::VertexHandle& anchor_vh,
                    const double& target_vol){

    auto current_vol = compute_total_volume(mesh);

    double scale_factor = std::pow(target_vol / current_vol, 1.0/3.0);

    auto anchor_p = mesh.vertex(anchor_vh);
    auto scaled_anchor_p = anchor_p * scale_factor;
    auto t = scaled_anchor_p - anchor_p;

    for(auto vh: mesh.vertices()){
        auto p = mesh.vertex(vh);

        mesh.set_vertex(vh, scale_factor * p - t);
    }
}


void scale_to_boundary_area(TetrahedralMesh &mesh,
                            const OpenVolumeMesh::VertexHandle &anchor_vh,
                            const double& target_boundary_area) {

    auto current_boundary_area = compute_boundary_surface_area(mesh);
    auto current_volume = compute_total_volume(mesh);

    scale_to_volume(mesh, anchor_vh, current_volume * std::pow(target_boundary_area / current_boundary_area, 3.0 / 2.0));
}



void print_boundary_quality_stats(const TetrahedralMesh& mesh){


    log(LogLevel::INFO)<<" ---------------------------------------------------------- "<<std::endl;
    log(LogLevel::INFO)<<" Checking quality of boundary triangulation..."<<std::endl;


    double min_len_ratio(std::numeric_limits<double>::max()), avg_len_ratio(0), max_len_ratio(0);
    double global_min_len(std::numeric_limits<double>::max()), global_max_len(0);
    double global_min_len_ratio(0);
    double min_angle_ratio(std::numeric_limits<double>::max()), avg_angle_ratio(0), max_angle_ratio(0);
    double global_min_angle_rad(std::numeric_limits<double>::max()),  global_avg_angle_rad(0),  global_max_angle_rad(0);
    double min_area(std::numeric_limits<double>::max()), avg_area(0), max_area(0);
    int zero_area_count(0), zero_len_count(0), zero_angle_count(0);
    double min_height(std::numeric_limits<double>::max()), avg_height(0), max_height(0);


    int boundary_face_count(0);
    //min-max-avg length, min-max-avg angle, min-max-avg area
    for(auto hf: mesh.halffaces()){
        if(mesh.is_boundary(hf)){
            auto f_vertices = mesh.get_halfface_vertices(hf);

            auto p0 = mesh.vertex(f_vertices[0]);
            auto p1 = mesh.vertex(f_vertices[1]);
            auto p2 = mesh.vertex(f_vertices[2]);

            Vec3d e1 = (p1 - p0);
            Vec3d e2 = (p2 - p0);
            Vec3d e12 = (p2 - p1);

            auto len0 =  e1.norm();
            auto len1 =  e2.norm();
            auto len2 =  e12.norm();

            if(len0 == 0 || len1 == 0 || len2 == 0){
                zero_len_count++;
            }

            auto min_len = std::min(len0, std::min(len1, len2));
            auto max_len = std::max(len0, std::max(len1, len2));
            auto len_ratio = min_len / max_len;

            global_min_len = std::min(global_min_len, min_len);
            global_max_len = std::max(global_max_len, max_len);

            min_len_ratio = std::min(min_len_ratio, len_ratio);
            max_len_ratio = std::max(max_len_ratio, len_ratio);

            avg_len_ratio += len_ratio;


            auto cross_prod = e1.cross(e2);
            auto area = cross_prod.norm();

            min_area = std::min(min_area, area);
            max_area = std::max(max_area, area);
            avg_area += area;

            if(area == 0.0){
                zero_area_count++;
            }

            auto h0 = 2 * area / len0;
            auto h1 = 2 * area / len1;
            auto h2 = 2 * area / len2;

            auto local_min_height = std::min(h0, std::min(h1, h2));
            auto local_max_height = std::max(h0, std::max(h1, h2));
            if(local_min_height < min_height){
                //log(LogLevel::INFO)<<" - updated min height to "<<local_min_height)<<" from face "<<hf<<":"<<mesh.get_halfface_vertices(hf)<<std::endl;;
            }

            min_height = std::min(min_height, local_min_height);
            avg_height += (h0 + h1 + h2);
            max_height = std::max(max_height, local_max_height);


            auto e1_dot_e2 = e1.dot(e2) / (len0 * len1);

            /*auto alpha0 = CGAL::approximate_angle(cgal_p2, cgal_p0, cgal_p1);
            auto alpha1 = CGAL::approximate_angle(cgal_p0, cgal_p1, cgal_p2);
            auto alpha2 = CGAL::approximate_angle(cgal_p1, cgal_p2, cgal_p0);*/
            double alpha0_rad = std::acos(e1.dot(e2) / (len0 * len1));
            double alpha1_rad = std::acos(-e1.dot(e12) / (len0 * len2));
            double alpha2_rad = std::acos(e2.dot(e12) / (len1 * len2));

            auto min_angle_rad = std::min(alpha0_rad, std::min(alpha1_rad, alpha2_rad));
            auto max_angle_rad = std::max(alpha0_rad, std::max(alpha1_rad, alpha2_rad));

            if(min_angle_rad == 0){
                zero_angle_count++;
            }
            auto angle_ratio = min_angle_rad / max_angle_rad;
            if(angle_ratio < min_angle_ratio){
                // log(LogLevel::INFO)<<" - updated min angle ratio to "<<min_angle<<" / "<<max_angle<<" = "<<angle_ratio<<" from face "<<hf<<":"<<mesh.get_halfface_vertices(hf)<<std::endl;;
            }
            min_angle_ratio = std::min(min_angle_ratio, angle_ratio);
            max_angle_ratio = std::max(max_angle_ratio, angle_ratio);
            avg_angle_ratio += angle_ratio;

            global_min_angle_rad = std::min(global_min_angle_rad, min_angle_rad);
            global_avg_angle_rad += (alpha0_rad + alpha1_rad + alpha2_rad);
            global_max_angle_rad = std::max(global_max_angle_rad, max_angle_rad);

            boundary_face_count++;
        }
    }


    avg_area             /= boundary_face_count;
    avg_len_ratio        /= (boundary_face_count * 3);
    avg_angle_ratio      /= (boundary_face_count * 3);
    global_avg_angle_rad /= (boundary_face_count * 3);
    avg_height           /= (boundary_face_count * 3);

    global_min_len_ratio = global_min_len / global_max_len;

    /*stats = {min_area), avg_area), max_area),
             min_len_ratio), avg_len_ratio), max_len_ratio),
             global_min_len_ratio),
             min_angle_ratio, avg_angle_ratio, max_angle_ratio,
             global_min_angle, global_avg_angle, global_max_angle,
             zero_area_count, zero_len_count, zero_angle_count,
             min_height), avg_height), max_height)};*/

    double global_min_angle_deg = global_min_angle_rad * 180.0 / M_PI;
    double global_avg_angle_deg = global_avg_angle_rad * 180.0 / M_PI;
    double global_max_angle_deg = global_max_angle_rad * 180.0 / M_PI;

    log(LogLevel::INFO)<<"                      min    |      avg      |       max"<<std::endl;
    log(LogLevel::INFO)<<"               area "<<std::setw(13)<<CGAL::to_double(min_area)<<" | "<<std::setw(13)<<CGAL::to_double(avg_area)<<" | "<<std::setw(10)<<CGAL::to_double(max_area)<<std::endl;
    log(LogLevel::INFO)<<"       length ratio "<<std::setw(13)<<CGAL::to_double(min_len_ratio)<<" | "<<std::setw(13)<<CGAL::to_double(avg_len_ratio)<<" | "<<std::setw(10)<<CGAL::to_double(max_len_ratio)<<std::endl;
    log(LogLevel::INFO)<<"        angle ratio "<<std::setw(13)<<CGAL::to_double(min_angle_ratio)<<" | "<<std::setw(13)<<avg_angle_ratio<<" | "<<std::setw(10)<<max_angle_ratio<<std::endl;
    log(LogLevel::INFO)<<"  global angle(deg) "<<std::setw(13)<<CGAL::to_double(global_min_angle_deg)<<" | "<<std::setw(13)<<global_avg_angle_deg<<" | "<<std::setw(10)<<global_max_angle_deg<<std::endl;
    log(LogLevel::INFO)<<"             height "<<std::setw(13)<<CGAL::to_double(min_height)<<" | "<<std::setw(13)<<CGAL::to_double(avg_height)<<" | "<<std::setw(10)<<CGAL::to_double(max_height)<<std::endl;
    log(LogLevel::INFO)<<" ---------------"<<std::endl;
    log(LogLevel::INFO)<<" global min length ratio: "<<(global_min_len / global_max_len)<<std::endl;
    log(LogLevel::INFO)<<" global min  angle ratio: "<<(global_min_angle_deg / global_max_angle_deg)<<std::endl;
    log(LogLevel::INFO)<<" ---------------------------------------------------------- "<<std::endl;
}


int count_fully_boundary_elements(const TetrahedralMesh& mesh){
    int count(0);
    for(auto c: mesh.cells()){
        int bv_count(0);
        for(auto v: mesh.get_cell_vertices(c)){
            bv_count += mesh.is_boundary(v);
        }

        count += bv_count == 4;
    }

    return count;
}


int split_fully_boundary_elements(TetrahedralMesh& mesh){

    log(LogLevel::INFO)<<" ============================================================ "<<std::endl;
    log(LogLevel::INFO)<<" splitting boundary-only elements..."<<std::endl;

    int split_count(0), edge_split_count(0), face_split_count(0);

    std::vector<OVM::EdgeHandle> edges_to_split;
    std::vector<OVM::FaceHandle> faces_to_split;

    for(auto c: mesh.cells()){
        int bv_count(0);
        for(auto v: mesh.get_cell_vertices(c)){
            bv_count += mesh.is_boundary(v);
        }

        if(bv_count == 4){
            split_count++;

            OVM::EdgeHandle to_split(-1);
            for(auto ce_it = mesh.ce_iter(c); ce_it.valid(); ce_it++){
                if(!mesh.is_boundary(*ce_it)){
                    to_split = *ce_it;
                }
            }

            if(to_split.is_valid()){
                edges_to_split.push_back(to_split);
            }else{
                OVM::FaceHandle f_to_split;
                for(auto cf_it = mesh.cf_iter(c); cf_it.valid(); cf_it++){
                    if(!mesh.is_boundary(*cf_it)){
                        f_to_split = *cf_it;
                    }
                }
                if(!f_to_split.is_valid()) {
                    log(LogLevel::ERROR) << " ERROR - Tet " << c << " is boundary only but has no interior face" << std::endl;
                    exit(EXIT_FAILURE);
                }
                faces_to_split.push_back(f_to_split);
            }
        }
    }

    for(auto to_split: faces_to_split){
        if(!mesh.is_deleted(to_split)){
            mesh.split_face(to_split);
            face_split_count++;
        }
    }

    for(auto to_split: edges_to_split){
        if(!mesh.is_deleted(to_split)){
            mesh.split_edge(to_split);
            edge_split_count++;
        }
    }

    log(LogLevel::INFO)<<" - split "<<split_count<<" boundary-only elements, by doing "<<std::endl;
    log(LogLevel::INFO)<<"     "<<edge_split_count<<" edge splits and "<<std::endl;
    log(LogLevel::INFO)<<"     "<<face_split_count<<" face splits"<<std::endl;

    log(LogLevel::INFO)<<" - updated mesh size (V, E, F, C): ("<<mesh.n_logical_vertices()<<", "<<mesh.n_logical_edges()<<", "<<mesh.n_logical_faces()<<", "<<mesh.n_logical_cells()<<")"<<std::endl;

    return split_count;
}


int split_cord_edges_and_faces(TetrahedralMesh& mesh){

    log(LogLevel::INFO)<<" ============================================================ "<<std::endl;
    log(LogLevel::INFO)<<" splitting chord edge and faces..."<<std::endl;

    int split_count(0), edge_split_count(0), face_split_count(0);

    std::vector<OVM::EdgeHandle> edges_to_split;
    std::vector<OVM::FaceHandle> faces_to_split;

    for(auto e: mesh.edges()){
        if(!mesh.is_boundary(e) &&
        mesh.is_boundary(mesh.edge(e).from_vertex()) &&
        mesh.is_boundary(mesh.edge(e).to_vertex())) {
            edges_to_split.push_back(e);
        }
    }


    for(auto to_split: edges_to_split){
        if(!mesh.is_deleted(to_split)){
            mesh.split_edge(to_split);
            edge_split_count++;
        }
    }


    for(auto f: mesh.faces()){
        if(!mesh.is_boundary(f)){
            bool all_boundary(true);
            for(auto v: mesh.get_halfface_vertices(mesh.halfface_handle(f, 0))){
                if(!mesh.is_boundary(v)){
                    all_boundary = false;
                    break;
                }
            }
            if(all_boundary){
                faces_to_split.push_back(f);
            }
        }
    }
    for(auto to_split: faces_to_split){
        if(!mesh.is_deleted(to_split)){
            mesh.split_face(to_split);
            face_split_count++;
        }
    }

    log(LogLevel::INFO)<<" - split "<<split_count<<" boundary-only elements, by doing "<<std::endl;
    log(LogLevel::INFO)<<"     "<<edge_split_count<<" edge splits and "<<std::endl;
    log(LogLevel::INFO)<<"     "<<face_split_count<<" face splits"<<std::endl;

    log(LogLevel::INFO)<<" - updated mesh size (V, E, F, C): ("<<mesh.n_logical_vertices()<<", "<<mesh.n_logical_edges()<<", "<<mesh.n_logical_faces()<<", "<<mesh.n_logical_cells()<<")"<<std::endl;

    return split_count;


}


bool collapse_would_create_cord(const TetrahedralMesh &mesh,
                                const OpenVolumeMesh::EdgeHandle &eh) {

    if (mesh.is_boundary(eh)) {
        return false;
    }

    auto from_v = mesh.edge(eh).from_vertex();
    auto to_v = mesh.edge(eh).to_vertex();

    if (!mesh.is_boundary(from_v) && !mesh.is_boundary(to_v)) {
        return false;
    }

    if (mesh.is_boundary(from_v) && mesh.is_boundary(to_v)) {
        std::cerr<<" ERROR - edge is already a chord"<<std::endl;
        exit(EXIT_FAILURE);
    }

    auto boundary_v = mesh.is_boundary(from_v) ? from_v : to_v;
    auto interior_v = mesh.is_boundary(from_v) ? to_v : from_v;

    //std::cerr<<" - boundary v: "<<boundary_v<<std::endl;
    //std::cerr<<" - interior v: "<<interior_v<<std::endl;

    std::unordered_set<OVM::VertexHandle> b_v_neighbors;
    for (auto vv_it = mesh.vv_iter(boundary_v); vv_it.valid(); vv_it++) {
        if (mesh.is_boundary(*vv_it)) {
            b_v_neighbors.insert(*vv_it);
            //std::cerr<<" -- boundary neighbor to boundary vertex: "<<(*vv_it)<<std::endl;
        }
    }

    bool found_cord(false);
    for (auto vv_it = mesh.vv_iter(interior_v); vv_it.valid(); vv_it++) {
        //std::cerr<<" -- checking neighbor "<<(*vv_it)<<", is boundary: "<<mesh.is_boundary(*vv_it)<<std::endl;
        if (*vv_it != boundary_v && mesh.is_boundary(*vv_it) && !b_v_neighbors.contains(*vv_it)) {
            //std::cerr<<" -- boundary neighbor "<<(*vv_it)<<" to interior vertex is not a neighbor of the boundary vertex "<<std::endl;
            found_cord = true;
            //break;
        }
    }

    return found_cord;
}


double halfface_SDE(const TetrahedralMesh& domain_mesh,
                    const TetrahedralMesh& codomain_mesh,
                    const OpenVolumeMesh::HalfFaceHandle& hfh,
                    const double ref_scaling){

    auto hf_verts = domain_mesh.get_halfface_vertices(hfh);
    //log(LogLevel::INFO)<<" ---- checking boundary face "<<hf_verts[0]<<", "<<hf_verts[1]<<", "<<hf_verts[2]<<std::endl;

    using ElemT = SDE_TRI_3D_AD;
    ElemT elem;

    ElemT::VecV vx;
    ElemT::VecC vc;

    for(int i(0); i<hf_verts.size(); i++){
        auto p     = codomain_mesh.vertex(hf_verts[i]);
        auto ref_p = ref_scaling * domain_mesh.vertex(hf_verts[i]);

        for(int j(0); j<p.size(); j++){
            vx[0 + 3 * i + j] = p[j];
            vc[1 + 3 * i + j] = ref_p[j];
        }
    }

    vc[0] = 1.0;

    double SDE = elem.eval_f(vx, vc);

    return SDE;
}



double cell_conformal_distortion(const TetrahedralMesh& mesh,
                                 const TetrahedralMesh& ref_mesh,
                                 const OpenVolumeMesh::CellHandle& ch){

    auto c_verts = mesh.get_cell_vertices(ch);
    //log(LogLevel::INFO)<<" ---- checking boundary face "<<hf_verts[0]<<", "<<hf_verts[1]<<", "<<hf_verts[2]<<std::endl;

    using ElemT = FOF_NO_EPS_3D_AD;
    ElemT elem;

    ElemT::VecV vx;
    ElemT::VecC vc;
    std::vector<Eigen::Vector3d> P_ref(4);

    for(int i(0); i < 4; i++){
        P_ref[i] = vec2vec(ref_mesh.vertex(c_verts[i]));
    }

    //const static int NC = 11; // [conformal w, volumetric w, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z]
    for(int i(0); i<c_verts.size(); i++){
        auto p = mesh.vertex(c_verts[i]);
        auto ref_p = P_ref[i];

        for(int j(0); j < p.size(); j++){
            vx[0 + 3 * i + j] = p[j];
        }
    }

    ElemT::compute_constants({1,0}, P_ref, vc);

    double conf_dist = elem.eval_f(vx, vc);

    return conf_dist;
}



double cell_SDE(const TetrahedralMesh& mesh,
                const TetrahedralMesh& ref_mesh,
                const OpenVolumeMesh::CellHandle& ch,
                const double ref_scaling){

    auto c_verts = mesh.get_cell_vertices(ch);

    using ElemT = SD_3D_AD;
    ElemT elem;

    ElemT::VecV vx;
    ElemT::VecC vc;
    std::vector<Eigen::Vector3d> P_ref(4);
    for(int i(0); i < 4; i++){
        P_ref[i] = vec2vec(ref_mesh.vertex(c_verts[i]));
    }

    for(int i(0); i<P_ref.size(); i++){
        P_ref[i] *= ref_scaling;
    }

    //const static int NC = 12; // [weight, size_scale, epsilon, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z]
    for(int i(0); i<c_verts.size(); i++){
        auto p = mesh.vertex(c_verts[i]);
        auto ref_p = P_ref[i];

        ElemT::compute_constants({1.0, 1.0, 0.0}, P_ref, vc);

        for(int j(0); j < p.size(); j++){
            vx[0 + 3 * i + j] = p[j];
        }
    }

    double SDE = elem.eval_f(vx, vc);

    return SDE;
}



double cell_SDE(const TetrahedralMesh& mesh,
                const OpenVolumeMesh::CellHandle& ch,
                const double ref_scaling){

    auto c_verts = mesh.get_cell_vertices(ch);
    //log(LogLevel::INFO)<<" ---- checking boundary face "<<hf_verts[0]<<", "<<hf_verts[1]<<", "<<hf_verts[2]<<std::endl;

    using ElemT = SD_B_3D_AD;
    ElemT elem;

    ElemT::VecV vx;
    ElemT::VecC vc;
    std::vector<Eigen::Vector3d> P_ref(4);
    P_ref[0] <<  1, 1, 1;
    P_ref[1] << -1, 1,-1;
    P_ref[2] <<  1,-1,-1;
    P_ref[3] << -1,-1, 1;

    Eigen::Matrix3d E_ref;
    for(int i(0); i<3; i++){
        E_ref.col(i) = P_ref[i+1] - P_ref[0];
    }
    double reg_volume = E_ref.determinant() / 6.0;
    double scale_factor = std::pow(1.0 / reg_volume, 1.0 / 3.0);
    for(int i(0); i<P_ref.size(); i++){
        P_ref[i] *= ref_scaling * scale_factor;
    }

    vc[0] = 1.0;
    vc[1] = -1.0;
    vc[2] = 1e20;
    for(int i(0); i<c_verts.size(); i++){
        auto p = mesh.vertex(c_verts[i]);
        auto ref_p = P_ref[i];

        ElemT::compute_constants({1,0,1}, P_ref, vc);

        for(int j(0); j < p.size(); j++){
            vx[0 + 3 * i + j] = p[j];
            //vc[3 + 3 * i + j] = ref_p[j];
        }
    }

    double SDE = elem.eval_f(vx, vc);

    return SDE;
}


int check_boundary_holonomy(const TetrahedralMesh& mesh) {

    log(LogLevel::DEBUG)<<" ------------------------------------------------------------------ "<<std::endl;
    log(LogLevel::DEBUG)<<" -> checking boundary edges holonomy..."<<std::endl;

    double min_dihedral_sum(std::numeric_limits<double>::max());
    double max_dihedral_sum(0.0);
    double avg_dihedral_sum(0.0);


    int more_than_one_turn_count(0);
    double eps(1e-7);

    int boundary_edges_count(0);
    for (auto eh: mesh.edges()) {
        if (mesh.is_boundary(eh)) {

            double sum = compute_dihedral_angle_sum_deg(mesh, eh);

            min_dihedral_sum = std::min(min_dihedral_sum, sum);
            max_dihedral_sum = std::max(max_dihedral_sum, sum);
            avg_dihedral_sum += sum;

            more_than_one_turn_count += (sum >= 360.0 - eps);

            if (sum >= 360.0 - eps) {
                //log(LogLevel::INFO)<<" - dihedral sum for edge "<<edge_to_str(mesh, eh)<<" = "<<sum<<std::endl;

                TetrahedralMesh one_ring_mesh;
                EdgeHandle local_eh;
                setup_mesh_as_one_ring_submesh(mesh, eh, one_ring_mesh, local_eh);

                OVM::IO::FileManager fm;
                fm.writeFile("wrong_holonomy_edge_"+std::to_string(eh.idx())+ ".ovm", one_ring_mesh);
            }


            boundary_edges_count++;
        }
    }

    avg_dihedral_sum /= boundary_edges_count;

    log(LogLevel::INFO)<<" - min dihedral sum = "<<min_dihedral_sum<<std::endl;
    log(LogLevel::INFO)<<" - max dihedral sum = "<<max_dihedral_sum<<std::endl;
    log(LogLevel::INFO)<<" - avg dihedral sum = "<<avg_dihedral_sum<<std::endl;
    log(LogLevel::INFO)<<" --> # of boundary edges with more than one turn: "<<more_than_one_turn_count<<" / "<<boundary_edges_count<<std::endl;
    log(LogLevel::INFO)<<" ------------------------------------------------------------------ "<<std::endl;

    return more_than_one_turn_count;
}


bool valid_holonomy(const TetrahedralMesh& mesh,
                            const OpenVolumeMesh::EdgeHandle& eh) {

    double eps(1e-7);

    double sum = compute_dihedral_angle_sum_deg(mesh, eh);
    return (sum < 360.0 + eps);

}

void check_boundary_SDE(const TetrahedralMesh& ref_mesh,
                        const TetrahedralMesh& mesh,
                        double& min_SDE,
                        double& avg_SDE,
                        double& max_SDE){



    int b_face_count(0);
    double min_area(std::numeric_limits<double>::max()), max_area(0);
    double ref_min_area(std::numeric_limits<double>::max()), ref_max_area(0);
    double total_boundary_area(0);
    double ref_total_boundary_area(0);

    for(auto hfh: mesh.halffaces()) {

        if (mesh.is_boundary(hfh)) {

            auto hf_verts = mesh.get_halfface_vertices(hfh);

            auto p0 = mesh.vertex(hf_verts[0]);
            auto p1 = mesh.vertex(hf_verts[1]);
            auto p2 = mesh.vertex(hf_verts[2]);

            auto ref_p0 = ref_mesh.vertex(hf_verts[0]);
            auto ref_p1 = ref_mesh.vertex(hf_verts[1]);
            auto ref_p2 = ref_mesh.vertex(hf_verts[2]);

            Vec3d e1 = p1 - p0;
            Vec3d e2 = p2 - p0;

            Vec3d ref_e1 = ref_p1 - ref_p0;
            Vec3d ref_e2 = ref_p2 - ref_p0;

            double area = e1.cross(e2).norm();
            double ref_area = ref_e1.cross(ref_e2).norm();

            total_boundary_area += area;
            ref_total_boundary_area += ref_area;

            min_area = std::min(min_area, area);
            max_area = std::max(max_area, area);

            ref_min_area = std::min(ref_min_area, ref_area);
            ref_max_area = std::max(ref_max_area, ref_area);

            b_face_count++;
        }
    }

    log(LogLevel::INFO)<<" ---------------------- domain ∂ area: "<<std::endl;
    log(LogLevel::INFO)<<" -        min: "<<ref_min_area<<std::endl;
    log(LogLevel::INFO)<<" -        avg: "<<(ref_total_boundary_area / b_face_count)<<std::endl;
    log(LogLevel::INFO)<<" -        max: "<<ref_max_area<<std::endl;
    log(LogLevel::INFO)<<" -      total: "<<ref_total_boundary_area<<std::endl;
    log(LogLevel::INFO)<<" ------------------------------------"<<std::endl;

    log(LogLevel::INFO)<<" ---------------------- codomain ∂ area: "<<std::endl;
    log(LogLevel::INFO)<<" -        min: "<<min_area<<std::endl;
    log(LogLevel::INFO)<<" -        avg: "<<(total_boundary_area / b_face_count)<<std::endl;
    log(LogLevel::INFO)<<" -        max: "<<max_area<<std::endl;
    log(LogLevel::INFO)<<" -      total: "<<total_boundary_area<<std::endl;
    log(LogLevel::INFO)<<" ------------------------------------"<<std::endl;

    double ref_to_mesh_area_ratio = ref_total_boundary_area / total_boundary_area;
    double ref_scaling = std::sqrt(ref_to_mesh_area_ratio);

    min_SDE = std::numeric_limits<double>::max();
    max_SDE = 0;
    avg_SDE = 0;
    for(auto hfh: mesh.halffaces()){

        if(mesh.is_boundary(hfh)) {

            //log(LogLevel::INFO)<<" -------------------- checking boundary halfface "<<hfh<<": "<<mesh.get_halfface_vertices(hfh)<<std::endl;

            double SDE = halfface_SDE(ref_mesh, mesh, hfh, ref_scaling);
            //log(LogLevel::INFO)<<" --> SDE = "<<SDE<<std::endl;

            //SDE *= (ref_area / total_boundary_area);
            //SDE *= mesh_to_ref_area_ratio;
            //log(LogLevel::INFO)<<" --> weighted SDE = "<<SDE<<std::endl;


            avg_SDE += SDE;
            min_SDE = std::min(min_SDE, SDE);
            max_SDE = std::max(max_SDE, SDE);

            //log(LogLevel::INFO)<<" stopping here for debugging"
            //return;
        }
    }

    avg_SDE /= b_face_count;

}


void check_full_SDE(TetrahedralMesh& mesh,
                    bool scale_ref_to_element_volume,
                    double& min_SDE,
                    double& avg_SDE,
                    double& max_SDE){

    double min_volume(std::numeric_limits<double>::max()), max_volume(0);
    double total_volume(0);

    auto cell_vol_prop = mesh.request_cell_property<double>();

    for(auto ch: mesh.cells()){
        auto c_verts = mesh.get_cell_vertices(ch);

        Eigen::Matrix3d E;

        for(int i(0); i<3; i++) {
            auto e = mesh.vertex(c_verts[i+1]) - mesh.vertex(c_verts[0]);
            //log(LogLevel::INFO)<<" - e"<<i<<": "<<e<<std::endl;
            for(auto j(0); j<3; j++){
                E(j,i) = e[j];
            }
        }

        double vol = E.determinant() / 6.0;

        cell_vol_prop[ch] = vol;

        min_volume = std::min(min_volume, vol);
        max_volume = std::max(max_volume, vol);
        total_volume += vol;
    }
    double avg_volume = total_volume / mesh.n_logical_cells();;

    log(LogLevel::INFO)<<" ---------------------- mesh volume: "<<std::endl;
    log(LogLevel::INFO)<<" -        min: "<<min_volume<<std::endl;
    log(LogLevel::INFO)<<" -        avg: "<<(avg_volume)<<std::endl;
    log(LogLevel::INFO)<<" -        max: "<<max_volume<<std::endl;
    log(LogLevel::INFO)<<" -      total: "<<total_volume<<std::endl;
    log(LogLevel::INFO)<<" ------------------------------------"<<std::endl;

    double scaling_factor = std::pow(avg_volume, 1.0 / 3.0);

    log(LogLevel::INFO)<<" -> scaling factor = "<<scaling_factor<<std::endl;

    min_SDE = std::numeric_limits<double>::max();
    max_SDE = 0;
    avg_SDE = 0;
    for(auto ch: mesh.cells()){

            //log(LogLevel::INFO)<<" -------------------- checking cell "<<ch<<": "<<mesh.get_cell_vertices(ch)<<std::endl;
            if(scale_ref_to_element_volume){
                scaling_factor = std::pow(cell_vol_prop[ch], 1.0 / 3.0);
                //log(LogLevel::INFO)<<" -> element scaling factor = "<<scaling_factor<<std::endl;
            }
            double SDE = cell_SDE(mesh, ch, scaling_factor);

        //log(LogLevel::INFO)<<" -> SDE = "<<SDE<<std::endl;

            avg_SDE += SDE;
            min_SDE = std::min(min_SDE, SDE);
            max_SDE = std::max(max_SDE, SDE);
    }

    avg_SDE /= mesh.n_logical_cells();

}


void check_full_SDE(const TetrahedralMesh& mesh,
                    const TetrahedralMesh& ref_mesh,
                    double& min_SDE,
                    double& avg_SDE,
                    double& max_SDE){

    double tot_vol = compute_total_volume(mesh);
    double ref_tot_vol = compute_total_volume(ref_mesh);
    double ref_to_mesh_scaling_factor = std::pow(tot_vol / ref_tot_vol, 1.0 / 3.0);

    log(LogLevel::INFO)<<" -> scaling factor = "<<ref_to_mesh_scaling_factor<<std::endl;

    min_SDE = std::numeric_limits<double>::max();
    max_SDE = 0;
    avg_SDE = 0;
    for(auto ch: mesh.cells()){

        double SDE = cell_SDE(mesh, ref_mesh, ch, ref_to_mesh_scaling_factor);

        avg_SDE += SDE;
        min_SDE = std::min(min_SDE, SDE);
        max_SDE = std::max(max_SDE, SDE);
    }

    avg_SDE /= mesh.n_logical_cells();

}

void check_full_conformal_distortion(const TetrahedralMesh& mesh,
                                     const TetrahedralMesh& ref_mesh,
                                     double& min_dist,
                                     double& avg_dist,
                                     double& max_dist){

    min_dist = std::numeric_limits<double>::max();
    avg_dist = 0;
    max_dist = 0;
    for(auto ch: mesh.cells()){

        double conf_dist = cell_conformal_distortion(mesh, ref_mesh, ch);

        avg_dist += conf_dist;
        min_dist = std::min(min_dist, conf_dist);
        max_dist = std::max(max_dist, conf_dist);
    }

    avg_dist /= mesh.n_logical_cells();
}



int euler_characteristic(const TetrahedralMesh& mesh){

    return (int)mesh.n_logical_vertices() - (int)mesh.n_logical_edges() + (int)mesh.n_logical_faces() - (int)mesh.n_logical_cells();
}


int is_surface_mesh_boundary(const TetrahedralMesh& boundary_mesh,
                                    const OpenVolumeMesh::VertexHandle& vh){

    int incident_boundary_edges_count(0);
    for(auto out_he: boundary_mesh.outgoing_halfedges(vh)){
        auto eh = boundary_mesh.edge_handle(out_he);
        int e_val = boundary_mesh.valence(eh);
        if(e_val < 1 || e_val > 2){
            log(LogLevel::ERROR)<<" ERROR - edge ("<<boundary_mesh.edge(eh).from_vertex()<<", "<<boundary_mesh.edge(eh).to_vertex()<<") has valence "<<e_val<<", which shouldn't happen on the boundary"<<std::endl;
            return -1;
        }
        incident_boundary_edges_count += e_val == 1;
    }

    if(incident_boundary_edges_count == 1 ||
       incident_boundary_edges_count > 2){
        log(LogLevel::ERROR)<<" ERROR - vertex "<<vh<<" is incident to "<<incident_boundary_edges_count<<" boundary edges, which shouldn't happen on the boundary"<<std::endl;
        return -1;
    }

    return incident_boundary_edges_count == 2;
}

double compute_min_IVL_energy(const TetrahedralMesh& mesh,
                              const OpenVolumeMesh::CellHandle ch){

    int n_vars = 4 * 3;
    COMISO::FiniteElementProblem fe_problem(n_vars);
    COMISO::FiniteElementSet<IVL_B_3D_AD> fe_ivl("inverted volume-length");
    std::vector<COMISO::LinearConstraint> boundary_eq_constraints;
    boundary_eq_constraints.reserve(9);

    IVL_B_3D_AD::VecI vi;
    IVL_B_3D_AD::VecC vc;
    vc << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    auto verts = mesh.get_cell_vertices(ch);
    for(int i(0); i<4; i++){

        auto pos = mesh.vertex(verts[i]);
        for(int j(0); j<3; j++){
            int gid = 3 * i + j;

            fe_problem.x()[gid] = pos[j];
            vi[gid] = 3 * i + j;
        }

        if(mesh.is_boundary(verts[i])){

            for (int j(0); j<3; j++) {

                int gid = 3 * i + j;

                COMISO::LinearConstraint::SVectorNC coeffs;
                coeffs.resize(fe_problem.n_unknowns());

                coeffs.coeffRef(gid) = 1.0;
                boundary_eq_constraints.emplace_back(coeffs, -pos[j]);
            }
        }
    }

    //if there are less than 2 vertices constrained then it should always be able to get to 0 anyway
    if(boundary_eq_constraints.size() < 6){
        return 0.0;
    }


    fe_ivl.instances().add_element(vi, vc);
    fe_problem.add_set(&fe_ivl);

    std::vector<COMISO::NConstraintInterface *> constraints_pointers;
    constraints_pointers.reserve(boundary_eq_constraints.size());
    for (auto &c: boundary_eq_constraints) {
        constraints_pointers.push_back(&c);
    }

    COMISO::TruncatedNewtonPCG tn(1e-8);
    tn.always_update_preconditioner() = true;
    tn.max_iters() = 1000;
    tn.max_pcg_iters() = 100;
    tn.allow_warmstart() = false;
    tn.tolerance_gdx() = 1e-3; //previously 1e-1
    tn.set_silent(true);
    tn.set_eps_constraints_violation(1e-9); //1e-9
    tn.solve_projected_normal_equation(&fe_problem, constraints_pointers);

    return fe_problem.eval_f(fe_problem.x().data());

}

Vec3d compute_visible_point(const TetrahedralMesh &mesh,
                            const OpenVolumeMesh::VertexHandle &vh,
                            const std::vector<OpenVolumeMesh::HalfFaceHandle>& incident_hfs) {

    Vec3d normal_acc = {0,0,0};
    double vec_norm_acc = 0;

    for (auto hf: incident_hfs) {

        Vec3d c, n;
        compute_centroid_and_normal(mesh,  hf,  c,  n);

        normal_acc += n;
        vec_norm_acc += c.norm();
    }

    return normal_acc.normalized() * vec_norm_acc / incident_hfs.size();


}

Vec3d compute_rotational_mid_point(const TetrahedralMesh& mesh,
                                   const OpenVolumeMesh::EdgeHandle& eh,
                                   const std::pair<OpenVolumeMesh::HalfFaceHandle, OpenVolumeMesh::HalfFaceHandle>& incident_hfs){

    auto log_level = LogLevel::DEBUG;

    log(log_level)<<" ------------------------------------------ "<<std::endl;
    //log(LogLevel::WARNING)<<" WARNING - still using centroids for norm "<<std::endl;


    auto left_hf = incident_hfs.first;
    auto right_hf = incident_hfs.second;

    Vec3d left_c, left_n, right_c, right_n;
    compute_centroid_and_normal(mesh,  left_hf,  left_c,  left_n);
    compute_centroid_and_normal(mesh, right_hf, right_c, right_n);

    log(log_level)<<" -  left c / n: "<<left_c<<" / " <<left_n<<std::endl;
    log(log_level)<<" - right c / n: "<<right_c<<" / " <<right_n<<std::endl;

    auto  left_verts = mesh.get_halfface_vertices( left_hf);
    auto right_verts = mesh.get_halfface_vertices(right_hf);


    auto from_v = mesh.edge(eh).from_vertex();
    auto   to_v = mesh.edge(eh).to_vertex();

    //do some rotation magic
    auto from_pos = mesh.vertex(from_v);
    auto   to_pos = mesh.vertex(to_v);
    auto  mid_pos = 0.5 * (from_pos + to_pos);

    auto  left_v = mesh.to_vertex_handle(mesh.next_halfedge_in_halfface(mesh.halfedge_handle(eh, 1),  left_hf));
    auto right_v = mesh.to_vertex_handle(mesh.next_halfedge_in_halfface(mesh.halfedge_handle(eh, 0), right_hf));

    if (left_v == from_v || left_v == to_v ||
        right_v == from_v || right_v == to_v) {

        left_v = mesh.to_vertex_handle(mesh.next_halfedge_in_halfface(mesh.halfedge_handle(eh, 0),  left_hf));
        right_v = mesh.to_vertex_handle(mesh.next_halfedge_in_halfface(mesh.halfedge_handle(eh, 1), right_hf));

        if (left_v == from_v || left_v == to_v ||
            right_v == from_v || right_v == to_v) {
            log(LogLevel::ERROR)<<" ERROR - either the left vertex "<<left_v<<" or the right vertex "<<right_v<<" is actually either of the end vertices "<<from_v<<" or "<<to_v<<std::endl;
            exit(EXIT_FAILURE);
        }
        }
    log(log_level)<<" -  left v: "<<left_v<<std::endl;
    log(log_level)<<" - right v: "<<right_v<<std::endl;

    auto  left_p = mesh.vertex( left_v) - mid_pos;
    auto right_p = mesh.vertex(right_v) - mid_pos;
    log(log_level)<<" -  left p: "<< left_v<<", norm = "<< left_p.norm()<<std::endl;
    log(log_level)<<" - right p: "<<right_v<<", norm = "<<right_p.norm()<<std::endl;

    double dot_prod = left_n.normalized().dot(right_n.normalized());
    double theta = std::acos(std::min(1.0, std::max(-1.0, dot_prod)));
    log(log_level)<<" - theta = "<<rad_to_deg(theta)<<"°"<<std::endl;

    //double thresh_deg = 1.0;
    log(log_level)<<" - using average normal as direction"<<std::endl;

    auto d = (left_n + right_n).normalized();
    auto norm = (left_p.norm() + right_p.norm()) * 0.5;
    log(log_level)<<" -    d: "<<d<<std::endl;
    log(log_level)<<" - norm: "<<norm<<std::endl;

    Vec3d x = mid_pos + d * norm;

    /*}else if (std::abs(theta) <= thresh_deg * (M_PI / 180.0)) {
        log(log_level)<<" - faces have the same normal"<<std::endl;
        x = (to_pos + from_pos) * 0.5 + left_n  * (right_v.norm() + left_v.norm()) * 0.5;

    }else {
        auto rot_axis = left_n.cross(right_n).normalized();

        auto rot_theta_deg = theta * 180 / M_PI;

        log(log_level)<<" - rotation  axis: "<<rot_axis<<std::endl;
        log(log_level)<<" - rotation angle: "<<rot_theta_deg<<std::endl;

        //NOTE: + π / 2.0 to shift from on-face vector to normal
        Eigen::Matrix3d rot_mat = Eigen::AngleAxisd(0.5 * theta + M_PI * 0.5, vec2vec(rot_axis)).toRotationMatrix();
        Eigen::Vector3d rot_left_v = rot_mat * vec2vec(left_v);
        x = from_pos + vec2vec(rot_left_v);
    }*/

    log(log_level)<<" ->  x: "<<x<<std::endl;

    return x;
}



Vec3d compute_optimal_virtual_vertex_pos(const TetrahedralMesh &mesh,
                                         const OpenVolumeMesh::EdgeHandle &eh,
                                         const std::vector<OpenVolumeMesh::HalfFaceHandle>& incident_hfs) {
    auto log_level = LogLevel::DEBUG;

    log(log_level)<<" ------------------------------------------ "<<std::endl;
    log(log_level)<<" computing optimal point for edge "<<edge_to_str(mesh, eh)<<std::endl;

    using Elem = FOF_B_3D_AD;

    /*auto left_hf = incident_hfs.first;
    auto right_hf = incident_hfs.second;

    auto  left_verts = mesh.get_halfface_vertices( left_hf);
    auto right_verts = mesh.get_halfface_vertices(right_hf);*/

    int n_vars = 5 * 3;

    COMISO::FiniteElementProblem fe_problem(n_vars);
    COMISO::FiniteElementSet<Elem> fe_virtual_elements;

    std::vector<COMISO::LinearConstraint> eq_constraints;
    eq_constraints.reserve(n_vars);

    typename Elem::VecC vc_virtual;


    std::vector<Eigen::Vector3d> P_ref(4);
    P_ref[0] <<  1, 1, 1;
    P_ref[1] << -1, 1,-1;
    P_ref[2] <<  1,-1,-1;
    P_ref[3] << -1,-1, 1;

    /*const double& epsilon        = _c[0];
    const double& conf_w         = _c[1];
    const double& vol_w          = _c[2];
    const double& barrier_w      = _c[3];
    const double& barrier_thresh = _c[4];*/

    std::vector parameters = {0.0, 1.0, 0.0, 0.0, 0.0};

    Elem::compute_constants(parameters, P_ref, vc_virtual);

    auto x_init = compute_rotational_mid_point(mesh, eh, {incident_hfs[0], incident_hfs[1]});
    log(log_level) << " - x_init (rot. mid-point):" <<x_init << std::endl;

    std::map<VertexHandle, int> var_id_map;

    int virtual_v_var_id(4);
    int var_id_count(0);

    for (auto hf: incident_hfs) {
        auto verts = mesh.get_halfface_vertices(hf);

        typename Elem::VecV vx_virtual;
        typename Elem::VecI vi_virtual;

        std::vector<OVM::VertexHandle> local_verts;
        for (int i(0); i < verts.size(); i++) {
            auto v = verts[i];
            auto pos = mesh.vertex(v);

            bool already_added_constraints(false);
            int var_id = -1;
            if (var_id_map.contains(v)) {
                var_id = var_id_map[v];
                already_added_constraints = true;
            } else {
                var_id = var_id_count++;
                var_id_map.insert({v, var_id});
                log(log_level) << " -- v"<<v<<" at : "<<pos<<" has index "<<var_id << std::endl;
            }


            for (int j(0); j < 3; j++) {
                int local_idx = 3 * i + j;
                int global_idx = 3 * var_id + j;

                vi_virtual[local_idx] = global_idx;
                fe_problem.x()[global_idx] = pos[j];
                vx_virtual[local_idx] = fe_problem.x()[global_idx];

                if (!already_added_constraints) {
                    COMISO::LinearConstraint::SVectorNC coeffs;
                    coeffs.resize(n_vars);
                    coeffs.coeffRef(global_idx) = 1;
                    eq_constraints.push_back(COMISO::LinearConstraint(coeffs, -pos[j]));
                }
            }
        }


        for (int j(0); j < 3; j++) {
            int local_idx = 3 * 3 + j;
            int global_idx = 3 * virtual_v_var_id + j;

            vi_virtual[local_idx] = global_idx;
            fe_problem.x()[global_idx] = x_init[j];
            vx_virtual[local_idx] = fe_problem.x()[global_idx];
        }

        log(log_level)<<" - vx: "<<vx_virtual.transpose()<<std::endl;
        log(log_level)<<" - vi: "<<vi_virtual.transpose()<<std::endl;

        Elem elem;
        auto e = elem.eval_f(vx_virtual, vc_virtual);
        log(log_level) << " --> energy = " << e << std::endl;

        if (!std::isfinite(e)) {
            log(LogLevel::WARNING) << " WARNING - infinite energy for virtual element incident to halfface "<<halfface_to_str(mesh, hf) << std::endl;
            return x_init;
        }
        fe_virtual_elements.instances().add_element(vi_virtual, vc_virtual);
    }

    fe_problem.add_set(&fe_virtual_elements);


    std::vector<COMISO::NConstraintInterface*> constraints_pointers;
    constraints_pointers.reserve(constraints_pointers.size());
    for( auto& c : eq_constraints) {
        constraints_pointers.push_back(&c);
    }
    //log(LogLevel::INFO) << " energies before optimization: " << std::endl;
    //fe_problem.print_objectives();

    double initial_energy = fe_problem.eval_f(fe_problem.x().data());
    if(initial_energy == std::numeric_limits<double>::infinity()){
        log(LogLevel::WARNING) << " WARNING - infinite energy for virtual element incident to edge "<<edge_to_str(mesh, eh) << std::endl;
        return x_init;
    }

    COMISO::TruncatedNewtonPCG tn(1e-5);
    tn.always_update_preconditioner() = true;
    tn.max_iters() = 100;
    tn.max_pcg_iters() = 100;
    tn.allow_warmstart() = false;
    tn.tolerance_gdx() = 1e-3;
    tn.set_silent(true);
    tn.set_eps_constraints_violation(1e-9);
    tn.set_silent(true);

    tn.solve_projected_normal_equation(&fe_problem, constraints_pointers);

    //log(LogLevel::INFO) << " energies after optimization: " << std::endl;
    //fe_problem.print_objectives();

    Vec3d x = {
        fe_problem.x()[3 * virtual_v_var_id + 0],
        fe_problem.x()[3 * virtual_v_var_id + 1],
        fe_problem.x()[3 * virtual_v_var_id + 2]};
    log(log_level) << " --> optimal x: "<<x << std::endl;

    return x;
}


void compute_centroid_and_normal(const TetrahedralMesh &mesh,
                                 const OpenVolumeMesh::HalfFaceHandle hfh,
                                 Vec3d &centroid,
                                 Vec3d &normal) {

    auto verts = mesh.get_halfface_vertices(hfh);

    std::vector p = {mesh.vertex(verts[0]), mesh.vertex(verts[1]), mesh.vertex(verts[2])};

    centroid = (p[0] + p[1] + p[2]) / 3.0;
    normal = (p[1] - p[0]).cross(p[2] - p[0]).normalized();
}


double compute_dihedral_angle_sum_deg(const TetrahedralMesh& mesh,
                                      const OVM::EdgeHandle& eh) {

    log(LogLevel::DEBUG)<<" ------------------------------------ checking boundary edge "<<edge_to_str(mesh, eh)<<std::endl;
    double interior_angle_sum_deg(0);
    for (auto ec_it = mesh.ec_iter(eh); ec_it.valid(); ec_it++) {
        log(LogLevel::DEBUG)<<" ------ checking cell "<<cell_to_str(mesh, *ec_it)<<std::endl;
        for (auto vh: mesh.get_cell_vertices(*ec_it)) {
            log(LogLevel::DEBUG)<<" -- vertex "<<vh<<" at "<<mesh.vertex(vh)<<std::endl;
        }

        auto angle_deg = compute_dihedral_angle_deg(mesh, *ec_it, eh);
        interior_angle_sum_deg += angle_deg;
        log(LogLevel::DEBUG)<<" --- angle = "<<angle_deg<<"°"<<std::endl;
    }
    log(LogLevel::DEBUG)<<" ------------------------------------ total = "<<interior_angle_sum_deg<<std::endl;

    return interior_angle_sum_deg;
}

double compute_dihedral_angle_deg(const TetrahedralMesh& mesh,
                                  const OVM::CellHandle ch,
                                  const OVM::EdgeHandle& eh){

    OVM::HalfFaceHandle hf1(-1);
    OVM::HalfFaceHandle hf2(-1);
    for(auto ehf_it = mesh.ehf_iter(eh); ehf_it.valid(); ehf_it++){
        if(mesh.incident_cell(*ehf_it) == ch){
            if(!hf1.is_valid()){
                hf1 = *ehf_it;
            }else{
                hf2 = *ehf_it;
                break;
            }
        }
    }

    if(!hf1.is_valid() || ! hf2.is_valid()){
        log(LogLevel::ERROR)<<" ERROR - couldn't find two incident halffaces incident to cell "<<ch<<" and edge "<<(eh)<<std::endl;
        return -1;
    }

    std::vector<Vec3d> normals;

    for(auto hf: std::vector<OVM::HalfFaceHandle>({hf1, hf2})){

        auto face_vertices = mesh.get_halfface_vertices(hf);

        std::vector<Vec3d> positions;
        for(auto v: face_vertices){
            auto incident_c = mesh.incident_cell(hf);
            auto incident_c_vertices = mesh.get_cell_vertices(incident_c);
            for(int i(0); i < incident_c_vertices.size(); i++){
                if(incident_c_vertices[i] == v){
                    positions.push_back(mesh.vertex(incident_c_vertices[i]));
                }
            }
        }

        auto c_evi  = positions[1] - positions[0];
        auto c_evi1 = positions[2] - positions[0];


        auto normal = (c_evi.cross(c_evi1));
        normal.normalize();

        normals.push_back(normal);
    }

    //auto dihedral_angle = (M_PI - acos(normals[0].dot(normals[1]))) * 180 / M_PI;
    auto dihedral_angle = (M_PI - acos(normals[0].dot(normals[1])));

    return 180.0 * dihedral_angle / M_PI;
}


double compute_boundary_dihedral_angle_deg(const TetrahedralMesh &mesh,
                                           const OpenVolumeMesh::EdgeHandle eh) {
    std::vector<OVM::HalfFaceHandle> boundary_halffaces;
    for (auto ehf_it = mesh.ehf_iter(eh); ehf_it.valid(); ehf_it++) {
        if (mesh.is_boundary(*ehf_it)) {
            boundary_halffaces.push_back(*ehf_it);
        }
    }

    auto l_verts = mesh.get_halfface_vertices(boundary_halffaces[1]);
    auto r_verts = mesh.get_halfface_vertices(boundary_halffaces[0]);

    auto from_v = mesh.edge(eh).from_vertex();
    auto to_v = mesh.edge(eh).to_vertex();

    auto from_v_pos = mesh.vertex(from_v);
    auto to_v_pos = mesh.vertex(to_v);

    OVM::VertexHandle l_v, r_v;
    for (int i(0); i < 3; i++) {
        if (from_v != l_verts[i] && to_v != l_verts[i]) {
            l_v = l_verts[i];
        }
        if (from_v != r_verts[i] && to_v != r_verts[i]) {
            r_v = r_verts[i];
        }
    }

    auto l_pos = mesh.vertex(l_v);
    auto r_pos = mesh.vertex(r_v);

    auto e = (to_v_pos - from_v_pos).normalized();

    auto l_normal = (to_v_pos - from_v_pos).cross(l_pos - from_v_pos).normalized();
    auto r_normal = (to_v_pos - from_v_pos).cross(r_pos - from_v_pos).normalized();

    auto s = (l_normal.cross(r_normal)).dot(e);
    auto c = (l_normal.dot(r_normal));

    double theta_deg = 180.0 * std::atan2(s, c) / M_PI;
    if (theta_deg < 0) { theta_deg += 360; } //to map to {0,2pi}


    return theta_deg;
}

double compute_interior_angle_deg(const TetrahedralMesh& mesh,
                                  const OpenVolumeMesh::FaceHandle fh,
                                  const OpenVolumeMesh::VertexHandle vh) {

    auto hf_vertices = mesh.get_halfface_vertices(mesh.halfface_handle(fh, 0));

    int vid = hf_vertices[0] == vh ? 0 : (hf_vertices[1] == vh ? 1 : 2);
    int lid = (3 + vid - 1) % 3;//-1 to get the previous one and +3 to avoid the modulo of a negative number
    int rid = (3 + vid + 1) % 3;

    OVM::VertexHandle lvh(hf_vertices[lid]), rvh(hf_vertices[rid]);

    auto l = vec2vec(mesh.vertex(hf_vertices[lid]) - mesh.vertex(hf_vertices[vid]));
    auto r = vec2vec(mesh.vertex(hf_vertices[rid]) - mesh.vertex(hf_vertices[vid]));

    double interior_angle = std::acos(l.normalized().dot(r.normalized()));

    return 180.0 * interior_angle / M_PI;
}


double compute_volume(const TetrahedralMesh& mesh,
                      const OpenVolumeMesh::CellHandle ch) {

    Eigen::Matrix3d E;

    auto c_verts = mesh.get_cell_vertices(ch);

    for (int i(0); i < 3; i++) {
        auto e = mesh.vertex(c_verts[i + 1]) - mesh.vertex(c_verts[0]);
        for (auto j(0); j < 3; j++) {
            E(j, i) = e[j];
        }
    }

    return E.determinant() / 6.0;
}


double compute_area(const TetrahedralMesh& mesh,
                    const OpenVolumeMesh::HalfFaceHandle  hfh){

    auto hf_verts = mesh.get_halfface_vertices(hfh);
    std::vector<Vec3d> P;
    for(auto v: hf_verts){
        P.push_back(mesh.vertex(v));
    }

    double area = 0.5 * (P[1] - P[0]).cross(P[2] - P[0]).norm();
    return area;
}

double compute_total_volume(const TetrahedralMesh& mesh){

    double vol(0);
    for(auto ch: mesh.cells()){
        vol += compute_volume(mesh, ch);
    }
    return vol;
}

double compute_min_edge_length(const TetrahedralMesh& mesh) {

    double min_len(std::numeric_limits<double>::max());
    for (auto eh: mesh.edges()) {
        auto len = (mesh.vertex(mesh.edge(eh).from_vertex()) - mesh.vertex(mesh.edge(eh).to_vertex())).norm();
        min_len = std::min(min_len, len);
    }
    return min_len;
}


double compute_total_edge_length(const TetrahedralMesh& mesh) {

    double len(0);
    for (auto eh: mesh.edges()) {
        len += (mesh.vertex(mesh.edge(eh).from_vertex()) + mesh.vertex(mesh.edge(eh).to_vertex())).norm();
    }
    return len;
}

int count_boundary_faces(const TetrahedralMesh& mesh){

    int count(0);
    for(auto fh: mesh.faces()){
        count += mesh.is_boundary(fh);
    }
    return count;
}

int count_boundary_vertices(const TetrahedralMesh& mesh){

    int count(0);
    for(auto vh: mesh.vertices()){
        count += mesh.is_boundary(vh);
    }
    return count;
}


int count_boundary_edges(const TetrahedralMesh& mesh){

    int count(0);
    for(auto eh: mesh.edges()){
        count += mesh.is_boundary(eh);
    }
    return count;
}


int count_boundary_vertices(const TetrahedralMesh& mesh,
                            const OpenVolumeMesh::CellHandle& ch){
    auto verts = mesh.get_cell_vertices(ch);

    return mesh.is_boundary(verts[0]) +
            mesh.is_boundary(verts[1]) +
            mesh.is_boundary(verts[2]) +
            mesh.is_boundary(verts[3]);
}

double compute_boundary_surface_area(const TetrahedralMesh& mesh){

    double area(0);
    for(auto fh: mesh.faces()){
        if(mesh.is_boundary(fh)) {
            area += compute_area(mesh, mesh.halfface_handle(fh, 0));
        }
    }
    return area;
}

double compute_average_boundary_surface_area(const TetrahedralMesh& mesh){
    return compute_boundary_surface_area(mesh) / count_boundary_faces(mesh);
}


void compute_boundary_edge_lengths(const TetrahedralMesh& mesh, double& min_length, double& max_length, double& avg_length){


    int boundary_edge_count(0);

    min_length = std::numeric_limits<double>::max();
    max_length = 0.0;
    avg_length = 0.0;

    for (auto eh: mesh.edges()) {
        if (mesh.is_boundary(eh)) {
            double length = edge_length(mesh, eh);
            avg_length += length;

            min_length = std::min(min_length, length);
            max_length = std::max(max_length, length);

            boundary_edge_count++;
        }
    }

    avg_length /= boundary_edge_count;
}

double compute_optimal_number_of_elements(const TetrahedralMesh& mesh){

    log(LogLevel::INFO)<<" ---------------------------------------- "<<std::endl;
    log(LogLevel::INFO)<<" computing optimal number of elements..."<<std::endl;

    double average_boundary_area = compute_average_boundary_surface_area(mesh);
    log(LogLevel::INFO)<<" - average boundary area = "<<average_boundary_area<<std::endl;

    double optimal_edge_length = std::sqrt(average_boundary_area * 4.0 / std::sqrt(3.0));
    log(LogLevel::INFO)<<" - optimal edge length = "<<optimal_edge_length<<std::endl;

    double optimal_element_volume = std::pow(optimal_edge_length, 3.0) / (6.0 * std::sqrt(2.0));
    log(LogLevel::INFO)<<" - optimal element volume = "<<optimal_element_volume<<std::endl;

    double total_volume = compute_total_volume(mesh);
    log(LogLevel::INFO)<<" - total volume = "<<total_volume<<std::endl;

    double optimal_number_of_elements = total_volume / optimal_element_volume;
    log(LogLevel::INFO)<<" --> optimal number of elements = "<<optimal_number_of_elements<<std::endl;

    log(LogLevel::INFO)<<" ---------------------------------------- "<<std::endl;

    return optimal_number_of_elements;
}

Eigen::Vector3d compute_barycenter(const TetrahedralMesh& mesh){
    int boundary_v_count(0);
    Vec3d center = {0,0,0};
    for (auto vh: mesh.vertices()) {
        if (mesh.is_boundary(vh)) {
            auto pos = mesh.vertex(vh);
            center += pos;
            boundary_v_count++;
        }
    }
    center /= boundary_v_count;

    return vec2vec(center);
}

void translate_barycenter_to_origin(TetrahedralMesh& mesh) {

    log(LogLevel::INFO)<<" - translating barycenter to origin..."<<std::endl;

    auto center = vec2vec(compute_barycenter(mesh));
    log(LogLevel::INFO)<<" -> barycenter at "<<center<<std::endl;

    for (auto vh: mesh.vertices()) {
        auto pos = mesh.vertex(vh);
        mesh.set_vertex(vh, pos - center);
    }
}



void scale_to_unit_bounding_box_diagonal(TetrahedralMesh& mesh) {

    log(LogLevel::INFO)<<" - scaling to unit diagonal bounding box..."<<std::endl;

    Eigen::Vector3d min_corner, max_corner;
    compute_bounding_box(mesh, min_corner, max_corner);
    log(LogLevel::INFO)<<" - diagonal norm before: "<<(max_corner - min_corner).norm()<<std::endl;


    double scaling_factor = 1.0 / (max_corner - min_corner).norm();
    for (auto vh: mesh.vertices()) {
        auto pos = mesh.vertex(vh);
        mesh.set_vertex(vh, pos * scaling_factor);
    }

    compute_bounding_box(mesh, min_corner, max_corner);
    log(LogLevel::INFO)<<" - diagonal norm after: "<<(max_corner - min_corner).norm()<<std::endl;
}


void compute_bounding_box(const TetrahedralMesh &mesh,
                          Eigen::Vector3d &min_corner,
                          Eigen::Vector3d &max_corner) {

    min_corner = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    max_corner = {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};

    for (auto vh: mesh.vertices()) {
        if (mesh.is_boundary(vh)) {
            auto pos = mesh.vertex(vh);
            for (int i(0); i<3; i++) {
                min_corner[i] = std::min(min_corner[i], pos[i]);
                max_corner[i] = std::max(max_corner[i], pos[i]);
            }
        }
    }

}

Eigen::Vector3d compute_element_bounding_box_diagonal(const TetrahedralMesh& mesh,
                                                      const OpenVolumeMesh::CellHandle ch){

    Eigen::Vector3d min_corner = Eigen::Vector3d::Identity() * std::numeric_limits<double>::max();
    Eigen::Vector3d max_corner = Eigen::Vector3d::Identity() * std::numeric_limits<double>::lowest();

    auto verts = mesh.get_cell_vertices(ch);
    for(auto v: verts){
        auto pos = mesh.vertex(v);
        for(int j(0); j<3; j++){
            min_corner[j] = std::min(min_corner[j], pos[j]);
            max_corner[j] = std::max(max_corner[j], pos[j]);
        }
    }

    //log(LogLevel::INFO)<<" - diagonal: "<<(max_corner - min_corner).transpose()<<std::endl;

    return max_corner - min_corner;
}


Eigen::Vector3d compute_global_bounding_box_diagonal(const TetrahedralMesh& mesh){

    Eigen::Vector3d min_corner = Eigen::Vector3d::Identity() * std::numeric_limits<double>::max();
    Eigen::Vector3d max_corner = Eigen::Vector3d::Identity() * std::numeric_limits<double>::lowest();

    for(auto v: mesh.vertices()){
        if(mesh.is_boundary(v)) {
            auto pos = mesh.vertex(v);
            for (int j(0); j < 3; j++) {
                min_corner[j] = std::min(min_corner[j], pos[j]);
                max_corner[j] = std::max(max_corner[j], pos[j]);
            }
        }
    }

    return max_corner - min_corner;
}


Eigen::Vector3d vec2vec(const OpenVolumeMesh::Vec3d& vec){
    return {vec[0], vec[1], vec[2]};
}

Eigen::Vector2d vec2vec(const OpenVolumeMesh::Vec2d &vec){
    return {vec[0], vec[1]};
}

OpenVolumeMesh::Vec3d vec2vec(const Eigen::Vector3d& vec){
    return {vec[0], vec[1], vec[2]};
}

OpenVolumeMesh::Vec2d vec2vec(const Eigen::Vector2d& vec){
    return {vec[0], vec[1]};
}


Eigen::Vector3d eigen_conversion(const Eigen::Vector2d& vec){
    Eigen::Vector3d v;
    v << vec[0], vec[1], vec[2];
    return v;
}

Eigen::Vector3d eigen_conversion(const Eigen::Vector3d& vec){
    return vec;
}

void setup_boundary_adapted_ref_positions(const TetrahedralMesh& domain_mesh,
                                          const TetrahedralMesh& codomain_mesh,
                                          const OpenVolumeMesh::CellHandle& ch,
                                          double global_scaling_factor,
                                          std::vector<Eigen::Vector3d>& P_ref) {

    auto verts = domain_mesh.get_cell_vertices(ch);
    std::array<bool, 4> constraints_mask;

    for(int i(0); i<4; i++){
        constraints_mask[i] = domain_mesh.is_boundary(verts[i]);
    }

    setup_constraints_adapted_ref_positions(domain_mesh, codomain_mesh, ch, global_scaling_factor, constraints_mask, P_ref);

}


void setup_constraints_adapted_ref_positions(const TetrahedralMesh& domain_mesh,
                                             const TetrahedralMesh& codomain_mesh,
                                             const OpenVolumeMesh::CellHandle& ch,
                                             double global_scaling_factor,
                                             const std::array<bool, 4>& constraints_mask,
                                             std::vector<Eigen::Vector3d>& P_ref){

    bool debug(false);

    using Vec = Eigen::Vector3d;

    if(debug) {
        log(LogLevel::INFO) << " ---------- setting up ref positions with boundary scaling for element " << ch << std::endl;
        log(LogLevel::INFO) << " - global scaling factor: " << global_scaling_factor << std::endl;
    }

    auto verts = domain_mesh.get_cell_vertices(ch);
    int constrained_vertices_count(0);

    std::vector<int> local_idx_prop(domain_mesh.n_cells(), -1);
    std::vector<Eigen::Vector3d> P_domain, P_codomain;
    //log(LogLevel::INFO)<<" - vertices: "<<std::endl;
    for(int i(0); i<4; i++){
        auto vh = verts[i];
        P_domain.push_back(global_scaling_factor * vec2vec(domain_mesh.vertex(vh)));
        P_codomain.push_back(vec2vec(codomain_mesh.vertex(vh)));
        local_idx_prop[vh.idx()] = i;
        //constrained_vertices_count += domain_mesh.is_boundary(vh);
        constrained_vertices_count += constraints_mask[i];
        //log(LogLevel::INFO)<<" --- "<<i<<": "<<(domain_mesh.is_boundary(vh) ? "boundary" : "interior")<<" vertex "<<vh<<" (scaled) domain pos at "<<P_domain.back().transpose()<<" and codomain pos at "<<P_codomain.back().transpose()<<std::endl;
    }
    auto domain_volume = compute_volume(domain_mesh, ch) * std::pow(global_scaling_factor, 3.0);

    if(debug) {
        log(LogLevel::INFO)<<" - scaled domain volume = "<<domain_volume<<std::endl;
    }

    //used as basis
    P_ref = P_domain;

    if(constrained_vertices_count == 3){

        auto unconstrained_v = !constraints_mask[0] ? verts[0] : (!constraints_mask[1] ? verts[1] : (!constraints_mask[2] ? verts[2] : verts[3]));

        OVM::HalfFaceHandle hfh = domain_mesh.vertex_opposite_halfface(ch, unconstrained_v);
        /*for(auto chf_it = domain_mesh.chf_iter(ch); chf_it.valid(); chf_it++){
            //if(domain_mesh.is_boundary(domain_mesh.opposite_halfface_handle(*chf_it))){
            if()
                hfh = *chf_it;
                break;
            }
        }*/
        if(!hfh.is_valid()){
            log(LogLevel::ERROR)<<" ERROR - failed to find opposite face for cell "<<ch<<std::endl;
            log(LogLevel::INFO)<<" - vertices: "<<std::endl;
            for(int i(0); i<4; i++){
                auto vh = verts[i];
                log(LogLevel::INFO)<<" -- v"<<vh<<" is constrained: "<<constraints_mask[i]<<std::endl;
            }

            /*log(LogLevel::INFO)<<" - edges: "<<std::endl;
            for(auto ce_it = domain_mesh.ce_iter(ch); ce_it.valid(); ce_it++){
                auto from_v = domain_mesh.edge(*ce_it).from_vertex();
                auto   to_v = domain_mesh.edge(*ce_it).to_vertex();

                log(LogLevel::INFO)<<" -- e"<<(*ce_it)<<": ("<<from_v<<", "<<to_v<<") is constrained: "<<domain_mesh.is_boundary(*ce_it)<<std::endl;
            }*/

            exit(EXIT_FAILURE);
        }

        double domain_constrained_area = compute_area(domain_mesh, hfh) * global_scaling_factor * global_scaling_factor;
        double codomain_constrained_area = compute_area(codomain_mesh, hfh);

        double area_ratio = codomain_constrained_area / domain_constrained_area;
        double constrained_scaling_factor = std::sqrt(area_ratio);

        double domain_height = domain_volume / domain_constrained_area;

        if(debug) {
            log(LogLevel::INFO)<<" -   scaled domain constrained area = "<<domain_constrained_area<<std::endl;
            log(LogLevel::INFO)<<" -        codomain constrained area = "<<codomain_constrained_area<<std::endl;
            log(LogLevel::INFO)<<" -     area ratio = "<<area_ratio<<std::endl;
            log(LogLevel::INFO)<<" - scaling factor = "<<constrained_scaling_factor<<std::endl;
            log(LogLevel::INFO)<<" - domain height = "<<domain_height<<std::endl;
        }

        int interior_v_idx(-1);
        for(int i(0); i<4; i++){
            const auto& vh = verts[i];
            //if(domain_mesh.is_boundary(vh)){
            if(constraints_mask[i]){
                P_ref[i] = P_codomain[i];
            }else{
                interior_v_idx = i;
            }
        }

        //log(LogLevel::INFO)<<" - interior vertex has id "<<interior_v_idx<<std::endl;
        //log(LogLevel::INFO)<<" - boundary vertices positions: "<<std::endl;
        auto hf_verts = domain_mesh.get_halfface_vertices(hfh);
        std::vector<Eigen::Vector3d> domain_hf_verts_pos, codomain_hf_verts_pos;
        for(auto v: hf_verts){
            int local_idx = local_idx_prop[v.idx()];
            domain_hf_verts_pos.push_back(P_domain[local_idx]);
            codomain_hf_verts_pos.push_back(P_codomain[local_idx]);
            //log(LogLevel::INFO)<<" -- v"<<v<<" with idx "<<local_idx<<" domain: "<<P_domain[local_idx].transpose()<<", codomain: "<<P_codomain[local_idx].transpose()<<std::endl;
        }

        Vec domain_e1 = (domain_hf_verts_pos[1] - domain_hf_verts_pos[0]).normalized();
        Vec domain_e3 = domain_e1.cross(domain_hf_verts_pos[2] - domain_hf_verts_pos[0]).normalized();
        Vec domain_e2 = domain_e3.cross(domain_e1).normalized();
        Vec int_p = P_domain[interior_v_idx] - domain_hf_verts_pos[0];
        /*log(LogLevel::INFO)<<" - domain coordinate system: "<<std::endl;
        log(LogLevel::INFO)<<" -- e1: "<<domain_e1.transpose()<<std::endl;
        log(LogLevel::INFO)<<" -- e2: "<<domain_e2.transpose()<<std::endl;
        log(LogLevel::INFO)<<" -- e3: "<<domain_e3.transpose()<<std::endl;*/

        Vec domain_int_p = {int_p.dot(domain_e1),
                            int_p.dot(domain_e2),
                            int_p.dot(domain_e3)};
        //log(LogLevel::INFO)<<" - interior vertex relative position: "<<int_p.transpose()<<std::endl;
        //log(LogLevel::INFO)<<" -             in domain coordinates: "<<domain_int_p.transpose()<<std::endl;

        Vec codomain_e1 = (codomain_hf_verts_pos[1] - codomain_hf_verts_pos[0]).normalized();
        Vec codomain_e3 = codomain_e1.cross(codomain_hf_verts_pos[2] - codomain_hf_verts_pos[0]).normalized();
        Vec codomain_e2 = codomain_e3.cross(codomain_e1).normalized();
        /*log(LogLevel::INFO)<<" - codomain coordinate system: "<<std::endl;
        log(LogLevel::INFO)<<" -- e1: "<<codomain_e1.transpose()<<std::endl;
        log(LogLevel::INFO)<<" -- e2: "<<codomain_e2.transpose()<<std::endl;
        log(LogLevel::INFO)<<" -- e3: "<<codomain_e3.transpose()<<std::endl;*/

        Vec codomain_int_p = domain_int_p[0] * codomain_e1
                           + domain_int_p[1] * codomain_e2
                           + domain_int_p[2] * codomain_e3;
        //log(LogLevel::INFO)<<" - interior vertex in codomain coordinates: "<<domain_int_p.transpose()<<std::endl;

        P_ref = P_codomain;
        P_ref[interior_v_idx] = codomain_hf_verts_pos[0] + codomain_int_p;
        //log(LogLevel::INFO)<<" - final interior vertex position: "<<P_ref[interior_v_idx].transpose()<<std::endl;

        Eigen::Matrix3d E_ref;
        for(int i(0); i<3; i++){
            E_ref.col(i) = P_ref[i+1] - P_ref[0];
        }
        double ref_volume = E_ref.determinant() / 6.0;
        double ref_height = ref_volume / codomain_constrained_area;

        if(debug) {
            log(LogLevel::INFO)<<" - updated ref volume = "<<ref_volume<<std::endl;
            log(LogLevel::INFO)<<" - updated ref height = "<<ref_height<<std::endl;
        }

        if(std::abs(ref_height - domain_height) > 1e-3){
            log(LogLevel::ERROR)<<" ERROR - ref height is "<<ref_height<<" but expected to be same as domain height "<<domain_height<<std::endl;
            exit(EXIT_FAILURE);
        }

    }else if(constrained_vertices_count == 2){

        std::vector<std::pair<OVM::VertexHandle, int>> end_vertices;
        for(int i(0); i<4; i++){
            auto vh = verts[i];
            //if(domain_mesh.is_boundary(vh)){
            if(constraints_mask[i]){
                end_vertices.push_back({vh, i});
            }
        }

        double   domain_length = (  domain_mesh.vertex(end_vertices[0].first) -   domain_mesh.vertex(end_vertices[1].first)).norm() * global_scaling_factor;
        double codomain_length = (codomain_mesh.vertex(end_vertices[0].first) - codomain_mesh.vertex(end_vertices[1].first)).norm();
        //log(LogLevel::INFO)<<" -   domain length = "<<domain_length<<std::endl;
        //log(LogLevel::INFO)<<" - codomain length = "<<codomain_length<<std::endl;

        double length_scaling_factor = codomain_length / domain_length;
        //log(LogLevel::INFO)<<" - scaling factor = "<<length_scaling_factor<<std::endl;

        auto domain_mid_point = 0.5 * global_scaling_factor * (domain_mesh.vertex(end_vertices[0].first) + domain_mesh.vertex(end_vertices[1].first));
        //log(LogLevel::INFO)<<" - mid-point at "<<domain_mid_point<<std::endl;

        for(auto end_v: end_vertices){
            P_ref[end_v.second] = vec2vec(domain_mid_point) + length_scaling_factor * (P_domain[end_v.second] - vec2vec(domain_mid_point));
            //log(LogLevel::INFO)<<" -- moved end vertex "<<end_v.first<<" from "<<P_domain[end_v.second].transpose()<<" to "<<P_ref[end_v.second].transpose()<<std::endl;
        }
        //log(LogLevel::INFO)<<" - scaled length = "<<(P_ref[end_vertices[0].second] - P_ref[end_vertices[1].second]).norm()<<std::endl;


    }else{
        //log(LogLevel::INFO)<<" -> no boundary vertices "<<std::endl;
    }
    if(debug) {
        log(LogLevel::INFO)<<" - P_ref: "<<std::endl;
        for(int i(0); i<4; i++) {
            log(LogLevel::INFO) << " -- v" << verts[i] << ": " << P_ref[i].transpose() << std::endl;
        }
        log(LogLevel::INFO)<<" ---------------- "<<std::endl;
    }
}


int setup_mesh_as_one_ring_submesh(const TetrahedralMesh& mesh,
                                    const OpenVolumeMesh::EdgeHandle& eh,
                                    TetrahedralMesh& submesh,
                                    OpenVolumeMesh::EdgeHandle& submesh_eh) {

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

                auto pos = mesh.vertex(v);

                auto local_v = submesh.add_vertex(pos);
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

    submesh_eh = submesh.edge_handle(local_to_split_heh);

    return 0;
}


double deg_to_rad(const double& theta_deg) {
    return (theta_deg / 180) * M_PI;
}


double rad_to_deg(const double& theta_rad) {
    return (theta_rad / M_PI) * 180.0 ;
}


double boundary_size_variation(const TetrahedralMesh& ref_mesh,
                               const TetrahedralMesh& mesh,
                               const OpenVolumeMesh::CellHandle& ch,
                               double global_scaling_factor){


    double boundary_size_variation = 0.0;
    auto verts = ref_mesh.get_cell_vertices(ch);

    std::vector<OVM::VertexHandle> boundary_vertices;
    for (auto v: verts) {
        if(mesh.is_boundary(v)){
            boundary_vertices.push_back(v);
            //log(LogLevel::INFO)<<" -> found boundary vertex "<<v<<" at domain position "<<ref_mesh.vertex(v)<<" and codomain position "<<mesh.vertex(v)<<std::endl;
        }
    }

    if(boundary_vertices.size() == 3){
        //log(LogLevel::INFO)<<" -> element incident to boundary face, checking face area variation"<<std::endl;

        for(auto chf_it = ref_mesh.chf_iter(ch); chf_it.valid(); chf_it++){
            if(ref_mesh.is_boundary(ref_mesh.opposite_halfface_handle(*chf_it))){

                auto hf_verts = ref_mesh.get_halfface_vertices(*chf_it);

                std::vector<Vec3d> ref_P;
                std::vector<Vec3d> P;
                for(auto v: hf_verts){
                    ref_P.push_back(global_scaling_factor * ref_mesh.vertex(v));
                    P.push_back(mesh.vertex(v));
                }

                double ref_mesh_area = 0.5 * (ref_P[1] - ref_P[0]).cross(ref_P[2] - ref_P[0]).norm();
                double     mesh_area = 0.5 * (P[1] - P[0]).cross(P[2] - P[0]).norm();

                //log(LogLevel::INFO)<<" -> ref mesh area = "<<ref_mesh_area<<std::endl;
                //log(LogLevel::INFO)<<" ->     mesh area = "<<mesh_area<<std::endl;

                boundary_size_variation = mesh_area / ref_mesh_area;
            }
        }

    }else if(boundary_vertices.size() == 2){
        //log(LogLevel::INFO)<<" -> element incident to boundary edge, checking edge length variation"<<std::endl;
        for(auto ce_it = ref_mesh.ce_iter(ch); ce_it.valid(); ce_it++) {
            if (ref_mesh.is_boundary(*ce_it)) {
                auto from_v = ref_mesh.edge(*ce_it).from_vertex();
                auto   to_v = ref_mesh.edge(*ce_it).to_vertex();

                double ref_mesh_length = global_scaling_factor * (ref_mesh.vertex(from_v) - ref_mesh.vertex(to_v)).norm();
                double     mesh_length = (mesh.vertex(from_v) - mesh.vertex(to_v)).norm();

                //log(LogLevel::INFO)<<" -> ref mesh length = "<<ref_mesh_length<<std::endl;
                //log(LogLevel::INFO)<<" ->     mesh length = "<<mesh_length<<std::endl;

                boundary_size_variation = mesh_length / ref_mesh_length;

            }
        }
    }else{
        //log(LogLevel::INFO)<<" -> less than two vertices on the boundary, variation is 0"<<std::endl;
    }

    return boundary_size_variation;

}




bool Tutte_embedding(TetrahedralMesh& mesh){

    using namespace Eigen;

    //stores the free vertices' index in the vector above
    auto variable_index_property = mesh.request_vertex_property<int>("variable index");

    int v_idx(0);
    std::vector<OVM::VertexHandle> free_vertices;
    for(auto v: mesh.vertices()){
        if(!mesh.is_boundary(v)){
            variable_index_property[v] = v_idx;
            v_idx++;
            free_vertices.push_back(v);
        }
    }
    log(LogLevel::INFO)<<" - "<<free_vertices.size()<<" / "<<mesh.n_vertices()<<" vertices are free"<<std::endl;


    //setting up K and b
    const size_t variable_count(free_vertices.size());
    SparseMatrix<double> K(variable_count, variable_count);
    VectorXd bx(variable_count);
    VectorXd by(variable_count);
    VectorXd bz(variable_count);

    for(size_t i(0); i<variable_count; i++){

        auto v = free_vertices[i];
        int degree(0);
        double coeff_x(0);
        double coeff_y(0);
        double coeff_z(0);

        for(auto out_he: mesh.outgoing_halfedges(v)){
            auto neighbor = mesh.to_vertex_handle(out_he);

            if(mesh.is_boundary(neighbor)){
                coeff_x += mesh.vertex(neighbor)[0];
                coeff_y += mesh.vertex(neighbor)[1];
                coeff_z += mesh.vertex(neighbor)[2];

            }else /*if(cube_face_property_[neighbor] == face_index)*/{
                //log(LogLevel::INFO)<<"  added it at index "<<variable_index_property[neighbor]<<std::endl;
                K.insert(i, variable_index_property[neighbor]) = -1;
            }
            degree++;

        }
        //log(LogLevel::INFO)<<"--------"<<std::endl;

        K.insert(i, i) = degree;
        bx[i] = coeff_x;
        by[i] = coeff_y;
        bz[i] = coeff_z;
    }

    //NOTE: I have no idea why the rhs is not a matrix as well...
    // But no time to fix this stuff right now and it's not time critical in this context
    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(K);
    if(solver.info()!=Success) {
        log(LogLevel::INFO)<<"decomposition failed !"<<std::endl;
        // decomposition failed
        return false;
    }

    auto x = solver.solve(bx);
    if(solver.info()!=Success) {
        log(LogLevel::INFO)<<"Solving x failed ! "<<std::endl;
        return false;
    }

    auto y = solver.solve(by);
    if(solver.info()!=Success) {
        log(LogLevel::INFO)<<"Solving y failed ! "<<std::endl;
        return false;
    }

    auto z = solver.solve(bz);
    if(solver.info()!=Success) {
        log(LogLevel::INFO)<<"Solving z failed ! "<<std::endl;
        return false;
    }


    //move face vertices at the right locations
    for(size_t i(0); i<variable_count; i++){
        Vec3d new_pos;
        //log(LogLevel::INFO)<<" pos = "<<x[i]<<", "<<y[i]<<std::endl;
        new_pos[0] = x[i];
        new_pos[1] = y[i];
        new_pos[2] = z[i];
        mesh.set_vertex(free_vertices[i], new_pos);
    }

    //log(LogLevel::INFO)<<" done !"<<std::endl;
    return true;
}


void export_elements_with_extreme_angles(const TetrahedralMesh& mesh,
                                         const std::string& output_path,
                                         double angle_thresh_deg){

    log(LogLevel::INFO)<<" ------------------------------------ "<<std::endl;
    log(LogLevel::INFO)<<" - exporting all tets with dihedral angle <"<<angle_thresh_deg<<" or >"<<(180.0 - angle_thresh_deg)<<" to "<<output_path<<std::endl;

    TetrahedralMesh copy_mesh = mesh;

    for(auto ch: mesh.cells()){
        double min_angle = std::numeric_limits<double>::max();
        double max_angle = std::numeric_limits<double>::lowest();

        for(auto ce_it = mesh.ce_iter(ch); ce_it.valid(); ce_it++){
            auto angle = compute_dihedral_angle_deg(mesh, ch, *ce_it);
            min_angle = std::min(min_angle, angle);
            max_angle = std::max(max_angle, angle);
        }

        if(min_angle > angle_thresh_deg && max_angle < (180.0 - angle_thresh_deg)){
            copy_mesh.delete_cell(ch);
        }
    }

    log(LogLevel::INFO)<<" --> exporting "<<copy_mesh.n_logical_cells()<<" elements with extreme angles over "<<copy_mesh.n_cells()<<" in total"<<std::endl;
    copy_mesh.collect_garbage();

    OpenVolumeMesh::IO::FileManager fm;

    fm.writeFile(output_path, copy_mesh);


    log(LogLevel::INFO)<<" ------------------------------------ "<<std::endl;
}


void export_boundary_as_triangle_mesh(const TetrahedralMesh& mesh,
                                      const std::string& output_path){

    log(LogLevel::INFO)<<" - exporting mesh boundary to "<<output_path<<std::endl;

    std::ofstream outfile(output_path);

    TetrahedralMesh copy = mesh;

    //translate_barycenter_to_origin(copy);
    //scale_to_unit_bounding_box_diagonal(copy);

    // Check if the file was opened successfully
    if (outfile.is_open()) {
        outfile << "# vertices" << std::endl;
        for (auto v: mesh.vertices()) {
            outfile << "v";
            for (auto x: copy.vertex(v)) {
                outfile << " " << x;
            }
            outfile << std::endl;
        }
        //no normals

        outfile << "# faces" << std::endl;
        for (auto fh: mesh.faces()) {
            if(mesh.is_boundary(fh)) {
                auto f_verts = mesh.get_halfface_vertices(mesh.halfface_handle(fh, 0));

                outfile << "f";
                for (auto v: f_verts) {
                    outfile << " " << (v.idx() + 1);
                }
                outfile << std::endl;
            }
        }

        outfile.close();
    }else{
        log(LogLevel::ERROR)<<" ERROR - couldn't open file "<<output_path<<" to export boundary mesh"<<
                 std::endl;
        exit(EXIT_FAILURE);
    }
}


void export_one_ring_neighborhood(const TetrahedralMesh& mesh,
                                  const OpenVolumeMesh::VertexHandle& vh,
                                  const std::string& output_path){

    TetrahedralMesh copy = mesh;

    std::unordered_set<OpenVolumeMesh::VertexHandle> one_ring;
    one_ring.insert(vh);

    for(auto vv_it = mesh.vv_iter(vh); vv_it.valid(); vv_it++){
        one_ring.insert(*vv_it);
    }

    for(auto v: mesh.vertices()){
        if(!one_ring.contains(v)){
            copy.delete_vertex(v);
        }
    }



    copy.collect_garbage();

    OVM::IO::FileManager fm;
    fm.writeFile(output_path, copy);

}


void fill_layer_grouping_prop(TetrahedralMesh& mesh,
                            TetrahedralMesh& codomain_mesh, 
                            int axis,
                            int cycle_count,
                            int steps_per_cycle){

    axis = std::min(axis, 2);
    axis = std::max(axis, 0);

    auto color_prop = mesh.request_cell_property<int>("group");
    auto codomain_color_prop = codomain_mesh.request_cell_property<int>("group");
    mesh.set_persistent(color_prop);
    codomain_mesh.set_persistent(codomain_color_prop);

    Vec3d min(1,1,1), max(1,1,1);
    min *= std::numeric_limits<double>::max();
    max *= std::numeric_limits<double>::lowest();

    for(auto vh: mesh.vertices()){
        auto pos = mesh.vertex(vh);
        for(int i(0); i<3; i++){
            min[i] = std::min(min[i], pos[i]);
            max[i] = std::max(max[i], pos[i]);
        }
    }

    for(auto ch: mesh.cells()){
        auto verts = mesh.get_cell_vertices(ch);
        Vec3d centroid(0,0,0);
        for(auto vh: verts){
            centroid += mesh.vertex(vh);
        }
        centroid /= 4.0;

        double t = (centroid[axis] - min[axis]) / (max[axis] - min[axis]);
        //double color = 0.5 + 0.5 * std::sin(cycle_count * t * 2.0 * M_PI);
        double cycle = t * cycle_count;
        double remainder = cycle - std::floor(cycle);
        int color = (int) (remainder * steps_per_cycle);
        color_prop[ch] = color;
        codomain_color_prop[ch] = color_prop[ch];
        log(LogLevel::INFO)<<" - t = "<<t<<", cycle = "<<cycle<<", remainder = "<<remainder<<", color = "<<color<<" -> group: "<< color_prop[ch]<<std::endl;
    }

}

std::string edge_to_str(const TetrahedralMesh& mesh,
                        const OpenVolumeMesh::EdgeHandle& eh){

    std::string str = std::to_string(eh.idx()) + " (" + std::to_string(mesh.edge(eh).from_vertex().idx()) + ", " + std::to_string(mesh.edge(eh).to_vertex().idx()) + ")";
    return str;
}

std::string halfface_to_str(const TetrahedralMesh& mesh,
                            const OpenVolumeMesh::HalfFaceHandle& hfh){

    std::string str;
    auto verts = mesh.get_halfface_vertices(hfh);
    str = std::to_string(hfh.idx())+": ("+std::to_string(verts[0].idx())+ ", "+std::to_string(verts[1].idx())+ ", "+std::to_string(verts[2].idx())+ ")";
    return str;
}

std::string cell_to_str(const TetrahedralMesh& mesh,
                        const OpenVolumeMesh::CellHandle& hfh){

    std::string str;
    auto verts = mesh.get_cell_vertices(hfh);
    str = std::to_string(hfh.idx())+": (" +
            std::to_string(verts[0].idx())+ ", " +
            std::to_string(verts[1].idx())+ ", " +
            std::to_string(verts[2].idx()) + ", " +
            std::to_string(verts[3].idx())+ ")";

    return str;
}

void manual_break(std::string msg){
    log(LogLevel::INFO)<<msg<<std::endl;
    std::string ans;
    std::cin>>ans;
}


void manual_stop(std::string msg){
    log(LogLevel::INFO)<<msg<<std::endl;
    exit(EXIT_FAILURE);
}

bool face_is_blocking(const TetrahedralMesh& mesh,
                      const OpenVolumeMesh::HalfFaceHandle& hfh,
                      const OpenVolumeMesh::VertexHandle& vh,
                      const Vec3d& ray) {

    // Möller-Trumbore algorith, code adapted from wikipedia: https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm

    constexpr double epsilon = std::numeric_limits<double>::epsilon();

    auto p = mesh.vertex(vh);

    auto verts = mesh.get_halfface_vertices(hfh);
    auto a = mesh.vertex(verts[0]);
    auto b = mesh.vertex(verts[1]);
    auto c = mesh.vertex(verts[2]);

    Vec3d edge1 = b - a;
    Vec3d edge2 = c - a;
    Vec3d ray_cross_e2 = ray.cross(edge2);
    double det = edge1.dot(ray_cross_e2);

    if (det > -epsilon && det < epsilon) {
        return false;    // This ray is parallel to this triangle.
    }

    double inv_det = 1.0 / det;
    Vec3d s = p - a;
    double u = inv_det * s.dot(ray_cross_e2);

    if ((u < 0 && abs(u) > epsilon) || (u > 1 && abs(u-1) > epsilon)) {
        return false;
    }

    Vec3d s_cross_e1 = s.cross(edge1);
    double v = inv_det * ray.dot(s_cross_e1);

    if ((v < 0 && abs(v) > epsilon) || (u + v > 1 && abs(u + v - 1) > epsilon)){
        return false;
    }

    // At this stage we can compute t to find out where the intersection point is on the line.
    double t = inv_det * edge2.dot(s_cross_e1);

    if (t > epsilon){ // ray intersection
        return true;
    }else {
        // This means that there is a line intersection but not a ray intersection.
        return false;
    }
}

std::vector<OpenVolumeMesh::EdgeHandle> find_edges_blocking_movement(const TetrahedralMesh &mesh,
                                                                     const OpenVolumeMesh::HalfFaceHandle &hfh,
                                                                     const OpenVolumeMesh::VertexHandle &vh,
                                                                     const Vec3d &d,
                                                                     double epsilon) {

    auto p = mesh.vertex(vh);

    Vec3d c,n;
    compute_centroid_and_normal(mesh, hfh, c, n);

    double t = compute_t_factor(c, n, p, d);

    auto pd = p + t * d;

    auto verts = mesh.get_halfface_vertices(hfh);


    auto incident_ch = mesh.incident_cell(hfh);

    std::vector<OpenVolumeMesh::EdgeHandle> blocking_edges;

    for (int i(0); i<3; i++) {
        int i_1 = (i+1)%3;

        auto vi = mesh.vertex(verts[i]);
        auto vi_1 = mesh.vertex(verts[i_1]);
        auto pi = pd - vi;
        auto ei = vi_1 - vi;

        auto ei_pi = ei.cross(pi);
        if (ei_pi.dot(n) <= epsilon) {
            auto pivot_eh = mesh.edge_handle(mesh.find_halfedge(verts[i], verts[i_1]));
            //blocking_edges.push_back(opposite_edge_in_cell(mesh, incident_ch, pivot_eh));
            blocking_edges.push_back(pivot_eh);
        }
    }

    /*auto v0 = mesh.vertex(verts[0]);
    auto v1 = mesh.vertex(verts[1]);
    auto v2 = mesh.vertex(verts[2]);

    auto p0 = pd - v0;
    auto p1 = pd - v1;
    auto p2 = pd - v2;

    auto e0 = v1 - v0;
    auto e1 = v2 - v1;
    auto e2 = v0 - v2;*/


    return blocking_edges;

}


OpenVolumeMesh::EdgeHandle opposite_edge_in_cell(const TetrahedralMesh &mesh,
                                                 const OpenVolumeMesh::CellHandle &ch,
                                                 const OpenVolumeMesh::EdgeHandle &eh) {

    auto from_v = mesh.edge(eh).from_vertex();
    auto   to_v = mesh.edge(eh).to_vertex();

    auto verts = mesh.get_cell_vertices(ch);

    OVM::VertexHandle op_from_v(-1), op_to_v(-1);
    for(auto v: verts){
        if(v != from_v && v != to_v){
            if(!op_from_v.is_valid()){
                op_from_v = v;
            }else{
                op_to_v = v;
            }
        }
    }
    if(!op_from_v.is_valid() || !op_to_v.is_valid()){
        log(LogLevel::ERROR)<< " ERROR - couldn't find both opposite vertices to edge ("<<from_v<<", "<<to_v<<")"<<std::endl;
        exit(EXIT_FAILURE);
    }

    return mesh.edge_handle(mesh.find_halfedge(op_from_v, op_to_v));
}


std::string resize_str(const std::string str, const int size){

    std::stringstream ss;

    ss << std::setw(size)<<std::left<<str;

    return ss.str();
}