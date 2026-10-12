//
// Created by Valentin Nigolian on 11.12.2024.
//

#define NON_CONFORMING_MESHES_C


#include "NonConformingMeshesT.hh"

#define ENABLE_ALL_CHECKS 0


using namespace OpenVolumeMesh;


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
void
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
initialize_props(double K_threshold){

    //log(LogLevel::INFO)<<" ---- INITIALIZING PROPS"<<std::endl;

    for(auto vh: this->vertices()){
        anchor_position_prop_[vh] = vec2vec(this->vertex(vh));
        initial_anchor_position_prop_[vh] = anchor_position_prop_[vh];
        original_vertex_prop_[vh] = true;
    }


    int projections_count(0);

    int local_cid(0);
    for(auto ch: elements()){

        //log(LogLevel::INFO)<<" --------------------------------------------------- handling element "<<element_to_str(ch)<<": "<<std::endl;

        VertexMatrix P;
        bool projected = setup_position_matrix(ch, K_threshold, P);
        projections_count += projected;

        if(projected){
            //log(LogLevel::INFO)<<" ================== projected cell "<<element_to_str(ch)<<std::endl;
        }

        std::vector<EigenVec> per_vertex_positions;
        for(int local_vid(0); local_vid < TopoDim; local_vid++){
            per_vertex_positions.push_back(P.col(local_vid));
        }

        //log(LogLevel::INFO)<<" - initial positions: "<<std::endl<<P.transpose()<<std::endl;

        // set the properties
        per_vertex_ref_pos_prop(ch) = per_vertex_positions;
        per_vertex_pos_prop(ch)     = per_vertex_ref_pos_prop(ch);

        //log(LogLevel::INFO)<<" -> initial volume = "<<non_conforming_element_volume(ch)<<std::endl;

        local_cid++;
    }
    //log(LogLevel::INFO)<<" **** projected "<<projections_count<<" elements"<<std::endl;

}

template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
void
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
print_mesh_size() {
    log(LogLevel::INFO) << " ------------------------------------------ "<<std::endl;
    log(LogLevel::INFO) << " - mesh size (logical / array size):" << std::endl;
    log(LogLevel::INFO) << " --   vertices: " <<std::setw(6)<<this->n_logical_vertices()<<" / "<<std::setw(6)<<this->n_vertices()<<std::endl;
    log(LogLevel::INFO) << " --      edges: " <<std::setw(6)<<this->n_logical_edges()<<" / "<<std::setw(6)<<this->n_edges()<<std::endl;
    //log(LogLevel::INFO) << " -- half edges: " <<std::setw(6)<<this->n_logical_halfedges()<<" / "<<std::setw(6)<<this->n_halfedges()<<std::endl;
    log(LogLevel::INFO) << " --      faces: " <<std::setw(6)<<this->n_logical_faces()<<" / "<<std::setw(6)<<this->n_faces()<<std::endl;
    //log(LogLevel::INFO) << " -- half faces: " <<std::setw(6)<<this->n_logical_halffaces()<<" / "<<std::setw(6)<<this->n_halffaces()<<std::endl;
    log(LogLevel::INFO) << " --      cells: " <<std::setw(6)<<this->n_logical_cells()<<" / "<<std::setw(6)<<this->n_cells()<<std::endl;
    log(LogLevel::INFO) << " ------------------------------------------ "<<std::endl;
}

template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>&
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
operator=(const NonConformingMeshT& other){


    BaseMesh::operator=(other);

    return *this;
}

template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>&
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
operator=(const BaseMesh& other){

    BaseMesh::operator=(other);

    return *this;
}



template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
void
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
set_current_duplicates_pos_as_ref(){

    log(LogLevel::INFO)<<" - setting current duplicates as reference"<<std::endl;
    for(auto ch: elements()){
        per_vertex_ref_pos_prop(ch) = per_vertex_pos_prop(ch);
    }
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
void
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
snap_boundary_duplicates(){

    log(LogLevel::INFO)<<" - snapping boundary duplicates to the boundary positions"<<std::endl;
    log(LogLevel::INFO)<<" - (Note that this does not apply to the reference positions)"<<std::endl;

    for(auto ch: elements()){
        auto verts = this->get_element_vertices(ch);
        for(int i(0); i<verts.size(); i++){
            const auto& v = verts[i];

            if(is_boundary(v)){
                //log(LogLevel::INFO)<<" - snapped duplicate "<<i<<" of boundary vertex "<<v<<" from "<< per_vertex_pos_prop(ch)[i].transpose()<<" to "<<this->vertex(v)<<std::endl;
                per_vertex_pos_prop(ch)[i] = vec2vec(this->vertex(v));
            }
        }
    }
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
void
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
scale_to_bounding_box(const EigenVec& min_corner,
                      const EigenVec& max_corner,
                      bool keep_aspect_ratio){

    log(LogLevel::INFO)<<" ------------------------------------------- scaling mesh to fit bounding box from "<<min_corner.transpose()<<" to "<<max_corner.transpose()<<std::endl;
    log(LogLevel::INFO)<<" keeping aspect ratio: "<<keep_aspect_ratio<<std::endl;

    EigenVec current_min_corner, current_max_corner;
    current_min_corner.setOnes();
    current_max_corner.setOnes();

    current_min_corner *= std::numeric_limits<double>::max();
    current_max_corner *= std::numeric_limits<double>::lowest();

    for(auto v: this->vertices()){
        const auto& pos = this->vertex(v);
        for(int i(0); i<GeoDim; i++){
            current_min_corner[i] = std::min(current_min_corner[i], pos[i]);
            current_max_corner[i] = std::max(current_max_corner[i], pos[i]);
        }
    }

    log(LogLevel::INFO)<<" - current bounding box is from "<<current_min_corner.transpose()<<" to "<<current_max_corner.transpose();

    EigenVec diagonal = max_corner - min_corner;
    EigenVec current_diagonal = current_max_corner - current_min_corner;

    double min_ratio(std::numeric_limits<double>::max());
    EigenVec ratios;
    for(int i(0); i<GeoDim; i++){
        ratios[i] = diagonal[i] / current_diagonal[i];
        if(ratios[i] < min_ratio){
            min_ratio = ratios[i];
        }
    }

    if(keep_aspect_ratio){
        ratios.setOnes();
        ratios *= min_ratio;
    }

    auto centroid = min_corner + diagonal * 0.5;
    auto current_centroid = current_min_corner + current_diagonal * 0.5;

    log(LogLevel::INFO)<<" - transforming positions as (p - ("<<current_centroid.transpose()<<")) * "<<min_ratio<<" + ("<<centroid.transpose()<<")"<<std::endl;


    for(auto v: this->vertices()){
        auto pos = vec2vec(this->vertex(v));

        pos -= current_centroid; // translate to origin
        pos *= min_ratio;        // scale to bounding box
        pos += centroid;         // and translate to target bounding box centroid

        this->set_vertex(v, vec2vec(pos));
    }

    log(LogLevel::INFO)<<" ------------------------------------------- "<<std::endl;
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
void
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
scale_duplicates_to_bounding_box(const EigenVec& min_corner,
                                 const EigenVec& max_corner,
                                 bool keep_aspect_ratio){

    log(LogLevel::INFO)<<" ------------------------------------------- scaling duplicates to fit bounding box from "<<min_corner.transpose()<<" to "<<max_corner.transpose()<<std::endl;
    log(LogLevel::INFO)<<" keeping aspect ratio: "<<keep_aspect_ratio<<std::endl;

    EigenVec current_min_corner, current_max_corner;
    current_min_corner.setOnes();
    current_max_corner.setOnes();

    current_min_corner *= std::numeric_limits<double>::max();
    current_max_corner *= std::numeric_limits<double>::lowest();

    for(auto e: this->elements()){
        for(int j(0); j<TopoDim; j++){
            const auto& pos = per_vertex_pos_prop(e)[j];
            for(int i(0); i<GeoDim; i++){
                current_min_corner[i] = std::min(current_min_corner[i], pos[i]);
                current_max_corner[i] = std::max(current_max_corner[i], pos[i]);
            }
        }
    }

    log(LogLevel::INFO)<<" - current bounding box is from "<<current_min_corner.transpose()<<" to "<<current_max_corner.transpose();

    EigenVec diagonal = max_corner - min_corner;
    EigenVec current_diagonal = current_max_corner - current_min_corner;

    double min_ratio(std::numeric_limits<double>::max());
    EigenVec ratios;
    for(int i(0); i<GeoDim; i++){
        ratios[i] = diagonal[i] / current_diagonal[i];
        if(ratios[i] < min_ratio){
            min_ratio = ratios[i];
        }
    }

    if(keep_aspect_ratio){
        ratios.setOnes();
        ratios *= min_ratio;
    }

    auto centroid = diagonal * 0.5;
    auto current_centroid = current_diagonal * 0.5;

    log(LogLevel::INFO)<<" - transforming positions as (p - ("<<current_centroid.transpose()<<")) * "<<min_ratio<<" + ("<<centroid.transpose()<<")"<<std::endl;

    for(auto e: this->elements()){
        for(int j(0); j<TopoDim; j++){
            auto& pos = per_vertex_pos_prop(e)[j];

            pos -= current_centroid; // translate to origin
            pos *= min_ratio;        // scale to bounding box
            pos += centroid;         // and translate to target bounding box centroid
        }
    }

    log(LogLevel::INFO)<<" ------------------------------------------- "<<std::endl;
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
void
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
set_protected(const OVM::EdgeHandle& eh){
    protected_edge_prop_[eh] = true;
}

template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
bool
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
is_protected(const OVM::EdgeHandle& eh) const{
    return protected_edge_prop_[eh];
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
OVM::HalfEdgeHandle
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
is_collapse_ok(const OVM::EdgeHandle& eh,
               EigenVec& surviving_v_pos,
               double chebyshev_radius_eps){


    HalfEdgeHandle to_collapse(this->collapse_preserves_boundary(eh));
    bool perform_collapse = false;

    if(to_collapse.is_valid()) {

        //log(LogLevel::INFO)<<" - collapse is valid in direction "<<this->from_vertex_handle(to_collapse)<< " -> "<<this->to_vertex_handle(to_collapse)<<", checking link condition..."<<std::endl;
        bool lc = link_condition(*this, eh);
        //log(LogLevel::INFO)<<" -> lc: "<<lc<<std::endl;
        if (lc) {

            bool boundary_edge = this->is_boundary(eh);
            if(boundary_edge){
                //log(LogLevel::INFO)<<" - boundary but collapsible edge"<<std::endl;
                //simulate the collapse in ref mesh
                auto real_from_v = this->from_vertex_handle(to_collapse);
                auto real_to_v   = this->to_vertex_handle(to_collapse);
                auto from_v_pos  = this->vertex(real_from_v);
                auto to_v_pos    = this->vertex(real_to_v);

                //log(LogLevel::INFO)<<" - original positions: "<<from_v_pos<<" - "<<to_v_pos<<std::endl;

                this->set_vertex(real_from_v, to_v_pos);

                std::vector<bool> incident_element_prop(this->n_elements(), false);
                std::vector<ElemHandle> neighbors;
                this->gather_neighbor_elements(eh, neighbors);
                for(auto neighbor: neighbors){
                    incident_element_prop[neighbor.idx()] = true;
                }

                //check the energy for the updated elements
                bool pos_volume(true);
                std::vector<VertexHandle> end_vertices = {real_from_v, real_to_v};
                for(auto end_v: end_vertices){

                    std::vector<ElemHandle> v_neighbors;
                    this->gather_neighbor_elements(end_v, v_neighbors);

                    for(auto el: v_neighbors) {
                        if(incident_element_prop[el.idx()]){
                            continue;
                        }
                        if(CGAL_wrapper::is_flipped(*this, el) ||
                           CGAL_wrapper::is_degenerate(*this, el)){
                            pos_volume = false;
                            break;
                        }

                        /*EdgeMatrix E;
                        //typename MeshT::ElemHandleT elh(vf_it->idx());
                        this->set_as_shape_matrix(el, E);
                        if(E.determinant() < volume_eps){
                            //log(LogLevel::INFO)<<" --> collapse simulation created inverted element "<<this->element_to_str(el)<<std::endl;
                            pos_volume = false;
                            break;
                        }*/
                    }
                }

                if(pos_volume){
                    //if(from_v == real_from_v) {
                    surviving_v_pos = vec2vec(to_v_pos);
                    /*}else{
                        log(LogLevel::INFO)<<" --> fliperoo!"<<std::endl;
                        surviving_v_pos = vec2vec(from_v_pos);
                    }*/
                    perform_collapse = true;
                }else{
                    this->set_vertex(real_from_v, from_v_pos);
                    //log(LogLevel::INFO)<<" -> moved from vertex back to "<<from_v_pos<<std::endl;
                }
            }else {
                //log(LogLevel::INFO)<<" - interior edge, checking Chebyshev centroid"<<std::endl;
                double cheb_r(0);
                surviving_v_pos = this->solve_edge_Chebyshev(eh, cheb_r);

                if (cheb_r > chebyshev_radius_eps) {
                    //log(LogLevel::INFO) << " - Chebyshev radius = " << cheb_r << ", edge is geometrically collapsible"
                    //          << std::endl;
                    perform_collapse = true;

                } else {
                    //log(LogLevel::INFO) << " - 1-ring of candidate edge is not star-shaped (or kernel is too small, with r = "
                    //          << cheb_r << ")" << std::endl;
                }
            }
        }
    }

    if(!perform_collapse){
        to_collapse = OVM::HalfEdgeHandle(-1);
    }

    return to_collapse;
}



template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
void
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
set_vertex_duplicates(const OVM::VertexHandle& vh, const EigenVec& pos, bool debug){

    std::vector<ElemHandle> neighbors;
    gather_neighbor_elements(vh, neighbors);

    for(auto neighbor: neighbors){
        auto verts = get_element_vertices(neighbor);
        for(int i(0); i<verts.size(); i++){
            if(verts[i] == vh){
                if(debug) {
                    log(LogLevel::INFO)<<" -- moved duplicate "<<i<<" of vertex "<<vh<<" in cell "<<element_to_str(neighbor)<<" from "<<per_vertex_pos_prop(neighbor)[i].transpose()<<" to "<<pos.transpose()<<std::endl;
                }
                per_vertex_pos_prop(neighbor)[i] = pos;
                break;
            }
        }
    }
}




template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
bool
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
check_boundary_duplicates() const {

    log(LogLevel::INFO)<<" ------------------------------------------------------------------ "<<std::endl;
    log(LogLevel::INFO)<<" - checking boundary duplicates..."<<std::endl;
    int wrong_boundary_duplicates_count(0);
    double max_diff_norm(0);
    OVM::VertexHandle max_diff_norm_vertex(-1);

    for(auto elh: this->elements()){

        auto verts = this->get_element_vertices(elh);
        for(int i(0); i<verts.size(); i++){
            auto vh = verts[i];
            if(this->is_boundary(vh)){
                EigenVec pos = vec2vec(this->vertex(vh));
                double diff_norm = (per_vertex_pos_prop(elh)[i] - pos).norm();
                //NOTE: doing exact comparison here because boundary duplicates should regularly be snapped to the
                //  correct positions (given by the boundary conditions of the map)
                if(diff_norm) {
                    log(LogLevel::INFO) << " --> duplicate " << i << " of boundary vertex " << vh << " in cell "<<element_to_str(elh)<<" is at "
                              << per_vertex_pos_prop(elh)[i].transpose() << " but vertex is at " << pos.transpose()
                              << ", diff norm = " << diff_norm << " while both positions should be the exact same"<<std::endl;

                    if(diff_norm > max_diff_norm){
                        max_diff_norm = diff_norm;
                        max_diff_norm_vertex = vh;
                    }
                    wrong_boundary_duplicates_count++;
                }
            }
        }
    }

    log(LogLevel::INFO)<<" --> found "<<wrong_boundary_duplicates_count<<" wrong boundary duplicates"<<std::endl;
    if(wrong_boundary_duplicates_count) {
        log(LogLevel::INFO) << " --> the worst one was " << max_diff_norm_vertex << ", which was off by " << max_diff_norm
                  << std::endl;
    }
    log(LogLevel::INFO)<<" ------------------------------------------------------------------ "<<std::endl;

    return wrong_boundary_duplicates_count == 0;
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
bool
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
contains_degenerate_or_flipped_elements(VERTEX_POSITION_SOURCE pos_source) const {

    for(auto elh: this->elements()){
        if(this->is_flipped_or_degenerate(elh, pos_source)){
            return true;
        }
    }
    return false;
}



template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
QualityVector
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
compute_quality_vector_in_neighborhood(const OVM::VertexHandle& vh,
                                       QUALITY_METRIC metric,
                                       ElementProperty<bool>& to_skip_prop,
                                       VERTEX_POSITION_SOURCE pos_source){

    std::vector<ElemHandle> neighbors;
    gather_neighbor_elements(vh, neighbors);

    QualityVector vec;
    for(auto neighbor: neighbors) {
        if(!to_skip_prop[neighbor]) {
            auto local_vector = compute_quality_vector(neighbor, metric, pos_source);
            vec.merge(local_vector);
        }
    }

    return vec;
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
QualityVector
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
compute_quality_vector_in_neighborhood(const OVM::EdgeHandle& eh,
                                       QUALITY_METRIC metric,
                                       ElementProperty<bool>& to_skip_prop,
                                       VERTEX_POSITION_SOURCE pos_source){

    auto visited_prop = BaseMesh:: template request_property<bool, Elem>("", false);

    std::vector<ElemHandle> neighbors;
    gather_neighbor_elements(this->edge(eh).from_vertex(), neighbors);
    gather_neighbor_elements(this->edge(eh).  to_vertex(), neighbors);

    QualityVector vec;
    for(auto neighbor: neighbors) {
        if(!to_skip_prop[neighbor] && !visited_prop[neighbor]) {
            auto neighbor_vec = compute_quality_vector(neighbor, metric, pos_source);
            vec.merge(neighbor_vec);
            visited_prop[neighbor] = true;
        }
    }

    return vec;
}

template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
QualityVector
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
compute_quality_vector_in_neighborhood(const OVM::FaceHandle& fh,
                                       QUALITY_METRIC metric,
                                       ElementProperty<bool>& to_skip_prop,
                                       VERTEX_POSITION_SOURCE pos_source){

    std::vector<CellHandle> incident_cells = {this->incident_cell(this->halfface_handle(fh, 0)),
                                              this->incident_cell(this->halfface_handle(fh, 1))};

    QualityVector vec;

    for(auto ch: incident_cells){
        if(ch.is_valid() && !to_skip_prop[ch]) {
            auto neighbor_vec = compute_quality_vector(ch, metric, pos_source);
            vec.merge(neighbor_vec);
        }
    }

    return vec;
}

template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
double
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
compute_min_volume_in_neighborhood(const OVM::EdgeHandle& eh){

    std::vector<ElemHandle> neighbors;
    gather_neighbor_elements(this->edge(eh).from_vertex(), neighbors);
    gather_neighbor_elements(this->edge(eh).  to_vertex(), neighbors);

    double min_vol(std::numeric_limits<double>::max());
    for(auto neighbor: neighbors) {
        min_vol = std::min(min_vol, compute_volume(*this, neighbor));
    }

    return min_vol;
}
template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
double
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
compute_min_volume_in_neighborhood(const OVM::EdgeHandle& eh,
                                              const ElementProperty<bool> &to_skip_prop){


    std::vector<ElemHandle> neighbors;
    gather_neighbor_elements(this->edge(eh).from_vertex(), neighbors);
    gather_neighbor_elements(this->edge(eh).  to_vertex(), neighbors);

    double min_vol(std::numeric_limits<double>::max());
    for(auto neighbor: neighbors) {
        if (!to_skip_prop[neighbor]) {
            min_vol = std::min(min_vol, compute_volume(*this, neighbor));
        }
    }

    return min_vol;
}



template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
double
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
compute_min_volume_in_neighborhood(const OVM::VertexHandle& vh){

    auto visited_prop = BaseMesh:: template request_property<bool, Elem>("", false);

    std::vector<ElemHandle> neighbors;
    gather_neighbor_elements(vh, neighbors);

    double min_vol(std::numeric_limits<double>::max());
    for(auto neighbor: neighbors) {
        min_vol = std::min(min_vol, compute_volume(*this, neighbor));
    }

    return min_vol;
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
double
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
compute_min_volume_in_neighborhood(const OVM::FaceHandle& fh){

    return std::min(compute_volume(*this, this->incident_cell(this->halfface_handle(fh, 0))),
                    compute_volume(*this, this->incident_cell(this->halfface_handle(fh, 1))));
}




template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
void
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
update_index_structure(){

    if(enable_anchor_vertices_){
        n_anchor_vertices_ = this->n_logical_vertices();
    }else{
        n_anchor_vertices_ = 0;
    }

    n_vars_ = n_vars_per_cell_ * (this->n_logical_elements()) + GeoDim * n_anchor_vertices_;
    anchor_vertices_index_start_ = n_vars_per_cell_*(this->n_logical_elements());
}



template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
void
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
set_as_non_conforming_shape_matrix(const ElemHandle& elh,
                    NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::EdgeMatrix& E) const {

    //log(LogLevel::INFO)<<" ------------------------------ setting shape matrix for non-conforming elements "<<element_to_str(elh)<<std::endl;

    E.resize(TopoDim-1, TopoDim-1);
    auto c_verts = this->get_element_vertices(elh);
    /*for(auto v: c_verts){
        log(LogLevel::INFO)<<" --- v"<<v<<" pos: "<<this->vertex(v)<<std::endl;
    }*/
    for(int i(0); i<TopoDim-1; i++) {
        EigenVec e = per_vertex_pos_prop(elh)[i + 1] - per_vertex_pos_prop(elh)[0];
        //auto e = this->vertex(c_verts[i+1]) - this->vertex(c_verts[0]);
        //log(LogLevel::INFO)<<" - e"<<i<<": "<<e<<std::endl;
        for(auto j(0); j<GeoDim; j++){
            E(j,i) = e[j];
        }
    }
}

template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
void
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
set_as_shape_matrix(const ElemHandle& ch,
                    NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::EdgeMatrix& E) const {

    //log(LogLevel::INFO)<<" ------------------------------ setting shape matrix for cell "<<ch<<std::endl;

    E.resize(TopoDim-1, TopoDim-1);
    auto c_verts = this->get_element_vertices(ch);
    /*for(auto v: c_verts){
        log(LogLevel::INFO)<<" --- v"<<v<<" pos: "<<this->vertex(v)<<std::endl;
    }*/
    for(int i(0); i<TopoDim-1; i++) {
        auto e = this->vertex(c_verts[i+1]) - this->vertex(c_verts[0]);
        //log(LogLevel::INFO)<<" - e"<<i<<": "<<e<<std::endl;
        for(auto j(0); j<GeoDim; j++){
            E(j,i) = e[j];
        }
    }
    //auto e2 = this->vertex(c_verts[2]) - this->vertex(c_verts[0]);
    //auto e3 = this->vertex(c_verts[3]) - this->vertex(c_verts[0]);

    /*E.col(0)<<e1[0],e1[1],e1[2];
    E.col(1)<<e2[0],e2[1],e2[2];
    E.col(2)<<e3[0],e3[1],e3[2];*/
}

template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
bool
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
setup_position_matrix(const ElemHandle& ch,
                      double K_threshold,
                      NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::VertexMatrix& P){



    EdgeMatrix E;
    set_as_shape_matrix(ch, E);

    //log(LogLevel::INFO)<<" - initial shape matrix: "<<std::endl<<E<<std::endl;
    //log(LogLevel::INFO)<<" - initial volume = "<<(E.determinant()/((TopoDim-1) * (TopoDim-2)))<<std::endl;

    auto c_verts = this->get_element_vertices(ch);
    VertexHandle v0(c_verts[0]);
    int i0(0);
    int interior_count(0);
    for(auto v: c_verts){
        if(!this->is_boundary(v)){
            interior_count++;
        }
    }
    //log(LogLevel::INFO)<<" - contains "<<interior_count<<" interior vertices"<<std::endl;
    // log(LogLevel::INFO)<<" using vertex n°"<<i0<<": "<<v0<<std::endl;

    EigenVec p0 = vec2vec(this->vertex(v0));

    Eigen::MatrixXd J = E * E_ref_inv_;
    //log(LogLevel::INFO)<<" - Er-1: "<<E_ref_inv_<<std::endl;
    //log(LogLevel::INFO)<<" - J: "<<std::endl<<J<<std::endl;

    auto K = condition_number(J);

    //log(LogLevel::INFO)<<" - K = "<<K<<std::endl;

    bool projected(false);

    if(std::isnan(K)){
        //log(LogLevel::INFO)<<" WARNING - condition number is NaN"<<std::endl;
    }

    if(K > K_threshold || std::isnan(K)){
        projected = true;

        //log(LogLevel::INFO)<<" - K = "<<K<<" > "<<K_threshold_<<", projecting..."<<std::endl;


        Eigen::MatrixXd U(GeoDim, GeoDim), V(GeoDim, GeoDim);
        Eigen::VectorXd S(GeoDim);
        //S = optimize_singular_values(S, K_threshold);
        signed_SVD(J, U,S,V);

        //if(this->is_boundary(ch)){
        //simple projection
        for(int i(1); i<GeoDim; i++){
            //log(LogLevel::INFO)<<" -- "<<i<<"-th ratio = "<<(S(0) / S(i))<<std::endl;
            auto ratio = S(0) / S(i);
            if(ratio > K_threshold || ratio < 0){
                S(i) = S(0) / K_threshold;
                //log(LogLevel::INFO)<<" --> updated to "<<S(i)<<std::endl;
            }
        }
        //temp
        S.setOnes();


        if(K == std::numeric_limits<double>::infinity()){
            //log(LogLevel::INFO)<<" ----> found flipped cell "<<std::endl;
        }

        //log(LogLevel::INFO)<<" - updated singular values: "<<S.transpose()<<std::endl;
        //log(LogLevel::INFO)<<" - updated K = "<<(S(0) / S(2))<<std::endl;

        J = U * S.asDiagonal() * V.transpose();
        E = J * E_ref_;

        //log(LogLevel::INFO)<<" - projected shape matrix: "<<std::endl<<E<<std::endl;
        //log(LogLevel::INFO)<<" - projected volume = "<<(E.determinant()/((TopoDim-1) * (TopoDim-2)))<<std::endl;
    }

    P.col(0) = p0;
    for(int i(1); i < TopoDim; i++){
        P.col(i) = E.col(i-1) + p0;
    }
    /*P.col(1) = E.col(0) + p0;
    P.col(2) = E.col(1) + p0;
    P.col(3) = E.col(2) + p0;*/

    return projected;
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
double
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
non_conforming_element_volume(const ElemHandle& ch) const {

    EdgeMatrix E;
    set_as_non_conforming_shape_matrix(ch, E);
    return E.determinant()/((TopoDim-1) * (TopoDim-2));
}

template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
double
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
element_volume(const ElemHandle& ch) const {

    EdgeMatrix E;
    set_as_shape_matrix(ch, E);
    return E.determinant()/((TopoDim-1) * (TopoDim-2));
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
double
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
min_element_volume() const {

    double min_vol(std::numeric_limits<double>::max());
    for(auto e: this->elements()){
        min_vol = std::min(min_vol, element_volume(e));
    }
    return min_vol;
}

template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
double
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
total_volume() const{

    double tot_vol(0);
    for(auto e: this->elements()){
        tot_vol += element_volume(e);
    }

    return tot_vol;
}

template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
Eigen::Vector<double, GeoDim>
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
nonconforming_element_bounding_box_diagonal(const ElemHandle& elh) const{


    EigenVec min_corner = EigenVec::Identity() * std::numeric_limits<double>::max();
    EigenVec max_corner = EigenVec::Identity() * std::numeric_limits<double>::lowest();

    for(int i(0); i<TopoDim; i++){
        auto pos = per_vertex_pos_prop(elh)[i];
        for(int j(0); j<GeoDim; j++){
            min_corner[j] = std::min(min_corner[j], pos[j]);
            max_corner[j] = std::max(max_corner[j], pos[j]);
        }
    }

    return max_corner - min_corner;
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
Eigen::Vector<double, GeoDim>
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
element_bounding_box_diagonal(const ElemHandle& elh) const{


    EigenVec min_corner = EigenVec::Identity() * std::numeric_limits<double>::max();
    EigenVec max_corner = EigenVec::Identity() * std::numeric_limits<double>::lowest();

    auto verts = this->get_element_vertices(elh);
    for(auto v: verts){
        auto pos = this->vertex(v);
        for(int j(0); j<GeoDim; j++){
            min_corner[j] = std::min(min_corner[j], pos[j]);
            max_corner[j] = std::max(max_corner[j], pos[j]);
        }
    }

    return max_corner - min_corner;
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
Eigen::Vector<double, GeoDim>
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
global_bounding_box_diagonal() const{

    EigenVec min_corner = EigenVec::Identity() * std::numeric_limits<double>::max();
    EigenVec max_corner = EigenVec::Identity() * std::numeric_limits<double>::lowest();

    for(auto v: this->vertices()){
        if(this->is_boundary(v)) {
            auto pos = this->vertex(v);
            for (int j(0); j < GeoDim; j++) {
                min_corner[j] = std::min(min_corner[j], pos[j]);
                max_corner[j] = std::max(max_corner[j], pos[j]);
            }
        }
    }

    return max_corner - min_corner;
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
double
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
condition_number(const Eigen::MatrixXd& M,
                 int ignore_n_smallest){

    using namespace Eigen;

    MatrixXd U,V;
    VectorXd S;
    signed_SVD(M, U,S,V);

    ignore_n_smallest = std::max(std::min(ignore_n_smallest, (int)S.size()), 0);

    int smallest_index = (int)S.size()-1-ignore_n_smallest;

    /*log(LogLevel::INFO)<<" --> smallest index = "<<smallest_index<<std::endl;
    log(LogLevel::INFO)<<" U: "<<U<<std::endl;
    log(LogLevel::INFO)<<" S: "<<S.transpose()<<std::endl;
    log(LogLevel::INFO)<<" V: "<<V<<std::endl;*/

    //assumes the sigmas are sorted in descending order
    if(S(smallest_index) < 0){
        return std::numeric_limits<double>::infinity();
    }else{
        return S(0)/S(smallest_index);
    }

    if(S.minCoeff() < 0){
        return std::numeric_limits<double>::infinity();
    }else{
        return S.maxCoeff() / S.minCoeff();
    }
}



template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
double
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
shape_matrix_condition_number(const ElemHandle& ch){

    EdgeMatrix E;
    set_as_shape_matrix(ch, E);

    //log(LogLevel::INFO)<<" - shape matrix: "<<std::endl<<E<<std::endl;

    return condition_number(E);

}



template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
double
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
compute_interior_angle_deg(const OVM::FaceHandle fh,
                           const OVM::VertexHandle vh,
                           VERTEX_POSITION_SOURCE pos_source,
                           const OVM::CellHandle* p_ch) const{

    if(pos_source != CONFORMING && !p_ch){
        log(LogLevel::ERROR)<<" ERROR - requesting non-conforming interior angles but cell handle pointer is null"<<std::endl;
        exit(EXIT_FAILURE);
    }

    auto hf_vertices = this->get_halfface_vertices(this->halfface_handle(fh, 0));

    int vid = hf_vertices[0] == vh ? 0 : (hf_vertices[1] == vh ? 1 : 2);
    int lid = (3 + vid - 1) % 3;//-1 to get the previous one and +3 to avoid the modulo of a negative number
    int rid = (3 + vid + 1) % 3;

    OVM::VertexHandle lvh(hf_vertices[lid]), rvh(hf_vertices[rid]);

    if(p_ch && pos_source != CONFORMING){
        lid = get_duplicate_local_id(*p_ch, lvh);
        rid = get_duplicate_local_id(*p_ch, rvh);
    }

    EigenVec l, r;
    switch(pos_source){
        case CONFORMING:{
            l = vec2vec(this->vertex(hf_vertices[lid]) - this->vertex(hf_vertices[vid]));
            r = vec2vec(this->vertex(hf_vertices[rid]) - this->vertex(hf_vertices[vid]));
            break;
        }
        case NON_CONFORMING:{
            l = this->per_vertex_pos_prop(*p_ch)[lid] - this->per_vertex_pos_prop(*p_ch)[vid];
            r = this->per_vertex_pos_prop(*p_ch)[rid] - this->per_vertex_pos_prop(*p_ch)[vid];
            break;
        }
        case REFERENCE:{
            l = this->per_vertex_ref_pos_prop(*p_ch)[lid] - this->per_vertex_ref_pos_prop(*p_ch)[vid];
            r = this->per_vertex_ref_pos_prop(*p_ch)[rid] - this->per_vertex_ref_pos_prop(*p_ch)[vid];
            break;
        }
    }

    return std::acos(l.normalized().dot(r.normalized()));
}



template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
void
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
setup_ref_shape_matrix(double target_ref_volume){


    VertexMatrix P_ref;
    this->setup_ref_shape_position_matrix(P_ref);

    log(LogLevel::DEBUG)<<" Pr: "<<P_ref<<std::endl;

    for(int i(0); i<GeoDim; i++) {
        E_ref_.col(i) = P_ref.col(i+1) - P_ref.col(0);
    }

    log(LogLevel::DEBUG)<<" Er: "<<E_ref_<<std::endl;

    double det_to_vol_ratio = ((TopoDim-1) * (TopoDim-2));
    log(LogLevel::DEBUG)<<" det to vol ratio: "<<det_to_vol_ratio<<std::endl;
    double reg_tet_vol = E_ref_.determinant() / det_to_vol_ratio;


    log(LogLevel::DEBUG)<<" - regular tet vol = "<<reg_tet_vol<<std::endl;

    double volume_ratio = target_ref_volume / reg_tet_vol;
    double scale_factor = std::pow(volume_ratio, 1./(double)(TopoDim-1));

    E_ref_ *= scale_factor;

    log(LogLevel::DEBUG)<<" - ref vol = " << (E_ref_.determinant() / det_to_vol_ratio)<<std::endl;

    E_ref_inv_ = E_ref_.inverse();

    log(LogLevel::DEBUG)<<" E_ref_inv: "<<std::endl<<E_ref_inv_<<std::endl;
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
void
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
print_duplicates_positions() const{
    int local_elid(0);
    for(auto el: elements()){
        log(LogLevel::INFO)<<" ---------------------------- checking element "<<el<<": ";

        auto c_verts = this->get_element_vertices(el);
        for(auto v: c_verts){
            log(LogLevel::INFO)<<" "<<v;
        }
        log(LogLevel::INFO)<<std::endl;

        int local_vid(0);
        for(auto v: c_verts){
            log(LogLevel::INFO) << " -- duplicate " << get_duplicate_gid(local_elid, local_vid) <<
                         " of vertex " << v <<
                         " is at " << per_vertex_pos_prop(el)[local_vid].transpose()<<
                         " and its reference position is at "<< per_vertex_ref_pos_prop(el)[local_vid].transpose() << std::endl;
            local_vid++;
        }
        local_elid++;
    }
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
int
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
get_duplicate_gid(int elid, int local_vid) const{
    return TopoDim * elid + local_vid;
}

template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
int
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
get_duplicate_gid(ElemHandleT elh, OVM::VertexHandle& vh) const{
    return get_duplicate_gid(elh.idx(), get_duplicate_local_id(elh, vh));
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
int
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
get_duplicate_local_id(ElemHandleT elh, OVM::VertexHandle vh) const{

    auto verts = get_element_vertices(elh);

    int vid = (vh == verts[0] ? 0 :
              (vh == verts[1] ? 1 :
              (vh == verts[2] ? 2 :
              (TopoDim == 3 ? -1 :
              (vh == verts[3] ? 3 : -1)))));

    return vid;
}

template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::EigenVec
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
get_duplicate_position(int duplicate_gid) const {
    ElemHandle elh(duplicate_gid / topo_dim);
    int local_vid = duplicate_gid % topo_dim;

    return per_vertex_pos_prop(elh)[local_vid];
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
void
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
set_duplicate_position(int duplicate_gid, EigenVec new_pos){

    ElemHandle elh(duplicate_gid / topo_dim);
    int local_vid = duplicate_gid % topo_dim;

    per_vertex_pos_prop(elh)[local_vid] = new_pos;
}

/*template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
const NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::EigenVec&
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
get_duplicate_position(int duplicate_gid) const{

    ElemHandle elh(duplicate_gid / topo_dim);
    int local_vid = duplicate_gid % topo_dim;

    return per_vertex_pos_prop(elh)[local_vid];
}*/

template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::EigenVec
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
get_duplicate_position(ElemHandleT elh, OVM::VertexHandle vh) const{

    return per_vertex_pos_prop(elh)[get_duplicate_local_id(elh, vh)];
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
int
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
get_anchor_gid(const VertexHandle& vh){
    return anchor_vertices_index_start_ + GeoDim * (vh.idx());
}

template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
bool
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
initialize_pos_props_from_string_prop() {

    if (!this-> template cell_property_exists<std::string>("str_per_vertex_pos_prop")) {
        return false;
    }

    auto str_per_vertex_pos_prop = this->template request_cell_property<std::string>("str_per_vertex_pos_prop");

    for (auto ch: this->cells()) {

        ElemHandleT elh(ch.idx());
        this->per_vertex_pos_prop(elh).clear();
        this->per_vertex_ref_pos_prop(elh).clear();

        const auto& str_pos = str_per_vertex_pos_prop[ch];
        //std::cout<<" str positions: "<<str_pos<<std::endl;

        std::stringstream ss(str_pos);

        double x, y, z;
        while (ss >> x >> y >> z) {

            EigenVec pos;
            pos << x, y, z;
            //std::cout<<" - duplicate of vertex "<<get_element_vertices(elh)[this->per_vertex_pos_prop(elh).size()]<<" is at "<<std::setprecision(12)<<pos.transpose()<<std::endl;
            this->per_vertex_pos_prop(elh).push_back(pos);
        }
    }

    return true;
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
void
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
serialize_duplicates_pos_prop() {

    log(LogLevel::INFO)<<" --------------------------------------------------------- "<<std::endl;
    log(LogLevel::INFO)<<" ------------------- serializing duplicates position prop"<<std::endl;
    auto str_per_vertex_pos_prop = this->template request_cell_property<std::string>("str_per_vertex_pos_prop");
    this->set_persistent(str_per_vertex_pos_prop);

    for (auto ch: this->cells()) {
        if (!this->is_deleted(ch)) {
            //std::cout<<" ----------------------------------- serializing cell "<<element_to_str(ElemHandleT(ch.idx()))<<std::endl;
            std::string str_pos;
            auto verts = this->get_cell_vertices(ch);
            for (int i(0); i<4; i++) {
                auto v = verts[i];
                auto pos = per_vertex_pos_prop(ElemHandleT(ch.idx()))[i];
                //std::cout<<" - duplicate of vertex "<<v<<" is at ";
                for (int j(0); j<3; j++) {

                    std::ostringstream oss;
                    oss << std::fixed << std::setprecision(100) << pos[j];
                    str_pos += oss.str() + " ";
                    //std::string str_x = std::to_string(pos[j]);
                    //std::cout<<str_x<<" ";
                    //str_pos += std::to_string(pos[j]) + " ";
                }
                //std::cout<<std::endl;
            }
            str_per_vertex_pos_prop[ch] = str_pos;
        }
    }
    log(LogLevel::INFO)<<" --------------------------------------------------------- "<<std::endl;
}

template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
void
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
export_non_conforming_mesh(const std::string& filename,
                       bool filter_based_on_glued_cell_prop,
                       std::vector<bool>* to_keep) {

    log(LogLevel::INFO)<<" --------------------------------------------------------- "<<std::endl;
    log(LogLevel::INFO)<<" ------------------- exporting duplicates mesh to "<<filename<<"..."<<std::endl;

    std::string extension = "";
    size_t dot_loc = filename.find_last_of(".");
    if (dot_loc == std::string::npos) {
        log(LogLevel::INFO)<<" WARNING - couldn't find extension in filename. Using .ovm as default"<<std::endl;
        extension = "ovm";
    }
    extension = filename.substr(dot_loc + 1);

    BaseMesh duplicates_mesh;

    auto cell_prop = duplicates_mesh.template request_vertex_property<int>("element_index");
    auto dup_prop = duplicates_mesh.template request_vertex_property<int>("duplicate_index");
    auto glued_cell_prop = duplicates_mesh.template request_cell_property<int>("last glued cell");
    auto glued_hf_prop = duplicates_mesh. template request_halfface_property<int>("last glued halfface");

    //auto original_face_prop = duplicates_mesh.template request_face_property<int>("original_face_index");

    duplicates_mesh.template set_persistent(cell_prop);
    duplicates_mesh.template set_persistent(dup_prop);
    duplicates_mesh.template set_persistent(glued_cell_prop);
    duplicates_mesh.template set_persistent(glued_hf_prop);
    //duplicates_mesh.template set_persistent(original_face_prop);


    auto last_glued_cell_prop = this-> template request_cell_property<int>("last glued cell");
    auto last_glued_hf_prop = this-> template request_halfface_property<int>("last glued halfface");

    //generic prop
    /*auto to_keep =  std::vector<bool>(n_elements(), false);
    for(auto c: this->elements()){
        if(filter_based_on_glued_cell_prop){

            CellHandle ch(c.idx());
            //only keep glued cells incident to the last glued ones
            if(last_glued_cell_prop[ch] > 0){
                to_keep[ch.idx()] = true;
                auto c_verts = this->get_cell_vertices(ch);
                for(auto v: c_verts){
                    for(auto vc_it = this->vc_iter(v); vc_it.valid(); vc_it++){
                        if(last_glued_cell_prop[*vc_it] >= 0){
                            to_keep[vc_it->idx()] = true;
                        }
                    }
                }
            }
        }else{
            to_keep[c.idx()] = true;
        }
    }*/
    if(to_keep && to_keep->size() != this->n_elements()){
        log(LogLevel::ERROR)<<" ERROR - to_keep property has size "<<to_keep->size()<<" but mesh has "<<this->n_elements()<<" elements"<<std::endl;
        exit(EXIT_FAILURE);
    }


    double min(TopoDim), max(0), avg(0);
    for(auto c: this->elements()){

        if(filter_based_on_glued_cell_prop && last_glued_cell_prop[CellHandle(c.idx())] < 0 ||
                (to_keep && !(*to_keep)[c.idx()])){
            //log(LogLevel::INFO)<<" -> skipping element "<<c<<std::endl;
            continue;
        }
        //log(LogLevel::INFO)<<" ---------------------- handling element "<<this->element_to_str(c)<<std::endl;
        //log(LogLevel::INFO)<<" - deleted: "<<this->is_deleted(c)<<std::endl;

        auto cverts = this->get_element_vertices(c);
        std::vector<VertexHandle> duplicate_vertices;

        //log(LogLevel::INFO)<<" - duplicate vertices: "<<std::endl;
        for(int i(0); i < cverts.size(); i++){
            auto v = cverts[i];
            auto vpos = per_vertex_pos_prop(c)[i];
            auto dup_v = duplicates_mesh.add_vertex(vec2vec(vpos));

            duplicate_vertices.push_back(dup_v);
            cell_prop[dup_v] = c.idx();
            dup_prop[dup_v] = v.idx();

            //log(LogLevel::INFO)<<" "<<dup_v<<" -> "<<v<<",";
        }
        //log(LogLevel::INFO)<<std::endl;

        //if(boundary_vertices_count){
        if(TopoDim == 4) {
            auto dup_c = duplicates_mesh.add_cell(duplicate_vertices);
            glued_cell_prop[dup_c] = last_glued_cell_prop[CellHandle(c.idx())];
            //log(LogLevel::INFO)<<" -- added cell "<<dup_c<<std::endl;
        }else if(TopoDim == 3){
            auto dup_c = duplicates_mesh.add_face(duplicate_vertices);
        }
    }


    //log(LogLevel::INFO)<<" - mesh #elements: "<<this->n_logical_elements()<<std::endl;
    int duplicate_elements_count(0);
    if(TopoDim == 4) {
        duplicate_elements_count = duplicates_mesh.n_logical_cells();
    }else if(TopoDim == 3){
        duplicate_elements_count = duplicates_mesh.n_logical_faces();
    }
    //log(LogLevel::INFO)<<" - dup mesh #elements: "<<duplicate_elements_count<<std::endl;

    /*if(this->n_logical_elements() != duplicate_elements_count){
        log(LogLevel::INFO)<<" ERROR - mesh has "<<this->n_logical_elements()<<" elements but duplicates mesh has "<<duplicate_elements_count<<std::endl;
        exit(EXIT_FAILURE);
    }*/

    /*avg /= this->n_logical_cells();
    log(LogLevel::INFO)<<" - per-cell #boundary vertices (min, max, avg): ("<<min<<", "<<max<<", "<<avg<<")"<<std::endl;*/

    if(extension == "ovm") {
        log(LogLevel::INFO) << " exporting mesh as a .ovm file (with the OVM FileManager)" << std::endl;

        OpenVolumeMesh::IO::FileManager fm;
        fm.writeFile(filename, duplicates_mesh);
    }else if(extension == "obj"){
        log(LogLevel::INFO) << " exporting mesh as a .obj file (with a quick-and-dirty exporter)" << std::endl;

        std::ofstream outfile(filename);

        // Check if the file was opened successfully
        if (outfile.is_open()) {
            outfile<<"# vertices"<<std::endl;
            for(auto v: duplicates_mesh.vertices()){
                outfile<<"v";
                for(auto x: duplicates_mesh.vertex(v)){
                    outfile<<" "<<x;
                }
                outfile<<std::endl;
            }
            //no normals

            outfile<<"# faces"<<std::endl;
            for(auto fh: duplicates_mesh.faces()){
                auto f_verts = duplicates_mesh.get_halfface_vertices(duplicates_mesh.halfface_handle(fh, 0));

                outfile<<"f";
                for(auto v: f_verts){
                    outfile<<" "<<(v.idx() + 1);
                }
                outfile<<std::endl;
            }

            outfile.close();
        }else{
            log(LogLevel::ERROR)<<" ERROR - couldn't open file "<<filename<<std::endl;
        }
    }else{
        log(LogLevel::WARNING)<<" WARNING - exportation to "<<extension<<" file is not supported"<<std::endl;
    }

    log(LogLevel::INFO)<<" --------------------------------------------------------- "<<std::endl;
    log(LogLevel::INFO)<<" --------------------------------------------------------- "<<std::endl;
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
void
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
export_conforming_mesh(const std::string &filename, bool center_and_scale) const {
    log(LogLevel::INFO) << " --------------------------------------------------------- " << std::endl;
    log(LogLevel::INFO) << " ------------------- exporting conforming mesh to " << filename << "..." << std::endl;

    BaseMesh temp_mesh = *this;
    temp_mesh.collect_garbage();

    if(center_and_scale) {
        translate_barycenter_to_origin(temp_mesh);
        scale_to_unit_bounding_box_diagonal(temp_mesh);
    }

    std::string extension = "";
    size_t dot_loc = filename.find_last_of(".");
    if (dot_loc == std::string::npos) {
        log(LogLevel::INFO) << " WARNING - couldn't find extension in filename. Using .ovm as default" << std::endl;
        extension = "ovm";
    }
    extension = filename.substr(dot_loc + 1);


    if (extension == "ovm") {
        log(LogLevel::INFO) << " exporting mesh as a .ovm file (with the OVM FileManager)" << std::endl;

        OpenVolumeMesh::IO::FileManager fm;
        fm.writeFile(filename, temp_mesh);
    } else if (extension == "obj") {
        log(LogLevel::INFO) << " exporting mesh as a .obj file (with a quick-and-dirty exporter)" << std::endl;

        std::ofstream outfile(filename);

        // Check if the file was opened successfully
        if (outfile.is_open()) {
            outfile << "# vertices" << std::endl;
            for (auto v: temp_mesh.vertices()) {
                outfile << "v";
                for (auto x: temp_mesh.vertex(v)) {
                    outfile << " " << x;
                }
                outfile << std::endl;
            }
            //no normals

            outfile << "# faces" << std::endl;
            for (auto fh: temp_mesh.faces()) {
                auto f_verts = temp_mesh.get_halfface_vertices(temp_mesh.halfface_handle(fh, 0));

                outfile << "f";
                for (auto v: f_verts) {
                    outfile << " " << (v.idx() + 1);
                }
                outfile << std::endl;
            }

            outfile.close();
        } else {
            log(LogLevel::ERROR) << " ERROR - couldn't open file " << filename << std::endl;
        }
    } else {
        log(LogLevel::WARNING) << " WARNING - exportation to " << extension << " file is not supported" << std::endl;
    }


    log(LogLevel::INFO)<<" --------------------------------------------------------- "<<std::endl;
    log(LogLevel::INFO)<<" --------------------------------------------------------- "<<std::endl;
}




template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
double
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
non_conforming_cell_condition_number(const ElemHandle& ch){
    auto vpos = per_vertex_pos_prop(ch);
    auto vpos_ref = per_vertex_ref_pos_prop(ch);

    EdgeMatrix E, E_ref;
    for(int i(0); i<GeoDim; i++){
        E.col(i)     = vpos[i+1] - vpos[0];
        E_ref.col(i) = vpos_ref[i+1] - vpos_ref[0];
    }

    EdgeMatrix E_ref_inv = E_ref.inverse();
    EdgeMatrix J = E * E_ref_inv;

    auto K = condition_number(J);
    return K;
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
void
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
print_condition_number_of_non_conforming_cells(){

    log(LogLevel::INFO)<<" '''''''''''''''''''''''''''''''''''''''''''''''''''''''' "<<std::endl;

    double min(std::numeric_limits<double>::max()), max(0), average(0);
    for(auto c: this->elements()){

        //auto K = non_conforming_cell_condition_number(c);
        double K = non_conforming_cell_condition_number(c);
        // log(LogLevel::INFO)<<" - K for cell "<<c<<" = "<<K<<std::endl;

        min = std::min(min, K);
        max = std::max(max, K);
        average += K;
    }

    average /= n_logical_elements();

    log(LogLevel::INFO)<<" - min / max / avg condition number: "<<min<<" / "<<max<<" / "<<average<<std::endl;
    log(LogLevel::INFO)<<" '''''''''''''''''''''''''''''''''''''''''''''''''''''''' "<<std::endl;
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
std::string
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
element_to_str(const ElemHandle& ch) const{
    std::stringstream ss;
    ss<<ch<<": (";
    for(auto v: this->get_element_vertices(ch)){
        ss<<v<<" ";
    }
    ss<<")";

    return ss.str();
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
std::string
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
halfface_to_str(const OVM::HalfFaceHandle& hfh) const{
    std::stringstream ss;
    if (hfh.is_valid()) {
        ss<<hfh<<": (";
        for(auto v: this->get_halfface_vertices(hfh)){
            ss<<v<<" ";
        }
        ss<<")";
    }else {
        ss<<"invalid halfface";
    }

    return ss.str();
}

/*
template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
bool
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
is_boundary(const OVM::HalfFaceHandle& hfh) const{
    return BaseMeshT::is_boundary(hfh);
}


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
bool
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
is_boundary(const OVM::FaceHandle& fh) const{
    return BaseMeshT::is_boundary(fh);
}*/


template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
int
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
glue_elements(double constr_eps, bool update_conforming_mesh){

    using namespace Eigen;


    log(LogLevel::INFO)<<" ------------------------ glueing elements back together "<<std::endl;
    if(update_conforming_mesh){
        log(LogLevel::INFO)<<" --> and updating conforming mesh"<<std::endl;
    }else{
        log(LogLevel::INFO)<<" --> but not updating conforming mesh"<<std::endl;
    }

    auto glued_vertex_prop = this-> template request_vertex_property<bool>();
    auto in_queue_prop     = this-> template request_vertex_property<bool>();
    std::vector<bool> glued_cell_prop(n_elements());

    std::queue<VertexHandle> to_glue_next;

    int interior_element_count(0);
    for(auto c: elements()){
        //log(LogLevel::INFO)<<" ------------------ glueing element "<<element_to_str(c)<<std::endl;
        auto c_verts = get_element_vertices(c);
        int i0(-1);
        int i(0);
        for(auto v: c_verts){
            if(this->is_boundary(v)){
                i0 = i;
                break;
            }
            i++;
        }

        if(i0 < 0){
            //log(LogLevel::INFO)<<" not a boundary tet, skipping"<<std::endl;
            interior_element_count++;
            continue;
        }
        //log(LogLevel::INFO)<<" - found boundary vertex "<<c_verts[i0]<<std::endl;

        auto vpos_init = this->vertex(c_verts[i0]);
        //log(LogLevel::INFO)<<" - which initially was at "<<vpos_init<<std::endl;

        EigenVec p0 = per_vertex_pos_prop(c)[i0];
        EigenVec p0p = vec2vec(vpos_init);//{vpos_init[0], vpos_init[1], vpos_init[2]};
        //log(LogLevel::INFO)<<" - and duplicate is at "<<p0.transpose()<<std::endl;

        EigenVec t = p0p - p0;

        //log(LogLevel::INFO)<<" - thus moving cell by "<<t.transpose()<<std::endl;

        i = 0;
        for(auto v: c_verts){
            auto& vpos = per_vertex_pos_prop(c)[i];
            vpos += t;
            //log(LogLevel::INFO)<<" - moved local duplicate n°"<<i<<" v"<<v<<" from "<<(vpos-t).transpose()<<" to "<<vpos.transpose()<<std::endl;
            if(update_conforming_mesh) {
                this->set_vertex(v, vec2vec(vpos));//{vpos[0], vpos[1], vpos[2]});
            }
            i++;
            if(!this->is_boundary(v) && !in_queue_prop[v]){
                to_glue_next.push(v);
                in_queue_prop[v] = true;
            }
            glued_vertex_prop[v] = true;
        }

        /*log(LogLevel::INFO)<<" - glued boundary tet "<<c<<": "<<cell_to_str(c)<<" to :"<<std::endl;
        for(auto v:c_verts){
            log(LogLevel::INFO)<<"  -- "<<v<<" at "<<mesh_.vertex(v)<<std::endl;
        }
        log(LogLevel::INFO)<<" -------------------- "<<std::endl;*/

        glued_cell_prop[c.idx()] = true;
    }

    //log(LogLevel::INFO)<<" ---- checking boundary glueing..."<<std::endl;
    int bad_boundary_count(0);
    int local_cid(0);
    for(auto c: elements()){

        auto c_verts = get_element_vertices(c);
        int i(0);
        for(auto v: c_verts){
            EigenVec diff = per_vertex_pos_prop(c)[i] - vec2vec(this->vertex(v));
            //Vector3d vpos = vec2vec(mesh_.vertex(v));
            /*log(LogLevel::INFO)<<" - dpos = "<<per_vertex_pos_prop_[c][i].transpose()<<std::endl;
            log(LogLevel::INFO)<<" - vpos = "<<vpos.transpose()<<std::endl;
            log(LogLevel::INFO)<<" - diff norm for duplicate "<<get_duplicate_gid(local_cid, i)<<" of vertex "<<v<<": "<<diff.norm()<<std::endl;*/
            if(this->is_boundary(v) && (diff).norm() > constr_eps){
                bad_boundary_count++;
                //log(LogLevel::INFO)<<" - duplicate "<<get_duplicate_gid(local_cid, i)<<" of vertex "<<v<<" in cell "<<c<<" is at "<<per_vertex_pos_prop_[c][i].transpose()<<" and mesh vertex is at "<<mesh_.vertex(v)<<", diff = "<<diff.norm()<<std::endl;
            }
            i++;
        }
        local_cid++;
    }


#if ENABLE_ALL_CHECKS
    for(auto c: mesh_.cells()){
        if(glued_cell_prop[c.idx()] && OVMtetToCGALtet(mesh_, c).orientation() != CGAL::POSITIVE){
            log(LogLevel::INFO)<<" ERROR - boundary tet "<<c<<": "<<cell_to_str(c)<<" is flipped"<<std::endl;
            return -1;
        }
    }
#endif

    //log(LogLevel::INFO)<<" interior elements count: "<<interior_element_count<<std::endl;
    if(bad_boundary_count){
        log(LogLevel::INFO)<<" WARNING: found "<<bad_boundary_count<<" bad boundary duplicates"<<std::endl;
        return 1;
    }else{
        log(LogLevel::INFO)<<" --> all boundary duplicates are positioned correctly"<<std::endl;
    }

    if(!interior_element_count){
        log(LogLevel::INFO)<<" and no interior element, all good"<<std::endl;
        return 0;
    }

    //then glue the rest

    //log(LogLevel::INFO)<<" glueing interior tets..."<<std::endl;
    //log(LogLevel::INFO)<<" - "<<to_glue_next.size()<<" vertices to glue at first"<<std::endl;
    int max_it(n_logical_elements());
    int it(0);
    int glued_count(0);
    while(!to_glue_next.empty() && it++ < max_it){
        auto to_glue = to_glue_next.front();
        to_glue_next.pop();

        //log(LogLevel::INFO)<<" ------------------- glueing cells around vertex "<<to_glue<<std::endl;
        //not optimal but the price to pay for genericity
        std::vector<ElemHandle> vertex_neighbors;
        gather_neighbor_elements(to_glue, vertex_neighbors);
        //log(LogLevel::INFO)<<" - found "<<vertex_neighbors.size()<<" neighboring elements"<<std::endl;

        for(auto c: vertex_neighbors){
            if(!glued_cell_prop[c.idx()]){
                //log(LogLevel::INFO)<<" -- glueing incident cell "<<element_to_str(c)<<std::endl;
                auto c_verts = get_element_vertices(c);
                EigenVec t;
                bool found_glued(false);
                for(int i(0); i<TopoDim; i++){
                    auto v = c_verts[i];
                    //log(LogLevel::INFO)<<" -- checking vertex "<<v<<", glued: "<<glued_vertex_prop[v]<<", in queue: "<<in_queue_prop[v]<<", found glued: "<<found_glued<<std::endl;
                    if(glued_vertex_prop[v] && !found_glued){
                        //log(LogLevel::INFO)<<" -- found already glued vertex "<<v<<std::endl;
                        t = vec2vec(this->vertex(v)) - per_vertex_pos_prop(c)[i];
                        //log(LogLevel::INFO)<<" -- t: "<<t.transpose()<<std::endl;
                        found_glued = true;
                    }else if(!in_queue_prop[v]){
                        to_glue_next.push(v);
                        in_queue_prop[v] = true;
                        //log(LogLevel::INFO)<<" -- added vertex "<<c_verts[i]<<" to queue"<<std::endl;
                    }
                }
                if(!found_glued){
                    log(LogLevel::ERROR)<<" ERROR - couldn't find glued vertex in element "<<element_to_str(c)<<std::endl;
                    return -1;
                }

                for(int i(0); i<TopoDim; i++){
                    auto& vpos = per_vertex_pos_prop(c)[i];
                    vpos += t;
                    //log(LogLevel::INFO)<<" - moved local duplicate n°"<<i<<" v"<<c_verts[i]<<" from "<<(vpos-t).transpose()<<" to "<<vpos.transpose()<<std::endl;
                    if(update_conforming_mesh) {
                        this->set_vertex(c_verts[i], vec2vec(vpos));
                    }
                    glued_vertex_prop[c_verts[i]] = true;
                }


#if ENABLE_ALL_CHECKS
                for(auto v: c_verts){
                    bool found_flipped(false);
                    for(auto vc_it = mesh_.vc_iter(v); vc_it.valid(); vc_it++){
                        auto cgal_tet = OVMtetToCGALtet(mesh_, *vc_it);
                        if(glued_cell_prop[*vc_it] && cgal_tet.orientation() != CGAL::POSITIVE){
                            log(LogLevel::INFO)<<" ERROR - moving vertex "<<v<<" to "<<mesh_.vertex(v)<<" flipped tet "<<cell_to_str(*vc_it)<<std::endl;
                            found_flipped = true;
                        }
                    }
                    if(found_flipped){
                        return -1;
                    }
                }

                if(OVMtetToCGALtet(mesh_, *vc_it).orientation() != CGAL::POSITIVE){
                    log(LogLevel::INFO)<<" ERROR - glued tet is flipped "<<std::endl;
                    return -1;
                }
#endif

                glued_cell_prop[c.idx()] = true;
                glued_count++;
            }
        }
    }
    log(LogLevel::INFO)<<" --> glued "<<glued_count<<" / "<<interior_element_count<<" interior tets"<<std::endl;

    if(glued_count != interior_element_count){
        log(LogLevel::ERROR)<<" ERROR - only glued "<<glued_count<<"/"<<interior_element_count<<" interior tets"<<std::endl;
    }
    for(auto c: elements()){
        if(!glued_cell_prop[c.idx()]){
            log(LogLevel::ERROR)<<" ERROR - element "<<c<<": "<<element_to_str(c)<<" was not glued"<<std::endl;
        }
    }

#if ENABLE_ALL_CHECKS
    if(BadTetFinder::meshContainsFlippedTets(mesh_)){
        log(LogLevel::INFO)<<" ERROR - mesh contains flipped tets after glueing tets"<<std::endl;
        auto bad_tets = BadTetFinder::findBadTets(mesh_);
        log(LogLevel::INFO)<<" - found "<<bad_tets.degenerate.size()<<" degenerate tets and "<<bad_tets.flipped.size()<<" flipped tets"<<std::endl;
        for(auto c: bad_tets.flipped){
            log(LogLevel::INFO)<<" - "<<c<<": "<<cell_to_str(c)<<" has volume "<<OVMtetToCGALtet(mesh_, c).volume()<<std::endl;
        }
    }
#endif

    return 0;

}




template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
void
NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>::
print_full_mesh() const{
    log(LogLevel::INFO)<<" ************************ vertices ("<<this->n_logical_vertices()<<"/"<<this->n_vertices()<<"):"<<std::endl;
    for(auto v: this->vertices()) {
        log(LogLevel::INFO) << " - v" << v << " at " << this->vertex(v) << std::endl;
    }

    log(LogLevel::INFO)<<" ************************ edges ("<<this->n_logical_edges()<<"/"<<this->n_edges()<<"):"<<std::endl;
    for(auto e: this->edges()) {
        log(LogLevel::INFO) << " - e" << e << ": " << this->edge(e).from_vertex()<<"-"<<this->edge(e).to_vertex() << std::endl;
    }

    log(LogLevel::INFO)<<" ************************ faces ("<<this->n_logical_faces()<<"/"<<this->n_faces()<<"):"<<std::endl;
    for(auto f: this->faces()) {
        log(LogLevel::INFO) << " - f" << f << ":";
        auto hf = this->halfface_handle(f, 0);
        for(auto v: this->get_halfface_vertices(hf)){
            log(LogLevel::INFO)<<" "<<v;
        }
        log(LogLevel::INFO)<<std::endl;
    }

    log(LogLevel::INFO)<<" ************************ cells ("<<this->n_logical_cells()<<"/"<<this->n_cells()<<"):"<<std::endl;
    for(auto c: this->cells()) {
        log(LogLevel::INFO) << " - c" << c << ":";
        for(auto v: this->get_cell_vertices(c)){
            log(LogLevel::INFO)<<" "<<v;
        }
        log(LogLevel::INFO)<<std::endl;
    }
}


#if 0
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
#endif
