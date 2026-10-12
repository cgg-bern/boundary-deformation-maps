//
// Created by Valentin Nigolian on 11.12.2024.
//

#pragma once

#include "Helpers.hh"
#include "TetMesh.hh"

namespace OVM = OpenVolumeMesh;

enum VERTEX_POSITION_SOURCE {CONFORMING, NON_CONFORMING, REFERENCE};


struct MeshAngles{
    double min_dihedral_angle_deg = std::numeric_limits<double>::max();
    double max_dihedral_angle_deg = std::numeric_limits<double>::lowest();

    double min_interior_angle_deg = std::numeric_limits<double>::max();
    double max_interior_angle_deg = std::numeric_limits<double>::lowest();
};


#warning TODO: rename "shape matrix" occurences for "edge matrix"

template<int GeoDim, int TopoDim, typename BaseMesh, typename ElemIter, typename ElemHandle, typename Elem>
class NonConformingMeshT : public BaseMesh {

public:

    static_assert(GeoDim == 2 || GeoDim == 3, "non-conforming mesh is only implemented to work with embeddings in 2 or 3-dimensional space");
    static_assert(TopoDim == 3 || TopoDim == 4, "non-conforming mesh is only implemented to work with triangle or tetrahedral meshes");

    using EigenVec = Eigen::Vector<double, GeoDim>;
    using EdgeMatrix = Eigen::Matrix<double, GeoDim, GeoDim>;
    using VertexMatrix = Eigen::Matrix<double, GeoDim, TopoDim>;
    using Jacobian = Eigen::Matrix<double, GeoDim, GeoDim>;
    using ElemHandleT = ElemHandle;
    using ElemT = Elem;
    template <typename T> using ElementProperty = OVM::PropertyPtr<T, Elem>;

    static constexpr int geo_dim  = GeoDim;
    static constexpr int topo_dim = TopoDim;

    using PointT = typename BaseMesh::PointT;
    using BaseMeshT = BaseMesh;

    using ElementMapT = std::pair<ElemHandle, std::map<OVM::VertexHandle, std::pair<EigenVec, EigenVec>>>;

    NonConformingMeshT(bool enable_anchor_vertices = false) :
            enable_anchor_vertices_(enable_anchor_vertices),
            per_vertex_pos_prop_(BaseMesh:: template request_cell_property<std::vector<EigenVec>>()),
    per_vertex_ref_pos_prop_(BaseMesh:: template request_cell_property<std::vector<EigenVec>>()),
    anchor_position_prop_(BaseMesh:: template request_vertex_property<EigenVec>()),
    initial_anchor_position_prop_(BaseMesh:: template request_vertex_property<EigenVec>()),
    original_vertex_prop_(BaseMesh:: template request_vertex_property<bool>()),
    protected_edge_prop_(BaseMesh:: template request_edge_property<bool>())
    {}


    NonConformingMeshT(const NonConformingMeshT<GeoDim, TopoDim, BaseMesh, ElemIter, ElemHandle, Elem>& other) :
            enable_anchor_vertices_(other.enable_anchor_vertices_),
            per_vertex_pos_prop_(BaseMesh:: template request_cell_property<std::vector<EigenVec>>()),
    per_vertex_ref_pos_prop_(BaseMesh:: template request_cell_property<std::vector<EigenVec>>()),
    anchor_position_prop_(BaseMesh:: template request_vertex_property<EigenVec>()),
    initial_anchor_position_prop_(BaseMesh:: template request_vertex_property<EigenVec>()),
    original_vertex_prop_(BaseMesh:: template request_vertex_property<bool>()),
    protected_edge_prop_(BaseMesh:: template request_edge_property<bool>())
    {
        *this = other;
    }


    /*NonConformingMeshT(const BaseMesh& other) : NonConformingMeshT(false)
    {
        *this = other;
        *this = BaseMesh::operator=(other);
    }*/

    NonConformingMeshT& operator=(const NonConformingMeshT& other);

    NonConformingMeshT& operator=(const BaseMesh& other);

    void update_index_structure();

    /** Does the following:
     *  1) project bad elements (with conditioning > K_thresh
     *  2) initialize the duplicates (potentially with projected positions)
     *  @return 0: ok, 1: invalid element found
  */
    void initialize_props(double K_thresh = std::numeric_limits<double>::max());

    void print_mesh_size();


    /*void print_condition_number_of_non_conforming_cells();


    int glue_elements_and_update_mesh(double constr_eps);

    int check_edge_glueing_constraints(double& max_diff_norm,
                                       bool only_check_boundary = false);*/

#warning TODO: rename those for "compute_..." to make it clear there's a cost to the function call
    double non_conforming_element_volume(const ElemHandle& ch) const;

    double element_volume(const ElemHandle& ch) const;

    double min_element_volume() const;

    double total_volume() const;

    EigenVec nonconforming_element_bounding_box_diagonal(const ElemHandle& elh) const;

    EigenVec element_bounding_box_diagonal(const ElemHandle& elh) const;

    EigenVec global_bounding_box_diagonal() const;

    double condition_number(const Eigen::MatrixXd& M,
                            int ignore_n_smallest = 0);


    double shape_matrix_condition_number(const ElemHandle& ch);


    double compute_interior_angle_deg(const OVM::FaceHandle fh,
                                      const OVM::VertexHandle vh,
                                      VERTEX_POSITION_SOURCE pos_source,
                                      const OVM::CellHandle* p_ch = nullptr) const;

    /* Projects the edge matrix if needed
     * @return whether or not projection was needed */
    bool setup_position_matrix(const ElemHandle& ch,
                               double K_thresh,
                               VertexMatrix& P);

    void set_as_non_conforming_shape_matrix(const ElemHandle& ch,
                                            EdgeMatrix& E) const;

    void set_as_shape_matrix(const ElemHandle& ch,
                             EdgeMatrix& E) const;


    void setup_ref_shape_matrix(double ref_volume);

    virtual bool initialize_pos_props_from_string_prop();

    virtual void serialize_duplicates_pos_prop();

    virtual void export_non_conforming_mesh(const std::string &filename,
                                            bool filter_based_on_glued_cell_prop = false,
                                            std::vector<bool> *to_keep = nullptr);

    virtual void export_conforming_mesh(const std::string& filename, bool center_and_scale = true) const;

    int n_vars() const{
        return n_vars_;
    }

    void print_duplicates_positions() const;

    int get_anchor_gid(const OVM::VertexHandle& vh);

    int get_duplicate_gid(ElemHandleT elh, OVM::VertexHandle& vh) const;

    int get_duplicate_gid(int elid, int local_vid) const;

    int get_duplicate_local_id(ElemHandleT elh, OVM::VertexHandle vh) const;

    EigenVec get_duplicate_position(int duplicate_gid) const;

    void set_duplicate_position(int duplicate_gid, EigenVec new_pos);

    //const EigenVec& get_duplicate_position(int duplicate_gid) const;

    EigenVec get_duplicate_position(ElemHandleT elh, OVM::VertexHandle vh) const;

    void print_condition_number_of_non_conforming_cells();

    double non_conforming_cell_condition_number(const ElemHandle& ch);

    void set_current_duplicates_pos_as_ref();

    void snap_boundary_duplicates();

    void scale_to_bounding_box(const EigenVec& min_corner,
                               const EigenVec& max_corner,
                               bool keep_aspect_ratio = true);


    void scale_duplicates_to_bounding_box(const EigenVec& min_corner,
                                          const EigenVec& max_corner,
                                          bool keep_aspect_ratio = true);

    void set_protected(const OVM::EdgeHandle& eh);
    bool is_protected(const OVM::EdgeHandle& eh) const;

    OVM::HalfEdgeHandle is_collapse_ok(const OVM::EdgeHandle& eh,
                                       EigenVec& surviving_v_pos,
                                       double chebyshev_radius_eps = 1e-6);

    void set_vertex_duplicates(const OVM::VertexHandle& vh, const EigenVec& pos, bool debug = false);

    bool check_boundary_duplicates() const;

    bool contains_degenerate_or_flipped_elements(VERTEX_POSITION_SOURCE pos_source) const;


    QualityVector compute_quality_vector_in_neighborhood(const OVM::VertexHandle &vh,
                                                         QUALITY_METRIC metric,
                                                         ElementProperty<bool> &to_skip_prop,
                                                         VERTEX_POSITION_SOURCE pos_source);

    /**
     * NOTE: checks all tets incident to either end-vertices of the edge.
     * Set the to_skip_prop accordingly if you only want to check the cells incident to the edge
     * */
    QualityVector compute_quality_vector_in_neighborhood(const OVM::EdgeHandle &eh,
                                                         QUALITY_METRIC metric,
                                                         ElementProperty<bool> &to_skip_prop,
                                                         VERTEX_POSITION_SOURCE pos_source);



    QualityVector compute_quality_vector_in_neighborhood(const OVM::FaceHandle& fh,
                                                         QUALITY_METRIC metric,
                                                         ElementProperty<bool>& to_skip_prop,
                                                         VERTEX_POSITION_SOURCE pos_source);

    double compute_min_volume_in_neighborhood(const OVM::VertexHandle &vh);

    double compute_min_volume_in_neighborhood(const OVM::EdgeHandle &eh);

    double compute_min_volume_in_neighborhood(const OVM::EdgeHandle &eh,
                                              const ElementProperty<bool> &to_skip_prop);


    double compute_min_volume_in_neighborhood(const OVM::FaceHandle &fh);


    ElemHandleT to_element_handle(const OVM::CellHandle& ch){
        return ElemHandleT(ch.idx());
    }
    ElemHandleT to_element_handle(const OVM::FaceHandle& fh){
        return ElemHandleT(fh.idx());
    }

    virtual void setup_boundary_adapted_init_positions(const ElemHandle& ch) = 0;

    virtual void set_regular_element_as_reference(double vol = 1.0) = 0;

    virtual int check_edge_glueing_constraints(double glueing_eps,
                                               double& max_diff_norm,
                                               bool only_check_boundary = false) const = 0;


    virtual int check_edge_holonomy(const OVM::EdgeHandle eh,
                                    double constr_eps,
                                    double& diff_to_two_pi,
                                    bool print_debug = false) = 0;

    /** glueing eps is used to check if a face has been glued or not */
    virtual int check_holonomy(double constr_eps) = 0;

    int glue_elements(double constr_eps, bool update_conforming_mesh);


    virtual OVM::VertexHandle split_edge_and_update_props(const OVM::EdgeHandle& to_split) = 0;


    // returns the valid collapse direction, or HalfEdgeHandle(-1) if there is none.
    // What a "valid" direction means depends on whether this is a 2D or 3D mesh
    virtual OVM::HalfEdgeHandle collapse_preserves_boundary(OVM::EdgeHandle eh) = 0;

    /* returns the surviving vertex */
    virtual OVM::VertexHandle collapse_edge_and_update_props(const OVM::HalfEdgeHandle& to_collapse) = 0;

    /* planar angles for 2D meshes, and solid angles for 3D ones
     * Warning: assumes that the property is zero-initialized */
    virtual void fill_angles_sum_property(OVM::VertexPropertyT<double>& angles_sum_prop) = 0;


    /** NOTE: uses the conforming mesh */
    virtual EigenVec solve_edge_Chebyshev(const OVM::EdgeHandle eh,
                                          double& radius,
                                          bool print_debug = 0) = 0;


    virtual bool is_flipped_or_degenerate(const ElemHandle& elh,
                                          VERTEX_POSITION_SOURCE pos_source) const = 0;


#if 0
    using ElementMapT = std::pair<typename MeshT::ElemHandleT,
            std::map<OVM::VertexHandle,
                    std::pair<typename MeshT::EigenVec, typename MeshT::EigenVec>>>;
     /** Moves mid-vertex at from_pos + t * (to_pos - from_pos) in each
     *  non-conforming cell and updates the props accordingly */
    void update_props_around_mid_vertex(const OVM::VertexHandle& from_v,
                                        const OVM::VertexHandle& mid_v,
                                        const OVM::VertexHandle& to_v,
                                        const double t,
                                        std::vector<ElementMapT>& original_cell_prop);

#endif

#warning TODO: rename for "gather_incident_elements"
    virtual void gather_neighbor_elements(const OVM::VertexHandle& vh,
                                          std::vector<ElemHandle>& neighbors) const = 0;

    virtual void gather_neighbor_elements(const OVM::EdgeHandle& eh,
                                          std::vector<ElemHandle>& neighbors) const = 0;


    virtual MeshAngles compute_angles_in_neighborhood(const OVM::VertexHandle& vh,
                                                      ElementProperty<bool>& to_skip_prop,
                                                      VERTEX_POSITION_SOURCE pos_source) const = 0;

    /* NOTE: check angles in union of the end vertices' neighborhoods */
    virtual MeshAngles compute_angles_in_neighborhood(const OVM::EdgeHandle& eh,
                                                      ElementProperty<bool>& to_skip_prop,
                                                      VERTEX_POSITION_SOURCE pos_source) const = 0;

    /** NOTE: only checks angles in incident cells */
    virtual MeshAngles compute_angles_in_neighborhood(const OVM::FaceHandle& fh,
                                                      VERTEX_POSITION_SOURCE pos_source) const = 0;

    virtual double compute_dihedral_angle_deg(const ElemHandle,
                                              const OVM::EdgeHandle& eh,
                                              VERTEX_POSITION_SOURCE pos_source) const = 0;



    virtual QualityVector compute_quality_vector(const ElemHandle& elh,
                                                       QUALITY_METRIC metric,
                                                       VERTEX_POSITION_SOURCE pos_source) const = 0;



    std::string element_to_str(const ElemHandle& ch) const;

    std::string halfface_to_str(const OVM::HalfFaceHandle& hfh) const;

    EigenVec& anchor_position_prop(const OVM::VertexHandle& vh){
        return anchor_position_prop_[vh];
    }
    const EigenVec& anchor_position_prop(const OVM::VertexHandle& vh) const{
        return anchor_position_prop_[vh];
    }

    EigenVec& initial_anchor_position_prop(const OVM::VertexHandle& vh){
        return initial_anchor_position_prop_[vh];
    }
    const EigenVec& initial_anchor_position_prop(const OVM::VertexHandle& vh) const{
        return initial_anchor_position_prop_[vh];
    }

    void print_full_mesh() const;

    bool original_vertex_prop(OVM::VertexHandle vh) const{
        return original_vertex_prop_[vh];
    }

    void set_as_original_vertex(OVM::VertexHandle vh, bool value){
        original_vertex_prop_[vh] = value;
    }



    bool is_boundary(const OVM::HalfFaceHandle& hfh) const{
        return BaseMesh::is_boundary(hfh);
    }


    bool is_boundary(const OVM::FaceHandle& fh) const{
        return BaseMesh::is_boundary(fh);
    }

    // ----------------------------------------------------------------------------------------- virtual methods

    virtual std::vector<EigenVec>& per_vertex_pos_prop(const ElemHandle& ch) = 0;
    virtual std::vector<EigenVec>  per_vertex_pos_prop(const ElemHandle& ch) const = 0;

    virtual std::vector<EigenVec>& per_vertex_ref_pos_prop(const ElemHandle& ch) = 0;
    virtual std::vector<EigenVec>  per_vertex_ref_pos_prop(const ElemHandle& ch) const = 0;

   /* virtual std::vector<bool>& glued_duplicates_prop(const ElemHandle& ch) = 0;
    virtual std::vector<bool>  glued_duplicates_prop(const ElemHandle& ch) const = 0;*/



    virtual void setup_ref_shape_position_matrix(VertexMatrix& P_ref) = 0;


    virtual std::vector<OVM::VertexHandle> get_element_vertices(ElemHandle ch) const = 0;

    virtual std::pair<ElemIter, ElemIter> elements() const = 0;

    virtual int n_logical_elements() const = 0;
    virtual int n_elements() const = 0;

    virtual bool is_boundary(const OVM::EdgeHandle& eh) const = 0;
    virtual bool is_boundary(const OVM::VertexHandle& eh) const = 0;

protected:

#warning TODO: probably replace vector with fixed-size array
    //OVM::VertexPropertyT<std::vector<int>> vertex_duplicates_prop_;
    OVM::CellPropertyT<std::vector<EigenVec>> per_vertex_pos_prop_;
    OVM::CellPropertyT<std::vector<EigenVec>> per_vertex_ref_pos_prop_;
    //OVM::CellPropertyT<std::vector<Eigen::Vector3d>> per_vertex_previous_pos_prop_;

    //OVM::CellPropertyT<std::vector<bool>> glued_duplicates_prop_;

    //OVM::VertexPropertyT<Eigen::Vector3d> anchor_previous_pos_prop_;
    OVM::VertexPropertyT<EigenVec> anchor_position_prop_;
    OVM::VertexPropertyT<EigenVec> initial_anchor_position_prop_;
    //OVM::FacePropertyT<double> last_normals_angles_prop_;

    EdgeMatrix E_ref_;
    EdgeMatrix E_ref_inv_;

    double ref_edge_length_ = 0;

    // centroid vertex is copied by k duplicates, with k = #incident tets
    constexpr static int n_vars_per_cell_ = GeoDim * TopoDim;
    int n_vars_;
    int n_anchor_vertices_;// last indices for "anchor" vertices
    int anchor_vertices_index_start_;//

    const bool enable_anchor_vertices_;

    OVM::VertexPropertyT<bool> original_vertex_prop_;

    OVM::EdgePropertyT<bool> protected_edge_prop_;
};



class NonConformingTetrahedralMesh : public NonConformingMeshT<3, 4, TetrahedralMesh, OVM::CellIter, OVM::CellHandle, OVM::Entity::Cell>{

public:

    using NonConformingMesh = NonConformingMeshT<3, 4, TetrahedralMesh, OVM::CellIter, OVM::CellHandle, OVM::Entity::Cell>;
    using NonConformingMesh::is_boundary;

    NonConformingTetrahedralMesh()
    : collapse_axis_prop_(this->request_vertex_property<std::pair<OVM::VertexHandle, OVM::VertexHandle>>())
    {}

    NonConformingTetrahedralMesh(const NonConformingTetrahedralMesh& other)
    : collapse_axis_prop_(this->request_vertex_property<std::pair<OVM::VertexHandle, OVM::VertexHandle>>())  {
        *this = other;
    }

    NonConformingTetrahedralMesh& operator=(const NonConformingTetrahedralMesh& other){

        NonConformingMeshT<3, 4, TetrahedralMesh, OVM::CellIter, OVM::CellHandle, OVM::Entity::Cell>::operator=(other);

        for(auto el: this->elements()){
            this->per_vertex_pos_prop(el) = other.per_vertex_pos_prop(el);
            this->per_vertex_ref_pos_prop(el) = other.per_vertex_ref_pos_prop(el);
        }
        return *this;
    }


    NonConformingTetrahedralMesh& operator=(const TetrahedralMesh& other){

        TetrahedralMesh::operator=(other);

        //log(LogLevel::INFO)<<" CALLED TET MESH OPERATOR="<<std::endl;
        /*for(auto el: this->elements()){
            this->per_vertex_pos_prop(el) = other.per_vertex_pos_prop(el);
            this->per_vertex_ref_pos_prop(el) = other.per_vertex_ref_pos_prop(el);
        }*/
        return *this;
    }


    NonConformingTetrahedralMesh& operator=(const tet_sweep::TetMesh& other){

        *this = other.to_OVM<NonConformingTetrahedralMesh, OVM::VertexHandle>();
        return *this; 
    }


        virtual int check_edge_glueing_constraints(double glueing_eps,
                                               double& max_diff_norm,
                                               bool only_check_boundary = false) const {


        log(LogLevel::INFO)<<" ======================================================== "<<std::endl;
        log(LogLevel::INFO)<<" - checking tetrahedral edge glueing constraints..."<<std::endl;

        max_diff_norm = 0;

        /*auto local_cid_prop =  this->template request_cell_property<int>();

        int local_cid(0);
        for(auto ch: this->faces()){
            local_cid_prop[ch] = local_cid;
            local_cid++;
        }*/

        int bad_interior_edge_duplicates_count(0);
        int bad_boundary_edge_duplicates_count(0);

        for(auto e: this->edges()){
            if(!only_check_boundary || this->is_boundary(e)){

                //log(LogLevel::INFO)<<" ------------------------------ checking edge "<<e<<": "<<this->edge(e)<<std::endl;
                auto from_v = this->edge(e).from_vertex();
                auto to_v   = this->edge(e).to_vertex();

                //log(LogLevel::INFO)<<" - boundary: "<<this->is_boundary(e)<<std::endl;

                auto edge_vec = this->vertex(to_v) - this->vertex(from_v);
                Eigen::Vector3d eigen_edge_vec = {edge_vec[0], edge_vec[1], edge_vec[2]};
                //log(LogLevel::INFO)<<" - mesh edge vector: "<<eigen_edge_vec.transpose()<<std::endl;

                for(auto ef_it = this->ef_iter(e); ef_it.valid(); ef_it++){

                    /*log(LogLevel::INFO)<<" ----------------- handling face ";
                    for(auto v: this->get_halfface_vertices(this->halfface_handle(*ef_it, 0))){
                        log(LogLevel::INFO)<<" "<<v;
                    }
                    log(LogLevel::INFO)<<std::endl;*/

                    auto left_hf = this->halfface_handle(*ef_it, 0);
                    auto right_hf = this->halfface_handle(*ef_it, 1);

                    auto left_c = this->incident_cell(left_hf);
                    auto right_c = this->incident_cell(right_hf);


                    int from_left_vid(-1);
                    int   to_left_vid(-1);

                    int from_right_vid(-1);
                    int   to_right_vid(-1);

                    Eigen::Vector3d left_e, right_e;

                    if(left_c.is_valid()){
                        const auto& left_c_verts = this->get_cell_vertices(left_c);

                        for(int local_vid(0); local_vid<4; local_vid++){
                            const auto& left_v = left_c_verts[local_vid];
                            if(left_v == from_v){
                                from_left_vid = local_vid;
                            }else if(left_v == to_v){
                                to_left_vid = local_vid;
                            }
                        }
                        //log(LogLevel::INFO)<<" - left from-to vids: "<<from_left_vid<<", "<<to_left_vid<<std::endl;
                        left_e = per_vertex_pos_prop_[left_c][to_left_vid] - per_vertex_pos_prop_[left_c][from_left_vid];
                        //log(LogLevel::INFO)<<" - left edge vector: "<<left_e.transpose()<<std::endl;

                    }

                    if(right_c.is_valid()){
                        const auto& right_c_verts = this->get_cell_vertices(right_c);

                        for(int local_vid(0); local_vid<4; local_vid++){
                            const auto& right_v = right_c_verts[local_vid];
                            if(right_v == from_v){
                                from_right_vid = local_vid;
                            }else if(right_v == to_v){
                                to_right_vid   = local_vid;
                            }
                        }
                        //log(LogLevel::INFO)<<" - right from-to vids: "<<from_right_vid<<", "<<to_right_vid<<std::endl;
                        right_e = per_vertex_pos_prop_[right_c][to_right_vid] - per_vertex_pos_prop_[right_c][from_right_vid];
                        //log(LogLevel::INFO)<<" - right edge vector: "<<right_e.transpose()<<std::endl;
                    }

                    if(left_c.is_valid() && right_c.is_valid()){

                        double diff_norm = (right_e - left_e).norm();
                        max_diff_norm = std::max(max_diff_norm, diff_norm);

                        if(diff_norm > glueing_eps){
                            //log(LogLevel::INFO)<<" --> ||left_e - right_e|| = "<<(right_e - left_e).norm()<<" > "<<glueing_eps<<" for "<<(this->is_boundary(e) ? "boundary":"interior")<<" edge ("<<from_v<<", "<<to_v<<")"<<std::endl;
                            bad_interior_edge_duplicates_count++;
                            //return 1;
                        }

                    }else{

                        if(right_c.is_valid()){
                            left_e = right_e;
                        }
                        if((eigen_edge_vec - left_e).norm() > glueing_eps){

                            //log(LogLevel::INFO)<<" --> ||edge_vec - ∂edge_vec|| = "<<(eigen_edge_vec - left_e).norm()<<" > "<<glueing_eps<<" for boundary edge ("<<from_v<<", "<<to_v<<")"<<std::endl;
                            bad_boundary_edge_duplicates_count++;
                            //return 1;
                        }

                    }

                }
            }
        }

        log(LogLevel::INFO)<<" - max misalignment = "<<max_diff_norm<<std::endl;

        if(bad_interior_edge_duplicates_count || bad_boundary_edge_duplicates_count){
            log(LogLevel::INFO)<<" ... done! Found "<<bad_interior_edge_duplicates_count<<" interior and "<<bad_boundary_edge_duplicates_count<<" boundary misaligned edge duplicates"<<std::endl;
        }else{
            log(LogLevel::INFO)<<" ... done! All good"<<std::endl;
        }

        log(LogLevel::INFO)<<" ======================================================== "<<std::endl;

        return bad_interior_edge_duplicates_count || bad_boundary_edge_duplicates_count;
    }


#warning TODO: this should be an enum
    /* return value:
     * 0: closed tet fan, valid holonomy
     * 1: closed tet fan, invalid holonomy
     * 2: open tet fan, valid holonomy (geometrically less than 2π
     * 3: open tet fan, invalid holonomy (geometrically more than 2π)
     * 4: multiple tet fans */
    int check_edge_holonomy(const OVM::EdgeHandle eh,
                            double glueing_eps,
                            double& angle_diff_to_360_deg,
                            bool print_debug = false){

        if(print_debug) {
            log(LogLevel::INFO) << " --------------------------------------- checking holonomy of "
                      << (this->is_boundary(eh) ? "boundary" : "interior") << " edge (" << edge(eh).from_vertex()
                      << ", " << edge(eh).to_vertex() << ")" << std::endl;
        }

        int non_glued_f_count(0);
        for (auto ef_it = this->ef_iter(eh); ef_it.valid(); ef_it++) {
            if(!this->is_boundary(*ef_it)) {
                auto left_hf = this->halfface_handle(*ef_it, 0);
                auto right_hf = this->opposite_halfface_handle(left_hf);

                auto hf_verts = this->get_halfface_vertices(left_hf);
                if (print_debug) {
                    log(LogLevel::INFO) << " ------------- checking face " << this->halfface_to_str(left_hf)
                              << std::endl;
                }

                auto left_ch = this->incident_cell(left_hf);
                auto right_ch = this->incident_cell(right_hf);

                if (print_debug) {
                    log(LogLevel::INFO) << "  left cell is " << this->element_to_str(left_ch) << std::endl;
                    log(LogLevel::INFO) << " right cell is " << this->element_to_str(right_ch) << std::endl;
                }

                if (left_hf.is_valid() && right_ch.is_valid()) {

                    auto left_ch_verts = this->get_cell_vertices(left_ch);
                    auto right_ch_verts = this->get_cell_vertices(right_ch);

                    int glued_dup_count(0);

                    for (auto hf_v: hf_verts) {

                        if (print_debug) {
                            log(LogLevel::INFO) << " -- checking that " << hf_v << " is glued" << std::endl;
                        }

                        int left_ch_idx = -1;
                        int right_ch_idx = -1;
                        for (int i(0); i < 4; i++) {
                            if (hf_v == left_ch_verts[i]) {
                                left_ch_idx = i;
                            }
                            if (hf_v == right_ch_verts[i]) {
                                right_ch_idx = i;
                            }
                        }

                        if (print_debug) {
                            log(LogLevel::INFO) << " -- l/r indices: " << left_ch_idx << "/" << right_ch_idx << std::endl;
                        }

                        auto diff_norm = (per_vertex_pos_prop(left_ch)[left_ch_idx] -
                                          per_vertex_pos_prop(right_ch)[right_ch_idx]).norm();

                        if (print_debug) {
                            log(LogLevel::INFO) << " -- duplicates distance: " << diff_norm << std::endl;
                        }
                        if (diff_norm < glueing_eps) {
                            glued_dup_count++;

                            if (print_debug) {
                                log(LogLevel::INFO) << " --> it's glued " << std::endl;
                            }
                        } else {

                            if (print_debug) {
                                log(LogLevel::INFO) << " --> it's not glued" << std::endl;
                            }
                        }
                    }

                    if (print_debug) {
                        log(LogLevel::INFO) << " -> face has " << glued_dup_count << " glued duplicates" << std::endl;
                    }
                    if (glued_dup_count == 3) {

                        if (print_debug) {
                            log(LogLevel::INFO) << " --> face is glued" << std::endl;
                        }
                    } else {
                        non_glued_f_count++;
                    }
                } else {
                    log(LogLevel::INFO) << " that's an interior face, which shouldn't happen..." << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        }

        if(print_debug) {
            log(LogLevel::INFO) << " -> edge is incident to " << non_glued_f_count << " non-glued faces" << std::endl;
        }

        if (non_glued_f_count <= 1) {

            if(print_debug) {
                if(!this->is_boundary(eh)) {
                    log(LogLevel::INFO) << " -> single tet fan, checking holonomy" << std::endl;
                }else{
                    log(LogLevel::INFO)<<" -> one or two tet fans but boundary edge, checking holonomy"<<std::endl;
                }
            }
            double angle_sum_deg(0);

            if(this->is_boundary(eh)){

                if (print_debug) {
                    log(LogLevel::INFO) << " - boundary edge, computing angle between incident boundary faces" << std::endl;
                }

                angle_sum_deg += compute_boundary_dihedral_angle_deg(*this, eh);
            }

            for (auto ec_it = this->ec_iter(eh); ec_it.valid(); ec_it++) {
                double theta = compute_dihedral_angle_deg(*ec_it, eh, NON_CONFORMING);

                if(print_debug) {
                    log(LogLevel::DEBUG)<<" - theta for cell "<<element_to_str(*ec_it)<<": "<<theta<<std::endl;
                }
                angle_sum_deg += theta;
            }


            double turn_count = angle_sum_deg / 360.0;

            if(print_debug) {
                log(LogLevel::INFO) << " - angle sum = " << angle_sum_deg << "° = " << turn_count << " turns" << std::endl;
            }

            double diff_to_turn = turn_count - 1;
            angle_diff_to_360_deg = diff_to_turn * 360.0;


            if(print_debug) {
                log(LogLevel::INFO) << " - diff to single full turn: " << diff_to_turn << " = " << angle_diff_to_360_deg << "°" << std::endl;
            }

            if (std::abs(diff_to_turn) > glueing_eps) {
                if (non_glued_f_count) {
                    if (angle_diff_to_360_deg > 0) {
                        if(print_debug) {
                            log(LogLevel::INFO) << " -> extra " << angle_diff_to_360_deg << "° of angular overlap" << std::endl;
                        }
                        return 3;
                    } else {
                        if(print_debug) {
                            log(LogLevel::INFO) << " -> missing " << (-angle_diff_to_360_deg) << " to close the fan" << std::endl;
                        }
                        return 2;
                    }
                } else {
                    if(print_debug) {
                        log(LogLevel::INFO) << " -> fan should be fully glued but there's an extra angle of " << angle_diff_to_360_deg
                                  << std::endl;
                    }
                    return 1;
                }
            }else{
                if(print_debug) {
                    log(LogLevel::INFO)<<" -> valid closed fan"<<std::endl;
                }
                return 0;
            }
        }else{
            return 4;
        }

        return -1;
    }


    /** NOTE: doesn't make sense before the mesh is fully glued, because wrong holonomy is to be expected for edges that are not fully glued
     */
    int check_holonomy(double eps){

        log(LogLevel::INFO)<<" ======================================================== "<<std::endl;
        log(LogLevel::INFO)<<" - checking holonomy..."<<std::endl;

        /** contains:
         * #  valid closed holonomy edges
         * #invalid closed holonomy edges
         * #  valid open holonomy edges
         * #invalid open holonomy edges
         * multiple fans edges
         * boundary edges
         * */
        std::vector<int> edge_types_count(6);

        for(auto eh: this->edges()){

            /*log(LogLevel::INFO) << " ------------- checking edge " << this->edge(eh).from_vertex() << " - "
                      << this->edge(eh).to_vertex() << std::endl;*/

            double angle(0);
            auto holonomy_type = check_edge_holonomy(eh, eps, angle);
            if(holonomy_type < 0 || holonomy_type >= 6){
                log(LogLevel::ERROR)<<" ERROR - unhandled holonomy type "<<holonomy_type<<std::endl;
                exit(EXIT_FAILURE);
            }
            edge_types_count[holonomy_type]++;
        }

        log(LogLevel::INFO)<<" ---------+----------+--------------"<<std::endl;
        log(LogLevel::INFO)<<" holonomy |    ok    |    not ok    "<<std::endl;
        log(LogLevel::INFO)<<" ---------+----------+--------------"<<std::endl;
        log(LogLevel::INFO)<<"  closed  | "<<std::setw(8)<<edge_types_count[0]<<" | "<<std::setw(8)<<edge_types_count[1]<<std::endl;
        log(LogLevel::INFO)<<" ---------+----------+--------------"<<std::endl;
        log(LogLevel::INFO)<<"   open   | "<<std::setw(8)<<edge_types_count[2]<<" | "<<std::setw(8)<<edge_types_count[3]<<std::endl;
        log(LogLevel::INFO)<<" ---------+----------+--------------"<<std::endl;
        log(LogLevel::INFO)<<" edges with more than one fan incident: "<<edge_types_count[4]<<std::endl;

        log(LogLevel::INFO)<<" ======================================================== "<<std::endl;
        return edge_types_count[1] | edge_types_count[3];
    }

    /* Warning: assumes that the property is zero-initialized */
    void fill_angles_sum_property(OVM::VertexPropertyT<double>& angles_sum_prop){
        log(LogLevel::ERROR)<<" ERROR - angle sum for non-conforming tet meshes is not implemented!"<<std::endl;
        exit(EXIT_FAILURE);
    }

    //tet
    EigenVec solve_edge_Chebyshev(OVM::EdgeHandle eh,
                                  double& radius,
                                  bool print_debug = false){
        using namespace OpenVolumeMesh;

        //print_debug = true;
        if(this->is_boundary(eh)){
            log(LogLevel::INFO)<<" WARNING - asking for edge Chebyshev center on a boundary edge"<<std::endl;
            return EigenVec();
        }
        std::vector<EigenVec> normals;
        std::vector<EigenVec> points;

        auto from_v = this->edge(eh).from_vertex();
        auto to_v   = this->edge(eh).to_vertex();

        //log(LogLevel::INFO)<<" ------- checking chebyshev centroid for 1-ring of edge "<<from_v<<" "<<to_v<<std::endl;

        std::vector<VertexHandle> end_verts = {from_v, to_v};

        for(auto v: end_verts){
            //log(LogLevel::INFO)<<" -- checking faces around end vertex "<<v<<std::endl;
            for(auto vc_it = this->vc_iter(v); vc_it.valid(); vc_it++){
                auto c_verts = get_element_vertices(*vc_it);
                //log(LogLevel::INFO)<<" --- checking cell "<<element_to_str(*vc_it)<<std::endl;
                HalfFaceHandle op_hf(-1);
                for(auto chf_it = this->chf_iter(*vc_it); chf_it.valid(); chf_it++){
                    if(this->halfface_opposite_vertex(*chf_it) == v){
                        op_hf = *chf_it;
                    }
                }
                if(!op_hf.is_valid()){
                    log(LogLevel::ERROR)<<" ERROR - couldn't find opposite halfface to vertex "<<v<<" in cell "<<this->element_to_str(*vc_it)<<std::endl;
                    exit(EXIT_FAILURE);
                }

                auto op_hf_verts = this->get_halfface_vertices(op_hf);
                /*log(LogLevel::INFO)<<" - opposite halfface: "<<op_hf_verts[0]<<", "<<op_hf_verts[1]<<", "<<op_hf_verts[2]<<std::endl;
                for(auto v: op_hf_verts){
                    log(LogLevel::INFO)<<" -- v "<<v<<" at "<<this->vertex(v)<<std::endl;
                }*/

                // faces incident to the edge are discarded
                auto left  = this->vertex(op_hf_verts[1]) - this->vertex(op_hf_verts[0]);
                auto right = this->vertex(op_hf_verts[2]) - this->vertex(op_hf_verts[0]);
                if(print_debug) {
                    log(LogLevel::INFO)<<" ----  left: "<<left<<std::endl;
                    log(LogLevel::INFO)<<" ---- right: "<<right<<std::endl;
                }

                EigenVec mid_point = vec2vec(this->vertex(op_hf_verts[0]) + this->vertex(op_hf_verts[1]) + this->vertex(op_hf_verts[2])) / 3.0;
                EigenVec normal = -vec2vec((left.cross(right)).normalized());

                if(print_debug) {
                    log(LogLevel::INFO)<<" ----           mid-point: "<<mid_point.transpose()<<std::endl;
                    log(LogLevel::INFO)<<" ----    normal direction: "<<(-vec2vec((left.cross(right)))).transpose()<<std::endl;
                    log(LogLevel::INFO)<<" ----              normal: "<<normal.transpose()<<std::endl;
                }
                normals.push_back(normal);
                points.push_back(mid_point);

            }
        }

        return solve_Chebyshev<EigenVec>(normals, points, radius, print_debug);
    }

    //tet
    bool is_flipped_or_degenerate(const OVM::CellHandle& ch,
                                  VERTEX_POSITION_SOURCE pos_source) const {

        using namespace CGAL_wrapper;

        auto cell_vertices = this->get_cell_vertices(ch);
        std::vector<Vec3d> verts;
        for(int i(0); i<4; i++){
            switch(pos_source) {
                case CONFORMING: {
                    verts.push_back(this->vertex(cell_vertices[i]));
                    break;
                }
                case NON_CONFORMING: {
                    verts.push_back(vec2vec(per_vertex_pos_prop(ch)[i]));
                    break;
                }
                case REFERENCE: {
                    verts.push_back(vec2vec(per_vertex_ref_pos_prop(ch)[i]));
                    break;
                }
            }
        }

        CGAL_Tetrahedron tet( ovm_to_cgal(verts[0]),
                              ovm_to_cgal(verts[1]),
                              ovm_to_cgal(verts[2]),
                              ovm_to_cgal(verts[3]));

        return tet.orientation() != CGAL::POSITIVE;
    }

    //tet
    OVM::VertexHandle split_edge_and_update_props(const OVM::EdgeHandle& to_split){

        //log(LogLevel::INFO)<<" ------------------------------- splitting edge "<<this->edge(to_split).from_vertex()<<", "<<this->edge(to_split).to_vertex()<<std::endl;
        //log(LogLevel::INFO)<<" - boundary: "<<this->is_boundary(to_split)<<std::endl;

        if(this->is_deleted(to_split)){
            log(LogLevel::INFO)<<" -> deleted edge, skipping"<<std::endl;
            return OVM::VertexHandle(-1);
        }

        using namespace Eigen;

        using CellMapT = std::pair<OVM::CellHandle, std::map<OVM::VertexHandle, std::pair<Eigen::Vector3d, Eigen::Vector3d>>>;

        // original cell and its vh->(pos, ref_pos) map
        auto original_cell_prop = this->request_cell_property<CellMapT>();

        auto from_v = this->edge(to_split).from_vertex();
        auto to_v   = this->edge(to_split).to_vertex();

        for(auto ec_it = this->ec_iter(to_split); ec_it.valid(); ec_it++){
            original_cell_prop[*ec_it].first = *ec_it;

            auto c_verts = this->get_cell_vertices(*ec_it);
            //log(LogLevel::INFO)<<" ---------- handling cell "<<element_to_str(*ec_it)<<std::endl;

            for(int local_vid(0); local_vid < 4; local_vid++){
                original_cell_prop[*ec_it].second.insert({c_verts[local_vid],
                                                          {per_vertex_pos_prop_[*ec_it][local_vid],
                                                           per_vertex_ref_pos_prop_[*ec_it][local_vid]}});
            }
        }


        auto mid_v = this->split_edge(to_split);

        if(original_vertex_prop(from_v)){
            if(original_vertex_prop(to_v)){
                collapse_axis_prop_[mid_v] = {from_v, to_v};
            }else{
                collapse_axis_prop_[mid_v] = collapse_axis_prop_[to_v];
            }
        }else if(original_vertex_prop(to_v)){
            collapse_axis_prop_[mid_v] = collapse_axis_prop_[from_v];
        }else{
            /*if(collapse_axis_prop_[from_v] != collapse_axis_prop_[to_v]){
                log(LogLevel::INFO)<<" ERROR - collapse axis for from-vertex "<<from_v
                         <<" is ("<<collapse_axis_prop_[from_v].first<<" - "<<collapse_axis_prop_[from_v].second<<")"
                         <<", but axis for to-vertex "<<to_v
                         <<" is ("<<collapse_axis_prop_[to_v].first<<" - "<<collapse_axis_prop_[to_v].second<<")"<<std::endl;
                exit(EXIT_FAILURE);
            }*/
            collapse_axis_prop_[mid_v] = {from_v, to_v};
        }
        //log(LogLevel::INFO)<<" -> set collapse axis for mid-vertex "<<mid_v<<" as ("<<collapse_axis_prop_[mid_v].first<<" - "<<collapse_axis_prop_[mid_v].second<<")"<<std::endl;

        //log(LogLevel::INFO)<<" - mid-vertex is "<<mid_v<<std::endl;

        double min_max_angle(181);
        double best_t(0.5);

        //skipping this stuff for now
#if 0
        double trials_count = json_data_["split_mid_point_trials_count"];
        for(double i(1); i <= trials_count; i++){

            auto t = i / (trials_count + 1);

            //log(LogLevel::INFO)<<" -- trying t = "<<t<<std::endl;

            update_props_around_mid_vertex(from_v, mid_v, to_v, t, original_cell_prop);

            double max_angle(0);
            for(auto vc_it = this->vc_iter(mid_v); vc_it.valid(); vc_it++){
                MeshT single_cell;

                auto c_verts = this->get_cell_vertices(*vc_it);
                //log(LogLevel::INFO)<<" --- checking cell "<<(*vc_it)<<":";
                //for(auto v: c_verts){log(LogLevel::INFO)<<" "<<v;}
                //log(LogLevel::INFO)<<std::endl;

                int local_vid(0);
                int local_from_vid(-1);
                int local_to_vid(-1);
                for(auto vpos: per_vertex_pos_prop_[*vc_it]){
                    single_cell.add_vertex({vpos[0], vpos[1], vpos[2]});
                    if(c_verts[local_vid] == from_v){
                        local_from_vid = local_vid;
                    }else if(c_verts[local_vid] == to_v){
                        local_to_vid = local_vid;
                    }
                    local_vid++;
                }

                if(local_from_vid == -1 && local_to_vid == -1){
                    log(LogLevel::INFO)<<" ---> couldn't find neither from- or to-vertex"<<std::endl;
                    exit(EXIT_FAILURE);
                }

                //log(LogLevel::INFO)<<" --- local from-to vertices are "<<local_from_vid<<", "<<local_to_vid<<std::endl;
                single_cell.add_cell({VertexHandle(0), VertexHandle(1), VertexHandle(2), VertexHandle(3)});

                double max_cell_angle(0);
                EdgeHandle worst_edge;
                for(auto ce_it = single_cell.ce_iter(CellHandle(0)); ce_it.valid(); ce_it++){
                    //log(LogLevel::INFO)<<" ---- checking edge "<<single_cell.edge(*ce_it)<<std::endl;
                    if(single_cell.edge(*ce_it).from_vertex() == VertexHandle(local_from_vid) ||
                       single_cell.edge(*ce_it).to_vertex()   == VertexHandle(local_from_vid) ||
                       single_cell.edge(*ce_it).from_vertex() == VertexHandle(local_to_vid) ||
                       single_cell.edge(*ce_it).to_vertex()   == VertexHandle(local_to_vid)){

                        //log(LogLevel::INFO)<<" ----- skipped edge "<<single_cell.edge(*ce_it)<<std::endl;
                        continue;
                    }
                    auto angle_rad = PEHelpers::compute_dihedral_angle(single_cell, *ce_it, CellHandle(0));
                    auto angle_deg = 180 * angle_rad/M_PI;

                    if(angle_deg > max_cell_angle){
                        worst_edge = *ce_it;
                    }
                    max_cell_angle = std::max(max_cell_angle, angle_deg);

                    //log(LogLevel::INFO)<<" ---- angle for edge = "<<angle_deg<<std::endl;
                }
                //log(LogLevel::INFO)<<" --- max angle for cell "<<(*vc_it)<<" = "<<max_cell_angle<<" at edge "<<single_cell.edge(worst_edge)<<std::endl;
                max_angle = std::max(max_angle, max_cell_angle);
            }
            //log(LogLevel::INFO)<<" -- max angle with t = "<<t<<": "<<max_angle<<std::endl;

            if(max_angle < min_max_angle){
                min_max_angle = max_angle;
                best_t = t;
                //log(LogLevel::INFO)<<" ---> updated best t to "<<t<<" with max angle "<<min_max_angle<<std::endl;
            }
        }
#endif

        //log(LogLevel::INFO)<<" - best t is "<<best_t<<", giving max angle = "<<min_max_angle<<std::endl;

        bool boundary_mid_v = this->is_boundary(mid_v);
        if(boundary_mid_v){
            log(LogLevel::INFO)<<" -> mid-vertex is boundary, adapting incident element to boundary"<<std::endl;
        }
        for(auto vc_it = this->vc_iter(mid_v); vc_it.valid(); vc_it++) {

            auto c_verts = this->get_cell_vertices(*vc_it);
            auto original_cell = original_cell_prop[*vc_it].first;
            //const auto &original_cell_vertices = original_cell_prop[*vc_it].second;

            //log(LogLevel::INFO) << " ---------- handling cell " << element_to_str(*vc_it) << std::endl;
            //log(LogLevel::INFO)<<" -- original cell: "<<original_cell<<std::endl;
            //log(LogLevel::INFO)<<" -- with vertices: "; for(auto v: original_cell_vertices)log(LogLevel::INFO)<<v<<" "; log(LogLevel::INFO)<<std::endl;

            auto &v_to_pos_map = original_cell_prop[original_cell].second;

            if (auto search = v_to_pos_map.find(from_v) == v_to_pos_map.end()) {
                log(LogLevel::ERROR) << " ERROR - couldnt' find from-vertex in position map" << std::endl;
                exit(EXIT_FAILURE);
                //return -1;
            }
            if (auto search = v_to_pos_map.find(to_v) == v_to_pos_map.end()) {
                log(LogLevel::ERROR) << " ERROR - couldnt' find to-vertex in position map" << std::endl;
                exit(EXIT_FAILURE);
                //return -1;
            }

            Vector3d mid_v_pos =     v_to_pos_map[from_v].first  + best_t * (v_to_pos_map[to_v].first  - v_to_pos_map[from_v].first);
            Vector3d ref_mid_v_pos = v_to_pos_map[from_v].second + best_t * (v_to_pos_map[to_v].second - v_to_pos_map[from_v].second);

            //log(LogLevel::INFO)<<" -     mid-v pos = "<<mid_v_pos.transpose()<<std::endl;
            //log(LogLevel::INFO)<<" - ref mid-v pos = "<<ref_mid_v_pos.transpose()<<std::endl;

            /*per_vertex_pos_prop_[*vc_it].resize(4);
            per_vertex_ref_pos_prop_[*vc_it].resize(4);*/
            //log(LogLevel::INFO)<<" -- per vertex pos prop size: "<<per_vertex_pos_prop_[*vc_it].size()<<std::endl;

            for (int local_vid(0); local_vid < 4; local_vid++) {
                if (c_verts[local_vid] == mid_v) {
                    if(this->is_boundary(mid_v)){
                        per_vertex_pos_prop_[*vc_it][local_vid]     = vec2vec(this->vertex(mid_v));
                        per_vertex_ref_pos_prop_[*vc_it][local_vid] = vec2vec(this->vertex(mid_v));
                    }else {
                        per_vertex_pos_prop_[*vc_it][local_vid] = mid_v_pos;
                        per_vertex_ref_pos_prop_[*vc_it][local_vid] = ref_mid_v_pos;
                    }
                    //auto global_id =
                    //vertex_duplicates_prop_[v].push_back(global_id);

                } else {
                    per_vertex_pos_prop_[*vc_it][local_vid] = v_to_pos_map[c_verts[local_vid]].first;
                    per_vertex_ref_pos_prop_[*vc_it][local_vid] = v_to_pos_map[c_verts[local_vid]].second;
                }
            }
            //if(boundary_mid_v){
            //    setup_boundary_adapted_init_positions(*vc_it);
            //}
        }

        //update the conforming mesh position
        PointT conforming_mid_v_pos =  this->vertex(from_v)  + best_t * (this->vertex(to_v) - this->vertex(from_v));

        return mid_v;
    }

    //tet
    OVM::HalfEdgeHandle collapse_preserves_boundary(OVM::EdgeHandle eh){


        bool enable_interior_to_boundary_collapse(false);

        OVM::HalfEdgeHandle invalid_heh(-1);
        if(is_protected(eh)){
            return invalid_heh;
        }

        auto from_v = this->edge(eh).from_vertex();
        auto to_v   = this->edge(eh).to_vertex();

        //log(LogLevel::INFO)<<" -------------------- checking whether collapse "<<from_v<<" -> "<<to_v<<" preserves the boundary"<<std::endl;


        auto axis_from = collapse_axis_prop_[from_v].first;
        auto axis_to   = collapse_axis_prop_[from_v].second;

        //log(LogLevel::INFO)<<" - from-vertex "<<from_v<<" on axis "<<axis_from<<" - "<<axis_to<<std::endl;

        auto heh0 = this->halfedge_handle(eh, 0);
        auto heh1 = this->halfedge_handle(eh, 1);

        if (!this->is_boundary(from_v) && !this->is_boundary(to_v)) { //if full interior, then ok
            //log(LogLevel::INFO)<<" -> full interior -> collapsible"<<std::endl;
            return heh0;
        }

        if(this->is_boundary(from_v) != this->is_boundary(to_v)){

            if(enable_interior_to_boundary_collapse) {
                //log(LogLevel::INFO) << " -> interior to boundary -> checking if it would create a cord..." << std::endl;


                OVM::VertexHandle boundary_vertex = this->is_boundary(from_v) ? from_v : to_v;
                OVM::VertexHandle interior_vertex = this->is_boundary(from_v) ? to_v : from_v;

                if (!collapse_would_create_cord(*this, eh)) {
                    return this->find_halfedge(interior_vertex, boundary_vertex);
                }
#if 0
                int boundary_vertices_in_edge_1_ring(0);
                for (auto ef_it = this->ef_iter(eh); ef_it.valid(); ef_it++) {
                    auto verts = this->get_halfface_vertices(this->halfface_handle(*ef_it, 0));
                    for (auto v: verts) {
                        if (v != from_v && v != to_v) {
                            boundary_vertices_in_edge_1_ring += this->is_boundary(v);
                        }
                    }
                }
                //log(LogLevel::INFO) << " - found " << boundary_vertices_in_edge_1_ring
                //          << " boundary vertices in edge 1-ring (excluding the boundary end vertex)" << std::endl;

                int boundary_neighbors(0);
                OVM::VertexHandle boundary_vertex = this->is_boundary(from_v) ? from_v : to_v;
                OVM::VertexHandle interior_vertex = this->is_boundary(from_v) ? to_v : from_v;

                for (auto out_he: this->outgoing_halfedges(interior_vertex)) {
                    boundary_neighbors += this->is_boundary(this->to_vertex_handle(out_he));
                }
                boundary_neighbors--; //to remove the boundary end vertex

                //log(LogLevel::INFO) << " - found " << boundary_neighbors << " boundary neighbors to interior vertex " << interior_vertex
                //          << " (excluding the boundary end vertex)" << std::endl;

                if (boundary_vertices_in_edge_1_ring > boundary_neighbors) {
                    log(LogLevel::ERROR)
                            << " ERROR - there should never be more boundary vertices in edge 1-ring than boundary vertices in interior vertex 1-ring"
                            << std::endl;
                    exit(EXIT_FAILURE);
                }
                if (boundary_vertices_in_edge_1_ring == boundary_neighbors) {
                    //log(LogLevel::INFO) << " --> same count, edge collapse would not create a cord" << std::endl;
                    //log(LogLevel::INFO) << " - now checking for cord face..." << std::endl;

                    auto heh = this->find_halfedge(interior_vertex, boundary_vertex);
                    std::vector <OVM::VertexHandle> edge_neighbors;
                    for (auto hehf_it = this->hehf_iter(heh); hehf_it.valid(); hehf_it++) {
                        auto neighbor = this->to_vertex_handle(this->next_halfedge_in_halfface(heh, *hehf_it));
                        edge_neighbors.push_back(neighbor);
                        //log(LogLevel::INFO) << " -- added neighbor " << neighbor << " to list" << std::endl;
                    }

                    if (edge_neighbors.size() < 3) {
                        log(LogLevel::ERROR) << " ERROR - less than 3 neighbors, which makes no sense" << std::endl;
                        exit(EXIT_FAILURE);
                    }

                    bool found_boundary_edge(false);
                    for (int i(0); i < edge_neighbors.size() - 3; i++) {
                        for (int j(i + 2); j < edge_neighbors.size() - 2; j++) {
                            auto from_neighbor = edge_neighbors[i];
                            auto to_neighbor = edge_neighbors[j];
                            //log(LogLevel::INFO) << " -- checking neighbors " << i << " and " << j << " -> " << from_neighbor << " and "
                            //          << to_neighbor << std::endl;

                            if (this->halfedge(from_neighbor, to_neighbor).is_valid()) {
                                //log(LogLevel::INFO) << " ---> edge exists, collapse would create a cord face" << std::endl;
                                found_boundary_edge = true;
                                //return invalid_heh; // commented out for now to test execution
                            }
                        }
                    }

                    if (found_boundary_edge) {
                        return invalid_heh;
                    }

                    //log(LogLevel::INFO) << " --> collapse in direction " << interior_vertex << " -> " << boundary_vertex << " is valid" << std::endl;
                    return heh;
                } else {
                    //log(LogLevel::INFO) << " --> boundary vertices in 1-ring of interior vertex but not of edge, this collapse would create a cord" << std::endl;
                }

#endif
            }else {
                //log(LogLevel::INFO) << " -> interior to boundary, not collapsing it" << std::endl;
            }

            return invalid_heh;
        }

        //if either vertex is boundary, then we have to check whether it's collapsing on an original axis

        //log(LogLevel::INFO)<<" - boundary edge, checking if preserving boundary"<<std::endl;
        //if it's collapsing to an original vertex, we gotta make sure the from-vertex is both not an original vertex AND on the right collapse axis
        if(original_vertex_prop(to_v)){
            if(original_vertex_prop(from_v)){
                //log(LogLevel::INFO)<<" -> original to original -> not collapsible"<<std::endl;
                return invalid_heh;

            }else if(axis_from == to_v || axis_to == to_v){
                //log(LogLevel::INFO)<<" -> from-vertex is not original, and on the right axis -> collapsible"<<std::endl;
                return heh0;

            }else{
                //log(LogLevel::INFO)<<" -> from-vertex is not original, but on the wrong axis -> not collapsible"<<std::endl;
                return invalid_heh;
            }
        }else if(original_vertex_prop(from_v)){//maybe we can switch the collapse direction

            axis_from = collapse_axis_prop_[to_v].first;
            axis_to   = collapse_axis_prop_[to_v].second;

            //log(LogLevel::INFO)<<" - updated axis to to-vertex' "<<to_v<<" on axis "<<axis_from<<" - "<<axis_to<<std::endl;

            if(axis_from == from_v || axis_to == from_v){
                //log(LogLevel::INFO)<<" -> to-vertex is not original, and on the right axis -> collapsible (flipped direction)"<<std::endl;
                return heh1;
            }else{
                //log(LogLevel::INFO)<<" -> to-vertex is not original, but on the wrong axis -> not collapsible (even with flip)"<<std::endl;
                return invalid_heh;
            }

        }else if(collapse_axis_prop_[from_v] == collapse_axis_prop_[to_v]){ //if both end vertices are not original, then we gotta make sure we're collapsing along the right axis
            //log(LogLevel::INFO)<<" -> both end-vertices on the same axis -> collapsible"<<std::endl;
            return heh0;
        }else{
            //log(LogLevel::INFO)<<" -> both end vertices are not original, but not on the same axis -> not collapsible"<<std::endl;
            return invalid_heh;
        }
    }

    //tet
    OVM::VertexHandle collapse_edge_and_update_props(const OVM::HalfEdgeHandle& to_collapse){

        auto from_v = this->from_vertex_handle(to_collapse);
        auto   to_v = this->to_vertex_handle(to_collapse);

        const int new_cells_idx_start(this->n_elements());

        //log(LogLevel::INFO)<<" ------------------------------- collapsing "<<(this->is_boundary(this->edge_handle(to_collapse)) ? "boundary" : "interior")<<" edge "<<from_v<<", "<<to_v<<std::endl;
        if(this->is_deleted(to_collapse)){
            //log(LogLevel::INFO)<<" -> deleted edge, skipping"<<std::endl;
            return OVM::VertexHandle(-1);
        }
        if(this->is_boundary(from_v) &&
           !this->is_boundary(to_v)){
            log(LogLevel::ERROR)<<" ERROR - collapsing inward, which should never happen"<<std::endl;
            exit(EXIT_FAILURE);
        }
        if(original_vertex_prop(from_v) &&
           this->is_boundary(from_v)){
            log(LogLevel::ERROR)<<" ERROR - collapsing original boundary vertex, which should never happen"<<std::endl;
            exit(EXIT_FAILURE);
        }


        auto deleted_cell_id_prop = this->request_cell_property<int>("", -1);
        int deleted_cell_id(0);
        auto duplicate_map_prop = this->request_cell_property<std::map<OVM::VertexHandle, EigenVec>>();
        for(auto vc_it = this->vc_iter(from_v); vc_it.valid(); vc_it++) {
            auto verts = this->get_element_vertices(*vc_it);
            //log(LogLevel::INFO)<<" - map for "<<deleted_cell_id<<"-th neighbor cell "<<this->element_to_str(*vc_it)<<std::endl;
            for(int i(0); i<verts.size(); i++){
                auto v = verts[i];
                duplicate_map_prop[*vc_it].insert({v, per_vertex_pos_prop(*vc_it)[i]});
                //log(LogLevel::INFO)<<" -- "<<v<<" -> "<<per_vertex_pos_prop(*vc_it)[i].transpose()<<std::endl;
            }
            deleted_cell_id_prop[*vc_it] = deleted_cell_id++;

            if(is_flipped_or_degenerate(*vc_it, NON_CONFORMING)){
                double vol = this->non_conforming_element_volume(*vc_it);
                log(LogLevel::ERROR)<<" ERROR - somehow cell "<<element_to_str(*vc_it)<<" got flipped or degenerate, with volume = "<<vol<<std::endl;
                log(LogLevel::ERROR)<<" comparing conforming positions with non-conforming ones:"<<std::endl;
                auto verts = this->get_cell_vertices(*vc_it);
                for(int i(0); i<verts.size(); i++){
                    log(LogLevel::ERROR)<<" -- v"<<verts[i]<<" conforming position: "<<this->vertex(verts[i])<<", non-conforming position: "<<per_vertex_pos_prop(*vc_it)[i].transpose()<<std::endl;
                }
                exit(EXIT_FAILURE);
            }
        }

        auto surviving_v = this->collapse_edge(to_collapse);

        //log(LogLevel::INFO)<<" - initializing new cells' duplicates..."<<std::endl;
        for(auto vc_it = this->vc_iter(surviving_v); vc_it.valid(); vc_it++) {
            int deleted_id = deleted_cell_id_prop[*vc_it];
            //log(LogLevel::INFO) << " -- cell " << element_to_str(*vc_it) << " was the " << deleted_id
            //          << "-th neighbor of the collapsed vertex" << std::endl;
            if(deleted_id >= 0) {
                //log(LogLevel::INFO)<<" --- fixing its duplicate positions"<<std::endl;
                auto verts = get_element_vertices(*vc_it);
                for(int i(0); i < verts.size(); i++){
                    auto v = verts[i];
                    auto key_v = v;
                    //log(LogLevel::INFO)<<" ----- copying duplicate of vertex "<<v<<std::endl;
                    if(v == surviving_v){
                        //log(LogLevel::INFO)<<" --> surviving vertex, using deleted one as key instead "<<std::endl;
                        key_v = from_v;
                    }
                    per_vertex_pos_prop(*vc_it)[i] = duplicate_map_prop[*vc_it].at(key_v);
                    //log(LogLevel::INFO)<<" ------ position for vertex "<<v<<": "<<per_vertex_pos_prop(*vc_it)[i].transpose()<<std::endl;
                }
                //log(LogLevel::INFO)<<" -> volume after correction = "<<non_conforming_element_volume(*vc_it)<<std::endl;
            }
        }

        /*log(LogLevel::INFO)<<" ------------------------------- cells incident to surviving vertex "<<surviving_v<<std::endl;
        for(auto vc_it(this->vc_iter(surviving_v)); vc_it.valid(); vc_it++){
            log(LogLevel::INFO)<<" -- "<<this->element_to_str(*vc_it)<<std::endl;
            auto verts = this->get_element_vertices(*vc_it);
            for(int i(0); i<verts.size(); i++) {
                log(LogLevel::INFO) << " ---- v" << verts[i] << " duplicate at " << this->per_vertex_pos_prop(*vc_it)[i].transpose()
                          << std::endl;
            }
        }*/

        if(this->is_boundary(surviving_v)){
            //log(LogLevel::INFO)<<" - checking for cord creation"<<std::endl;
            for(auto out_he: this->outgoing_halfedges(surviving_v)){
                if(!this->is_boundary(this->edge_handle(out_he))){
                    auto neighbor = this->to_vertex_handle(out_he);
                    if(this->is_boundary(neighbor)){
                        log(LogLevel::ERROR)<<" ERROR - collapsing edge ("<<from_v<<", "<<to_v<<") created a cord ("<<surviving_v<<" - "<<neighbor<<")"<<std::endl;
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }

        if(this->is_boundary(surviving_v)) {
            auto surviving_v_pos = this->vertex(surviving_v);
            //log(LogLevel::INFO)<<" -> surviving vertex is on the boundary, updating its duplicates"<<std::endl;
            for(int el_idx(new_cells_idx_start); el_idx < this->n_elements(); el_idx++){
                OVM::CellHandle ch(el_idx);
                //log(LogLevel::INFO)<<" --- updating new cell "<<this->element_to_str(ch)<<std::endl;
                auto verts = this->get_element_vertices(ch);
                EigenVec surviving_v_initial_pos;

                int surviving_v_local_idx(-1);
                for(int i(0); i<verts.size(); i++){
                    auto vh = verts[i];
                    if(vh == surviving_v){
                        surviving_v_local_idx = i;
                        surviving_v_initial_pos = per_vertex_pos_prop(ch)[i];

                        //log(LogLevel::INFO)<<" -- moved surviving vertex from "<<surviving_v_initial_pos.transpose()<<" to "<<surviving_v_pos<<std::endl;
                        per_vertex_pos_prop(ch)[i] = vec2vec(surviving_v_pos);
                        break;
                    }
                }
                if(surviving_v_local_idx == -1){
                    log(LogLevel::ERROR)<<" ERROR - couldn't find surviving vertex in new cell, which is SUPER WEIRD"<<std::endl;
                    exit(EXIT_FAILURE);
                }

                double updated_vol = non_conforming_element_volume(ch);
                //log(LogLevel::INFO)<<" -- updated volume = "<<updated_vol<<std::endl;
                if(updated_vol <= DBL_EPSILON){
                    //log(LogLevel::INFO)<<" --> updating boundary duplicate flipped element, using boundary init instead"<<std::endl;
                    per_vertex_pos_prop(ch)[surviving_v_local_idx] = surviving_v_initial_pos;

                    setup_boundary_adapted_init_positions(ch);

                }else{
                    //log(LogLevel::INFO)<<" -- duplicate update was enough"<<std::endl;
                }
            }
        }

        return surviving_v;
    }



    //tet
    void setup_boundary_adapted_init_positions(const OpenVolumeMesh::CellHandle& ch) {

        using Vec = Eigen::Vector3d;

        const bool debug(false);

        if (debug) {
            log(LogLevel::ERROR) << " ----------------- setting up codomain pre-opti init positions for element " << this->element_to_str(ch) << std::endl;
        }

        const auto& verts = this->get_cell_vertices(ch);
        auto per_vertex_pos = this->per_vertex_pos_prop(ch);
        int boundary_vertices_count(0);

        //log(LogLevel::INFO)<<" - current duplicate positions:"<<std::endl;
        std::vector<Eigen::Vector3d> P_init;
        for (int i(0); i < verts.size(); i++) {
            auto vh = verts[i];
            boundary_vertices_count += this->is_boundary(vh);
            //log(LogLevel::INFO)<<" -- "<<i<<" -> "<<vh<<" at "<<per_vertex_pos[i].transpose()<<std::endl;
        }

        if (boundary_vertices_count == 1 || boundary_vertices_count == 2){
            //translation
            //log(LogLevel::INFO) << " -> one or two boundary vertices, doing a translation" << std::endl;

            std::vector<std::pair<OVM::VertexHandle, int>> boundary_vertices;
            for (int local_vid(0); local_vid < verts.size(); local_vid++) {
                const auto &vh = verts[local_vid];
                if (this->is_boundary(vh)) {
                    boundary_vertices.push_back({vh, local_vid});
                    //log(LogLevel::INFO)<<" -- added boundary vertex "<<vh<<" with local id "<<local_vid<<" to list"<<std::endl;
                }
            }

            auto    anchor_v = boundary_vertices[0].first;
            auto anchor_v_id = boundary_vertices[0].second;
            Vec anchor_pos = vec2vec(this->vertex(anchor_v));
            auto current_anchor_pos = per_vertex_pos[anchor_v_id];
            Vec t = anchor_pos - current_anchor_pos;
            //log(LogLevel::INFO) << " -- using boundary vertex " << anchor_v << " as anchor, initially located at "<<current_anchor_pos.transpose() << std::endl;
            //log(LogLevel::INFO) << " -- t: " << t.transpose() << std::endl;

            for (int local_vid(0); local_vid < verts.size(); local_vid++) {
                //log(LogLevel::INFO)<<" -- translated vertex "<<verts[local_vid]<<" from "<<per_vertex_pos[local_vid].transpose()<<" to "<<(per_vertex_pos[local_vid] + t).transpose()<<std::endl;
                per_vertex_pos[local_vid] += t;
            }

            if (boundary_vertices_count == 2) {
                //log(LogLevel::INFO)<<" - adding scaling and rotation to match the boundary edge"<<std::endl;

                auto rotating_v    = boundary_vertices[1].first;
                auto rotating_v_id = boundary_vertices[1].second;

                //log(LogLevel::INFO)<<" - rotating vertex is "<<rotating_v<<std::endl;

                double current_length = (per_vertex_pos[rotating_v_id] - per_vertex_pos[anchor_v_id]).norm();
                //log(LogLevel::INFO)<<" - current length = "<<current_length<<std::endl;

                Vec target_edge = vec2vec(this->vertex(rotating_v) - this->vertex(anchor_v));
                double target_length = target_edge.norm();
                //log(LogLevel::INFO)<<" - target length = "<<target_length<<std::endl;
                //log(LogLevel::INFO)<<" - target edge: "<<target_edge.transpose()<<std::endl;

                double scaling = target_length / current_length;
                //log(LogLevel::INFO)<<" -> scaling = "<<scaling<<std::endl;

                for(int i(0); i<verts.size(); i++){
                    Vec new_pos = anchor_pos + scaling * (per_vertex_pos[i] - anchor_pos);
                    //log(LogLevel::INFO)<<" -- scaled vertex "<<verts[i]<<" from "<<per_vertex_pos[i].transpose()<<" to "<<new_pos.transpose()<<std::endl;
                    per_vertex_pos[i] = new_pos;
                }

                Vec current_edge = (per_vertex_pos[rotating_v_id] - per_vertex_pos[anchor_v_id]);
                //log(LogLevel::INFO)<<" - current edge: "<<current_edge.transpose()<<std::endl;
                //log(LogLevel::INFO)<<" - scaled norm = "<<current_edge.norm()<<std::endl;

                double angle = std::acos(target_edge.normalized().dot(current_edge.normalized()));
                //log(LogLevel::INFO)<<" - rotation angle = "<<angle<<std::endl;
                Vec axis = target_edge.cross(current_edge).normalized();
                //log(LogLevel::INFO)<<" - around axis: "<<axis.transpose()<<std::endl;
                Eigen::AngleAxis<double> rot(-angle, axis);
                Eigen::Matrix3d rot_mat = rot.toRotationMatrix();

                for(int i(0); i<verts.size(); i++){
                    Vec new_pos = anchor_pos + rot_mat * (per_vertex_pos[i] - anchor_pos);
                    //log(LogLevel::INFO)<<" -- rotating vertex "<<verts[i]<<" at "<<per_vertex_pos[i].transpose()<<std::endl;
                    //log(LogLevel::INFO)<<" --         edge vector: "<<(per_vertex_pos[i] - anchor_pos).transpose()<<std::endl;
                    //log(LogLevel::INFO)<<" -- rotated edge vector: "<<(rot_mat * (per_vertex_pos[i] - anchor_pos)).transpose()<<std::endl;
                    per_vertex_pos[i] = new_pos;

                    if(this->is_boundary(verts[i])) {
                        if ((per_vertex_pos[i] - vec2vec(this->vertex(verts[i]))).norm() > 1e-7) {
                            log(LogLevel::ERROR) << " ERROR - boundary vertex " << verts[i] << " is at "
                                      << per_vertex_pos[i].transpose() << " but should be at " << this->vertex(verts[i])
                                      << std::endl;
                            exit(EXIT_FAILURE);
                        }
                        //snapping to get the positions exactly right
                        per_vertex_pos[i] = vec2vec(this->vertex(verts[i]));
                    }
                }
            }
        }else if(boundary_vertices_count == 3) {
            //log(LogLevel::INFO) << " -> three boundary vertices, using some arbitrary height for now" << std::endl;
            OVM::HalfFaceHandle boundary_hfh(-1);
            for(auto chf_it = this->chf_iter(ch); chf_it.valid(); chf_it++){
                if(TetrahedralMesh::is_boundary(opposite_halfface_handle(*chf_it))){
                    boundary_hfh = *chf_it;
                    break;
                }
            }

            auto hf_verts = this->get_halfface_vertices(boundary_hfh);
            Vec3d normal = (this->vertex(hf_verts[1]) - this->vertex(hf_verts[0])).cross(this->vertex(hf_verts[2]) - this->vertex(hf_verts[0])).normalized();
            Vec3d centroid = 1.0 / 3.0 * (this->vertex(hf_verts[0]) + this->vertex(hf_verts[1]) + this->vertex(hf_verts[2]));

            //log(LogLevel::INFO)<<" - normal is "<<normal<<std::endl;
            //log(LogLevel::INFO)<<" - centroid at "<<centroid<<std::endl;

            double vol = this->non_conforming_element_volume(ch);
            double boundary_area = compute_area(*this, boundary_hfh);
            double target_height = vol / boundary_area;

            //log(LogLevel::INFO)<<" - non-conforming element volume is "<<vol<<" and target boundary area is "<<boundary_area<<", giving a target height of "<<target_height<<std::endl;

            Vec3d interior_v_pos = centroid + target_height * normal;
            //log(LogLevel::INFO)<<" - interior vertex position is "<<interior_v_pos<<std::endl;

            for(int i(0); i<verts.size(); i++){
                if(this->is_boundary(verts[i])){
                    per_vertex_pos[i] = vec2vec(this->vertex(verts[i]));
                }else {
                    per_vertex_pos[i] = vec2vec(interior_v_pos);
                }
            }

        }else {
            //log(LogLevel::INFO) << " interior element, no need to do any extra init stuff" << std::endl;
        }

        this->per_vertex_pos_prop(ch) = per_vertex_pos;
        //log(LogLevel::INFO)<<" - volume after setup: "<<this->non_conforming_element_volume(ch)<<std::endl;

        //log(LogLevel::INFO)<<" ---------- "<<std::endl;
    }



    void set_regular_element_as_reference(double vol = 1.0){
        log(LogLevel::ERROR)<<" ERROR - not implemented"<<std::endl;
        exit(EXIT_FAILURE);
    }


    virtual void gather_neighbor_elements(const OVM::VertexHandle& vh,
                                          std::vector<OVM::CellHandle>& neighbors) const {
        for(auto vc_it = this->vc_iter(vh); vc_it.valid(); vc_it++){
            neighbors.push_back(*vc_it);
        }
    }

    virtual void gather_neighbor_elements(const OVM::EdgeHandle& eh,
                                          std::vector<OVM::CellHandle>& neighbors) const {

        for(auto ec_it = this->ec_iter(eh); ec_it.valid(); ec_it++){
            neighbors.push_back(*ec_it);
        }
    }

    //tet
    virtual MeshAngles compute_angles_in_neighborhood(const OVM::VertexHandle& vh,
                                                      ElementProperty<bool>& to_skip_prop,
                                                      VERTEX_POSITION_SOURCE pos_source) const {


        MeshAngles angles;

        for(auto vc_it = this->vc_iter(vh); vc_it.valid(); vc_it++){
            if(!to_skip_prop[*vc_it]){

                if(is_flipped_or_degenerate(*vc_it, pos_source)){
                    log(LogLevel::INFO)<<" --> cannot compute dihedral angles in inverted or degenerate element "<<this->element_to_str(*vc_it)<<" with volume "<<CGAL::to_double(CGAL_wrapper::exact_volume(*this, *vc_it))<<std::endl;
                    log(LogLevel::INFO)<<" - position source: "<<pos_source<<std::endl;
                    exit(EXIT_FAILURE);
                }

                const OVM::CellHandle* p_ch = (pos_source == CONFORMING ? nullptr :  &(*vc_it));

                for(auto ce_it = this->ce_iter(*vc_it); ce_it.valid(); ce_it++){

                    auto angle_deg = this->compute_dihedral_angle_deg(*vc_it, *ce_it, pos_source);

                    angles.min_dihedral_angle_deg = std::min(angles.min_dihedral_angle_deg, angle_deg);
                    angles.max_dihedral_angle_deg = std::max(angles.max_dihedral_angle_deg, angle_deg);
                }

                //NOTE: not very efficient since each face will be dealt with twice but that's ok
                for(auto cf_it = this->cf_iter(*vc_it); cf_it.valid(); cf_it++){
                    for(auto fv_it = this->fv_iter(*cf_it); fv_it.valid(); fv_it++){
                        auto angle_deg = this->compute_interior_angle_deg(*cf_it, *fv_it, pos_source, p_ch);

                        angles.min_interior_angle_deg = std::min(angles.min_interior_angle_deg, angle_deg);
                        angles.max_interior_angle_deg = std::max(angles.max_interior_angle_deg, angle_deg);
                    }
                }
            }
        }

        return angles;
    }

    //tet
    virtual MeshAngles compute_angles_in_neighborhood(const OVM::EdgeHandle& eh,
                                                      ElementProperty<bool>& to_skip_prop,
                                                      VERTEX_POSITION_SOURCE pos_source) const {


        MeshAngles angles;

        auto from_v = this->edge(eh).from_vertex();
        auto   to_v = this->edge(eh).to_vertex();

        auto from_v_angles = compute_angles_in_neighborhood(from_v, to_skip_prop, pos_source);
        auto   to_v_angles = compute_angles_in_neighborhood(  to_v, to_skip_prop, pos_source);


        angles.min_dihedral_angle_deg = std::min(from_v_angles.min_dihedral_angle_deg, to_v_angles.min_dihedral_angle_deg);
        angles.max_dihedral_angle_deg = std::max(from_v_angles.max_dihedral_angle_deg, to_v_angles.max_dihedral_angle_deg);

        angles.min_interior_angle_deg = std::min(from_v_angles.min_interior_angle_deg, to_v_angles.min_interior_angle_deg);
        angles.max_interior_angle_deg = std::max(from_v_angles.max_interior_angle_deg, to_v_angles.max_interior_angle_deg);


        return angles;
    }


    virtual MeshAngles compute_angles_in_neighborhood(const OVM::FaceHandle& fh,
                                                      VERTEX_POSITION_SOURCE pos_source) const {

        MeshAngles angles;

        std::vector<OVM::CellHandle> incident_cells = {this->incident_cell(this->halfface_handle(fh, 0)),
                                                  this->incident_cell(this->halfface_handle(fh, 1))};

        for(auto ch: incident_cells) {
            if (!ch.is_valid()) {
                if (!this->is_boundary(fh)) {
                    log(LogLevel::ERROR) << " ERROR - couldn't find one of the incident cells and it's not a boundary face"
                              << std::endl;
                    exit(EXIT_FAILURE);
                }
            } else {

                if (is_flipped_or_degenerate(ch, pos_source)) {
                    log(LogLevel::INFO) << " --> cannot compute dihedral angles in inverted or degenerate element "
                              << this->element_to_str(ch) << std::endl;
                    exit(EXIT_FAILURE);
                }

                for (auto ce_it = this->ce_iter(ch); ce_it.valid(); ce_it++) {

                    auto angle_deg = this->compute_dihedral_angle_deg(ch, *ce_it, pos_source);
                    angles.min_dihedral_angle_deg = std::min(angles.min_dihedral_angle_deg, angle_deg);
                    angles.max_dihedral_angle_deg = std::max(angles.max_dihedral_angle_deg, angle_deg);
                }


                const OVM::CellHandle* p_ch = (pos_source == CONFORMING ? nullptr :  &ch);

                //NOTE: not very efficient since each face will be dealt with twice but that's ok
                for (auto cf_it = this->cf_iter(ch); cf_it.valid(); cf_it++) {
                    for (auto fv_it = this->fv_iter(*cf_it); fv_it.valid(); fv_it++) {
                        auto angle_deg = this->compute_interior_angle_deg(*cf_it, *fv_it, pos_source, p_ch);

                        angles.min_interior_angle_deg = std::min(angles.min_interior_angle_deg, angle_deg);
                        angles.max_interior_angle_deg = std::max(angles.max_interior_angle_deg, angle_deg);
                    }
                }
            }
        }


        return angles;
    }



    virtual double compute_dihedral_angle_deg(const OVM::CellHandle ch,
                                              const OVM::EdgeHandle& eh,
                                              VERTEX_POSITION_SOURCE pos_source) const{

        OVM::HalfFaceHandle hf1(-1);
        OVM::HalfFaceHandle hf2(-1);
        for(auto ehf_it = this->ehf_iter(eh); ehf_it.valid(); ehf_it++){
            if(this->incident_cell(*ehf_it) == ch){
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

            auto face_vertices = this->get_halfface_vertices(hf);

            std::vector<Vec3d> positions;
            for(auto v: face_vertices){
                auto incident_c = this->incident_cell(hf);
                auto incident_c_vertices = this->get_cell_vertices(incident_c);
                for(int i(0); i < incident_c_vertices.size(); i++){
                    if(incident_c_vertices[i] == v){

                        switch(pos_source){
                            case CONFORMING:{
                                positions.push_back(this->vertex(v));
                                break;
                            }
                            case NON_CONFORMING:{
                                positions.push_back(vec2vec(per_vertex_pos_prop(incident_c)[i]));
                                break;
                            }
                            case REFERENCE:{
                                positions.push_back(vec2vec(per_vertex_pos_prop(incident_c)[i]));
                                break;
                            }
                        }
                    }
                }
            }

            Vec3d c_evi  = positions[1] - positions[0];
            Vec3d c_evi1 = positions[2] - positions[0];

            /*log(LogLevel::INFO)<<" - hf positions: "<<std::endl;
            for(auto v: face_vertices){
                log(LogLevel::INFO)<<"  -- "<<v<<" at "<<this->vertex(v)<<std::endl;
            }*/

            Vec3d normal = (c_evi.cross(c_evi1));
            normal.normalize();

            normals.push_back(normal);

            //log(LogLevel::INFO)<<" - normal for face "<<this->get_halfface_vertices(hf)<<" = "<<normal<<std::endl;
            //log(LogLevel::INFO)<<" ---------"<<std::endl;
        }

        auto dihedral_angle = (M_PI - acos(normals[0].dot(normals[1])));

        return 180.0 * dihedral_angle / M_PI;
    }


    //tet
    virtual QualityVector compute_quality_vector(const OpenVolumeMesh::CellHandle& ch,
                                               QUALITY_METRIC metric,
                                               VERTEX_POSITION_SOURCE pos_source) const{

        //UNKNOWN_METRIC=-1, MIN_DIHEDRAL_ANGLE, MAX_DIHEDRAL_ANGLE_DIFF_TO_PI, MIN_DIHEDRAL_ANGLE_DIFF, VOLUME_LENGTH, WICKE

        QualityVector vec;
        switch(metric){

            case DIHEDRAL_ANGLE:{

                for (auto ce_it = this->ce_iter(ch); ce_it.valid(); ce_it++) {
                    auto angle_deg = this->compute_dihedral_angle_deg(ch, *ce_it, pos_source);
                    vec.insert_sort(angle_deg);
                }
                break;
            }

            case NEGATIVE_DIHEDRAL_ANGLE: {

                for (auto ce_it = this->ce_iter(ch); ce_it.valid(); ce_it++) {
                    auto angle_deg = this->compute_dihedral_angle_deg(ch, *ce_it, pos_source);
                    vec.insert_sort(-angle_deg);
                }
                break;
            }

            case DIHEDRAL_ANGLE_DIFF: {

                for (auto ce_it = this->ce_iter(ch); ce_it.valid(); ce_it++) {
                    auto angle_deg = this->compute_dihedral_angle_deg(ch, *ce_it, pos_source);
                    vec.insert_sort(std::min(180.0 - angle_deg, angle_deg));
                }
                break;
            }

            case VOLUME_LENGTH: {

                double vol(0);
                if(pos_source == CONFORMING){
                    vol = this->element_volume(ch);
                }else if(pos_source == NON_CONFORMING){
                    vol = this->non_conforming_element_volume(ch);
                }else{
                    log(LogLevel::ERROR)<<" ERROR - volume length quality is not implemented using position source "<<pos_source<<std::endl;
                    exit(EXIT_FAILURE);
                }

                double sum_of_squares(0);
                double edge_val(0);
                for (auto ce_it = this->ce_iter(ch); ce_it.valid(); ce_it++) {
                    auto from_v = this->edge(*ce_it).from_vertex();
                    auto   to_v = this->edge(*ce_it).to_vertex();

                    double edge_len(0);
                    if(pos_source == CONFORMING){
                        edge_len = (this->vertex(from_v) - this->vertex(to_v)).norm();
                    }else{
                        auto from_v_id = get_duplicate_local_id(ch, from_v);
                        auto   to_v_id = get_duplicate_local_id(ch,   to_v);
                        edge_len = (this->per_vertex_pos_prop(ch)[from_v_id] - this->per_vertex_pos_prop(ch)[to_v_id]).norm();
                    }

                    sum_of_squares += edge_len * edge_len;

                    edge_val++;
                }
                double l_RMS = std::sqrt(sum_of_squares /= edge_val);

                double vl = 6.0 * std::sqrt(2.0) * vol / (l_RMS * l_RMS * l_RMS);

                vec.insert_sort(vl);

                break;
            }

            case CONSTANT_METRIC:{
                vec.insert_sort(0.0);
                break;
            }

            default:{
                log(LogLevel::ERROR)<<" ERROR - asking for unhandled metric "<<metric<<std::endl;
                exit(EXIT_FAILURE);
            }
        }

        return vec;
    }

    virtual bool is_boundary(const OVM::EdgeHandle& eh) const{
        return TetrahedralMesh::is_boundary(eh);
    }

    virtual bool is_boundary(const OVM::VertexHandle& vh) const{
        return TetrahedralMesh::is_boundary(vh);
    }

    virtual void setup_ref_shape_position_matrix(VertexMatrix& P_ref){


        P_ref<< std::sqrt(8.)/3., -std::sqrt(2.)/3, -std::sqrt(2.)/3.,         0.0,
                0.0,  std::sqrt(2./3.),-std::sqrt(2./3.),         0.0,
                0.0,               0.0,              0.0, 1.0+1.0/3.0;
    }

    virtual std::vector<EigenVec>& per_vertex_pos_prop(const OVM::CellHandle& ch){
        //log(LogLevel::INFO)<<" tet getter"<<std::endl;
        return per_vertex_pos_prop_[ch];
    }
    virtual std::vector<EigenVec> per_vertex_pos_prop(const OVM::CellHandle& ch) const {
        //log(LogLevel::INFO)<<" tet getter"<<std::endl;
        return per_vertex_pos_prop_[ch];
    }


    virtual std::vector<EigenVec>& per_vertex_ref_pos_prop(const OVM::CellHandle& ch){
        return per_vertex_ref_pos_prop_[ch];
    }
    virtual std::vector<EigenVec> per_vertex_ref_pos_prop(const OVM::CellHandle& ch) const {
        return per_vertex_ref_pos_prop_[ch];
    }


    /*virtual std::vector<bool>& glued_duplicates_prop(const OVM::CellHandle& ch){
        return glued_duplicates_prop_[ch];
    }
    virtual std::vector<bool> glued_duplicates_prop(const OVM::CellHandle& ch) const{
        return glued_duplicates_prop_[ch];
    }*/

    virtual std::vector<OVM::VertexHandle> get_element_vertices(OVM::CellHandle ch) const {
        return this->get_cell_vertices(ch);
    }


    virtual std::pair<OVM::CellIter, OVM::CellIter> elements() const {
        return this->cells();
    }


    virtual int n_logical_elements() const {
        return this->n_logical_cells();
    }

    virtual int n_elements() const {
        return this->n_cells();
    }


private:


    OVM::VertexPropertyT<std::pair<OVM::VertexHandle, OVM::VertexHandle>> collapse_axis_prop_;

};




template < class VecT, class MatT>
static void SVD(const MatT& M,
                MatT& U,
                VecT& S,
                MatT& V){

    using namespace Eigen;

    /*JacobiSVD<MatT> svd;
    svd.compute(M, ComputeThinU | ComputeThinV);*/

    JacobiSVD<MatT, ComputeThinU | ComputeThinV> svd;
    svd.compute(M);
    /*log(LogLevel::INFO)<<" ------------------------------------"<<std::endl;
    log(LogLevel::INFO)<<"  M  = "<<std::endl<<M<<std::endl;
    log(LogLevel::INFO)<<" |M| = "<<M.determinant()<<std::endl;
    log(LogLevel::INFO)<<"  U  = "<<std::endl<<svd.matrixU()<<std::endl;
    log(LogLevel::INFO)<<"  S  = "<<std::endl<<svd.singularValues()<<std::endl;
    log(LogLevel::INFO)<<"  V  = "<<std::endl<<svd.matrixV()<<std::endl;*/

    //U.resize(M.rows())

    U = svd.matrixU();
    S = svd.singularValues();
    V = svd.matrixV();

    /*log(LogLevel::INFO)<<" |U| = "<<U.determinant()<<std::endl;
    log(LogLevel::INFO)<<" |S| = "<<S.prod()<<std::endl;
    log(LogLevel::INFO)<<" |V| = "<<V.determinant()<<std::endl;
    log(LogLevel::INFO)<<" ------------------------------------"<<std::endl;*/
}

template < class VecT, class MatT>
void signed_SVD(const MatT& M,
                MatT& U,
                VecT& S,
                MatT& V){

    SVD(M, U,S,V);

    auto det_U = U.determinant();
    auto det_V = V.determinant();

    auto det_M = det_U * S.prod() * det_V;

    int smallest_index((int)S.size()-1);

    //log(LogLevel::INFO)<<" det M = "<<det_M<<std::endl;
    //log(LogLevel::INFO)<<" smallest index = "<<smallest_index<<std::endl;

    if(det_M > 0){
        //probably better to check sign rather than compare to -1.0
        if(det_U < 0 && det_V < 0){
            U.col(0) *= -1;
            V.col(0) *= -1;
        }//else, all good
    }else if (det_M < 0){
        //log(LogLevel::INFO)<<" - negative det"<<std::endl;
        if(det_U < 0){
            S(smallest_index) *= -1;
            U.col(smallest_index) *= -1;
        }else if (det_V < 0){
            /*log(LogLevel::INFO)<<" negative det V"<<std::endl;
            log(LogLevel::INFO)<<" S size = "<<S.size()<<std::endl;
            log(LogLevel::INFO)<<" V size = "<<V.cols()<<std::endl;
            log(LogLevel::INFO)<<" col: "<<V.col(smallest_index).transpose()<<std::endl;*/
            S(smallest_index) *= -1;
            V.col(smallest_index) *= -1;

        }else{
            log(LogLevel::INFO)<<" - unhandled case with |U| = |S| = |V| = -1"<<std::endl;
            std::string ans;
            std::cin>>ans;
        }
    }else{
        for(int i(0); i<S.size(); i++){
            //log(LogLevel::INFO)<<" S_"<<i<<" = "<<S(i)<<std::endl;
            if(S(i) == 0){
                if(U.determinant() < 0){
                    U.col(i) *= -1;
                }
                if(V.determinant() < 0){
                    V.col(i) *= -1;
                }
            }
        }
    }

    /*log(LogLevel::INFO)<<" - updated determinants: "<<std::endl;
    log(LogLevel::INFO)<<" |U| = "<<U.determinant()<<std::endl;
    log(LogLevel::INFO)<<" |S| = "<<S.prod()<<std::endl;
    log(LogLevel::INFO)<<" |V| = "<<V.determinant()<<std::endl;
    log(LogLevel::INFO)<<" ------------------------------------"<<std::endl;*/
}



#if !defined(NON_CONFORMING_MESHES_C)
#define NON_CONFORMING_MESHES_TEMPLATES
#include "NonConformingMeshesT.cc"
#endif
