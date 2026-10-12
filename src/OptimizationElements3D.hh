/*===========================================================================*\
 *                                                                           *
 *                               Optimization Elements                       *
 *      Copyright (C) 2021 by Computer Graphics Group, University of Bern    *
 *                           http://cgg.unibe.ch                             *
 *                                                                           *
\*===========================================================================*/


#pragma once

//== INCLUDES =================================================================

#include <fstream>

#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<Eigen/StdVector>
#include<cmath>

#include <CoMISo/NSolver/FiniteElementTinyAD.hh>
#include <CoMISo/NSolver/FiniteElementHessianProjection.hh>
#include <CoMISo/Utils/Polynomials.hh>



//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================


class EdgeLengthBarrier3DElement
{
public:

    // define dimensions
    const static int NV = 6; // [ax,ay,az,bx,by,bz] = two 3D end points of an edge
    const static int NC = 3; // [weight, bound, target]

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    static std::string element_name(){ return "Edge Length Bound"; }
    static std::string element_code(){ return "ELB"; }


    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {

        using VecT  = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        using MatT  = Eigen::Matrix<ScalarT,3,3,Eigen::ColMajor>;
        using MatdT = Eigen::Matrix<double,3,3,Eigen::RowMajor>;

        // get constants
        const double& w = _c[0];
        const double& len_bound = _c[1];
        const double& len_target = _c[2];


        // get input points
        Eigen::Map<VecT> A( (ScalarT*) &(_x[0]) );
        Eigen::Map<VecT> B( (ScalarT*) &(_x[3]) );

        /*std::cerr<<" ----------- "<<std::endl;
        std::cerr<<"         A: "<<A.transpose()<<std::endl;
        std::cerr<<"         B: "<<B.transpose()<<std::endl;
        std::cerr<<"      diff: "<<(B - A).transpose()<<std::endl;
        std::cerr<<"   length = "<<(B - A).norm()<<std::endl;
        std::cerr<<" ----------- "<<std::endl;*/

        ScalarT len = (B-A).norm();

        if (w > 0.0) {
            ScalarT E_bound = len - len_bound;

            if (len <= len_bound) {
                return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
            } else {
                return w * ((len_target - len_bound) / (len-len_bound) + (len - len_bound) / (len_target - len_bound) - 2);
            }
        }else {
            return 0.0;
        }
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
    {
        double t=2.0;
        for(int i=0; i<50; ++i)
        {
            VecV xn = _x+t*_v;
            if(std::isfinite(eval_f(xn, _c)))
                return t;
            else
                t *= 0.7;
        }

        return 0;
    }
};

using ELEN_B_3D_AD = COMISO::FiniteElementTinyAD<EdgeLengthBarrier3DElement>;
//using ELEN_B_3D_PH_AD = COMISO::FiniteElementHessianProjection< ELEN_B_3D_AD >;






class EdgeVec3DElement
{
public:

    // define dimensions
    const static int NV = 6; // [ax,ay,az,bx,by,bz] = two 3D end points of an edge
    const static int NC = 4; // [weight, target edge vector]

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    static std::string element_name(){ return "Edge Vec Diff"; }


    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {

        using VecT  = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        using MatT  = Eigen::Matrix<ScalarT,3,3,Eigen::ColMajor>;
        using MatdT = Eigen::Matrix<double,3,3,Eigen::RowMajor>;

        // get constants
        const double& w = _c[0];


        Eigen::Vector3d target_vec;
        target_vec <<_c[1], _c[2], _c[3];

        // get input points
        Eigen::Map<VecT> A( (ScalarT*) &(_x[0]) );
        Eigen::Map<VecT> B( (ScalarT*) &(_x[3]) );

        /*std::cerr<<" ----------- "<<std::endl;
        std::cerr<<" A: "<<A.transpose()<<std::endl;
        std::cerr<<" B: "<<B.transpose()<<std::endl;
        std::cerr<<" target: "<<target_vec.transpose()<<std::endl;
        std::cerr<<"   diff: "<<((B - A) - target_vec).transpose()<<std::endl;
        std::cerr<<" ----------- "<<std::endl;*/

        return w * ((B - A) - target_vec).squaredNorm();
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
    {
        return DBL_MAX;
    }
};

using EV3D_AD = COMISO::FiniteElementTinyAD<EdgeVec3DElement>;
//using EV3D_PH_AD = COMISO::FiniteElementHessianProjection< EV3D_AD >;



//== BASE FOLDOVER-FREE ELEMENT ===============================================================

class FoldoverFree_F
{
public:

    // define dimensions
    const static int NV = 12; // [ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz] = 3D points A,B,C,D of a tetrahedron
    const static int NC = 4; // [weight_conformal, weight_det, epsilon, average volume]

    static std::string element_name(){ return "Foldover-Free"; }
    static std::string element_code(){ return "FOF"; }


    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {

        Eigen::Vector3d v0, v1, v2, v3;
        v0 << 1.0, 1.0, 1.0;
        v1 << 1.0, -1.0, -1.0;
        v2 << -1.0, -1.0, 1.0;
        v3 << -1.0, 1.0, -1.0;
        Eigen::Matrix3d Shape;
        Shape.col(0) = v1 - v0; Shape.col(1) = v2 - v0; Shape.col(2) = v3 - v0;

        using VecT  = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        using MatT  = Eigen::Matrix<ScalarT,3,3,Eigen::ColMajor>;
        using MatdT = Eigen::Matrix<double,3,3,Eigen::RowMajor>;

        // get constants
        const double& w_conformal = _c[0];
        const double& w_det       = _c[1];
        const double& epsilon     = _c[2];
        const double& average_vol = _c[3];

        // get input points
        Eigen::Map<VecT> A( (ScalarT*) &(_x[0]) );
        Eigen::Map<VecT> B( (ScalarT*) &(_x[3]) );
        Eigen::Map<VecT> C( (ScalarT*) &(_x[6]) );
        Eigen::Map<VecT> D( (ScalarT*) &(_x[9]) );

        // 3./8. to normalize shape matrix to unit tet
        Shape = Shape * pow(average_vol, 1./3.) * 3./8.;

        /*std::cout<<" -----------------------------"<<std::endl;
        std::cout<<" - shape = "<<Shape<<std::endl;
        std::cout<<" - shape det = "<<Shape.determinant()<<std::endl;*/

        // calc Jacobi Matrix
        MatT J;
        J.col(0) = B-A;
        J.col(1) = C-A;
        J.col(2) = D-A;
        //std::cout<<" - Jacobian = "<<J<<std::endl;
        //std::cout<<" - Jacobian det = "<<J.determinant()<<std::endl;

        J = J * Shape.inverse();
        // calculate modified determinant first
        ScalarT det = J.determinant();

        /*std::cout<<" - mapped J = "<<J<<std::endl;
        std::cout<<" - mapped J det = "<<J.determinant()<<std::endl;
        std::cout<<" -----------------------------"<<std::endl;*/

        if(epsilon > 0.0)
        {
            ScalarT det_reg = 0.5*(det + sqrt( pow(det,2) + epsilon*epsilon));
            return w_conformal*(J.transpose()*J).trace()/pow(det_reg,2./3.) + w_det*((pow(det,2)+1.0)/det_reg);
        }
        else
        {
            if (det <= 0.0)
                return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
            else
                return w_conformal*(J.transpose()*J).trace()/pow(det,2./3.) + w_det*(det + 1.0/det);
        }
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
    {
        return DBL_MAX;
    }
};

using FF_AD = COMISO::FiniteElementTinyAD<FoldoverFree_F>;
using FF_PH_AD = COMISO::FiniteElementHessianProjection< FF_AD >;




class WickeEnergyElement
{
public:
    const static int NV = 12; // positions of 4 vertices
    const static int NC = 1; // weight

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x,
                          const Eigen::Matrix<double, NC, 1> &_c) const
    {

        using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        using Mat3T    = Eigen::Matrix<ScalarT,3,3,Eigen::ColMajor>;

        Vec3T p0 =  {_x[0], _x[1], _x[2]};
        Vec3T p1 =  {_x[3], _x[4], _x[5]};
        Vec3T p2 =  {_x[6], _x[7], _x[8]};
        Vec3T p3 =  {_x[9], _x[10], _x[11]};

        ScalarT w = _c[0];

        /*std::cout<<" - p0: "<<p0.transpose()<<std::endl;
        std::cout<<" - p1: "<<p1.transpose()<<std::endl;
        std::cout<<" - p2: "<<p2.transpose()<<std::endl;
        std::cout<<" - p3: "<<p3.transpose()<<std::endl;*/

        Mat3T E;
        E.col(0) = p1-p0;
        E.col(1) = p2-p0;
        E.col(2) = p3-p0;

        //std::cout<<" E: "<<E<<std::endl;

        //order and sign don't matter since we only care about the lengths
        Mat3T E2;
        E2.col(0) = p1-p2;
        E2.col(1) = p2-p3;
        E2.col(2) = p3-p1;

        //std::cout<<" E2: "<<E2<<std::endl;

        ScalarT vol = E.determinant() / 6.0;
        //std::cout<<" vol = "<<vol<<std::endl;

        if(vol <= 0){
            return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
        }

        ScalarT l_harm(0);
        ScalarT l_rms(0);

        for(int i(0); i<6; i++){

            ScalarT len(0);

            //not super clean but hey...
            if(i<3){
                len = E.col(i).norm();
            }else{
                len = E2.col(i-3).norm();
            }
            //std::cout<<" len "<<i<<" = "<<len<<std::endl;

            if(len == 0){
                return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
            }

            l_harm += 1.0/len;
            l_rms += len * len;
        }

        l_harm = 6.0 / l_harm;
        l_rms = sqrt(l_rms/6.0);

        // actually 1/wicke to penalize way more the elements near zero
        auto wicke = (6.0 * sqrt(2.0) * vol * l_harm) / (l_rms * l_rms * l_rms * l_rms) ;
        //std::cout<<" wicke: "<<wicke<<std::endl;
        auto e = - w * log(wicke);
        //std::cout<<" - e = "<<e<<std::endl;
        return e;
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        // HACK (add exact computation or better line search)
        double t=2.0;
        for(int i=0; i<50; ++i)
        {
            VecV xn = _x+t*_v;
            if(std::isfinite(eval_f(xn, _c)))
                return t;
            else
                t *= 0.7;
        }

        return 0;
    }

};

using WICKE_AD = COMISO::FiniteElementTinyAD<WickeEnergyElement>;
using WICKE_PH_AD = COMISO::FiniteElementHessianProjection< WICKE_AD >;


#if 0
class MinConditionNumberElement
{
public:
    // define dimensions
    const static int NV = 12; // positions of 4 vertices
    const static int NC = 18; // reference shape (edge matrix) and inverse of reference shape

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x,
                          const Eigen::Matrix<double, NC, 1> &_c) const
    {

        using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        using Mat3T    = Eigen::Matrix<ScalarT,3,3,Eigen::ColMajor>;

        Mat3T E_ref = {_c[0], _c[1], _c[2],
                       _c[3], _c[4], _c[5],
                       _c[6], _c[7], _c[8]};


        Mat3T E_ref_inv = {_c[9], _c[10], _c[11],
                           _c[12], _c[13], _c[14],
                           _c[15], _c[16], _c[7]};


        Vec3T p0 =  {_x[0], _x[1], _x[2]};
        Vec3T p1 =  {_x[3], _x[4], _x[5]};
        Vec3T p2 =  {_x[6], _x[7], _x[8]};
        Vec3T p3 =  {_x[9], _x[10], _x[11]};

        Mat3T E;
        E.col(0) = p1-p0;
        E.col(1) = p2-p0;
        E.col(2) = p3-p0;


        Mat3T J = E * E_ref_inv;

        Mat3T U, V;
        Vec3T S;
        //ShrinkAndExpand::PEHelpers::signed_SVD(J, U,S,V);



        return 0;
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        return std::numeric_limits<double>::max();
    }

  };

  using MIN_K_AD = COMISO::FiniteElementTinyAD<MinConditionNumberElement>;
  using MIN_K_PH_AD = COMISO::FiniteElementHessianProjection< MIN_K_AD >;

#endif


class MinDistanceTranslationElement_F
{
public:
    // define dimensions
    const static int NV = 15; // positions of 4 anchor points, displacement vector
    const static int NC = 12; // positions of 4 vertices

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x,
                          const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;

        Eigen::Map<Vec3T> t((ScalarT*) &(_x[12]) );

        //std::cout<<" -------------- "<<std::endl;
        //std::cout<<" - translation: "<<t.transpose()<<std::endl;


        ScalarT e(0);
        for(int i(0); i<4; i++){
            Eigen::Map<Vec3T> a((ScalarT*) &(_x[3*i]) );
            Eigen::Vector3d p;
            p<<_c[3 * i + 0], _c[3 * i + 1], _c[3 * i + 2];


            auto d = p + t - a;

            /*std::cout<<"------"<<std::endl;
            std::cout<<" - p"<<i<<": "<<p.transpose()<<std::endl;
            std::cout<<" - a"<<i<<": "<<a.transpose()<<std::endl;
            std::cout<<" - d"<<i<<": "<<d.transpose()<<std::endl;
            std::cout<<" - e"<<i<<" = "<<(d.dot(d))<<std::endl;*/

            e += d.dot(d);
        }

        return e;

    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        return std::numeric_limits<double>::max();
    }

};

using MDT_AD = COMISO::FiniteElementTinyAD<MinDistanceTranslationElement_F>;
using MDT_PH_AD = COMISO::FiniteElementHessianProjection< MDT_AD >;




//== CONSTRAINED VECTOR DIFF DEFINITION =========================================================
class ConstrainedVectorDiffElement_F
{
public:
    // define dimensions
    const static int NV = 3; // three singular values to optimize
    const static int NC = 5; // three reference singular values, minimum and penalty mu

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x,
                          const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;

        Eigen::Map<Vec3T> sigma_var((ScalarT*) &(_x[0]) );

        Eigen::Vector3d sigma_ref;
        sigma_ref<<_c[0], _c[1], _c[2];

        //std::cout<<" - sigma: "<<sigma_ref.transpose()<<std::endl;
        //std::cout<<" - sigma': "<<sigma_var.transpose()<<std::endl;

        //std::cout<<" -******************************************* mu = "<<_c[4]<<std::endl;

        auto diff = sigma_var - sigma_ref;

        auto diff_e = diff.dot(diff);

        ScalarT barrier_e = 0;
        for(int i(0); i<3; i++){
            for(int j(0); j<3; j++){
                if(i!=j){
                    barrier_e += -log(_c[3] * _x[i] - _x[j]);
                }
            }
        }
        //std::cout<<" - arg = "<<(_c[3] * _x[2] - _x[0])<<std::endl;
        //barrier_e += -log(_c[3] * _x[2] - _x[0]);
        barrier_e *= _c[4];

        //diff_e *= _c[4];

        //std::cout<<" -    diff energy = "<<diff_e<<std::endl;
        //std::cout<<" - barrier energy = "<<barrier_e<<std::endl;

        return diff_e + barrier_e;
        return diff.dot(diff) + _x.minCoeff() -_c[4] * log(-_x.maxCoeff() + _c[3] * _x.minCoeff());


    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        return std::numeric_limits<double>::max();
    }

};

using CVECDIFF_AD = COMISO::FiniteElementTinyAD<ConstrainedVectorDiffElement_F>;
using CVECDIFF_PH_AD = COMISO::FiniteElementHessianProjection< CVECDIFF_AD >;





//== TET HEIGHT DEFINITION =========================================================

class TetHeightElement_F
{
public:

    // define dimensions
    const static int NV = 12; // [] = 3D points A, B, C and D of tet
    const static int NC = 1; // [nothing really]

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        using Mat3T    = Eigen::Matrix<ScalarT,3,3,Eigen::ColMajor>;

        Eigen::Map<Vec3T> p0( (ScalarT*) &(_x[0]) );
        Eigen::Map<Vec3T> p1( (ScalarT*) &(_x[3]) );
        Eigen::Map<Vec3T> p2( (ScalarT*) &(_x[6]) );
        Eigen::Map<Vec3T> p3( (ScalarT*) &(_x[9]) );

        //std::cout<<" - element x = "<<_x.transpose()<<std::endl;

        std::ofstream solver_out;
        solver_out.open("solver_out.txt", std::fstream::app);
        solver_out<<" ----------------------"<<std::endl;
        solver_out << " - element x = "<<_x.transpose()<<std::endl;


        //compute edges
        Vec3T e1 = p1-p0;
        Vec3T e2 = p2-p0;
        Vec3T e3 = p3-p0;

        auto A1 = 0.5 * e1.cross(e2).norm();
        auto A2 = 0.5 * e1.cross(e3).norm();
        auto A3 = 0.5 * e2.cross(e3).norm();
        auto A4 = 0.5 * (p2-p1).cross(p3-p1).norm();

        Mat3T M;
        M<< e1[0], e2[0], e3[0],
                e1[1], e2[1], e3[1],
                e1[2], e2[2], e3[2];

        auto vol = M.determinant() / 6.0;

        /*solver_out << " - A1 = "<<A1<<std::endl;
        solver_out << " - A2 = "<<A2<<std::endl;
        solver_out << " - A3 = "<<A3<<std::endl;
        solver_out << " - A4 = "<<A4<<std::endl;

        solver_out << " - vol = "<<vol<<std::endl;
        solver_out << " - dbl eps "<<DBL_EPSILON<<std::endl;*/

        if(vol <= 0.0){
            return 1e15;
            return std::numeric_limits<ScalarT>::max();
        }else{

            vol *= 6;

            return (A1 + A2 + A3 + A4) / (vol);
            return (A1 * A1 + A2 * A2 + A3 * A3 + A4 * A4) / (vol * vol);
            return (A1 * A1 * A1 + A2 * A2 * A2 + A3 * A3 * A3 + A4 * A4 * A4) / (vol * vol * vol);

        }
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        return std::numeric_limits<double>::max();
    }

};

using TETH_3D_AD = COMISO::FiniteElementTinyAD<TetHeightElement_F>;
using TETH_PH_AD = COMISO::FiniteElementHessianProjection< TETH_3D_AD >;





//== CONFORMITY INITIALIZATION DEFINITION =========================================================


class InitializationElement3D_F
{
public:

    // define dimensions
    const static int NV = 2; // t_i, t_j
    const static int NC = 12; // [c_i, n_i, c_j, n_j] all in 3D

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {

        using Vec3d    = Eigen::Matrix<double,3,1,Eigen::ColMajor>;

        Eigen::Map<Vec3d> ci( (double*) &(_c[0]) );
        Eigen::Map<Vec3d> ni( (double*) &(_c[3]) );
        Eigen::Map<Vec3d> cj( (double*) &(_c[6]) );
        Eigen::Map<Vec3d> nj( (double*) &(_c[9]) );

        auto pi = ci + _x[0] * ni;
        auto pj = cj + _x[1] * nj;


        return (pi-pj).dot(pi-pj);
        return _x[0];
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        return std::numeric_limits<double>::max();
    }

};

using INIT_AD = COMISO::FiniteElementTinyAD<InitializationElement3D_F>;
using INIT_PH_AD = COMISO::FiniteElementHessianProjection< INIT_AD >;




//== LOG BARRIER DEFINITION =========================================================


class LogBarrierTerm_F
{
public:

    // define dimensions
    const static int NV = 1; // t_i
    const static int NC = 1; // mu

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        return -_c[0] * log(_x[0]);
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        return std::numeric_limits<double>::max();
    }

};

using LB_AD = COMISO::FiniteElementTinyAD<LogBarrierTerm_F>;
using LB_PH_AD = COMISO::FiniteElementHessianProjection< LB_AD >;


//== HALFFACE GLUEING DEFINITION =========================================================


//For 4 3D points, as follows, the energy is the dihedral angle at the axis (a,b).
// The epsilon parameter allows to set a barrier angle at which the energy is infinite.
//    c
//   / \
//  /   \
// a-----b --->
// \    /
//  \ /
//   c'
class HalfFaceGlueingElement_F
{
public:


    const static int NV = 4 * 3; // 4 3D points (a, b, c, c')
    const static int NC = 3; // weight, angle epsilon, barrier weight

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        Eigen::Map<Vec3T> a( (ScalarT*) &(_x[0]) );
        Eigen::Map<Vec3T> b( (ScalarT*) &(_x[3]) );
        Eigen::Map<Vec3T> c( (ScalarT*) &(_x[6]) );
        Eigen::Map<Vec3T> cp( (ScalarT*) &(_x[9]) );

        auto barrier_mu = _c[2];

        //std::cout<<" - barrier w = "<<barrier_mu<<std::endl;

        auto eps = _c[1];
        auto cos_eps = std::cos(eps);
        auto sin_eps = std::sin(eps);

        Vec3T edge_vec = (b - a).normalized();

        Eigen::Matrix<ScalarT, 3,3> rot_mat = Eigen::AngleAxis<ScalarT>(-eps, edge_vec).toRotationMatrix();

        Vec3T nl = edge_vec.cross(c - a).normalized();
        nl = rot_mat * nl;


        //NOTE: opposite normal so we try to get the angle between them to zero
        Vec3T nr = -(cp - a).cross(edge_vec).normalized();

        /*std::cout<<" -      eps = "<<eps<<std::endl;
        std::cout<<" - cos(eps) = "<<cos_eps<<std::endl;
        std::cout<<" - sin(eps) = "<<sin_eps<<std::endl;
        std::cout<<" -        a: "<<a.transpose()<<std::endl;
        std::cout<<" -        b: "<<b.transpose()<<std::endl;
        std::cout<<" -        c: "<<c.transpose()<<std::endl;
        std::cout<<" -       c': "<<cp.transpose()<<std::endl;
        std::cout<<" -         nr: "<<nr.transpose()<<std::endl;
        std::cout<<" -         nl: "<<edge_vec.cross(c - a).normalized().transpose()<<std::endl;
        std::cout<<" - rotated nl: "<<nl.transpose()<<std::endl;
        std::cout<<" -          e: "<<edge_vec.transpose()<<std::endl;*/


        auto nl_cross_nr = nl.cross(nr);
        auto sin_theta = nl_cross_nr.norm();
        nl_cross_nr.normalize();

        auto nl_cross_nr_dot_e = nl_cross_nr.dot(edge_vec);
        auto nl_dot_nr = nl.dot(nr);

        /*std::cout<<" - nl x nr ^: "<<nl_cross_nr.transpose()<<std::endl;
        std::cout<<" - (nl x nr ^) . e: "<<nl_cross_nr_dot_e<<std::endl;
        std::cout<<" - nl . nr = cos(theta): "<<nl_dot_nr<<std::endl;
        std::cout<<" - sin(theta): "<<sin_theta<<std::endl;*/

        if(nl_cross_nr_dot_e == cos_eps || nl_dot_nr == 1 || nl_dot_nr == -1){
            return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
        }

        bool use_angles_version(true);

        ScalarT fx = 1; //to move the minimum to 1
        if(use_angles_version) {
            fx = 0;
        }

        if(use_angles_version) {

            ScalarT angle = acos(nl_dot_nr);
            //std::cout<<" - base angle: "<<angle<<std::endl;

            if (nl_cross_nr_dot_e < 0) { //lower half
                //std::cout<<" - lower half"<<std::endl;
                if (angle < 0) {
                    angle = angle;

                } else {
                    angle = (2 * M_PI - angle);
                    //std::cout << " - lower right sub-epsilon quadrant" << std::endl;
                }

            } else {
                angle = angle;
            }

            //first term for -epsilon and second one for M_PI, where the derivative of acos is not defined
            ScalarT barrier = barrier_mu * sqrt((1.0/ ((1-nl_dot_nr) * (1-nl_dot_nr)))
                    /* + 1.0/ ((1 + nl_dot_nr) * (1 + nl_dot_nr))*/);

            //std::cout<<" - barrier = "<<barrier<<std::endl;

            angle -= eps;
            if(angle < 0){
                angle *= -1;
            }

            //std::cout<<" - final angle = "<<angle<<std::endl;

            fx = angle * angle + barrier;
            //fx = angle + barrier;

        }else{
            if(nl_cross_nr_dot_e < 0){ //lower half
                //std::cout<<" - lower half"<<std::endl;

                if (nl_dot_nr < cos_eps) {
                    fx += (2 + nl_dot_nr);

                } else {
                    fx += -nl_dot_nr;
                    //std::cout << " lower right sub-epsilon quadrant" << std::endl;
                }

                fx += 1e-5 / ((cos_eps - nl_dot_nr) * (cos_eps - nl_dot_nr));
            }else { //upper half
                //std::cout<<" - upper half"<<std::endl;

                fx += -nl_dot_nr;
            }
        }


        //std::cout<<" - base f(x) = "<<fx<<std::endl;


        //std::cout<<" - with barrier = "<<fx<<std::endl;

        return _c[0] * fx;
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        return std::numeric_limits<double>::max();
    }

};



//For 4 3D points, as follows, the energy is the dihedral angle at the axis (a,b).
// The epsilon parameter allows to set a barrier angle at which the energy is infinite.
//    c
//   / \
//  /   \
// a-----b --->
// \    /
//  \ /
//   c'
class HalfFaceGlueingElementV2_F
{
public:


    const static int NV = 4 * 3; // 4 3D points (a, b, c, c')
    const static int NC = 3; // weight, angle epsilon, barrier weight

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        Eigen::Map<Vec3T> a( (ScalarT*) &(_x[0]) );
        Eigen::Map<Vec3T> b( (ScalarT*) &(_x[3]) );
        Eigen::Map<Vec3T> c( (ScalarT*) &(_x[6]) );
        Eigen::Map<Vec3T> cp( (ScalarT*) &(_x[9]) );

        auto barrier_mu = _c[2];

        //std::cout<<" - barrier w = "<<barrier_mu<<std::endl;

        auto eps = _c[1];
        auto cos_eps = std::cos(eps);
        auto sin_eps = std::sin(eps);

        Vec3T edge_vec = (b - a).normalized();

        //Eigen::Matrix<ScalarT, 3,3> rot_mat = Eigen::AngleAxis<ScalarT>(-eps, edge_vec).toRotationMatrix();

        Vec3T nl = edge_vec.cross(c - a).normalized();
        //nl = rot_mat * nl;


        //NOTE: opposite normal so we try to get the angle between them to zero
        Vec3T nr = -(cp - a).cross(edge_vec).normalized();

        /*std::cout<<" -      eps = "<<eps<<std::endl;
        std::cout<<" - cos(eps) = "<<cos_eps<<std::endl;
        std::cout<<" - sin(eps) = "<<sin_eps<<std::endl;
        std::cout<<" -        a: "<<a.transpose()<<std::endl;
        std::cout<<" -        b: "<<b.transpose()<<std::endl;
        std::cout<<" -        c: "<<c.transpose()<<std::endl;
        std::cout<<" -       c': "<<cp.transpose()<<std::endl;
        std::cout<<" -         nr: "<<nr.transpose()<<std::endl;
        std::cout<<" -         nl: "<<edge_vec.cross(c - a).normalized().transpose()<<std::endl;
        std::cout<<" - rotated nl: "<<nl.transpose()<<std::endl;
        std::cout<<" -          e: "<<edge_vec.transpose()<<std::endl;*/


        auto nl_cross_nr = nl.cross(nr);
        //auto sin_theta = nl_cross_nr.norm();
        //nl_cross_nr.normalize();

        auto nl_cross_nr_dot_e = nl_cross_nr.dot(edge_vec);
        auto nl_dot_nr = nl.dot(nr);

        auto& cos_theta = nl_dot_nr;
        auto& sin_theta = nl_cross_nr_dot_e;

        /*std::cout<<" - nl x nr ^: "<<nl_cross_nr.transpose()<<std::endl;
        std::cout<<" - (nl x nr ^) . e: "<<nl_cross_nr_dot_e<<std::endl;
        std::cout<<" - nl . nr = cos(theta): "<<nl_dot_nr<<std::endl;
        std::cout<<" - sin(theta): "<<sin_theta<<std::endl;*/


        ScalarT fx = 0;
        if(nl_cross_nr_dot_e < 0 && nl_dot_nr > 0){
            fx = 3 + nl_cross_nr_dot_e;
        }else {
            fx = 1 - nl_cross_nr_dot_e;
        }

        ScalarT barrier = 1 / (1 - nl_dot_nr);
        //std::cout<<" - base f(x) = "<<fx<<std::endl;


        //std::cout<<" - with barrier = "<<fx<<std::endl;

        return _c[0] * fx + barrier_mu * barrier;
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        // HACK (add exact computation or better line search)
        double t=2.0;
        for(int i=0; i<50; ++i)
        {
            VecV xn = _x+t*_v;
            if(std::isfinite(eval_f(xn, _c)))
                return t;
            else
                t *= 0.7;
        }

        return 0;
    }

};



//For 4 3D points, as follows, the energy is the difference between the dihedral angle at the axis (a,b) and a target angle given by the epsilon constant.
// the left face is (a,b,c) and the right face is (a,c',b).
// The normal for the right face is based on the (a,b,c') orientation since we want left and right faces to have the same normal (so zero angle between them)
//    c
//   / \
//  /   \
// a-----b --->
// \    /
//  \ /
//   c'
class HalfFaceGlueingElementV3_F
{
public:


    const static int NV = 4 * 3; // 4 3D points (a, b, c, c')
    const static int NC = 2; // weight, angle epsilon

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        Eigen::Map<Vec3T> a( (ScalarT*) &(_x[0]) );
        Eigen::Map<Vec3T> b( (ScalarT*) &(_x[3]) );
        Eigen::Map<Vec3T> c( (ScalarT*) &(_x[6]) );
        Eigen::Map<Vec3T> cp( (ScalarT*) &(_x[9]) );

        auto barrier_mu = 0;

        //std::cout<<" - barrier w = "<<barrier_mu<<std::endl;

        auto eps = _c[1];
        auto cos_eps = std::cos(eps);
        auto sin_eps = std::sin(eps);

        Vec3T edge_vec = (b - a).normalized();

        Vec3T nl = edge_vec.cross(c - a).normalized();
        //NOTE: opposite normal so we try to get the angle between them to zero
        Vec3T nr = edge_vec.cross(cp - a).normalized();

        /*std::cout<<" -      eps = "<<eps<<std::endl;
        std::cout<<" - cos(eps) = "<<cos_eps<<std::endl;
        std::cout<<" - sin(eps) = "<<sin_eps<<std::endl;
        std::cout<<" -        a: "<<a.transpose()<<std::endl;
        std::cout<<" -        b: "<<b.transpose()<<std::endl;
        std::cout<<" -        c: "<<c.transpose()<<std::endl;
        std::cout<<" -       c': "<<cp.transpose()<<std::endl;
        std::cout<<" -       nr: "<<nr.transpose()<<std::endl;
        std::cout<<" -       nl: "<<nl.transpose()<<std::endl;
        std::cout<<" -        e: "<<edge_vec.transpose()<<std::endl;*/


        auto nl_cross_nr = nl.cross(nr);
        auto nl_cross_nr_norm = nl_cross_nr.norm();
        //nl_cross_nr.normalize();


        ScalarT nl_cross_nr_dot_e = nl_cross_nr.dot(edge_vec);
        ScalarT nl_dot_nr = nl.dot(nr);

        auto& cos_theta = nl_dot_nr;
        auto& sin_theta = nl_cross_nr_dot_e;

        double dbl_cos_theta = TinyAD::to_passive(cos_theta);
        double dbl_sin_theta = TinyAD::to_passive(sin_theta);

        /*std::cout<<" - cos theta = "<<dbl_cos_theta<<std::endl;
        std::cout<<" - sin theta = "<<dbl_sin_theta<<std::endl;
        std::cout<<" - theta from cos = "<<std::acos(dbl_cos_theta)<<std::endl;
        std::cout<<" - theta from sin = "<<std::asin(dbl_sin_theta)<<std::endl;
        std::cout<<" ---"<<std::endl;
        std::cout<<" - cos diff = "<<(dbl_cos_theta - cos_eps)<<std::endl;
        std::cout<<" - sin diff = "<<(dbl_sin_theta - sin_eps)<<std::endl;
        std::cout<<" - cos diff squared = "<<((dbl_cos_theta - cos_eps) * (dbl_cos_theta - cos_eps))<<std::endl;
        std::cout<<" - sin diff squared = "<<((dbl_sin_theta - sin_eps) * (dbl_sin_theta - sin_eps))<<std::endl;
        std::cout<<" -----------------------"<<std::endl;

        //(n0.n1-cos(alpha))^2 + ((n0 x n1).e.normalized()-sin(alpha))^2

        std::cout<<" - nl x nr ^: "<<nl_cross_nr.transpose()<<std::endl;
        std::cout<<" - (nl x nr ^) . e: "<<nl_cross_nr_dot_e<<std::endl;
        std::cout<<" - nl . nr = cos(theta): "<<nl_dot_nr<<std::endl;
        std::cout<<" - sin(theta): "<<sin_theta<<std::endl;*/


        ScalarT fx = (cos_theta - cos_eps) * (cos_theta - cos_eps) + (sin_theta - sin_eps) * (sin_theta - sin_eps);
        //(n0.n1-cos(alpha))^2 + ((n0 x n1).e.normalized()-sin(alpha))^2

        return _c[0] * fx;
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        // HACK (add exact computation or better line search)
        double t=2.0;
        for(int i=0; i<50; ++i)
        {
            VecV xn = _x+t*_v;
            if(std::isfinite(eval_f(xn, _c)))
                return t;
            else
                t *= 0.7;
        }

        return 0;
    }

};



using HF_GLUEING_AD = COMISO::FiniteElementTinyAD<HalfFaceGlueingElementV3_F>;
//using HF_GLUEING_PH_AD = COMISO::FiniteElementHessianProjection< HF_GLUEING_AD >;



//For 4 3D points, as shown below, the energy is the difference between the dihedral angle theta at the axis (a,b) and a target angle theta_t
// the left face is (a,b,c) and the right face is (a,c',b).
// The normal for the right face is based on the (a,b,c') orientation since we want left and right faces to have the same normal (so zero angle between them)
// (a,b,c) are constants, and the target angle theta_t is
//     0 if |theta| < pi/3
//     (3/pi) * (|theta| - pi/3)^2 if pi/3 <= |theta| <= pi/2
//     |theta| - pi/2 + pi/12 if pi / 2 <= theta
//So basically going from constant between -pi/3 and pi/3, then quadratically transitioning to a linear function to +-pi/2
//
// there's also barrier terms maintaining theta_t between theta < pi/2 and theta
//
//    c
//   / \
//  /   \
// a-----b --->
// \    /
//  \ /
//   c'
class HalfFaceGlueingConstantTargetFaceElement_F
{
public:


    const static int NV = 3; // one 3D points c'
    const static int NC = 1 + 1 + 3 * 3; // weight, theta_t (rad), three 3D points a,b and c

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        using Vec3C    = Eigen::Matrix<double,3,1, Eigen::ColMajor>;
        Vec3C a = {_c[2], _c[3], _c[4]};
        Vec3C b = {_c[5], _c[6], _c[7]};
        Vec3C c = {_c[8], _c[9], _c[10]};

        Eigen::Map<Vec3T> cp( (ScalarT*) &(_x[0]) );

        auto theta_t = _c[1];
        auto cos_t = cos(theta_t);
        auto sin_t = sin(theta_t);

        //std::cout<<" - barrier w = "<<barrier_mu<<std::endl;


        Vec3T edge_vec = (b - a).normalized();

        Vec3T nl = edge_vec.cross(c - a).normalized();
        //NOTE: opposite normal so we try to get the angle between them to zero
        Vec3T nr = edge_vec.cross(cp - a).normalized();

        /*std::cout<<" -  theta_t = "<<theta_t<<std::endl;
        std::cout<<" - cos(theta_t) = "<<cos_t<<std::endl;
        std::cout<<" - sin(theta_t) = "<<sin_t<<std::endl;
        std::cout<<" -             a: "<<a.transpose()<<std::endl;
        std::cout<<" -             b: "<<b.transpose()<<std::endl;
        std::cout<<" -             c: "<<c.transpose()<<std::endl;
        std::cout<<" -            c': "<<cp.transpose()<<std::endl;
        std::cout<<" -            nr: "<<nr.transpose()<<std::endl;
        std::cout<<" -            nl: "<<nl.transpose()<<std::endl;
        std::cout<<" -             e: "<<edge_vec.transpose()<<std::endl;*/


        auto nl_cross_nr = nl.cross(nr);
        auto nl_cross_nr_norm = nl_cross_nr.norm();
        //nl_cross_nr.normalize();


        ScalarT nl_cross_nr_dot_e = nl_cross_nr.dot(edge_vec);
        ScalarT nl_dot_nr = nl.dot(nr);

        auto& cos_theta = nl_dot_nr;
        auto& sin_theta = nl_cross_nr_dot_e;

        double dbl_cos_theta = TinyAD::to_passive(cos_theta);
        double dbl_sin_theta = TinyAD::to_passive(sin_theta);

        /*std::cout<<" - cos theta = "<<dbl_cos_theta<<std::endl;
        std::cout<<" - sin theta = "<<dbl_sin_theta<<std::endl;
        std::cout<<" - theta from cos = "<<std::acos(dbl_cos_theta)<<std::endl;
        std::cout<<" - theta from sin = "<<std::asin(dbl_sin_theta)<<std::endl;
        std::cout<<" ---"<<std::endl;
        std::cout<<" - cos diff = "<<(dbl_cos_theta - cos_t)<<std::endl;
        std::cout<<" - sin diff = "<<(dbl_sin_theta - sin_t)<<std::endl;
        std::cout<<" - cos diff squared = "<<((dbl_cos_theta - cos_t) * (dbl_cos_theta - cos_t))<<std::endl;
        std::cout<<" - sin diff squared = "<<((dbl_sin_theta - sin_t) * (dbl_sin_theta - sin_t))<<std::endl;
        std::cout<<" -----------------------"<<std::endl;

        //(n0.n1-cos(alpha))^2 + ((n0 x n1).e.normalized()-sin(alpha))^2

        std::cout<<" - nl x nr ^: "<<nl_cross_nr.transpose()<<std::endl;
        std::cout<<" - (nl x nr ^) . e: "<<nl_cross_nr_dot_e<<std::endl;
        std::cout<<" - nl . nr = cos(theta): "<<nl_dot_nr<<std::endl;
        std::cout<<" - sin(theta): "<<sin_theta<<std::endl;*/



        ScalarT fx = (cos_theta - cos_t) * (cos_theta - cos_t) + (sin_theta - sin_t) * (sin_theta - sin_t);
        //(n0.n1-cos(alpha))^2 + ((n0 x n1).e.normalized()-sin(alpha))^2

        return _c[0] * fx;
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        // HACK (add exact computation or better line search)
        double t=2.0;
        for(int i=0; i<50; ++i)
        {
            VecV xn = _x+t*_v;
            if(std::isfinite(eval_f(xn, _c)))
                return t;
            else
                t *= 0.7;
        }

        return 0;
    }

};



using HF_GLUEING_CONSTANT_TARGET_AD = COMISO::FiniteElementTinyAD<HalfFaceGlueingConstantTargetFaceElement_F>;
//using HF_GLUEING_PH_AD = COMISO::FiniteElementHessianProjection< HF_GLUEING_AD >;


//For 4 3D points, as follows, the energy is the dihedral angle at the axis (a,b).
// The epsilon parameter allows to set a barrier angle at which the energy is infinite.
//    c
//   / \
//  /   \
// a-----b --->
// \    /
//  \ /
//   c'
class HalfFaceGlueingBarrierElement_F
{
public:


    const static int NV = 4 * 3; // 4 3D points (a, b, c, c')
    const static int NC = 3; // weight, angle epsilon, barrier weight

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        Eigen::Map<Vec3T> a( (ScalarT*) &(_x[0]) );
        Eigen::Map<Vec3T> b( (ScalarT*) &(_x[3]) );
        Eigen::Map<Vec3T> c( (ScalarT*) &(_x[6]) );
        Eigen::Map<Vec3T> cp( (ScalarT*) &(_x[9]) );

        auto barrier_mu = _c[2];

        //std::cout<<" - barrier w = "<<barrier_mu<<std::endl;

        auto eps = _c[1];
        auto cos_eps = std::cos(eps);
        auto sin_eps = std::sin(eps);

        Vec3T edge_vec = (b - a).normalized();

        //Eigen::Matrix<ScalarT, 3,3> rot_mat = Eigen::AngleAxis<ScalarT>(-eps, edge_vec).toRotationMatrix();

        Vec3T nl = edge_vec.cross(c - a).normalized();
        //nl = rot_mat * nl;


        //NOTE: opposite normal so we try to get the angle between them to zero
        Vec3T nr = -(cp - a).cross(edge_vec).normalized();

        /*std::cout<<" -      eps = "<<eps<<std::endl;
        std::cout<<" - cos(eps) = "<<cos_eps<<std::endl;
        std::cout<<" - sin(eps) = "<<sin_eps<<std::endl;
        std::cout<<" -        a: "<<a.transpose()<<std::endl;
        std::cout<<" -        b: "<<b.transpose()<<std::endl;
        std::cout<<" -        c: "<<c.transpose()<<std::endl;
        std::cout<<" -       c': "<<cp.transpose()<<std::endl;
        std::cout<<" -         nr: "<<nr.transpose()<<std::endl;
        std::cout<<" -         nl: "<<edge_vec.cross(c - a).normalized().transpose()<<std::endl;
        std::cout<<" - rotated nl: "<<nl.transpose()<<std::endl;
        std::cout<<" -          e: "<<edge_vec.transpose()<<std::endl;*/


        auto nl_cross_nr = nl.cross(nr);
        //nl_cross_nr.normalize();

        auto nl_cross_nr_dot_e = nl_cross_nr.dot(edge_vec);
        auto nl_dot_nr = nl.dot(nr);

        auto& cos_theta = nl_dot_nr;
        auto& sin_theta = nl_cross_nr_dot_e;

        /*std::cout<<" - nl x nr ^: "<<nl_cross_nr.transpose()<<std::endl;
        std::cout<<" - (nl x nr ^) . e: "<<nl_cross_nr_dot_e<<std::endl;
        std::cout<<" - nl . nr = cos(theta): "<<nl_dot_nr<<std::endl;
        std::cout<<" - sin(theta): "<<sin_theta<<std::endl;*/

        //lower  quadrant: out of domain

        ScalarT fx = 1 - nl_dot_nr;
        ScalarT barrier = 0;

        //upper quadrants: barrier at pi
        if(sin_theta > 0){
            //std::cout<<" upper quadrants, barrier at π"<<std::endl;
            barrier = 1.0 / (1 + cos_theta);
        }else{
            //lower right quadrant: barrier at -pi/2
            if(cos_theta > 0){
                //std::cout<<" lower right quadrant, barrier at -π/2"<<std::endl;
                barrier = 1.0 / (1 + sin_theta);
            }else{ // lower left quadrant: infinity
                //std::cout<<" lower left quadrant, infinity"<<std::endl;
                return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
            }
        }

        //ScalarT barrier = 1 / (1 - cos_theta); + 1 / (1 + sin_theta);
        //std::cout<<" - base f(x) = "<<fx<<std::endl;


        //std::cout<<" - with barrier = "<<barrier<<std::endl;

        return _c[0] * fx + barrier_mu * barrier;
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        // HACK (add exact computation or better line search)
        double t=2.0;
        for(int i=0; i<50; ++i)
        {
            VecV xn = _x+t*_v;
            if(std::isfinite(eval_f(xn, _c)))
                return t;
            else
                t *= 0.7;
        }

        return 0;
    }

};

using HF_GLUEING_BARRIER_AD = COMISO::FiniteElementTinyAD<HalfFaceGlueingBarrierElement_F>;




//For 4 3D points, as follows, the energy is the squared norm of the difference between the normals
// of the faces (a, b, c) and (a, c', b), with a and c being constant
//    c
//   / \
//  /   \
// a-----b --->
// \    /
//  \ /
//   c'
class HalfFaceGlueingElementFixedHinge_F
{
public:

    const static int NV = 2 * 3; // 2 3D points c and c'
    const static int NC = 1 + 2 * 3; // a weight and two fixed 3D points a and b,

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec3T = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        using Vec3C = Eigen::Matrix<double, 3,1, Eigen::ColMajor>;

        Eigen::Map<Vec3C> a( (double*) &(_c[1]) );
        Eigen::Map<Vec3C> b( (double*) &(_c[4]) );
        Eigen::Map<Vec3T> c( (ScalarT*) &(_x[0]) );
        Eigen::Map<Vec3T> cp( (ScalarT*) &(_x[3]) );

        /*auto eps = _c[1];
        auto cos_eps = std::cos(eps);
        auto sin_eps = std::sin(eps);*/

        Vec3T edge_vec = (b - a).normalized();

        Vec3T nl = edge_vec.cross(c - a).normalized();
        //NOTE: opposite normal so we try to get the angle between them to zero
        Vec3T nr = edge_vec.cross(cp - a).normalized();

        ScalarT fx = (nl-nr).squaredNorm();


        /*auto nl_cross_nr = nl.cross(nr);
        auto nl_cross_nr_norm = nl_cross_nr.norm();
        //nl_cross_nr.normalize();


        ScalarT nl_cross_nr_dot_e = nl_cross_nr.dot(edge_vec);
        ScalarT nl_dot_nr = nl.dot(nr);

        auto& cos_theta = nl_dot_nr;
        auto& sin_theta = nl_cross_nr_dot_e;

        double dbl_cos_theta = TinyAD::to_passive(cos_theta);
        double dbl_sin_theta = TinyAD::to_passive(sin_theta);*/

        /*std::cout<<" - cos theta = "<<dbl_cos_theta<<std::endl;
        std::cout<<" - sin theta = "<<dbl_sin_theta<<std::endl;
        std::cout<<" - theta from cos = "<<std::acos(dbl_cos_theta)<<std::endl;
        std::cout<<" - theta from sin = "<<std::asin(dbl_sin_theta)<<std::endl;
        std::cout<<" ---"<<std::endl;
        std::cout<<" - cos diff = "<<(dbl_cos_theta - cos_eps)<<std::endl;
        std::cout<<" - sin diff = "<<(dbl_sin_theta - sin_eps)<<std::endl;
        std::cout<<" - cos diff squared = "<<((dbl_cos_theta - cos_eps) * (dbl_cos_theta - cos_eps))<<std::endl;
        std::cout<<" - sin diff squared = "<<((dbl_sin_theta - sin_eps) * (dbl_sin_theta - sin_eps))<<std::endl;
        std::cout<<" -----------------------"<<std::endl;

        //(n0.n1-cos(alpha))^2 + ((n0 x n1).e.normalized()-sin(alpha))^2

        std::cout<<" - nl x nr ^: "<<nl_cross_nr.transpose()<<std::endl;
        std::cout<<" - (nl x nr ^) . e: "<<nl_cross_nr_dot_e<<std::endl;
        std::cout<<" - nl . nr = cos(theta): "<<nl_dot_nr<<std::endl;
        std::cout<<" - sin(theta): "<<sin_theta<<std::endl;*/


        //ScalarT fx = (cos_theta - cos_eps) * (cos_theta - cos_eps) + (sin_theta - sin_eps) * (sin_theta - sin_eps);
        //(n0.n1-cos(alpha))^2 + ((n0 x n1).e.normalized()-sin(alpha))^2

        return _c[0] * fx;
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        // HACK (add exact computation or better line search)
        double t=2.0;
        for(int i=0; i<50; ++i)
        {
            VecV xn = _x+t*_v;
            if(std::isfinite(eval_f(xn, _c)))
                return t;
            else
                t *= 0.7;
        }

        return 0;
    }

};



using HF_GLUEING_FIXED_HINGE_AD = COMISO::FiniteElementTinyAD<HalfFaceGlueingElementFixedHinge_F>;
using HF_GLUEING_FIXED_HINGE_PH_AD = COMISO::FiniteElementHessianProjection<HF_GLUEING_FIXED_HINGE_AD>;




//== DISTANCE DEFINITION =========================================================


class DistanceBarrierTerm_F
{
public:

    // define dimensions
    const static int NV = 6; // 2 3d points
    const static int NC = 1; // mu

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        Eigen::Map<Vec3T> A( (ScalarT*) &(_x[0]) );
        Eigen::Map<Vec3T> B( (ScalarT*) &(_x[3]) );

        auto diff = A-B;


        /*std::cout<<" ------- "<<std::endl;
        std::cout<<" A: "<<A.transpose()<<std::endl;
        std::cout<<" B: "<<B.transpose()<<std::endl;
        std::cout<<" diff: "<<diff.transpose()<<std::endl;
        std::cout<<" e = "<<(_c[0] * diff.dot(diff))<<std::endl;*/
        //squared distance
        return _c[0] * diff.dot(diff);
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        // HACK (add exact computation or better line search)
        double t=2.0;
        for(int i=0; i<50; ++i)
        {
            VecV xn = _x+t*_v;
            if(std::isfinite(eval_f(xn, _c)))
                return t;
            else
                t *= 0.7;
        }

        return 0;
    }

};

using DB_AD = COMISO::FiniteElementTinyAD<DistanceBarrierTerm_F>;
using DB_PH_AD = COMISO::FiniteElementHessianProjection< DB_AD >;


//== Max DISTANCE LOG BARRIER DEFINITION =========================================================

class MaxDistanceBarrierTerm_F
{
public:

    // define dimensions
    const static int NV = 6; // 2 3d points
    const static int NC = 2; // mu, r

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        Eigen::Map<Vec3T> A( (ScalarT*) &(_x[0]) );
        Eigen::Map<Vec3T> B( (ScalarT*) &(_x[3]) );

        auto diff = A-B;
        auto r = (A-B).norm();
        auto r_diff = _c[1] - r;

        /*std::cout<<" -    diff: "<<diff.transpose()<<std::endl;
        std::cout<<" -       r: "<<r<<std::endl;
        std::cout<<" -  r_diff: "<<r_diff<<std::endl;*/

        if(r_diff <= 0){
            return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
        }else{

            //return -log(10*log(3*r_diff) + 7);
            return -_c[0] * log(1.*r_diff);
        }
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        // HACK (add exact computation or better line search)
        double t=2.0;
        for(int i=0; i<50; ++i)
        {
            VecV xn = _x+t*_v;
            if(std::isfinite(eval_f(xn, _c)))
                return t;
            else
                t *= 0.7;
        }

        return 0;
    }

};

using MaxDLB_AD = COMISO::FiniteElementTinyAD<MaxDistanceBarrierTerm_F>;
//using MinDLB_PH_AD = COMISO::FiniteElementHessianProjection< MinDLB_AD >;





//== GLUEING DEFINITION =========================================================


class EdgeGlueingElement3D_F
{
public:

    // define dimensions
    const static int NV = 12; // [] = two 3D points A_left, B_left, A_right, B_right, of two edges to glue
    const static int NC = 1; // [weight]

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;

        // get constants
        const double& w         = _c[0];

        Eigen::Map<Vec3T> Al( (ScalarT*) &(_x[0]) );
        Eigen::Map<Vec3T> Bl( (ScalarT*) &(_x[3]) );
        Eigen::Map<Vec3T> Ar( (ScalarT*) &(_x[6]) );
        Eigen::Map<Vec3T> Br( (ScalarT*) &(_x[9]) );

        /*std::cout<<" Al: "<<Al.transpose()<<std::endl;
        std::cout<<" Bl: "<<Bl.transpose()<<std::endl;
        std::cout<<" Ar: "<<Ar.transpose()<<std::endl;
        std::cout<<" Br: "<<Br.transpose()<<std::endl;*/

        Vec3T left = Bl - Al;
        Vec3T right = Br - Ar;

        //std::cout<<" - left edge: "<<left.transpose()<<std::endl;
        //std::cout<<" - right edge: "<<right.transpose()<<std::endl;

        Vec3T diff = left - right;

        //std::cout<<" - diff = "<<diff.transpose()<<std::endl;
        //std::cout<<" - e = "<<(left - right).squaredNorm()<<std::endl;


        //return w * (left - right).template lpNorm<1>();
        return w * (left - right).squaredNorm();

    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        return std::numeric_limits<double>::max();
    }

};

using EGLUE3D_AD = COMISO::FiniteElementTinyAD<EdgeGlueingElement3D_F>;
using EGLUE3D_PH_AD = COMISO::FiniteElementHessianProjection< EGLUE3D_AD >;




class GlueingElement3D_F
{
public:

    // define dimensions
    const static int NV = 18; // [] = 3D points A_left, B_left, C_left,  A_right, B_right, C_right of two halffaces to glue
    const static int NC = 2; // [mu, weight]

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;

        // get constants
        const double& mu         = _c[0];
        const double& w          = _c[1];

        Eigen::Map<Vec3T> Al( (ScalarT*) &(_x[0]) );
        Eigen::Map<Vec3T> Bl( (ScalarT*) &(_x[3]) );
        Eigen::Map<Vec3T> Cl( (ScalarT*) &(_x[6]) );
        Eigen::Map<Vec3T> Ar( (ScalarT*) &(_x[9]) );
        Eigen::Map<Vec3T> Br( (ScalarT*) &(_x[12]) );
        Eigen::Map<Vec3T> Cr( (ScalarT*) &(_x[15]) );

        return mu * w * ((Al - Ar).dot(Al - Ar) + (Bl - Br).dot(Bl - Br) + (Cl - Cr).dot(Cl - Cr));

    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        return std::numeric_limits<double>::max();
    }

};

using GLUE3D_AD = COMISO::FiniteElementTinyAD<GlueingElement3D_F>;
using GLUE3D_PH_AD = COMISO::FiniteElementHessianProjection< GLUE3D_AD >;




//== SYMMETRIC DIRICHLET DEFINITION =========================================================

class SymmetricDirichletElement3D_F
{
public:

    static std::string element_name(){ return "Symmetric Dirichlet"; }
    static std::string element_code(){ return "SD"; }

    // define dimensions
    const static int NV = 12; // [ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz] = 3D points A, B, C, D of a tetrahedron
    const static int NC = 12; // [weight, size_scale, epsilon, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z]

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    static void scale_weights(double scale_factor, VecC& vc) {
        vc[0] *= scale_factor;
    }

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        using Vec3dT   = Eigen::Matrix<double ,3,1,Eigen::ColMajor>;
        using Mat3x3T  = Eigen::Matrix<ScalarT,3,3,Eigen::ColMajor>;
        using Mat3x4T  = Eigen::Matrix<ScalarT,3,4,Eigen::ColMajor>;
        using Mat3x3dT = Eigen::Matrix<double ,3,3,Eigen::ColMajor>;

        // get constants
        const double& w          = _c[0];
        const double& size_scale = _c[1];
        const double& epsilon    = _c[2];

        if (epsilon > 0.0) {
            std::cerr<<" ERROR - epsilon SDE is not working anymore"<<std::endl;
            exit(EXIT_FAILURE);
        }

        // get gradient computation matrix
        Eigen::Map<Mat3x3dT> G( (double*) &(_c[3]));

        // get points of tetrahedron
        Eigen::Map<Mat3x4T> P((ScalarT*) &(_x[0]));

        // get edge vectors
        Mat3x3T E;
        E.col(0) = P.col(1)-P.col(0);
        E.col(1) = P.col(2)-P.col(0);
        E.col(2) = P.col(3)-P.col(0);
//    ScalarT V = 1.0/6.0*E.determinant();
        /*if(G.determinant() <= 0.0){
            G = Mat3x3dT::Identity();
        }*/

        // calc (scaled) Jacobi Matrix of map
        Mat3x3T J = (size_scale*E)*G.transpose();

        // calculate determinant
        ScalarT d = J.determinant();

        if(false) {
            auto e = w*(J.squaredNorm() + J.inverse().squaredNorm() - 6.0);;
            std::cout << " - det E = " << E.determinant() << std::endl;
            std::cout << " - det G = " << G.determinant() << std::endl;
            std::cout << " - G = " << G << std::endl;
            std::cout << " - d = " << d << std::endl;
            std::cout << " - eps = " << epsilon << std::endl;
            std::cout << " - J = " << J << std::endl;
            std::cout<<" - SDE = "<<e<<std::endl;
            std::cout << " ----------" << std::endl;
        }

        if(epsilon == 0.0)
        {
            if (d <= 0.0)
            {
                return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
            }
            else
            {
                auto e = w*(J.squaredNorm() + J.inverse().squaredNorm() - 6.0);;
                //std::cout<<" - e = "<<e<<std::endl;
                return e;
                //return w*(J.squaredNorm() + J.inverse().squaredNorm() - 6.0);
                //        return w*(1.0+d)*(J.squaredNorm() + J.inverse().squaredNorm() - 6.0);
                //        return w*V*(J.squaredNorm() + J.inverse().squaredNorm() - 6.0);
                //        return w*V*(1.0+d)*(J.squaredNorm() + J.inverse().squaredNorm() - 6.0);
            }
        }
        else
        {
            // regularized determinant
            ScalarT d_reg = d+epsilon;

            if (d_reg <= 0.0)
                return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
            else
            {
                ScalarT di = 1.0/d_reg;

                // compute regularized inverse
                Mat3x3T Ji;
                for(int i=0; i<3; ++i)
                    for(int j=0; j<3; ++j)
                    {
                        int i1 = (i+1)%3;
                        int j1 = (j+1)%3;
                        int i2 = (i+2)%3;
                        int j2 = (j+2)%3;

                        Ji(j,i) = di*(J(i1,j1)*J(i2,j2)-J(i1,j2)*J(i2,j1));
                    }

                return w*(J.squaredNorm() + Ji.squaredNorm() - 6.0);
            }
        }
    }

    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        return dirichlet_max_feasible_step(_x, _v, _c);

    }

    //double t_factor_for

    inline double volume_max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){


        return std::numeric_limits<double>::max();

        using Mat3x3T = Eigen::Matrix<double, 3,3,Eigen::ColMajor>;
        using Mat3x4T = Eigen::Matrix<double, 3,4,Eigen::ColMajor>;

        // get constants
        const double& size_scale = _c[1];
        const double& epsilon    = _c[2];

        // get points of tetrahedron
        Eigen::Map<Mat3x4T> P((double*) &(_x[0]));

        // get edge vectors
        Mat3x3T E;
        E.col(0) = P.col(1)-P.col(0);
        E.col(1) = P.col(2)-P.col(0);
        E.col(2) = P.col(3)-P.col(0);

        // get points of tetrahedron
        Eigen::Map<Mat3x4T> dP((double*) &(_v[0]));

        // get edge vectors
        Mat3x3T dE;
        dE.col(0) = dP.col(1)-dP.col(0);
        dE.col(1) = dP.col(2)-dP.col(0);
        dE.col(2) = dP.col(3)-dP.col(0);

    }




    inline double dirichlet_max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
    {
        using Mat3x3T = Eigen::Matrix<double, 3,3,Eigen::ColMajor>;
        using Mat3x4T = Eigen::Matrix<double, 3,4,Eigen::ColMajor>;

        // get constants
        //const double& size_scale = _c[1];
        //const double& epsilon    = _c[2];
        const double& size_scale = 1;
        const double& epsilon    = 0;

        // get gradient computation matrix
        Eigen::Map<Mat3x3T> G( (double*) &(_c[3]));

        // get points of tetrahedron
        Eigen::Map<Mat3x4T> P((double*) &(_x[0]));

        // get edge vectors
        Mat3x3T E;
        E.col(0) = P.col(1)-P.col(0);
        E.col(1) = P.col(2)-P.col(0);
        E.col(2) = P.col(3)-P.col(0);

        // get points of tetrahedron
        Eigen::Map<Mat3x4T> dP((double*) &(_v[0]));

        // get edge vectors
        Mat3x3T dE;
        dE.col(0) = dP.col(1)-dP.col(0);
        dE.col(1) = dP.col(2)-dP.col(0);
        dE.col(2) = dP.col(3)-dP.col(0);

        // coefficients computed via Mathematica (+index shift since Mathematica is 1 based instead of 0 based)
        double c0 =-(E(0,2)*E(1,1)*E(2,0)) + E(0,1)*E(1,2)*E(2,0) + E(0,2)*E(1,0)*E(2,1) -
                   E(0,0)*E(1,2)*E(2,1) - E(0,1)*E(1,0)*E(2,2) + E(0,0)*E(1,1)*E(2,2);

        double c1 = -(dE(2,2)*E(0,1)*E(1,0)) + dE(2,1)*E(0,2)*E(1,0) + dE(2,2)*E(0,0)*E(1,1) -
                    dE(2,0)*E(0,2)*E(1,1) - dE(2,1)*E(0,0)*E(1,2) + dE(2,0)*E(0,1)*E(1,2) +
                    dE(1,2)*E(0,1)*E(2,0) - dE(1,1)*E(0,2)*E(2,0) - dE(0,2)*E(1,1)*E(2,0) +
                    dE(0,1)*E(1,2)*E(2,0) - dE(1,2)*E(0,0)*E(2,1) + dE(1,0)*E(0,2)*E(2,1) +
                    dE(0,2)*E(1,0)*E(2,1) - dE(0,0)*E(1,2)*E(2,1) + dE(1,1)*E(0,0)*E(2,2) -
                    dE(1,0)*E(0,1)*E(2,2) - dE(0,1)*E(1,0)*E(2,2) + dE(0,0)*E(1,1)*E(2,2);

        double c2 = -(dE(1,2)*dE(2,1)*E(0,0)) + dE(1,1)*dE(2,2)*E(0,0) + dE(1,2)*dE(2,0)*E(0,1) -
                    dE(1,0)*dE(2,2)*E(0,1) - dE(1,1)*dE(2,0)*E(0,2) + dE(1,0)*dE(2,1)*E(0,2) +
                    dE(0,2)*dE(2,1)*E(1,0) - dE(0,1)*dE(2,2)*E(1,0) - dE(0,2)*dE(2,0)*E(1,1) +
                    dE(0,0)*dE(2,2)*E(1,1) + dE(0,1)*dE(2,0)*E(1,2) - dE(0,0)*dE(2,1)*E(1,2) -
                    dE(0,2)*dE(1,1)*E(2,0) + dE(0,1)*dE(1,2)*E(2,0) + dE(0,2)*dE(1,0)*E(2,1) -
                    dE(0,0)*dE(1,2)*E(2,1) - dE(0,1)*dE(1,0)*E(2,2) + dE(0,0)*dE(1,1)*E(2,2);

        double c3 = -(dE(0,2)*dE(1,1)*dE(2,0)) + dE(0,1)*dE(1,2)*dE(2,0) + dE(0,2)*dE(1,0)*dE(2,1) -
                    dE(0,0)*dE(1,2)*dE(2,1) - dE(0,1)*dE(1,0)*dE(2,2) + dE(0,0)*dE(1,1)*dE(2,2);

        // adjust coefficients
        double s = std::pow(size_scale,3)*G.determinant();
        c0 *= s;
        c1 *= s;
        c2 *= s;
        c3 *= s;
        // shift w.r.t. epsilon
        c0 += epsilon;

        COMISO::Monomial<3> f;
        f.coeffs() << c0, c1, c2, c3;
        double t;
        if(COMISO::Polynomials::first_root_in_interval(f, 0.0, 10.0, t))
        {
//      // check result
//      std::cout << "determinant at estimated root = " << (size_scale*(E+t*dE)*G.transpose()).determinant() + epsilon << std::endl;

            // reduce a little bit to avoid infeasibility of max_step
            return 0.99*t;
        }
        else
            return DBL_MAX;
        return DBL_MAX;
    }





    // NC = 12; // [weight, size_scale, epsilon, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z]
    // calculate constants from reference tetrahedron
//epsilon = 0 if all tets are initially valid
    static int compute_constants(const std::vector<double> parameters,
                                 const std::vector<Eigen::Vector3d>& _p,
                                 VecC& _c)
    {
        _c[0] = parameters[0];
        _c[1] = parameters[1];
        _c[2] = parameters[2];

        // compute gradients of pwl basis functions g1, g2, g3

        // get edge vectors
        Eigen::Matrix3d E;
        /*E.col(0) = _p1-_p0;
        E.col(1) = _p2-_p0;
        E.col(2) = _p3-_p0;*/
        E.col(0) = _p[1]-_p[0];
        E.col(1) = _p[2]-_p[0];
        E.col(2) = _p[3]-_p[0];

        double d  = E.determinant();
        double di = 1.0/d;
        if(!std::isfinite(di) || d < 0.0)
        {
            std::cout << "ERROR: determinant of edge vectors numerically not valid ---> " << d << std::endl;
            /*std::cout<<" p0 = "<<_p0.transpose()<<std::endl;
            std::cout<<" p1 = "<<_p1.transpose()<<std::endl;
            std::cout<<" p2 = "<<_p2.transpose()<<std::endl;
            std::cout<<" p3 = "<<_p3.transpose()<<std::endl;
            std::cout<<" edge vectors matrix: "<<std::endl<<E<<std::endl;
            std::cout<<" d = "<<d<<std::endl;
            std::cout<<" di = "<<di<<std::endl;
            std::cout<<" -----------------------------"<<std::endl;*/
            return 1;

            di = 0.0;
        }

        _c.segment<3>(3) = di*(E.col(1)).cross(E.col(2));
        _c.segment<3>(6) = di*(E.col(2)).cross(E.col(0));
        _c.segment<3>(9) = di*(E.col(0)).cross(E.col(1));

        // -------------------
        // test and debug code
        // -------------------
        if(0)
        {
            Eigen::Matrix3d J;
            J =     _c.segment<3>(3)*E.col(0).transpose() +
                    _c.segment<3>(6)*E.col(1).transpose() +
                    _c.segment<3>(9)*E.col(2).transpose();

            J -= Eigen::Matrix3d::Identity();

            double Jn = J.norm();
            if(Jn > 1e-6)
            {
                std::cout << "ERROR: Jacobi matrix was not computed correctly!" << std::endl;

                J =     _c.segment<3>(3)*E.col(0).transpose() +
                        _c.segment<3>(6)*E.col(1).transpose() +
                        _c.segment<3>(9)*E.col(2).transpose();
                std::cout << J << std::endl;
            }
            else
                std::cout << "Jacobi matrix within tolerance ---> " << Jn << std::endl;
            std::cout<<" debug J = "<<std::endl<<J<<std::endl;
        }

        return 0;
    }


    static double determine_suitable_epsilon(const Eigen::Matrix<double, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c, const double _target_penalty)
    {
        using Mat3x3dT = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;
        using Mat3x4dT = Eigen::Matrix<double, 3, 4, Eigen::ColMajor>;

        // get constants
        const double &w = _c[0];
        const double &size_scale = _c[1];
        const double &epsilon = _c[2];

        // get gradient computation matrix
        Eigen::Map<Mat3x3dT> G((double *) &(_c[3]));

        // get points of tetrahedron
        Eigen::Map<Mat3x4dT> P((double *) &(_x[0]));

        // get edge vectors
        Mat3x3dT E;
        E.col(0) = P.col(1) - P.col(0);
        E.col(1) = P.col(2) - P.col(0);
        E.col(2) = P.col(3) - P.col(0);

        // calc (scaled) Jacobi Matrix of map
        Mat3x3dT J = size_scale * E * G.transpose();

        // calculate determinant
        double d = J.determinant();

        if (d > 0.0)
            return 0.0;
        else
            return 1.0 / _target_penalty - d;
    }

};

using SD_3D_AD = COMISO::FiniteElementTinyAD<SymmetricDirichletElement3D_F>;
//using SD_3D_PH_AD = COMISO::FiniteElementHessianProjection< SD_3D_AD >;



//== AUGMENTED SYMMETRIC DIRICHLET DEFINITION =========================================================

class AugmentedSymmetricDirichletElement3D_F : public SymmetricDirichletElement3D_F
{
public:

    static std::string element_name(){ return "Augmented Symmetric Dirichlet"; }
    static std::string element_code(){ return "AUG_SD"; }

    // define dimensions
    const static int NV = 12; // [ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz] = 3D points A, B, C, D of a tetrahedron
    const static int NC = 12; // [weight, size_scale, epsilon, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z]

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    static void scale_weights(double scale_factor, VecC& vc) {
        vc[0] *= scale_factor;
    }

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        using Vec3dT   = Eigen::Matrix<double ,3,1,Eigen::ColMajor>;
        using Mat3x3T  = Eigen::Matrix<ScalarT,3,3,Eigen::ColMajor>;
        using Mat3x4T  = Eigen::Matrix<ScalarT,3,4,Eigen::ColMajor>;
        using Mat3x3dT = Eigen::Matrix<double ,3,3,Eigen::ColMajor>;

        // get constants
        const double& w          = _c[0];
        const double& size_scale = _c[1];
        const double& sqrt_reg_epsilon    = _c[2];

        // get gradient computation matrix
        Eigen::Map<Mat3x3dT> G( (double*) &(_c[3]));

        // get points of tetrahedron
        Eigen::Map<Mat3x4T> P((ScalarT*) &(_x[0]));

        // get edge vectors
        Mat3x3T E;
        E.col(0) = P.col(1)-P.col(0);
        E.col(1) = P.col(2)-P.col(0);
        E.col(2) = P.col(3)-P.col(0);
//    ScalarT V = 1.0/6.0*E.determinant();
        /*if(G.determinant() <= 0.0){
            G = Mat3x3dT::Identity();
        }*/

        // calc (scaled) Jacobi Matrix of map
        Mat3x3T J = (size_scale*E)*G.transpose();

        // calculate determinant
        ScalarT d = J.determinant();

        /*std::cout<<" - det E = "<<E.determinant()<<std::endl;
        std::cout<<" - det G = "<<G.determinant()<<std::endl;
        std::cout<<" - G = "<<G<<std::endl;
        std::cout<<" - d = "<<d<<std::endl;
        std::cout<<" - eps = "<<epsilon<<std::endl;
        std::cout<<" - J = "<<J<<std::endl;
        std::cout<<" ----------"<<std::endl;*/

        if (d <= 0.0)
        {
            return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
        }
        else
        {
            auto SDE = J.squaredNorm() + J.inverse().squaredNorm() - 6.0;
            //std::cout<<" - e = "<<e<<std::endl;
            return w * (SDE + SDE*SDE + sqrt(SDE + sqrt_reg_epsilon) - std::sqrt(sqrt_reg_epsilon));
        }
    }


    static int compute_constants(const std::vector<double> parameters,
                                 const std::vector<Eigen::Vector3d>& _p,
                                 VecC& _c){
        SymmetricDirichletElement3D_F::compute_constants(parameters, _p, _c);

        //replacing the SDE epsilon
        _c[2] = parameters[2];

        return 0;
    }

};

using AUG_SD_3D_AD = COMISO::FiniteElementTinyAD<AugmentedSymmetricDirichletElement3D_F>;
//using SD_3D_PH_AD = COMISO::FiniteElementHessianProjection< SD_3D_AD >;


//== SYMMETRIC DIRICHLET WITH BOUND DEFINITION =========================================================

class SymmetricDirichletWithBoundElement3D_F
{
public:

    // define dimensions
    const static int NV = 12; // [ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz] = 3D points A, B, C, D of a tetrahedron
    const static int NC = 12; // [weight, barrier_w, barrier_c, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z]

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    static std::string element_name(){ return "Symmetric Dirichlet with Bound"; }
    static std::string element_code(){ return "SD_B"; }

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        using Vec3dT   = Eigen::Matrix<double ,3,1,Eigen::ColMajor>;
        using Mat3x3T  = Eigen::Matrix<ScalarT,3,3,Eigen::ColMajor>;
        using Mat3x4T  = Eigen::Matrix<ScalarT,3,4,Eigen::ColMajor>;
        using Mat3x3dT = Eigen::Matrix<double ,3,3,Eigen::ColMajor>;

        // get constants
        const double& w         = _c[0];
        const double& barrier_w = _c[1];
        const double& barrier_min = _c[2];

        // get gradient computation matrix
        Eigen::Map<Mat3x3dT> G( (double*) &(_c[3]));

        // get points of tetrahedron
        Eigen::Map<Mat3x4T> P((ScalarT*) &(_x[0]));

        // get edge vectors
        Mat3x3T E;
        E.col(0) = P.col(1)-P.col(0);
        E.col(1) = P.col(2)-P.col(0);
        E.col(2) = P.col(3)-P.col(0);

        //*std::cout<<" - E: "<<std::endl<<E<<std::endl;

//    ScalarT V = 1.0/6.0*E.determinant();
        /*if(G.determinant() <= 0.0){
            G = Mat3x3dT::Identity();
        }*/

        // calc (scaled) Jacobi Matrix of map
        Mat3x3T J = E*G.transpose();

        // calculate determinant
        ScalarT d = J.determinant();

        if(barrier_w < 0) {
            std::cout<<" - E: "<<std::endl<<E<<std::endl;
            std::cout<<" - G: "<<std::endl<<G<<std::endl;
            std::cout<<" - J: "<<std::endl<<J<<std::endl;
            std::cout<<" - det = "<<d<<std::endl;
            std::cout<<" - w = "<<w<<std::endl;
            std::cout<<" - barrier w = "<<barrier_w<<std::endl;
            std::cout<<" - barrier min = "<<barrier_min<<std::endl;
            auto E_sd = J.squaredNorm() + J.inverse().squaredNorm() - 6.0;
            std::cout<<" - SDE = "<<E_sd<<std::endl;

            ScalarT E_sd_bound = w * E_sd + barrier_w / (barrier_min - E_sd);
            std::cout<<" SDE with bound = "<<E_sd_bound<<std::endl;
        }

        if (d <= 0.0)
        {
            return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
        }
        else
        {
            auto E_sd = J.squaredNorm() + J.inverse().squaredNorm() - 6.0;

            if(barrier_w > 0){
                if(barrier_min - E_sd <= 0.0) {
                    return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
                }else{
                    return w * (E_sd
                                + barrier_w / (barrier_min - E_sd)
                                - barrier_w / barrier_min);//last term to have actual zero energy
                }
            }else{
                return w * E_sd;
            }

        }

    }

    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){

        if(_c[1] <= 0){
            //using the standard feasible step if the barrier term is disabled
            SymmetricDirichletElement3D_F elem;

            VecC c = _c;
            c[1] = 1;
            c[2] = 0;

            return elem.dirichlet_max_feasible_step(_x, _v, _c);

        }else{
            // HACK (add exact computation or better line search)
            double t=2.0;
            for(int i=0; i<50; ++i)
            {
                VecV xn = _x+t*_v;
                if(std::isfinite(eval_f(xn, _c)))
                    return t;
                else
                    t *= 0.7;
            }
        }

        return 0;

    }


    // calculate constants from reference tetrahedron
    // const double _w, const double _barrier_w, const double _barrier_c,
    static int compute_constants(const std::vector<double> parameters,
                                 const std::vector<Eigen::Vector3d>& _p,
            //const Eigen::Vector3d& _p0, const Eigen::Vector3d& _p1, const Eigen::Vector3d& _p2, const Eigen::Vector3d& _p3,
                                 VecC& _c)
    {
        if(parameters.size() != 3){
            std::cout<<" ERROR - Symmetric Dirichlet Energy element requires 3 parameters."<<std::endl;
            exit(EXIT_FAILURE);
        }

        _c[0] = parameters[0];
        _c[1] = parameters[1];
        _c[2] = parameters[2];

        // compute gradients of pwl basis functions g1, g2, g3

        // get edge vectors
        Eigen::Matrix3d E;
        /*E.col(0) = _p1-_p0;
        E.col(1) = _p2-_p0;
        E.col(2) = _p3-_p0;*/
        E.col(0) = _p[1]-_p[0];
        E.col(1) = _p[2]-_p[0];
        E.col(2) = _p[3]-_p[0];

        double d  = E.determinant();
        double di = 1.0/d;
        if(!std::isfinite(di) || d < 0.0)
        {
            std::cout << "ERROR: determinant of edge vectors numerically not valid ---> " << d << std::endl;

            return 1;

            di = 0.0;
        }

        _c.segment<3>(3) = di*(E.col(1)).cross(E.col(2));
        _c.segment<3>(6) = di*(E.col(2)).cross(E.col(0));
        _c.segment<3>(9) = di*(E.col(0)).cross(E.col(1));

        return 0;
    }
};

using SD_B_3D_AD = COMISO::FiniteElementTinyAD<SymmetricDirichletWithBoundElement3D_F>;
//using SD_B_3D_PH_AD = COMISO::FiniteElementHessianProjection< SD_B_3D_AD >;




// ============================================================================================
// ============================= SYMMETRIC DIRICHLET WITH BOUND AND VARIABLE REFERENCE
// ============================================================================================



class SymmetricDirichletWithBarrierAndVariableReferenceElement3D {
public:

    // define dimensions
    const static int NV = 24; // 4 3D points of a tet and its reference shape
    //const static int NC = 12; // [weight, barrier_w, barrier_min, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z]
    const static int NC = 3; // [weight, barrier_w, barrier_min]

    using VecV = Eigen::Matrix<double, NV, 1>;
    using VecC = Eigen::Matrix<double, NC, 1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const {
        using Vec3T = Eigen::Matrix<ScalarT, 3, 1, Eigen::ColMajor>;
        using Mat3x3T = Eigen::Matrix<ScalarT, 3, 3, Eigen::ColMajor>;
        using Mat3x4T = Eigen::Matrix<ScalarT, 3, 4, Eigen::ColMajor>;
        using Mat4x3T = Eigen::Matrix<ScalarT, 4, 3, Eigen::RowMajor>;
        using Vec3dT = Eigen::Matrix<double, 3, 1, Eigen::ColMajor>;
        using Mat3x3dT = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;
        using Mat3x4dT = Eigen::Matrix<double, 3, 4, Eigen::RowMajor>;


        // get constants
        const double &w = _c[0];
        const double &barrier_w = _c[1];
        const double &barrier_min = _c[2];


        // get points of triangle
        //Eigen::Map <Mat2x3T> P((ScalarT *) &(_x[0]));
        Mat4x3T P, P_ref;
        P <<_x[0], _x[1],_x[2],
                _x[3],_x[4],_x[5],
                _x[6],_x[7],_x[8],
                _x[9],_x[10],_x[11];

        P_ref <<_x[12], _x[13],_x[14],
                _x[15],_x[16],_x[17],
                _x[18],_x[19],_x[20],
                _x[21],_x[22],_x[23];


        // get edge vectors
        Mat3x3T E, E_ref;
        E.col(0) = P.row(1) - P.row(0);
        E.col(1) = P.row(2) - P.row(0);
        E.col(2) = P.row(3) - P.row(0);

        E_ref.col(0) = P_ref.row(1) - P_ref.row(0);
        E_ref.col(1) = P_ref.row(2) - P_ref.row(0);
        E_ref.col(2) = P_ref.row(3) - P_ref.row(0);


        // calc  Jacobi Matrix of map
        Mat3x3T J = E * E_ref.inverse();

        //regular ref element
        //Mat2x2T J = jacobian(_x);


        // calculate determinant
        ScalarT d = J.determinant();


        if(_c[1] < 0) {
            std::cout<<" P: "<<std::endl<<P<<std::endl;
            std::cout<<" P_ref: "<<std::endl<<P_ref<<std::endl;
            std::cout<<" E: "<<std::endl<<E<<std::endl;
            std::cout<<" E_ref: "<<std::endl<<E_ref<<std::endl;

            std::cout<<" J: "<<std::endl<<J<<std::endl;
            std::cout<<" d = "<<d<<std::endl;

            ScalarT E_sd = (J.squaredNorm() + J.inverse().squaredNorm() - 6.0);
            std::cout<<" SDE = "<<E_sd<<std::endl;

            ScalarT E_sd_barrier = w * E_sd + barrier_w / (barrier_min - E_sd);
            std::cout<<" SDE with barrier = "<<E_sd_barrier<<std::endl;

        }


        if (d <= 0.0) {
            return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
        } else {
            //auto E_sd = (1.0+d)*(J.squaredNorm() + J.inverse().squaredNorm() - 4.0);
            ScalarT E_sd = (J.squaredNorm() + J.inverse().squaredNorm() - 6.0);
            //std::cout<<" E_sd = "<<E_sd<<std::endl;

            //return w * E_sd;

            if (barrier_w && barrier_min - E_sd <= 0.0) {
                //std::cout<<" inf barrier term"<<std::endl;
                return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
            }
            //return w * E_sd;

            return w * E_sd + barrier_w / (barrier_min - E_sd);
        }
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c) {

        // HACK (add exact computation or better line search)
        double t = 2.0;
        for (int i = 0; i < 50; ++i) {
            VecV xn = _x + t * _v;
            //std::cout<<" ------------------------- t"<<i<<" = "<<t<<std::endl;
            if (std::isfinite(eval_f(xn, _c)))
                //if(eval_f(xn, _c) < _c[2]) //alternative, simply not letting the energy reach the maximum
                return t;
            else
                t *= 0.7;
        }

        return 0;
    }


    static int compute_constants(const double _w, const double _barrier_w, const double _barrier_c,
                                 VecC &_c){

        _c[0] = _w;
        _c[1] = _barrier_w;
        _c[2] = _barrier_c;
        return 0;

    }
};

using SDEB_VR_3D_AD = COMISO::FiniteElementTinyAD<SymmetricDirichletWithBarrierAndVariableReferenceElement3D>;
//using SDEB_VR_2D_PH_AD = COMISO::FiniteElementHessianProjection< SDEB_VR_2D_AD >;



//== FOLDOVER-FREE WITHOUT EPSILON DEFINITION =========================================================
class FoldoverFreeWithoutEpsilonElement3D
{
public:

    // define dimensions
    const static int NV = 12; // [ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz] = 3D points A, B, C, D of a tetrahedron
    const static int NC = 11; // [conformal w, volumetric w, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z]

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    static std::string element_name(){ return "Foldover-Free without epsilon"; }
    static std::string element_code(){ return "FOF_NO_EPS"; }

    static void scale_weights(double scale_factor, VecC& vc) {
        vc[0] *= scale_factor;
        vc[1] *= scale_factor;
    }

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        //using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        //using Vec3dT   = Eigen::Matrix<double ,3,1,Eigen::ColMajor>;
        using Mat3x3T  = Eigen::Matrix<ScalarT,3,3,Eigen::ColMajor>;
        using Mat3x4T  = Eigen::Matrix<ScalarT,3,4,Eigen::ColMajor>;
        using Mat3x3dT = Eigen::Matrix<double ,3,3,Eigen::ColMajor>;

        // get constants
        const double& conf_w         = _c[0];
        const double& vol_w          = _c[1];

        // get gradient computation matrix
        Eigen::Map<Mat3x3dT> G( (double*) &(_c[2]));

        // get points of tetrahedron
        Eigen::Map<Mat3x4T> P((ScalarT*) &(_x[0]));

        // get edge vectors
        Mat3x3T E;
        E.col(0) = P.col(1)-P.col(0);
        E.col(1) = P.col(2)-P.col(0);
        E.col(2) = P.col(3)-P.col(0);


        // calc (scaled) Jacobi Matrix of map
        Mat3x3T J = E*G.transpose();

        // calculate determinant
        ScalarT d = J.determinant();

        if(conf_w < 0) {
            std::cout<<" - E: "<<std::endl<<E<<std::endl;
            std::cout<<" - G: "<<std::endl<<G<<std::endl;
            std::cout<<" - J: "<<std::endl<<J<<std::endl;
            std::cout<<" - det = "<<d<<std::endl;
            std::cout<<" - conformal w = "<<conf_w<<std::endl;
            std::cout<<" -    volume w = "<<vol_w<<std::endl;
            auto E_sd = J.squaredNorm() + J.inverse().squaredNorm() - 6.0;
            std::cout<<" - SDE = "<<E_sd<<std::endl;
            auto E_conf = (J.transpose() * J).trace() / pow(d, 2.0 / 3.0) - 3.0;
            std::cout<<" - E_conf = "<<E_conf<<std::endl;
        }


        if (d <= 0.0)
            return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
        else
            return conf_w * ((J.transpose() * J).trace() / pow(d, 2.0 / 3.0) - 3.0) + vol_w * (d + 1.0 / d - 2.0);
        //      return w_conformal*((J.transpose()*J).trace()/pow(d,2.0/3.0)) + w_det*(d/target_det + target_det/d);

    }

    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c) {
        using Mat3x3T = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;
        using Mat3x4T = Eigen::Matrix<double, 3, 4, Eigen::ColMajor>;

        // get constants
        const double &size_scale = 1.0; //_c[1];
        const double &epsilon = 0.0; //_c[2];

        // get gradient computation matrix
        Eigen::Map<Mat3x3T> G((double *) &(_c[2]));

        // get points of tetrahedron
        Eigen::Map<Mat3x4T> P((double *) &(_x[0]));

        // get edge vectors
        Mat3x3T E;
        E.col(0) = P.col(1) - P.col(0);
        E.col(1) = P.col(2) - P.col(0);
        E.col(2) = P.col(3) - P.col(0);

        // get points of tetrahedron
        Eigen::Map<Mat3x4T> dP((double *) &(_v[0]));

        // get edge vectors
        Mat3x3T dE;
        dE.col(0) = dP.col(1) - dP.col(0);
        dE.col(1) = dP.col(2) - dP.col(0);
        dE.col(2) = dP.col(3) - dP.col(0);

        // coefficients computed via Mathematica (+index shift since Mathematica is 1 based instead of 0 based)
        double c0 = -(E(0, 2) * E(1, 1) * E(2, 0)) + E(0, 1) * E(1, 2) * E(2, 0) + E(0, 2) * E(1, 0) * E(2, 1) -
                    E(0, 0) * E(1, 2) * E(2, 1) - E(0, 1) * E(1, 0) * E(2, 2) + E(0, 0) * E(1, 1) * E(2, 2);

        double c1 = -(dE(2, 2) * E(0, 1) * E(1, 0)) + dE(2, 1) * E(0, 2) * E(1, 0) + dE(2, 2) * E(0, 0) * E(1, 1) -
                    dE(2, 0) * E(0, 2) * E(1, 1) - dE(2, 1) * E(0, 0) * E(1, 2) + dE(2, 0) * E(0, 1) * E(1, 2) +
                    dE(1, 2) * E(0, 1) * E(2, 0) - dE(1, 1) * E(0, 2) * E(2, 0) - dE(0, 2) * E(1, 1) * E(2, 0) +
                    dE(0, 1) * E(1, 2) * E(2, 0) - dE(1, 2) * E(0, 0) * E(2, 1) + dE(1, 0) * E(0, 2) * E(2, 1) +
                    dE(0, 2) * E(1, 0) * E(2, 1) - dE(0, 0) * E(1, 2) * E(2, 1) + dE(1, 1) * E(0, 0) * E(2, 2) -
                    dE(1, 0) * E(0, 1) * E(2, 2) - dE(0, 1) * E(1, 0) * E(2, 2) + dE(0, 0) * E(1, 1) * E(2, 2);

        double c2 = -(dE(1, 2) * dE(2, 1) * E(0, 0)) + dE(1, 1) * dE(2, 2) * E(0, 0) + dE(1, 2) * dE(2, 0) * E(0, 1) -
                    dE(1, 0) * dE(2, 2) * E(0, 1) - dE(1, 1) * dE(2, 0) * E(0, 2) + dE(1, 0) * dE(2, 1) * E(0, 2) +
                    dE(0, 2) * dE(2, 1) * E(1, 0) - dE(0, 1) * dE(2, 2) * E(1, 0) - dE(0, 2) * dE(2, 0) * E(1, 1) +
                    dE(0, 0) * dE(2, 2) * E(1, 1) + dE(0, 1) * dE(2, 0) * E(1, 2) - dE(0, 0) * dE(2, 1) * E(1, 2) -
                    dE(0, 2) * dE(1, 1) * E(2, 0) + dE(0, 1) * dE(1, 2) * E(2, 0) + dE(0, 2) * dE(1, 0) * E(2, 1) -
                    dE(0, 0) * dE(1, 2) * E(2, 1) - dE(0, 1) * dE(1, 0) * E(2, 2) + dE(0, 0) * dE(1, 1) * E(2, 2);

        double c3 = -(dE(0, 2) * dE(1, 1) * dE(2, 0)) + dE(0, 1) * dE(1, 2) * dE(2, 0) + dE(0, 2) * dE(1, 0) * dE(2, 1)
                    -
                    dE(0, 0) * dE(1, 2) * dE(2, 1) - dE(0, 1) * dE(1, 0) * dE(2, 2) + dE(0, 0) * dE(1, 1) * dE(2, 2);

        // adjust coefficients
        double s = std::pow(size_scale, 3) * G.determinant();
        c0 *= s;
        c1 *= s;
        c2 *= s;
        c3 *= s;
        // shift w.r.t. epsilon
        c0 += epsilon;

        COMISO::Monomial<3> f;
        f.coeffs() << c0, c1, c2, c3;
        double t;
        if (COMISO::Polynomials::first_root_in_interval(f, 0.0, 10.0, t)) {
            //      // check result
            //      std::cerr << "determinant at estimated root = " << (size_scale*(E+t*dE)*G.transpose()).determinant() + epsilon << std::endl;

            // reduce a little bit to avoid infeasibility of max_step
            return 0.99 * t;
        } else
            return std::numeric_limits<double>::max();

        return std::numeric_limits<double>::max();
    }


    // calculate constants from reference tetrahedron
    // parameters: (_epsilon, _weight_conformal, _weight_volume, _weight_barrier, _barrier_thresh)
    static int compute_constants(const std::vector<double> parameters,
                                 const std::vector<Eigen::Vector3d>& _p,
            //const Eigen::Vector3d& _p0, const Eigen::Vector3d& _p1, const Eigen::Vector3d& _p2, const Eigen::Vector3d& _p3,
                                 VecC& _c)
    {
        if(parameters.size() != 2){
            std::cout<<" ERROR - Foldover-free element requires 2 parameters."<<std::endl;
            exit(EXIT_FAILURE);
        }

        _c[0] = parameters[0];
        _c[1] = parameters[1];

        // get edge vectors
        Eigen::Matrix3d E;
        E.col(0) = _p[1]-_p[0];
        E.col(1) = _p[2]-_p[0];
        E.col(2) = _p[3]-_p[0];

        double d  = E.determinant();
        double di = 1.0/d;
        if(!std::isfinite(di) || d < 0.0)
        {
            std::cout << "ERROR: determinant of edge vectors numerically not valid ---> " << d << std::endl;

            return 1;

            di = 0.0;
        }

        _c.segment<3>(2) = di*(E.col(1)).cross(E.col(2));
        _c.segment<3>(5) = di*(E.col(2)).cross(E.col(0));
        _c.segment<3>(8) = di*(E.col(0)).cross(E.col(1));

        return 0;
    }
};

using FOF_NO_EPS_3D_AD = COMISO::FiniteElementTinyAD<FoldoverFreeWithoutEpsilonElement3D>;
//using FOF_NO_EPS_3D_PH_AD = COMISO::FiniteElementHessianProjection< FOF_NO_EPS_3D_AD >;



//== AUGMENTED FOLDOVER-FREE WITHOUT EPSILON DEFINITION =========================================================
class AugmentedFoldoverFreeWithoutEpsilonElement3D : public FoldoverFreeWithoutEpsilonElement3D
{
public:

    static constexpr double eps_norm_regularizer = 1e-4;

    // define dimensions
    const static int NV = 12; // [ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz] = 3D points A, B, C, D of a tetrahedron
    const static int NC = 11; // [conformal w, volumetric w, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z]

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    static std::string element_name(){ return "Augmented Foldover-Free without epsilon"; }
    static std::string element_code(){ return "AUG_FOF_NO_EPS"; }

    static void scale_weights(double scale_factor, VecC& vc) {
        vc[0] *= scale_factor;
        vc[1] *= scale_factor;
    }

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        //using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        //using Vec3dT   = Eigen::Matrix<double ,3,1,Eigen::ColMajor>;
        using Mat3x3T  = Eigen::Matrix<ScalarT,3,3,Eigen::ColMajor>;
        using Mat3x4T  = Eigen::Matrix<ScalarT,3,4,Eigen::ColMajor>;
        using Mat3x3dT = Eigen::Matrix<double ,3,3,Eigen::ColMajor>;

        // get constants
        const double& conf_w         = _c[0];
        const double& vol_w          = _c[1];

        // get gradient computation matrix
        Eigen::Map<Mat3x3dT> G( (double*) &(_c[2]));

        // get points of tetrahedron
        Eigen::Map<Mat3x4T> P((ScalarT*) &(_x[0]));

        // get edge vectors
        Mat3x3T E;
        E.col(0) = P.col(1)-P.col(0);
        E.col(1) = P.col(2)-P.col(0);
        E.col(2) = P.col(3)-P.col(0);

        //*std::cout<<" - E: "<<std::endl<<E<<std::endl;

//    ScalarT V = 1.0/6.0*E.determinant();
        /*if(G.determinant() <= 0.0){
            G = Mat3x3dT::Identity();
        }*/

        // calc (scaled) Jacobi Matrix of map
        Mat3x3T J = E*G.transpose();

        // calculate determinant
        ScalarT d = J.determinant();

        if(conf_w < 0) {
            std::cout<<" - E: "<<std::endl<<E<<std::endl;
            std::cout<<" - G: "<<std::endl<<G<<std::endl;
            std::cout<<" - J: "<<std::endl<<J<<std::endl;
            std::cout<<" - det = "<<d<<std::endl;
            std::cout<<" - conformal w = "<<conf_w<<std::endl;
            std::cout<<" -    volume w = "<<vol_w<<std::endl;
            auto E_sd = J.squaredNorm() + J.inverse().squaredNorm() - 6.0;
            std::cout<<" - SDE = "<<E_sd<<std::endl;
        }


        if (d <= 0.0) {
            return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
        }else {

            ScalarT E_conf = ((J.transpose() * J).trace() / pow(d, 2.0 / 3.0) - 3.0);
            ScalarT E_vol = (d + 1.0 / d - 2.0);

            return conf_w * (E_conf + E_conf*E_conf + sqrt(E_conf + eps_norm_regularizer) - sqrt(eps_norm_regularizer)) +
                    vol_w * (E_vol /* + E_vol*E_vol   + sqrt(E_vol  + eps_norm_regularizer) - sqrt(eps_norm_regularizer)*/);
        }
        //      return w_conformal*((J.transpose()*J).trace()/pow(d,2.0/3.0)) + w_det*(d/target_det + target_det/d);

    }

};

using AUG_FOF_NO_EPS_3D_AD = COMISO::FiniteElementTinyAD<AugmentedFoldoverFreeWithoutEpsilonElement3D>;


//== AUGMENTED NORMAL ANGLE DEFINITION =========================================================
class AugmentedDihedralAngleElement3D
{
public:

    static constexpr double eps_norm_regularizer = 1e-5;
    static constexpr double acos_eps = 1e-12;

    // define dimensions
    const static int NV = 12; // [ax,ay,az,bx,by,bz,cx,cy,cz,vx,vy,vz] = 3D points A, B, C, D of a tetrahedron
    const static int NC = 6; // [w, theta, sqrt_w, pow, bound_w, bound] =
    // weight, target dihedral angle at edge a,b, weight for the sqrt part, power to use on the higher-order term, bound weight, bound value

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    static std::string element_name(){ return "Augmented Dihedral Angle"; }
    static std::string element_code(){ return "AUG_DA"; }

    static void scale_weights(double scale_factor, VecC& vc) {
        vc[0] *= scale_factor;
    }

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        using Vec3dT   = Eigen::Matrix<double ,3,1,Eigen::ColMajor>;

        // get constants
        const double& w       = _c[0];
        const double& theta   = _c[1];
        const double& sqrt_w  = _c[2];
        const double& k       = _c[3];
        const double& bound_w = _c[4];
        const double& bound   = _c[5];


        // get points of tetrahedron
        Eigen::Map<Vec3T> a((ScalarT*) &(_x[0]));
        Eigen::Map<Vec3T> b((ScalarT*) &(_x[3]));
        Eigen::Map<Vec3T> c((ScalarT*) &(_x[6]));
        Eigen::Map<Vec3T> v((ScalarT*) &(_x[9]));

        Vec3T n_b = (b-a).cross(c-a).normalized();
        Vec3T n_v = (v-a).cross(b-a).normalized();

        ScalarT dot_prod = n_b.dot(-n_v);


        if(dot_prod <= -1.0){
            dot_prod = -1.0 + acos_eps;
        }
        if(dot_prod >= 1.0){
            dot_prod = 1.0 - acos_eps;
        }

        ScalarT angle_diff_deg = 180.0 * (acos(dot_prod) - theta) / M_PI;
        ScalarT E_angle = angle_diff_deg * angle_diff_deg;

        /*ScalarT pow_term = 1.0;
        for (int i(0); i<k; i++) {
            pow_term *= E_angle;
        }*/

        ScalarT E_aug = E_angle +
                        pow(E_angle, k) +
                        sqrt_w * (sqrt(E_angle + eps_norm_regularizer) - std::sqrt(eps_norm_regularizer));

        /*std::cout<<" -               weight: "<<w<<std::endl;
        std::cout<<" -                theta: "<<(180 * theta / M_PI)<<std::endl;
        std::cout<<" - boundary face normal: "<<n_b.transpose()<<std::endl;
        std::cout<<" -  virtual face normal: "<<n_v.transpose()<<std::endl;
        std::cout<<" -            dot prod = "<<n_b.dot(n_v)<<std::endl;
        std::cout<<" -               angle = "<<(180 * acos(dot_prod) / M_PI)<<std::endl;
        std::cout<<" -      diff to target = "<<angle_diff_deg<<std::endl;
        std::cout<<" -             bound w = "<<bound_w<<std::endl;
        std::cout<<" -          diff bound = "<<bound<<std::endl;
        std::cout<<" -         base energy = "<<E_angle<<std::endl;
        std::cout<<" -    augmented energy = "<<E_aug<<std::endl;*/

        if (bound_w > 0.0) {
            ScalarT E_bound = bound - sqrt(angle_diff_deg * angle_diff_deg + eps_norm_regularizer) + sqrt(eps_norm_regularizer);

            if (E_bound <= 0.0) {
                return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
            } else {
                return w * E_aug - bound_w * log(E_bound / bound);
            }
        }else {
            return w * E_aug;
        }

        //return w * E_aug;
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c) {

        if (_c[4] > 0) {
            double t=2.0;
            for(int i=0; i<50; ++i)
            {
                VecV xn = _x+t*_v;
                if(std::isfinite(eval_f(xn, _c)))
                    return t;
                else
                    t *= 0.7;
            }
            return 0;
        }else {
            return DBL_MAX;
        }
    }

    static int compute_constants(const std::vector<double> parameters,
                                const std::vector<Eigen::Vector3d>& _p,
                                VecC& _c) {

        _c[0] = parameters[0];
        _c[1] = parameters[1];
        _c[2] = parameters[2];
        _c[3] = parameters[3];
        _c[4] = parameters[4];
        _c[5] = parameters[5];
        return 0;
    }


    };

using AUG_ANGLE_3D_AD = COMISO::FiniteElementTinyAD<AugmentedDihedralAngleElement3D>;
using AUG_ANGLE_3D_PH_AD = COMISO::FiniteElementHessianProjection< AUG_ANGLE_3D_AD >;




//== FOLDOVER-FREE WITH BOUND DEFINITION =========================================================

class FoldoverFreeWithBoundElement3D
{
public:

    // define dimensions
    const static int NV = 12; // [ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz] = 3D points A, B, C, D of a tetrahedron
    const static int NC = 14; // [epsilon, conformal w, volumetric w, bound w, bound, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z]

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    static std::string element_name(){ return "Foldover-Free with Bound"; }
    static std::string element_code(){ return "FOF_B"; }

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        using Vec3dT   = Eigen::Matrix<double ,3,1,Eigen::ColMajor>;
        using Mat3x3T  = Eigen::Matrix<ScalarT,3,3,Eigen::ColMajor>;
        using Mat3x4T  = Eigen::Matrix<ScalarT,3,4,Eigen::ColMajor>;
        using Mat3x3dT = Eigen::Matrix<double ,3,3,Eigen::ColMajor>;

        // get constants
        const double& epsilon        = _c[0];
        const double& conf_w         = _c[1];
        const double& vol_w          = _c[2];
        const double& barrier_w      = _c[3];
        const double& barrier_thresh = _c[4];

        // get gradient computation matrix
        Eigen::Map<Mat3x3dT> G( (double*) &(_c[5]));

        // get points of tetrahedron
        Eigen::Map<Mat3x4T> P((ScalarT*) &(_x[0]));

        // get edge vectors
        Mat3x3T E;
        E.col(0) = P.col(1)-P.col(0);
        E.col(1) = P.col(2)-P.col(0);
        E.col(2) = P.col(3)-P.col(0);

        //*std::cout<<" - E: "<<std::endl<<E<<std::endl;

//    ScalarT V = 1.0/6.0*E.determinant();
        /*if(G.determinant() <= 0.0){
            G = Mat3x3dT::Identity();
        }*/

        // calc (scaled) Jacobi Matrix of map
        Mat3x3T J = E*G.transpose();

        // calculate determinant
        ScalarT d = J.determinant();

        if(barrier_w < 0) {
            std::cout<<" - E: "<<std::endl<<E<<std::endl;
            std::cout<<" - G: "<<std::endl<<G<<std::endl;
            std::cout<<" - J: "<<std::endl<<J<<std::endl;
            std::cout<<" - det = "<<d<<std::endl;
            std::cout<<" - conformal w = "<<conf_w<<std::endl;
            std::cout<<" -    volume w = "<<vol_w<<std::endl;
            std::cout<<" -   barrier w = "<<barrier_w<<std::endl;
            std::cout<<" - barrier thresh = "<<barrier_thresh<<std::endl;
            auto E_sd = J.squaredNorm() + J.inverse().squaredNorm() - 6.0;
            std::cout<<" - SDE = "<<E_sd<<std::endl;
        }

        ScalarT E_fof(0);
        if(epsilon > 0.0)
        {
            //ScalarT d_s   = d/target_det;
            ScalarT d_reg = 0.5 * (d + sqrt(pow(d, 2) + epsilon * epsilon));

            E_fof = conf_w * ((J.transpose()*J).trace() / pow(d_reg, 2.0/3.0)) + vol_w * ((pow(d,2) + 1.0) / d_reg);
        }
        else
        {
            if (d <= 0.0)
                return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
            else
                E_fof = conf_w * ((J.transpose()*J).trace()/pow(d, 2.0/3.0) - 3.0) + vol_w * (d + 1.0 / d - 2.0);
            //      return w_conformal*((J.transpose()*J).trace()/pow(d,2.0/3.0)) + w_det*(d/target_det + target_det/d);
        }


        if (barrier_w > 0.0) {
            if (barrier_thresh - E_fof <= 0.0){
                return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
            }

            return E_fof + barrier_w / (barrier_thresh - E_fof);
        }else {
            return E_fof;
        }


    }

    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){

        double t=2.0;
        for(int i=0; i<50; ++i)
        {
            VecV xn = _x+t*_v;
            if(std::isfinite(eval_f(xn, _c)))
                return t;
            else
                t *= 0.7;
        }

        return 0;
    }


    // calculate constants from reference tetrahedron
    // parameters: (_epsilon, _weight_conformal, _weight_volume, _weight_barrier, _barrier_thresh)
    static int compute_constants(const std::vector<double> parameters,
                                 const std::vector<Eigen::Vector3d>& _p,
            //const Eigen::Vector3d& _p0, const Eigen::Vector3d& _p1, const Eigen::Vector3d& _p2, const Eigen::Vector3d& _p3,
                                 VecC& _c)
    {
        if(parameters.size() != 5){
            std::cout<<" ERROR - Foldover-free element requires 5 parameters."<<std::endl;
            exit(EXIT_FAILURE);
        }

        _c[0] = parameters[0];
        _c[1] = parameters[1];
        _c[2] = parameters[2];
        _c[3] = parameters[3];
        _c[4] = parameters[4];

        // compute gradients of pwl basis functions g1, g2, g3

        // get edge vectors
        Eigen::Matrix3d E;
        /*E.col(0) = _p1-_p0;
        E.col(1) = _p2-_p0;
        E.col(2) = _p3-_p0;*/
        E.col(0) = _p[1]-_p[0];
        E.col(1) = _p[2]-_p[0];
        E.col(2) = _p[3]-_p[0];

        double d  = E.determinant();
        double di = 1.0/d;
        if(!std::isfinite(di) || d < 0.0)
        {
            std::cout << "ERROR: determinant of edge vectors numerically not valid ---> " << d << std::endl;

            return 1;

            di = 0.0;
        }

        _c.segment<3>(5) = di*(E.col(1)).cross(E.col(2));
        _c.segment<3>(8) = di*(E.col(2)).cross(E.col(0));
        _c.segment<3>(11) = di*(E.col(0)).cross(E.col(1));

        return 0;
    }
};

using FOF_B_3D_AD = COMISO::FiniteElementTinyAD<FoldoverFreeWithBoundElement3D>;
using FOF_B_3D_PH_AD = COMISO::FiniteElementHessianProjection< FOF_B_3D_AD >;




//== VOLUME LENGTH DEFINITION =========================================================

class InvertedVolumeLengthWithBoundElement3D
{
public:

    // define dimensions
    const static int NV = 12; // [ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz] = 3D points A, B, C, D of a tetrahedron
    const static int NC = 6; // [weight, bound weight, bound, volume bound weight, volume bound, optimal volume]
    //NOTE: the optimal volume can typically be the average volume, and the volume bound weight should be low enough that the IVL is still dominant

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    static std::string element_name(){ return "Inverted Volume-Length with Bound"; }
    static std::string element_code(){ return "IVL_B"; }

    static constexpr double eps_energy_regularizer = 1e-9; //to avoir numerics shenanigans in the energy computation

    static void scale_weights(double scale_factor, VecC& vc) {
        vc[0] *= scale_factor;
        vc[1] *= scale_factor;
        vc[3] *= scale_factor;
    }

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        using Vec3dT   = Eigen::Matrix<double ,3,1,Eigen::ColMajor>;
        using Mat3x3T  = Eigen::Matrix<ScalarT,3,3,Eigen::ColMajor>;
        using Mat3x4T  = Eigen::Matrix<ScalarT,3,4,Eigen::ColMajor>;
        using Mat3x3dT = Eigen::Matrix<double ,3,3,Eigen::ColMajor>;

        // get constants
        const double &w =             _c[0];
        const double &barrier_w =     _c[1];
        const double &barrier_bound = _c[2];
        const double &vol_w =         _c[3];
        const double &vol_bound =     _c[4];
        const double &vol_opti =      _c[5];

        //here you should weight the barrier with 1 / #cells

        // get points of tetrahedron
        Eigen::Map<Mat3x4T> P((ScalarT*) &(_x[0]));

        // get edge vectors
        Mat3x3T E;
        E.col(0) = P.col(1)-P.col(0);
        E.col(1) = P.col(2)-P.col(0);
        E.col(2) = P.col(3)-P.col(0);

        ScalarT d  = E.determinant();
        if (barrier_w < 0) {
            std::cerr<<" vc: "<<_c.transpose()<<std::endl;
            std::cerr<<" w = "<<w<<std::endl;
            std::cerr<<" E = "<<std::endl<<E<<std::endl;
            std::cerr<<" d = "<<d<<std::endl;
        }

        if(d <= 0.0) {
            return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
        }

        ScalarT vol = d / 6.0;

        ScalarT sum_of_squares =   (P.col(1) - P.col(0)).squaredNorm()
                                   + (P.col(2) - P.col(0)).squaredNorm()
                                   + (P.col(3) - P.col(0)).squaredNorm()
                                   + (P.col(2) - P.col(1)).squaredNorm()
                                   + (P.col(3) - P.col(1)).squaredNorm()
                                   + (P.col(3) - P.col(2)).squaredNorm();



        ScalarT E_ivl = 0.0;
        if (true) {
            //ScalarT l_RMS = sqrt(sum_of_squares / 6.0);
            //ScalarT E_ivl =  (l_RMS * l_RMS * l_RMS) / (6.0 * std::sqrt(2.0) * vol) - 1.0;
            E_ivl = pow(sum_of_squares / 6.0, 3.0 / 2.0) / (6.0 * std::sqrt(2.0) * vol) - 1.0 + eps_energy_regularizer;
        }else {
            E_ivl = (sum_of_squares / 6.0) / pow((6.0 * std::sqrt(2.0) * vol), 2.0 / 3.0) - 1.0 + eps_energy_regularizer;
        }

        if (barrier_w < 0) {
            std::cerr<<" sum of squares = "<<sum_of_squares<<std::endl;
            std::cerr<<" l_RMS = "<<sqrt(sum_of_squares / 6.0)<<std::endl;
            std::cerr<<" E_ivl = "<<E_ivl<<std::endl;
        }

        if (E_ivl < 0.0 || !isfinite(E_ivl)) {
            return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
        }

        if (barrier_w > 0.0) {
            ScalarT E_ivl_bound = barrier_bound - E_ivl;

            if (E_ivl_bound <= 0.0) {
                return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
            }else if (vol_w > 0.0 && vol_bound > 0.0) {
                ScalarT E_vol_bound = vol - vol_bound;

                if (E_vol_bound <= 0.0) {
                    return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
                } else {
                    ScalarT vol_reg = (vol_opti/vol + vol / vol_opti - 2.0);
                    ScalarT vol_log_barrier = 0.1 * (log(vol - vol_bound) - log(vol_opti - vol_bound));
                    return w * (E_ivl - barrier_w * log(E_ivl_bound / barrier_bound) + vol_w * (vol_reg - vol_log_barrier)) ;
                    //return 1e-2 * exp(-w * log(E_bound / bound));
                }
            }else {
                    return w * E_ivl - barrier_w * log(E_ivl_bound / barrier_bound);
            }
        }else {
            return w * E_ivl;
        }

        /*if (w > 0 && barrier_w > 0) {
            if (barrier_bound - E_ivl <= 0.0) {
                return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
            } else {
                return w * (E_ivl
                            + barrier_w / (barrier_bound - E_ivl)
                            - barrier_w / barrier_bound); //last term to have actual zero energy when E_ivl is zero
            }
        } else {
            return w * E_ivl;
        }*/
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){

        double t=2.0;
        for(int i=0; i<50; ++i)
        {
            VecV xn = _x+t*_v;
            if(std::isfinite(eval_f(xn, _c)))
                return t;
            else
                t *= 0.7;
        }

        return 0;
    }


    // calculate constants from reference tetrahedron
    static int compute_constants(const std::vector<double> parameters,
                                 const std::vector<Eigen::Vector3d>& _p,
                                 VecC& _c) {
        if(parameters.size() != NC){
            std::cout<<" ERROR - inverted volume-length element requires "<<NC<<" parameter."<<std::endl;
            exit(EXIT_FAILURE);
        }

        _c[0] = parameters[0];
        _c[1] = parameters[1];
        _c[2] = parameters[2];
        _c[3] = parameters[3];
        _c[4] = parameters[4];
        _c[5] = parameters[5];
        return 0;
    }
};

using IVL_B_3D_AD = COMISO::FiniteElementTinyAD<InvertedVolumeLengthWithBoundElement3D>;
//using IVL_3_3D_PH_AD = COMISO::FiniteElementHessianProjection< IVL_B_3D_AD >;



class InvertedVolumeLengthElementWithSlack3D
{
public:

    static constexpr double eps_norm_regularizer = 1e-4;
    // define dimensions
    const static int NV = 13; // [ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,s] = 3D points A, B, C, D of a tetrahedron and a slack variable
    const static int NC = 3; // [weight, barrier bound, slack weight]

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    static std::string element_name(){ return "Inverted Volume-Length with slack variable"; }
    static std::string element_code(){ return "IVL_S"; }

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
        using Vec3dT   = Eigen::Matrix<double ,3,1,Eigen::ColMajor>;
        using Mat3x3T  = Eigen::Matrix<ScalarT,3,3,Eigen::ColMajor>;
        using Mat3x4T  = Eigen::Matrix<ScalarT,3,4,Eigen::ColMajor>;
        using Mat3x3dT = Eigen::Matrix<double ,3,3,Eigen::ColMajor>;

        // get constants
        const double &w =           _c[0];
        const double &barrier_bound =   _c[1];
        const double &slack_w = _c[2];

        if (w <= 0.0) {
            return 0.0;
        }

        // get points of tetrahedron
        Eigen::Map<Mat3x4T> P((ScalarT*) &(_x[0]));
        const ScalarT&           s = _x[12];

        // get edge vectors
        Mat3x3T E;
        E.col(0) = P.col(1)-P.col(0);
        E.col(1) = P.col(2)-P.col(0);
        E.col(2) = P.col(3)-P.col(0);

        ScalarT d  = E.determinant();
        if(d <= 0.0) {
            return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
        }

        ScalarT vol = d / 6.0;

        ScalarT sum_of_squares =   (P.col(1) - P.col(0)).squaredNorm()
                                   + (P.col(2) - P.col(0)).squaredNorm()
                                   + (P.col(3) - P.col(0)).squaredNorm()
                                   + (P.col(2) - P.col(1)).squaredNorm()
                                   + (P.col(3) - P.col(1)).squaredNorm()
                                   + (P.col(3) - P.col(2)).squaredNorm();


        ScalarT l_RMS = sqrt(sum_of_squares / 6.0);

        ScalarT E_ivl =  (l_RMS * l_RMS * l_RMS) / (6.0 * std::sqrt(2.0) * vol) - 1.0;
        ScalarT E_ivl_bound = barrier_bound + s - E_ivl;

        if (E_ivl_bound <= 0.0 || E_ivl <= 0.0 || !isfinite(E_ivl)) {
            return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
        } else {
            return         w * (E_ivl
                                - log(E_ivl_bound / barrier_bound))
                           + slack_w * (pow(s, 2)
                                        + sqrt(pow(s, 2) + eps_norm_regularizer)
                                        - sqrt(eps_norm_regularizer));
        }

        /* const double ctd = std::sqrt(3.0) / 6.0;
        ScalarT tri_distortion_3d = ctd * E.squaredNorm() / (E.col(0).cross(E.col(1))).norm() - 1.0;
        ScalarT tri_distortion_3d_constraint = distortion_bound + s - tri_distortion_3d;

        if (tri_distortion_3d_constraint <= 0.0 || !isfinite(tri_distortion_3d_constraint))
        return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
        else
        return -w * log(tri_distortion_3d_constraint/distortion_bound) + w_slack_penalty*(pow(s,2) + sqrt(pow(s,2)+eps_norm_regularizer));*/
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){

        double t=2.0;
        for(int i=0; i<50; ++i)
        {
            VecV xn = _x+t*_v;
            if(std::isfinite(eval_f(xn, _c)))
                return t;
            else
                t *= 0.7;
        }

        return 0;
    }


    // calculate constants from reference tetrahedron
    static int compute_constants(const std::vector<double> parameters,
                                 const std::vector<Eigen::Vector3d>& _p,
                                 VecC& _c) {
        if(parameters.size() != 3){
            std::cout<<" ERROR - inverted volume-length element requires 3 parameter."<<std::endl;
            exit(EXIT_FAILURE);
        }

        _c[0] = parameters[0];
        _c[1] = parameters[1];
        _c[2] = parameters[2];
        return 0;
    }
};

using IVL_S_3D_AD = COMISO::FiniteElementTinyAD<InvertedVolumeLengthElementWithSlack3D>;
using IVL_S_3D_PH_AD = COMISO::FiniteElementHessianProjection< IVL_S_3D_AD >;
