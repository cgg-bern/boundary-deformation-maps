

#ifndef DEFORMATIONELEMENTS2D_HH
#define DEFORMATIONELEMENTS2D_HH


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


// ===========================================================

/** Code from Maksim Fomin */
class TriMeshTinyADFOFEnergy_F{


public:
    // define dimensions
    const static int NV = 6; // [ax,ay,bx,by,cx,cy] = 2D points A, B, C of a triangle,
                             // assuming A maps to (0,0), B to (1,0) and C to (0,1)
    const static int NC = 2; // [epsilon, lambda]

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;
    using Triplet = Eigen::Triplet<double>;

    TriMeshTinyADFOFEnergy_F(){
        S_inv_ << 2, 1,
                  0, std::sqrt(3.0);
        S_inv_ *= 0.5;
        det_S_ = S_inv_.determinant();

        S_inv_ = S_inv_.inverse().eval();
    }

    /**
     * @tparam ScalarT
     * @param _x 6d vector containing the 2d positions of all three vertices: (ax, ay, bx, by, cx, cy)
     * @param _c 2d vector containing (epsilon, lambda)
     * @return the energy of the triangle defined by the three points
     */
    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x,
                          const Eigen::Matrix<double, NC, 1> &_c) const{



        using VecT  = Eigen::Matrix<ScalarT,2,1,Eigen::ColMajor>;
        using MatT  = Eigen::Matrix<ScalarT,2,2,Eigen::ColMajor>;

        // get constants
        const double& epsilon = _c[0];
        const double& lambda  = _c[1];

        // get input points
        Eigen::Map<VecT> A( (ScalarT*) &(_x[0]) );
        Eigen::Map<VecT> B( (ScalarT*) &(_x[2]) );
        Eigen::Map<VecT> C( (ScalarT*) &(_x[4]) );

        /*std::cout<<" ---------------- "<<std::endl;
        std::cout<<" element with x = "<<A<<", "<<B<<", "<<C<<std::endl;
        std::cout<<" - eps/lambda = "<<epsilon<<"/"<<lambda<<std::endl;*/

        //TODO

        //std::cout<<" epsilon = "<<epsilon<<std::endl;
        //std::cout<<" lambda = "<<lambda<<std::endl;

        JacobianT<ScalarT> J = jacobian(_x);
        //std::cout<<" J = "<<J<<std::endl;
        ScalarT det = J.determinant();
        //std::cout<<" det = "<<det<<std::endl;

        ScalarT chi = det / 2 + sqrt(pow(epsilon,2) + pow(det,2))/2;
        //std::cout<<" chi = "<<chi<<std::endl;

        ScalarT chip = .5 + det/(2 * sqrt(pow(epsilon,2) + pow(det,2)));
        //std::cout<<" chip = "<<chip<<std::endl;

        ScalarT f = (J.transpose() * J).trace() / chi;
        //std::cout<<" f = "<<f<<std::endl;

        ScalarT g;
        if (lambda == 0)
            g = 0;
        else
            g = (pow(det,2) + 1.0)/chi;
        //std::cout<<" g = "<<g<<std::endl;



        return f + lambda * g;


        //temp stuff just to compute something
        /*VecT AB = B-A;
        VecT AC = C-A;

        return AB.norm() + AC.norm();*/

        //return AB.normalized().dot(AC.normalized());
    }


    //template<class ScalarT>
    inline void   eval_modified_hessian (const VecV& _x, const VecC& _c, std::vector<Triplet>& _triplets) const{

        using VecT  = Eigen::Matrix<double,2,1,Eigen::ColMajor>;
        using MatT  = Eigen::Matrix<double,2,2,Eigen::ColMajor>;

        // get constants
        const double& epsilon = _c[0];
        const double& lambda  = _c[1];

        JacobianT<double> J = jacobian(_x);

        double det = J.determinant();

        double chi = det / 2 + sqrt(pow(epsilon,2) + pow(det,2))/2;

        double chip = .5 + det/(2 * sqrt(pow(epsilon,2) + pow(det,2)));

        //TODO: set the modified hessian matrix here
        Eigen::Matrix<double,4,1,Eigen::ColMajor> a, b;
        a.block(0,0,2,1) = J.col(0);
        a.block(2,0,2,1) = J.col(1);
        //Hopefully that's the right way to initialize a.
        //turns out it wasn't. It should have been the columns, not the rows.

        b[0] =  a[3];
        b[1] = -a[2];
        b[2] = -a[1];
        b[3] =  a[0];

        Eigen::Matrix<double,3,2,Eigen::ColMajor> suppZ;

        suppZ << -1, -1,
                  1,  0,
                  0,  1;
        //suppZ.col(0) << -1, 1, 0;
        //suppZ.col(1) << -1, 0, 1;                                               //Supplementary matrix for Z of shape (-1 ... -1)
        Eigen::Matrix<double,3,2,Eigen::ColMajor> Z = suppZ * S_inv_;                                        // (    I    )


        Eigen::Matrix<double,4,5,Eigen::ColMajor> supp1M;    //Supplementary matrix at the beginning and at the end of the multiplication
        supp1M.setIdentity();
        supp1M.col(4) << b;

        double chi2  = chi*chi;
        double chi3  = chi2*chi;
        double chip2 = chip*chip;
        double det2  = det*det;

        auto I4 = Eigen::Matrix4d::Identity();

        auto mp1 = (2. / chi) * I4;
        auto mp23 = -2. * (chip / chi2) * a;
        auto mp4 = 2 * a.squaredNorm() * chip2/chi3 + lambda * (2./chi - 4.*det*chip/chi2 + 2 * (1 + det2) * chip2/chi3);

        Eigen::Matrix<double,5,5,Eigen::ColMajor> supp2M;    //Supplementary matrix in the middle of multiplication
        supp2M.block(0,0,4,4) = (2. / chi) * I4;
        supp2M.block(0,4, 4, 1) = -2. * (chip / chi2) * a;
        supp2M.block(4,0, 1, 4) = supp2M.block(0,4, 4, 1).transpose();
        supp2M(4,4) = 2 * a.squaredNorm() * chip2/chi3 + lambda * (2./chi - 4.*det*chip/chi2 + 2 * (1 + det2) * chip2/chi3);


        auto Mp = supp1M * supp2M * supp1M.transpose();                    //Fill out the M+ matrix!

        double det_S = S_inv_.inverse().determinant();
        MatT I3 = MatT::Identity();
        for (int j = 0; j <= 2; j++) {
            for (int i = 0; i <= 2; i++) {
                MatT Hplus = MatT::Zero();
                for (int m = 0; m < 2; m++) {
                    for (int l = 0; l < 2; l++) {
                        Hplus += (Z(j, m) * I3).transpose() * Mp.block(2 * m, 2 * l, 2, 2) /*MpQuadrants(m, l)*/ * (Z(i, l) * I3);
                    }
                }
                Hplus *= 0.5*det_S;

                //auto Hplus = 0.5 * det_S * sum;

                _triplets.push_back({2*j,2*i,Hplus(0,0)});
                _triplets.push_back({2*j,2*i+1,Hplus(0,1)});
                _triplets.push_back({2*j+1,2*i,Hplus(1,0)});
                _triplets.push_back({2*j+1,2*i+1,Hplus(1,1)});

            }
        }
    }


    //not sure how to use this yet
    inline double max_feasible_step(const VecV &_x,
                                    const VecV &_v,
                                    const VecC &_c) const{
        //code from QuadOptimizer version. Not sure we should change it

        return DBL_MAX;

#if 0
        //std::cout<<" checking max feasible step"<<std::endl;
        if (_c[2] >= 0.0)
            return DBL_MAX;
        else {
            using VecT  = Eigen::Matrix<double, 2, 1, Eigen::ColMajor>;
            using MatT  = Eigen::Matrix<double, 2, 2, Eigen::ColMajor>;

            // get input points
            Eigen::Map<VecT> A((double *) &(_x[0]));
            Eigen::Map<VecT> B((double *) &(_x[2]));
            Eigen::Map<VecT> C((double *) &(_x[4]));

            Eigen::Map<VecT> DA((double *) &(_v[0]));
            Eigen::Map<VecT> DB((double *) &(_v[2]));
            Eigen::Map<VecT> DC((double *) &(_v[4]));


            // calc Jacobi Matrix
            MatT J;
            J.col(0) = B - A;
            J.col(1) = C - A;

            MatT DJ;
            DJ.col(0) = DB - DA;
            DJ.col(1) = DC - DA;

            // check roots of polynomial a t^2 + bt + c
            // find first positive root if it exists
            double a = DJ.determinant();
            double b = J(0, 0) * DJ(1, 1) + J(1, 1) * DJ(0, 0) - J(0, 1) * DJ(1, 0) - J(1, 0) * DJ(0, 1);
            double c = J.determinant();

            double t0, t1;
            int n_roots = COMISO::Polynomials::positive_roots_of_quadratic(a, b, c, t0, t1);

            if (n_roots >= 1)
                return t0;
            else
                return DBL_MAX;
        }
#endif
    }


    //TODO: jacobian method
    template<class ScalarT>
    using JacobianT  = Eigen::Matrix<ScalarT,2,2,Eigen::ColMajor>;

    template<class ScalarT>
        static JacobianT<ScalarT> jacobian(const Eigen::Matrix<ScalarT, NV, 1> &_x) {

        using VecT  = Eigen::Matrix<ScalarT,2,1,Eigen::ColMajor>;
        using MatT  = Eigen::Matrix<ScalarT,2,2,Eigen::ColMajor>;


        Eigen::Map<VecT> A((ScalarT *) &(_x[0]));
        Eigen::Map<VecT> B((ScalarT *) &(_x[2]));
        Eigen::Map<VecT> C((ScalarT *) &(_x[4]));

        MatT J;
        J.col(0) = B - A;
        J.col(1) = C - A;

        TriMeshTinyADFOFEnergy_F energy;
        return J * energy.S_inv_;

    }

private:
    ///
    double det_S_;
    Eigen::Matrix2d S_inv_;
};

using FF2D_AD = COMISO::FiniteElementTinyAD<TriMeshTinyADFOFEnergy_F>;
using FF2D_PH_AD = COMISO::FiniteElementHessianProjection< FF2D_AD >;




class SquaredDistanceElement {
public :


    // define dimensions
    const static int NV = 2; // [vertex_x, vertex_y]
    const static int NC = 2; // [target_x, target_y]

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec2T   = Eigen::Matrix<ScalarT,2,1,Eigen::ColMajor>;

        //std::cout<<" - x: "<<_x.transpose()<<std::endl;
        //std::cout<<" - c: "<<_c.transpose()<<std::endl;

        Vec2T d = (_x - _c);

        //std::cout<<" - e = "<<(d.dot(d))<<std::endl;

        return d.dot(d);

    }

    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        return std::numeric_limits<double>::max();
    }

};

using QD_2D_AD = COMISO::FiniteElementTinyAD<SquaredDistanceElement>;
using QD_2D_PH_AD = COMISO::FiniteElementHessianProjection< QD_2D_AD >;


class PullingElement : public SquaredDistanceElement {
public :

    // define dimensions
    const static int NV = 2; // [vertex_x, vertex_y]
    const static int NC = 2; // [direction_x, direction_y]

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        VecC xc;
        xc << TinyAD::to_passive(_x[0]), TinyAD::to_passive(_x[1]);

        return SquaredDistanceElement::eval_f(_x, xc + _c);
    }
};


using PULL_2D_AD = COMISO::FiniteElementTinyAD<PullingElement>;
using PULL_2D_PH_AD = COMISO::FiniteElementHessianProjection< PULL_2D_AD >;




// ================================================================================
// ============================= BASE CLASSES
// ================================================================================


class RegularRefElement{

public:

    RegularRefElement(){
        S_inv_ << 2, 1,
                0, std::sqrt(3.0);
        S_inv_ *= 0.5;
        S_inv_ = S_inv_.inverse().eval();
    }


protected:

    Eigen::Matrix2d S_inv_;

};


class JacobianElement2D : public RegularRefElement {
public :

    JacobianElement2D() : RegularRefElement() {}

    template<class ScalarT>
    using JacobianT = Eigen::Matrix<ScalarT, 2, 2, Eigen::ColMajor>;


    template<class ScalarT>
    JacobianT<ScalarT> jacobian(const Eigen::Matrix<ScalarT, 6, 1> &_x) const {
        using VecT = Eigen::Matrix<ScalarT, 2, 1, Eigen::ColMajor>;
        using MatT = Eigen::Matrix<ScalarT, 2, 2, Eigen::ColMajor>;

        Eigen::Map <VecT> A((ScalarT *) &(_x[0]));
        Eigen::Map <VecT> B((ScalarT *) &(_x[2]));
        Eigen::Map <VecT> C((ScalarT *) &(_x[4]));

        //std::cout<<" - S-1: "<<S_inv_<<std::endl;

        MatT J;
        J.col(0) = B - A;
        J.col(1) = C - A;

        return J * S_inv_;
    }

    template<class ScalarT>
    static JacobianT<ScalarT> compute_jacobian(const Eigen::Matrix<ScalarT, 6, 1> &_x) {
        JacobianElement2D elem;
        return elem.jacobian(_x);
    }
};



// ================================================================================
// ============================= SYMMETRIC DIRICHLET WITH BARRIER
// ================================================================================


#define USE_DAVIDS_VERSION 0

class SymmetricDirichletWithBarrierElement2D : public JacobianElement2D {
public:

    // define dimensions
    const static int NV = 6; // [ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz] = 2D points A, B, C of a triangle
    //const static int NC = 12; // [weight, barrier_w, barrier_min, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z]
    const static int NC = 7; // [weight, barrier_w, barrier_min, E_ref00, E_ref01, E_ref10, E_ref11]

    using VecV = Eigen::Matrix<double, NV, 1>;
    using VecC = Eigen::Matrix<double, NC, 1>;

    static std::string element_name(){ return "Symmetric Dirichlet With Barrier 2D"; }

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const {
        using Vec2T = Eigen::Matrix<ScalarT, 2, 1, Eigen::ColMajor>;
        using Mat2x2T = Eigen::Matrix<ScalarT, 2, 2, Eigen::ColMajor>;
        using Mat2x3T = Eigen::Matrix<ScalarT, 2, 3, Eigen::ColMajor>;
        using Mat3x2T = Eigen::Matrix<ScalarT, 3, 2, Eigen::RowMajor>;
        using Vec2dT = Eigen::Matrix<double, 2, 1, Eigen::ColMajor>;
        using Mat2x2dT = Eigen::Matrix<double, 2, 2, Eigen::ColMajor>;
        using Mat2x3dT = Eigen::Matrix<double, 2, 3, Eigen::RowMajor>;


        // get constants
        const double &w = _c[0];
        const double &barrier_w = _c[1];
        const double &barrier_min = _c[2];


#if USE_DAVIDS_VERSION

        // get image positions
        Eigen::Map<Mat3x2T> U( (ScalarT*) &(_x[0]) );

        Eigen::Map<Mat2x3dT> G( (double*) &(_c[3]) );

        Mat2x2T J = U.transpose()*G.transpose();

        ScalarT d = J.determinant();

#else

        // get ref edge  matrix
        //Eigen::Map <Mat2x2dT> E_ref((double *) &(_c[3]));
        //Mat2x2dT temp;
        //temp << _c[3], _c[4], _c[5], _c[6];
        //Mat2x2dT E_ref = temp.transpose();
        Mat2x2dT E_ref;

        E_ref.col(0) = _c.segment<2>(3);
        E_ref.col(1) = _c.segment<2>(5);

        // get points of triangle
        //Eigen::Map <Mat2x3T> P((ScalarT *) &(_x[0]));
        Mat3x2T P;
        P <<_x[0], _x[1],_x[2],_x[3],_x[4],_x[5];

        //std::cout<<" P: "<<std::endl<<P<<std::endl;

        // get edge vectors
        Mat2x2T E;
        E.col(0) = P.row(1) - P.row(0);
        E.col(1) = P.row(2) - P.row(0);

        //std::cout<<" E: "<<std::endl<<E<<std::endl;

        // calc  Jacobi Matrix of map
        Mat2x2T J = E * E_ref.inverse();

        //regular ref element
        //Mat2x2T J = jacobian(_x);
        //std::cout<<" J: "<<std::endl<<J<<std::endl;

        // calculate determinant
        ScalarT d = J.determinant();
        //std::cout<<" d = "<<d<<std::endl;

#endif

        //David's version

        /*std::cout<<" - x: "<<_x.transpose()<<std::endl;
        std::cout<<" - E: "<<std::endl<<E<<std::endl;
        std::cout<<" - E_ref: "<<std::endl<<E_ref<<std::endl;
        std::cout<<" - det = "<<d<<std::endl;
        std::cout<<" - w = "<<w<<std::endl;
        std::cout<<" - barrier w = "<<barrier_w<<std::endl;
        std::cout<<" - barrier min = "<<barrier_min<<std::endl;*/

        if (d <= 0.0) {
            return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
        } else {
            //auto E_sd = (1.0+d)*(J.squaredNorm() + J.inverse().squaredNorm() - 4.0);
            ScalarT E_sd = (J.squaredNorm() + J.inverse().squaredNorm() - 4.0);
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

    // const double _w, const double _barrier_w, const double _barrier_c,
    static int compute_constants(const std::vector<double> parameters,
            //const Eigen::Vector3d& _p0, const Eigen::Vector3d& _p1, const Eigen::Vector3d& _p2, const Eigen::Vector3d& _p3,
                                 const std::vector <Eigen::Vector2d> &_p,
                                 VecC &_c){

#if USE_DAVIDS_VERSION
        return compute_davids_constants(_w, _barrier_w, _barrier_c, _p, _c);

#else

        _c[0] = parameters[0];
        _c[1] = parameters[1];
        _c[2] = parameters[2];

        Eigen::Matrix<double ,2,2,Eigen::ColMajor> E;
        E.col(0) = _p[1]-_p[0];
        E.col(1) = _p[2]-_p[0];

        double d  = E.determinant();
        double di = 1.0/d;
        if(!std::isfinite(di) || d < 0.0)
        {
            std::cerr << "ERROR: determinant of edge vectors numerically not valid ---> " << d << std::endl;
            di = 0.0;
            return -1;
        }

        _c.segment<2>(3) = E.col(0);
        _c.segment<2>(5) = E.col(1);

        return 0;
#endif

    }

    // calculate constants from reference tetrahedron
    static int compute_davids_constants( const double _w, const double _barrier_w, const double _barrier_c,
            //const Eigen::Vector3d& _p0, const Eigen::Vector3d& _p1, const Eigen::Vector3d& _p2, const Eigen::Vector3d& _p3,
                                         const std::vector<Eigen::Vector2d>& _p,
                                         VecC& _c)
    {
        _c[0] = _w;
        _c[1] = _barrier_w;
        _c[2] = _barrier_c;
        Eigen::Matrix<double ,2,3,Eigen::ColMajor> E;
        E.col(0) = _p[2]-_p[1];
        E.col(1) = _p[0]-_p[2];
        E.col(2) = _p[1]-_p[3];

        double d  = E.block<2,2>(0,0).determinant();
        double di = 1.0/d;
        if(!std::isfinite(di) || d < 0.0)
        {
            std::cerr << "ERROR: determinant of edge vectors numerically not valid ---> " << d << std::endl;
            di = 0.0;
        }

        auto rot90 = [](const Eigen::Vector2d& _v) { return Eigen::Vector2d(-_v[1],_v[0]); };

        Eigen::Matrix<double ,2,3,Eigen::ColMajor> G;
        G.col(0) = di*rot90(E.col(0));
        G.col(1) = di*rot90(E.col(1));
        G.col(2) = di*rot90(E.col(2));


        _c.segment<3>(3) = G.row(0);
        _c.segment<3>(6) = G.row(1);

        // -------------------
        // test and debug code
        // -------------------
#if 0
        {
            Eigen::Matrix<double,3,2> U;
            U.row(0) = _p0;
            U.row(1) = _p1;
            U.row(2) = _p2;

            Eigen::Map<Eigen::Matrix<double ,2,3,Eigen::RowMajor> > G2( (double*) &(_c[1]) );

            Eigen::Matrix2d J = U.transpose()*G2.transpose();

            J -= Eigen::Matrix2d::Identity();

            double Jn = J.norm();
            if(Jn > 1e-6)
            {
                std::cerr << "ERROR: Jacobi matrix was not computed correctly!" << std::endl;
                J = U.transpose()*G2.transpose();
                std::cerr << J << std::endl;
            }
            else
                std::cerr << "Jacobi matrix within tolerance ---> " << Jn << std::endl;
        }
#endif

        return 0;
    }


};

using SDEB2D_AD = COMISO::FiniteElementTinyAD<SymmetricDirichletWithBarrierElement2D>;
using SDEB2D_PH_AD = COMISO::FiniteElementHessianProjection< SDEB2D_AD >;



// ================================================================================
// ============================= SYMMETRIC DIRICHLET FOR 3D TRIANGLE
// ================================================================================


class SymmetricDirichlet3D {
public:

    // define dimensions
    const static int NV = 9; //  3D points A, B, C of a triangle
    const static int NC = 10; // [weight, 3D points of a reference triangle

    using VecI = Eigen::Matrix<double, NV, 1>;
    using VecV = Eigen::Matrix<double, NV, 1>;
    using VecC = Eigen::Matrix<double, NC, 1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const {
        using Vec2T = Eigen::Matrix<ScalarT, 2, 1, Eigen::ColMajor>;
        using Vec3T = Eigen::Matrix<ScalarT, 3, 1, Eigen::ColMajor>;
        using Mat2x2T = Eigen::Matrix<ScalarT, 2, 2, Eigen::ColMajor>;
        using Mat2x3T = Eigen::Matrix<ScalarT, 2, 3, Eigen::ColMajor>;
        using Mat3x3T = Eigen::Matrix<ScalarT, 3, 3, Eigen::RowMajor>;

        using Vec2dT = Eigen::Matrix<double, 2, 1, Eigen::ColMajor>;
        using Vec3dT = Eigen::Matrix<double, 3, 1, Eigen::ColMajor>;
        using Mat2x2dT = Eigen::Matrix<double, 2, 2, Eigen::ColMajor>;
        using Mat2x3dT = Eigen::Matrix<double, 2, 3, Eigen::RowMajor>;


        // get constants
        const double &w = _c[0];

        Vec3T p0, p1, p2;
        p0 << _x[0], _x[1],_x[2];
        p1 << _x[3], _x[4],_x[5];
        p2 << _x[6], _x[7],_x[8];

        Vec3dT ref_p0 = _c.segment<3>(1);
        Vec3dT ref_p1 = _c.segment<3>(4);
        Vec3dT ref_p2 = _c.segment<3>(7);

        Vec3T e1 = p1 - p0;
        Vec3T e2 = p2 - p0;
        //std::cout<<" - e1: "<<e1.transpose()<<std::endl;
        //std::cout<<" - e2: "<<e2.transpose()<<std::endl;

        Vec3dT ref_e1 = ref_p1 - ref_p0;
        Vec3dT ref_e2 = ref_p2 - ref_p0;
        //std::cout<<" - ref e1: "<<ref_e1.transpose()<<std::endl;
        //std::cout<<" - ref e2: "<<ref_e2.transpose()<<std::endl;

        Vec3T n = e1.cross(e2).normalized();
        Vec3dT ref_n = ref_e1.cross(ref_e2).normalized();
        //std::cout<<" -     n: "<<n.transpose()<<std::endl;
        //std::cout<<" - ref_n: "<<ref_n.transpose()<<std::endl;

        Vec3T x_axis = e1.normalized();
        Vec3dT ref_x_axis = ref_e1.normalized();

        Vec3T y_axis = n.cross(x_axis).normalized();
        Vec3dT ref_y_axis = ref_n.cross(ref_x_axis).normalized();
        //std::cout<<" -     y axis: "<<y_axis.transpose()<<std::endl;
        //std::cout<<" - ref y axis: "<<ref_y_axis.transpose()<<std::endl;

        Vec2T proj_e1 = {e1.dot(x_axis), 0};
        Vec2T proj_e2 = {e2.dot(x_axis), e2.dot(y_axis)};
        //std::cout<<" -     projected e1: "<<proj_e1.transpose()<<std::endl;
        //std::cout<<" -     projected e2: "<<proj_e2.transpose()<<std::endl;

        Vec2dT ref_proj_e1 = {ref_e1.dot(ref_x_axis), 0};
        Vec2dT ref_proj_e2 = {ref_e2.dot(ref_x_axis), ref_e2.dot(ref_y_axis)};
        //std::cout<<" - ref projected e1: "<<ref_proj_e1.transpose()<<std::endl;
        //std::cout<<" - ref projected e2: "<<ref_proj_e2.transpose()<<std::endl;

        Mat2x2T E;
        E << proj_e1[0], proj_e2[0], proj_e1[1], proj_e2[1];
        //std::cout<<" - E: "<<std::endl<<E<<std::endl;

        Mat2x2dT ref_E;
        ref_E << ref_proj_e1[0], ref_proj_e2[0], ref_proj_e1[1], ref_proj_e2[1];

        //std::cout<<" - E_ref: "<<std::endl<<ref_E<<std::endl;

        Mat2x2T J = ref_E.inverse() * E;
        //std::cout<<" - J: "<<std::endl<<J<<std::endl;

        ScalarT d = J.determinant();

        if (d <= 0.0) {
            return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
        } else {
            //auto E_sd = (1.0+d)*(J.squaredNorm() + J.inverse().squaredNorm() - 4.0);
            ScalarT E_sd = (J.squaredNorm() + J.inverse().squaredNorm() - 4.0);
            //std::cout<<" E_sd = "<<E_sd<<std::endl;

            return w * E_sd;
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

};

using SDE_TRI_3D_AD = COMISO::FiniteElementTinyAD<SymmetricDirichlet3D>;
using SDE_TRI_3D_PH_AD = COMISO::FiniteElementHessianProjection< SDE_TRI_3D_AD >;



// ============================================================================================
// ============================= SYMMETRIC DIRICHLET WITH BARRIER AND VARIABLE REFERENCE
// ============================================================================================


#define USE_DAVIDS_VERSION 0

class SymmetricDirichletWithBarrierAndVariableReferenceElement2D {
public:

    // define dimensions
    const static int NV = 12; // 3 2D points A, B, C of a triangle and 3 2D points A', B' and C' of a reference element
    //const static int NC = 12; // [weight, barrier_w, barrier_min, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z]
    const static int NC = 3; // [weight, barrier_w, barrier_min]

    using VecV = Eigen::Matrix<double, NV, 1>;
    using VecC = Eigen::Matrix<double, NC, 1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const {
        using Vec2T = Eigen::Matrix<ScalarT, 2, 1, Eigen::ColMajor>;
        using Mat2x2T = Eigen::Matrix<ScalarT, 2, 2, Eigen::ColMajor>;
        using Mat2x3T = Eigen::Matrix<ScalarT, 2, 3, Eigen::ColMajor>;
        using Mat3x2T = Eigen::Matrix<ScalarT, 3, 2, Eigen::RowMajor>;
        using Vec2dT = Eigen::Matrix<double, 2, 1, Eigen::ColMajor>;
        using Mat2x2dT = Eigen::Matrix<double, 2, 2, Eigen::ColMajor>;
        using Mat2x3dT = Eigen::Matrix<double, 2, 3, Eigen::RowMajor>;


        // get constants
        const double &w = _c[0];
        const double &barrier_w = _c[1];
        const double &barrier_min = _c[2];


        // get points of triangle
        //Eigen::Map <Mat2x3T> P((ScalarT *) &(_x[0]));
        Mat3x2T P, P_ref;
        P <<_x[0], _x[1],_x[2],_x[3],_x[4],_x[5];
        P_ref <<_x[6], _x[7],_x[8],_x[9],_x[10],_x[11];

        //std::cout<<" P: "<<std::endl<<P<<std::endl;
        //std::cout<<" P_ref: "<<std::endl<<P_ref<<std::endl;

        // get edge vectors
        Mat2x2T E, E_ref;
        E.col(0) = P.row(1) - P.row(0);
        E.col(1) = P.row(2) - P.row(0);

        E_ref.col(0) = P_ref.row(1) - P_ref.row(0);
        E_ref.col(1) = P_ref.row(2) - P_ref.row(0);

        //std::cout<<" E: "<<std::endl<<E<<std::endl;
        //std::cout<<" E_ref: "<<std::endl<<E_ref<<std::endl;

        // calc  Jacobi Matrix of map
        Mat2x2T J = E * E_ref.inverse();

        //regular ref element
        //Mat2x2T J = jacobian(_x);
        //std::cout<<" J: "<<std::endl<<J<<std::endl;

        // calculate determinant
        ScalarT d = J.determinant();
        //std::cout<<" d = "<<d<<std::endl;


        //David's version

        /*std::cout<<" - x: "<<_x.transpose()<<std::endl;
        std::cout<<" - E: "<<std::endl<<E<<std::endl;
        std::cout<<" - E_ref: "<<std::endl<<E_ref<<std::endl;
        std::cout<<" - det = "<<d<<std::endl;
        std::cout<<" - w = "<<w<<std::endl;
        std::cout<<" - barrier w = "<<barrier_w<<std::endl;
        std::cout<<" - barrier min = "<<barrier_min<<std::endl;*/

        if (d <= 0.0) {
            return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
        } else {
            //auto E_sd = (1.0+d)*(J.squaredNorm() + J.inverse().squaredNorm() - 4.0);
            ScalarT E_sd = (J.squaredNorm() + J.inverse().squaredNorm() - 4.0);
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

using SDEB_VR_2D_AD = COMISO::FiniteElementTinyAD<SymmetricDirichletWithBarrierAndVariableReferenceElement2D>;
//using SDEB_VR_2D_PH_AD = COMISO::FiniteElementHessianProjection< SDEB_VR_2D_AD >;


// ================================================================================
// ============================= EDGE GLUEING ELEMENT
// ================================================================================


//probably make this generic
class EdgeGlueingElement2D
{
public:

    // define dimensions
    const static int NV = 8; // [] = two 2D points A_left, B_left, A_right, B_right, of two edges to glue
    const static int NC = 1; // [weight]

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec2T    = Eigen::Matrix<ScalarT,2,1,Eigen::ColMajor>;

        // get constants
        const double& w         = _c[0];

        Eigen::Map<Vec2T> Al( (ScalarT*) &(_x[0]) );
        Eigen::Map<Vec2T> Bl( (ScalarT*) &(_x[2]) );
        Eigen::Map<Vec2T> Ar( (ScalarT*) &(_x[4]) );
        Eigen::Map<Vec2T> Br( (ScalarT*) &(_x[6]) );

        /*std::cout<<" Al: "<<Al.transpose()<<std::endl;
        std::cout<<" Bl: "<<Bl.transpose()<<std::endl;
        std::cout<<" Ar: "<<Ar.transpose()<<std::endl;
        std::cout<<" Br: "<<Br.transpose()<<std::endl;*/

        Vec2T left = Bl - Al;
        Vec2T right = Br - Ar;

        //std::cout<<" - left edge: "<<left.transpose()<<std::endl;
        //std::cout<<" - right edge: "<<right.transpose()<<std::endl;

        Vec2T diff = left - right;

        //std::cout<<" - diff = "<<diff.transpose()<<std::endl;
        //std::cout<<" - e = "<<(left - right).squaredNorm()<<std::endl;

        ScalarT norm_w = 1.0;//0.5 * (left.norm() + right.norm());

        //return w * (left - right).template lpNorm<1>();
        return (w * (left - right).squaredNorm()) / norm_w;

    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        return std::numeric_limits<double>::max();
    }

};

using EGLUE2D_AD = COMISO::FiniteElementTinyAD<EdgeGlueingElement2D>;
using EGLUE2D_PH_AD = COMISO::FiniteElementHessianProjection< EGLUE2D_AD >;



// ================================================================================
// ============================= TRANSLATION ELEMENT
// ================================================================================

template<int DIM>
class TranslationElement
{
public:

    // define dimensions
    const static int NV = 2 * DIM; // [] = two 2D translation vectors for the left and right triangles
    const static int NC = 4 * DIM; // 2D points A_left, B_left, A_right, B_right of two edges (connecting duplicates of vertices a and b)

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using VecCT    = Eigen::Matrix<double,DIM,1,Eigen::ColMajor>;
        using VecT    = Eigen::Matrix<ScalarT,DIM,1,Eigen::ColMajor>;

        // get constants
        Eigen::Map<VecCT> Al( (double*) &(_c[0 * DIM]) );
        Eigen::Map<VecCT> Bl( (double*) &(_c[1 * DIM]) );
        Eigen::Map<VecCT> Ar( (double*) &(_c[2 * DIM]) );
        Eigen::Map<VecCT> Br( (double*) &(_c[3 * DIM]) );

        Eigen::Map<VecT> tl( (ScalarT*) &(_x[0 * DIM]) );
        Eigen::Map<VecT> tr( (ScalarT*) &(_x[1 * DIM]) );

        /*std::cout<<" Al: "<<Al.transpose()<<std::endl;
        std::cout<<" Bl: "<<Bl.transpose()<<std::endl;
        std::cout<<" Ar: "<<Ar.transpose()<<std::endl;
        std::cout<<" Br: "<<Br.transpose()<<std::endl;*/

        VecT a_diff = (Al + tl) - (Ar + tr);
        VecT b_diff = (Bl + tl) - (Br + tr);

        /*std::cout<<" - A diff : "<<a_diff.transpose()<<std::endl;
        std::cout<<" - B diff : "<<b_diff.transpose()<<std::endl;*/

        ScalarT e = a_diff.squaredNorm() + b_diff.squaredNorm();

        //std::cout<<" - e = "<<e<<std::endl;

        return e;

    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        return std::numeric_limits<double>::max();
    }

};

template< int DIM>
using TRANSLATION_AD = COMISO::FiniteElementTinyAD<TranslationElement<DIM>>;

template< int DIM>
using TRANSLATION_PH_AD = TRANSLATION_AD<DIM>;





//minimizes the distance ||a + t - b || with a, b being two fix points and t being a translation vector
template<int DIM>
class SinglePointTranslationElement
{
public:

    // define dimensions
    const static int NV = DIM; // [] = translation vector
    const static int NC = 2 * DIM; // 2 points a and b

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using VecCT    = Eigen::Matrix<double,DIM,1,Eigen::ColMajor>;
        using VecT    = Eigen::Matrix<ScalarT,DIM,1,Eigen::ColMajor>;

        // get constants
        Eigen::Map<VecCT> A( (double*) &(_c[0 * DIM]) );
        Eigen::Map<VecCT> B( (double*) &(_c[1 * DIM]) );
        Eigen::Map<VecT> t( (ScalarT*) &(_x[0 * DIM]) );

        VecT diff = (A + t) - B;

        ScalarT e = diff.squaredNorm();

        return e;
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        return std::numeric_limits<double>::max();
    }

};

template< int DIM>
using SP_TRANSLATION_AD = COMISO::FiniteElementTinyAD<SinglePointTranslationElement<DIM>>;

template< int DIM>
using SP_TRANSLATION_PH_AD = SP_TRANSLATION_AD<DIM>;




// ================================================================================
// ============================= FOLDOVER-FREE ELEMENT
// ================================================================================


class FoldoverFreeElement2D : public JacobianElement2D
{
public:

    // define dimensions
    const static int NV = 6; // [ax,ay,bx,by,cx,cy] = 2D points A,B,C of a triangle
    const static int NC = 4; // [weight_conformal, weight_vol, epsilon, average area]


    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;
    using VecI = Eigen::Matrix<int,NV,1>;

    static std::string element_name(){ return "Foldover-Free 2D"; }

    FoldoverFreeElement2D() : JacobianElement2D() {}


    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {

        using VecT  = Eigen::Matrix<ScalarT,2,1,Eigen::ColMajor>;
        using MatT  = Eigen::Matrix<ScalarT,2,2,Eigen::ColMajor>;
        using MatdT = Eigen::Matrix<double,2,2,Eigen::RowMajor>;



        // get constants
        const double& w_conformal = _c[0];
        const double& w_vol       = _c[1];
        const double& epsilon     = _c[2];
        const double& average_vol = _c[3];


        auto J = jacobian(_x);
        auto det = J.determinant();

        if(epsilon > 0.0)
        {

            ScalarT det_reg = 0.5*(det + sqrt( pow(det,2) + epsilon*epsilon));
            ScalarT e =  w_conformal * (J.transpose() * J).trace()/det_reg + w_vol * ((pow(det,2)+1.0)/det_reg);
            //ScalarT max_e = _c[3];
            //return log10(e);
            return e;
        }
        else
        {
            if (det <= 0.0)
                return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
            else
                return w_conformal * (J.transpose() * J).trace()/pow(det,2./2.) + w_vol * (det + 1.0/det);
        }
    }

    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
    {
        double t = 2.0;
        for (int i = 0; i < 50; ++i) {
            VecV xn = _x + t * _v;
            if (std::isfinite(eval_f(xn, _c)))
                return t;
            else
                t *= 0.7;
        }
        return 0;
    }

    static int compute_constants(const double _w, const double _barrier_w, const double _barrier_c,
            //const Eigen::Vector3d& _p0, const Eigen::Vector3d& _p1, const Eigen::Vector3d& _p2, const Eigen::Vector3d& _p3,
                                 const std::vector <Eigen::Vector2d> &_p,
                                 VecC &_c){

        _c[0] = 1.0; //conformal
        _c[1] = 0.0; //volumetric
        _c[2] = 0.0; //epsilon
        _c[3] = 1.0; //avg vol
        return 0;
    }
};

using FF_2D_AD = COMISO::FiniteElementTinyAD<FoldoverFreeElement2D>;
using FF_2D_PH_AD = COMISO::FiniteElementHessianProjection< FF_2D_AD >;


///


#endif
