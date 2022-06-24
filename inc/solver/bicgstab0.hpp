#pragma once
#include <tuple>
#include <vector>
#include "math.hpp"
namespace solver{
    using namespace std;
    using namespace math;
    //  TODO :debug!! this cann't not be used anymore!!
    template<typename matrixT>
    class bicgstab0
    {
        public:

        bicgstab0(){ }
        bicgstab0(const matrixT & lhs_mat){init(lhs_mat);}

        virtual ~bicgstab0(){ }

        // ! init
        inline bool init(const matrixT & lhs_mat)
        {
            lhs_mat_ = lhs_mat;
            length_ = lhs_mat.row();

            xm_.resize(length_);
            r_hat_.resize(length_);
            p_.resize(length_);
            r_.resize(length_);
            rm_.resize(length_);
            s_.resize(length_);
            v_.resize(length_);
            t_.resize(length_);
            h_.resize(length_);

            return true;
        }

        void setTolerance(double tolerance) { zeta_ = tolerance;}

        std::pair<int, double> solve( vector<double> & rhs, vector<double> & x); 

        std::pair<int, double> operator ()(const std::vector<double> & rhs, std::vector<double> & x)
        { return solve(rhs, x); }


        private:
            //bicgstab variables
            int length_;

            double zeta_ = 1e-5;

            int iters_max_= 3000;

            double accurate_bicg_ = 1.0e-10;

            double norm0_, rho_ , alpha_, omeaga_, rhom_;

            vector<double> xm_;
            vector<double> r_hat_;
            vector<double> p_;
            vector<double> r_;
            vector<double> rm_;
            vector<double> s_;
            vector<double> v_;
            vector<double> t_;
            vector<double> h_;


            matrixT lhs_mat_;
            // --------------------------
            inline void 
            minus(std::vector<double> & a, std::vector<double> & b, vector<double> & r)
            {
                #pragma omp parallel
                {
                    #pragma omp for 
                    for (size_t i = 0; i < a.size()  ; ++i)
                        r[i] = a[i] - b[i];
                }
            }
            // --------------------------

    };



    // ! main
    template<typename matrixT>
    inline std::pair<int, double> bicgstab0<matrixT>::
    solve(vector<double> & rhs, vector<double> & x)
    {
        zero(x);
        zero(r_);

        lhs_mat_.multiply_omp(x, p_);

        minus(rhs, p_, r_);

        copy(r_, r_hat_);

        rho_ = 1.0;
        alpha_ = 1.0;
        omeaga_ = 1.0;

        zero(p_);
        zero(v_);
        zero(s_);
        zero(t_);
        zero(h_);
        zero(rm_);
        zero(xm_);

        norm0_ = dot(rhs, rhs);

        const auto tolerance = pow((accurate_bicg_* norm0_), 2);

        int iter= 1; double error;
        for (double beta_ = 0.0; iter < iters_max_ ;++iter){

            copy(x, xm_);

            rhom_ = rho_;

            rho_ = dot(r_hat_, r_);

            beta_ = (rho_ / rhom_) * (alpha_ / omeaga_);
            
            for (int i = 0 ; i < length_ ;++i )
            { p_[i] = r_[i] + beta_ * (p_[i] - omeaga_ * v_[i]); }

            lhs_mat_.multiply_omp(p_, v_);

            alpha_ = rho_ / dot(r_hat_, v_);

            for (int i = 0 ; i < length_ ;++i )
            { h_[i] =  xm_[i]+ alpha_ * p_[i]; }   

            for (int i = 0 ; i < length_ ;++i )
            { s_[i] = r_[i] - alpha_ * v_[i]; }


            lhs_mat_.multiply_omp(s_, t_);

            omeaga_ = dot(t_, s_) / dot(t_, t_);

            for (int i = 0 ; i < length_ ;++i )
                x[i] = h_[i] + omeaga_ * s_[i];


            error = getError(xm_, x);
            if (error < tolerance) { break;}
             
            for (int i = 0 ; i < length_ ;++i )
            { r_[i] = s_[i] - omeaga_ * t_[i]; }

        }

    return std::make_pair( iter, sqrt(error) / norm0_);
    }
}

