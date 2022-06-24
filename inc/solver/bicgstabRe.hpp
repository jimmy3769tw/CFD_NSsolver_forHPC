#pragma once
#include <vector>
#include <cmath>

namespace solver{
    using namespace std;

    template<typename matrixT>
    class bicgstabRe
    {
        public:

        bicgstabRe(){ }

        bicgstabRe(matrixT & lhs_mat){ init(lhs_mat);}

        virtual ~bicgstabRe(){ }

        // ! init
        inline bool init(const matrixT & lhs_mat)
        {
            lhs_mat_ = lhs_mat;
            length_ = lhs_mat.row();
            R0_.resize(length_);
            R1_.resize(length_);
            P1_.resize(length_);
            S1_.resize(length_);
            AP_.resize(length_);
            AS_.resize(length_);
            PA_.resize(length_);
            return true;
        }

        void setTolerance(double tolerance) { zeta_ = tolerance;}

        std::pair<int, double> solve(const vector<double> & rhs, vector<double> & x); 


        std::pair<int, double> operator ()(const std::vector<double> & rhs, std::vector<double> & x)
        { return solve(rhs, x); }



        private:
            //bicgstab variables
            int length_;
            double zeta_{1e-5};
            double tolerance_{1e-5};
            int iters_max_{3000};
            double  alpha_, beta_, norm0_, 
                    norm_, sum_, scal_, 
                    norm1_, norm2_,
                    omega_, rho1_, rho2_;

            vector<double> R0_;
            vector<double> R1_;
            vector<double> P1_;
            vector<double> S1_;
            vector<double> AP_;
            vector<double> AS_;
            vector<double> PA_;
            int pIterTotal_ = 0, 
                pRestarts_  = 0;

            matrixT lhs_mat_;

        inline double inner_product(const vector<double> & a, const vector<double> & b)
        {
            double r{0.0};

            #pragma omp simd reduction(+:r)
            for(int i=0; i<length_; ++i)
            {  r += a[i] * b[i]; }

            return r;
        }


    };


    // ! main
    template<typename matrixT>
    inline std::pair<int, double> 
    bicgstabRe<matrixT>::solve(
        const vector<double> & rhs, vector<double> & x)
    {
        length_ = lhs_mat_.row();

        auto normVector = [&](auto &a){ return std::sqrt(inner_product(a, a)); };

        double
            rr = 0, r0 = 0, a = 0, w = 0, b = 0, 
            norm0 = normVector(rhs);

        for (int i = 0; i < length_; ++i)
        {
            x[i] = 0.0;
            R0_[i] = rhs[i];
            R1_[i] = R0_[i];
            P1_[i] = R0_[i];
        }
        rr = inner_product(R1_,R0_);

        int pIter = 0;

        while (pIter < iters_max_)
        {
            lhs_mat_.multiply_omp(P1_, AP_);

            a  = rr / inner_product(AP_, R0_);

            for (int i=0; i<length_; ++i) { S1_[i] = R1_[i] - a*AP_[i]; } // 3.

            lhs_mat_.multiply_omp(S1_, AS_);
        
            w  = inner_product(AS_, S1_) / inner_product(AS_, AS_);  //2. 

            for (int i=0; i<length_; ++i)
            {
                x[i] += a*P1_[i] + w*S1_[i];
                R1_[i] = S1_[i] - w*AS_[i];
            }

            if ((normVector(R1_) / norm0) < tolerance_) { break; }

            r0 = rr;
            rr = inner_product(R1_, R0_);
            b  = (a / w) * (rr / r0);

            // Check rho for restart
            if (abs(rr) > tolerance_*0.6*tolerance_)
            {
                for (int i=0; i<length_; ++i) { P1_[i] = R1_[i] + b*(P1_[i] - w*AP_[i]); }
            }
            else
            {
                for (int i = 0; i < length_; ++i)
                {
                    R0_[i] = R1_[i];
                    P1_[i] = R1_[i];
                }
                ++pRestarts_;
            }
            ++pIter;
            ++pIterTotal_;
        }



        return std::make_pair( pIter, std::sqrt(abs(rr)) / norm0);
    }
}


