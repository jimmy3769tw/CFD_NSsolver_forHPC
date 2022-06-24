#ifndef __BICGSTAB_H__
#define	__BICGSTAB_H__

#include <yakutat/matrix_definition.hpp>
#include <yakutat/backend/Matrix_traits.hpp>

namespace yakutat
{
    using namespace yakutat::backend;

    template<typename MatrixType>
    class bicgstab
        : private MatrixType
    {
        // Member types
        using self_type = bicgstab<MatrixType>;

        

        public:

            // Member types
            using size_type = typename MatrixTraits<MatrixType>::size_type;
            using value_type = typename MatrixTraits<MatrixType>::value_type;
            using array_size_type = typename MatrixTraits<MatrixType>::array_size_type;
            using array_value_type = typename MatrixTraits<MatrixType>::array_value_type;
            using multi_output_type = std::tuple<size_type, value_type>;
            
            // --- Constructor & Destructor ---
            //---------------------------------------------//
            bicgstab();
            virtual ~bicgstab();
            //---------------------------------------------//


            bool initialize(const MatrixType & lhs_mat);
            bool setTolerance(value_type tolerance) noexcept;
            multi_output_type solve(const array_value_type & rhs, array_value_type & x); 



        private:
            //---------------------------------------------// Local variables
            //bicgstab variables
            size_type length_;
            value_type zeta_{1e-5};
            size_type iters_max_{3000};
            value_type alpha_,beta_,norm0_,norm_,sum_,scal_,norm1_,norm2_,omega_,rho1_,rho2_;
            array_value_type p_;
            array_value_type r_;
            array_value_type r2_;
            array_value_type v_;
            array_value_type ss_;
            array_value_type t_;
            MatrixType lhs_mat_;
            //---------------------------------------------//




            value_type inner_product(const array_value_type & a, const array_value_type & b);


        
    };


    //##################################################//
    //                Implementations                   //
    //##################################################//


    template<typename MatrixType>
	inline bicgstab<MatrixType>::bicgstab()
	{
		
	}

    template<typename MatrixType>
    inline bool bicgstab<MatrixType>::initialize(const MatrixType & lhs_mat)
    {
        lhs_mat_ = lhs_mat;
        length_ = lhs_mat.row();
        p_.resize(length_);
        r_.resize(length_);
        r2_.resize(length_);
        v_.resize(length_);
        ss_.resize(length_);
        t_.resize(length_);

        return true;
    }

    template<typename MatrixType>
    inline bool bicgstab<MatrixType>::setTolerance(value_type tolerance) noexcept
    {
        zeta_ = tolerance;
        return true;
    }

    template<typename MatrixType>
    inline typename bicgstab<MatrixType>::multi_output_type bicgstab<MatrixType>::solve(const array_value_type & rhs, array_value_type & x)
    {
        size_type iters{0};
        value_type norm{0};

        // ============== start calculate ============== //
        norm0_ = inner_product(rhs, rhs);

        norm0_ = sqrt(norm0_);

        p_ = lhs_mat_ * x;

        
        #pragma omp simd
        for(size_type i=0; i<length_; ++i)
            r_[i] = rhs[i] - p_[i];

        #pragma omp simd
        for(size_type i=0; i<length_; ++i) 
            r2_[i] = r_[i];

        rho1_  = 1; alpha_ = 1; omega_ = 1;

        #pragma omp simd
        for(size_type i=0; i<length_; ++i)
        {
            v_[i] = 0;
            p_[i] = 0;
        }

        norm = 0;

        norm = inner_product(r_, r_);
        norm = sqrt(norm) / norm0_;

        

        iters = 0;
        while(norm>zeta_ && iters<iters_max_)
        {
            ++ iters;

            rho2_ = inner_product(r2_, r_);
            beta_ = (rho2_/rho1_) * (alpha_/omega_);

            #pragma omp simd
            for(size_type i=0; i<length_; ++i)
                p_[i] = r_[i] + beta_ * (p_[i] - omega_ * v_[i]);


            v_ = lhs_mat_ * p_;

            alpha_ = rho2_ / inner_product(r2_, v_);

            #pragma omp simd
            for(size_type i=0; i<length_; ++i)
                ss_[i] = r_[i] - alpha_ * v_[i];


            t_ = lhs_mat_ * ss_;

            omega_ = inner_product(t_, ss_) / inner_product(t_, t_);

            #pragma omp simd
            for(size_type i=0; i<length_; ++i)
                x[i] += alpha_ * p_[i] + omega_ * ss_[i];

            #pragma omp simd
            for(size_type i=0; i<length_; ++i)
                r_[i] = ss_[i] - omega_ * t_[i];

            rho1_ = rho2_;

            norm = 0;

            norm = inner_product(r_, r_);
            norm = sqrt(norm) / norm0_;
        }

        
        return {iters, norm};
    }





    template<typename MatrixType>
    inline typename MatrixTraits<MatrixType>::value_type bicgstab<MatrixType>::inner_product(const array_value_type & a, const array_value_type & b)
    {
        value_type c{0.0};

        #pragma omp simd reduction(+:c)
        for(size_type i=0; i<length_; ++i)
        {
            c += a[i] * b[i];
        }
        return c;
    }


    
    template<typename MatrixType>
    inline bicgstab<MatrixType>::~bicgstab() {};


    
    


} // end namespace yakutat::mpi
#endif