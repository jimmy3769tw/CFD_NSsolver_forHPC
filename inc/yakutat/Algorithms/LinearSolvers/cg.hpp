#ifndef __CG_H__
#define	__CG_H__

#include <yakutat/matrix_definition.hpp>
#include <yakutat/backend/Matrix_traits.hpp>

namespace yakutat
{
    using namespace yakutat::backend;

    template<typename MatrixType>
    class cg
        : private MatrixType
    {
        // Member types
        using self_type = cg<MatrixType>;
        
        public:

            // Member types
            using size_type = typename MatrixTraits<MatrixType>::size_type;
            using value_type = typename MatrixTraits<MatrixType>::value_type;
            using array_size_type = typename MatrixTraits<MatrixType>::array_size_type;
            using array_value_type = typename MatrixTraits<MatrixType>::array_value_type;
            using multi_output_type = std::tuple<size_type, value_type>;

            // --- Constructor & Destructor ---
            //--------------------------------------------- //
            cg();
            virtual ~cg();
            //---------------------------------------------//



            bool setTolerance(value_type tolerance) noexcept;
            bool initialize(const MatrixType & lhs_mat);
            multi_output_type solve(const array_value_type & rhs, array_value_type & x); 


        private:
            //---------------------------------------------// Local variables
            //bicgstab variables
            size_type length_;
            value_type zeta_{1e-5};
            size_type iters_max_{3000};
            value_type alpha_, beta_, norm0_, z_, nu_, mu_;
            array_value_type p_;
            array_value_type q_;
            array_value_type r_;
            MatrixType lhs_mat_;
            //---------------------------------------------//

            value_type inner_product(const array_value_type & a, const array_value_type & b);
    };

    //##################################################//
    //                Implementations                   //
    //##################################################//


    template<typename MatrixType>
	inline cg<MatrixType>::cg()	{}


    template<typename MatrixType>
    inline bool cg<MatrixType>::setTolerance(value_type tolerance) noexcept
    {
        zeta_ = tolerance;
        return true;
    }

    template<typename MatrixType>
    inline bool cg<MatrixType>::initialize(const MatrixType & lhs_mat)
    {
        lhs_mat_ = lhs_mat;
        length_ = lhs_mat.row();
        p_.resize(length_);
        q_.resize(length_);
        r_.resize(length_);
        return true;
    }

    template<typename MatrixType>
    inline typename cg<MatrixType>::multi_output_type cg<MatrixType>::solve(const array_value_type & rhs, array_value_type & x)
    {
        size_type iters{0};
        value_type norm{0};
        
        
        // ============== start calculate ============== //
        norm0_ = inner_product(rhs, rhs);
        
        p_ = lhs_mat_ * x;
        
        #pragma omp simd
        for(size_type i=0; i<length_; ++i)
            r_[i] = rhs[i] - p_[i];
        
        #pragma omp simd
        for(size_type i=0; i<length_; ++i)
            p_[i] = r_[i];

        nu_ = inner_product(r_, r_);
        
        


        norm = 0.0;
        norm = inner_product(r_, r_);
        norm = sqrt(norm) / norm0_;
        

        iters = 0;
        while(norm>zeta_ && iters<iters_max_)
        {
            ++ iters;


            q_ = lhs_mat_ * p_;

            alpha_ =  nu_ / inner_product(p_, q_);

            #pragma omp simd 
            for(size_type i=0; i<length_; ++i)
                x[i] += alpha_ * p_[i];

            #pragma omp simd
            for(size_type i=0; i<length_; ++i)
                r_[i] -= alpha_ * q_[i];

            mu_ = inner_product(r_, r_);

            beta_ = mu_ / nu_;

            #pragma omp simd
            for(size_type i=0; i<length_; ++i)
                p_[i] = r_[i] +  beta_ * p_[i];

            nu_ = mu_;

            norm = 0.0;
            norm = inner_product(r_, r_);
            norm = sqrt(norm) / norm0_;

        }
        return {iters, norm};
    }




    template<typename MatrixType>
    inline typename MatrixTraits<MatrixType>::value_type cg<MatrixType>::inner_product(const array_value_type & a, const array_value_type & b)
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
    inline cg<MatrixType>::~cg() {};
} // end namespace yakutat
#endif