#ifndef __MPI_BICGSTAB_H__
#define	__MPI_BICGSTAB_H__

#include <yakutat/mpi/matrix_definition.hpp>
#include <yakutat/backend/Matrix_traits.hpp>
#include <yakutat/backend/mpi_tools.hpp>

namespace yakutat::mpi
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


            bool partition(size_type myid, size_type nproc);
            bool partition(mpi_tools & world);
            bool initialize(const MatrixType & lhs_mat);
            bool setTolerance(value_type tolerance) noexcept;
            multi_output_type solve(const array_value_type & rhs, array_value_type & x); 
            bool MPI_allgather(array_value_type & x, size_type itag);


        private:
            //---------------------------------------------// Local variables
            //MPI variables
            mpi_tools world_;

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
	inline bicgstab<MatrixType>::bicgstab() {}

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
    inline bool bicgstab<MatrixType>::partition(size_type myid, size_type nproc)
    {
        world_.MPI_division(length_, myid, nproc);
        lhs_mat_.partition(myid, nproc);
        return true;
    }

    template<typename MatrixType>
    inline bool bicgstab<MatrixType>::partition(mpi_tools & world)
    {
        if(world.length==0)
        {
            throw std::invalid_argument("mpi bicgbstab partition: Your world no division.");
        }

        if(length_ != world.length)
        {
            throw std::invalid_argument("mpi bicgbstab partition: wrong partition!");
        }

        world_ = world;
        lhs_mat_.partition(world_);
        return true;
    }

    template<typename MatrixType>
    inline bool bicgstab<MatrixType>::setTolerance(value_type tolerance) noexcept
    {
        zeta_ = tolerance;
        return true;
    }

    template<typename MatrixType>
    inline bool bicgstab<MatrixType>::MPI_allgather(array_value_type & x, size_type itag)
    {
        return world_.MPI_vector_collect(x, itag);
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
        for(size_type i=world_.start; i<world_.end+1; ++i)
            r_[i] = rhs[i] - p_[i];

        #pragma omp simd
        for(size_type i=world_.start; i<world_.end+1; ++i)   
            r2_[i] = r_[i];

        rho1_  = 1; alpha_ = 1; omega_ = 1;

        #pragma omp simd
        for(size_type i=world_.start; i<world_.end+1; ++i)
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
            for(size_type i=world_.start; i<world_.end+1; ++i)  
                p_[i] = r_[i] + beta_ * (p_[i] - omega_ * v_[i]);

            world_.MPI_vector_collect(p_, 100);

            v_ = lhs_mat_ * p_;

            alpha_ = rho2_ / inner_product(r2_, v_);

            #pragma omp simd
            for(size_type i=world_.start; i<world_.end+1; ++i)
                ss_[i] = r_[i] - alpha_ * v_[i];

            world_.MPI_vector_collect(ss_, 101);

            t_ = lhs_mat_ * ss_;

            omega_ = inner_product(t_, ss_) / inner_product(t_, t_);

            #pragma omp simd
            for(size_type i=world_.start; i<world_.end+1; ++i)
                x[i] += alpha_ * p_[i] + omega_ * ss_[i];

            #pragma omp simd
            for(size_type i=world_.start; i<world_.end+1; ++i)
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
        value_type c1{0.0};

        #pragma omp simd reduction(+:c)
        for(size_type i=world_.start; i<world_.end+1; ++i)
        {
            c += a[i] * b[i];
        }
        MPI_Allreduce((void *)&c, (void *)&c1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return c1;
    }


    
    template<typename MatrixType>
    inline bicgstab<MatrixType>::~bicgstab() {};


    
    


} // end namespace yakutat::mpi
#endif