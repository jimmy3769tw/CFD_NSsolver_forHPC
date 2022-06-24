#ifndef __MPI_CG_H__
#define	__MPI_CG_H__

// #include <yakutat/mpi/matrix_definition.hpp>
#include <yakutat/backend/Matrix_traits.hpp>
#include <yakutat/backend/mpi_tools.hpp>

namespace yakutat::mpi
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




            bool partition(size_type myid, size_type nproc);
            bool partition(mpi_tools & world);
            bool setTolerance(value_type tolerance) noexcept;
            bool initialize(const MatrixType & lhs_mat);
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
	inline cg<MatrixType>::cg()
	{

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
    inline bool cg<MatrixType>::partition(size_type myid, size_type nproc)
    {
        world_.MPI_division(length_, myid, nproc);
        lhs_mat_.partition(myid, nproc);
        return true;
    }

    template<typename MatrixType>
    inline bool cg<MatrixType>::partition(mpi_tools & world)
    {
        if(world.length==0)
        {
            throw std::invalid_argument("mpicg partition: Your world no division.");
        }

        if(length_ != world.length)
        {
            throw std::invalid_argument("mpi cg partition: wrong partition!");
        }

        world_ = world;
        lhs_mat_.partition(world_);
        return true;
    }

    template<typename MatrixType>
    inline bool cg<MatrixType>::setTolerance(value_type tolerance) noexcept
    {
        zeta_ = tolerance;
        return true;
    }

    template<typename MatrixType>
    inline bool cg<MatrixType>::MPI_allgather(array_value_type & x, size_type itag)
    {
        return world_.MPI_vector_collect(x, itag);
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
        for(size_type i=world_.start; i<world_.end+1; ++i)
            r_[i] = rhs[i] - p_[i];

        #pragma omp simd
        for(size_type i=world_.start; i<world_.end+1; ++i)
            p_[i] = r_[i];

        nu_ = inner_product(r_, r_);
        
        


        norm = 0.0;
        norm = inner_product(r_, r_);
        norm = sqrt(norm) / norm0_;

        iters = 0;
        while(norm>zeta_ && iters<iters_max_)
        {
            ++ iters;

            world_.MPI_vector_collect(p_, 100);

            q_ = lhs_mat_ * p_;

            alpha_ =  nu_ / inner_product(p_, q_);

            #pragma omp simd 
            for(size_type i=world_.start; i<world_.end+1; ++i)
                x[i] += alpha_ * p_[i];

            #pragma omp simd
            for(size_type i=world_.start; i<world_.end+1; ++i)
                r_[i] -= alpha_ * q_[i];

            mu_ = inner_product(r_, r_);

            beta_ = mu_ / nu_;

            #pragma omp simd
            for(size_type i=world_.start; i<world_.end+1; ++i)
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
    inline cg<MatrixType>::~cg() {};
} // end namespace yakutat::mpi
#endif
