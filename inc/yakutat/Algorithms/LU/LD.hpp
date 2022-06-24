#ifndef __LD_H__
#define	__LD_H__
#include <yakutat/matrix_definition.hpp>
#include <yakutat/mpi/matrix_definition.hpp>
#include <yakutat/backend/Matrix_operations.hpp>
#include <yakutat/backend/Matrix_traits.hpp>

namespace yakutat
{
    using namespace yakutat::backend;

    template<typename MatrixType>
    class LD
        : private MatrixType
    {
        // Member types
        using self_type = LD<MatrixType>;

        public:
            // Member types
            using size_type = typename MatrixTraits<MatrixType>::size_type;
            using value_type = typename MatrixTraits<MatrixType>::value_type;
            using array_size_type = typename MatrixTraits<MatrixType>::array_size_type;
            using array_value_type = typename MatrixTraits<MatrixType>::array_value_type;
            using multi_output_type = std::tuple<MatrixType &, MatrixType &>;
            using matops_type = MatOps<MatrixType>;

            LD();
            virtual ~LD();

            bool factorize(const MatrixType & matrix);
            MatrixType & matrixLD_inplace(void);
            MatrixType & matrixLD(void);
            bool solve(const array_value_type & rhs, array_value_type & x);
            bool solve(const MatrixType & rhs_mat, MatrixType & x_mat);


        private:
            MatrixType LD_inplace_, LD_;
            size_type row_{0}, col_{0};
            matops_type matops_;
    };

    //##################################################//
    //                Implementations                   //
    //##################################################//

    template<typename MatrixType>
	inline LD<MatrixType>::LD() {}

    template<typename MatrixType>
    inline LD<MatrixType>::~LD() {}

    template<typename MatrixType>
    inline bool LD<MatrixType>::factorize(const MatrixType & matrix)
    {
        row_ = matrix.row();
        col_ = matrix.col();

        value_type Lij{0.0};
        LD_inplace_ = matrix;
        for(size_type i=1; i< row_; ++i)
        {

            for (size_type j = i; j >= 1; --j) 
            {
                Lij=LD_inplace_(i,i-j);
                
                for (size_type k = i; k > j; k--)
                {
                    Lij -= LD_inplace_(i,i-k) * LD_inplace_(i-j,i-k) * LD_inplace_(i-k,i-k);
                }
                Lij /= LD_inplace_(i-j,i-j);
                LD_inplace_.set(i,i-j, Lij);
            }
            

            for (size_type k = i; k >= 1; k--)
            {
                LD_inplace_.set(i,i, LD_inplace_(i,i) - LD_inplace_(i-k,i-k)*LD_inplace_(i,i-k)*LD_inplace_(i,i-k) );
            }
        }
        return true;
    }

    template<typename MatrixType>
    inline MatrixType & 
    LD<MatrixType>::matrixLD(void)
    {
        LD_ = LD_inplace_;
        for(size_type i=0; i<row_; ++i)
        {
            for(size_type j=i+1; j<col_; ++j)
            {
                LD_.set(i,j,0);
            }
        }
        return LD_;
    }

    template<typename MatrixType>
    inline MatrixType & 
    LD<MatrixType>::matrixLD_inplace(void)
    {
        return LD_inplace_;
    }



    template<typename MatrixType>
    inline bool 
    LD<MatrixType>::solve(const array_value_type & rhs, array_value_type & x)
    {
        array_value_type temp(row_);

        /*   Solve for intermediate column vector solution  */
        for(size_type i=0; i<row_; ++i)
        {
            temp[i] = rhs[i];
            for (size_type k = i; k >= 1; --k) 
            {
                temp[i] -= LD_inplace_.at(i,i-k) * temp[i-k];
            }
        }

        /*   Solve for final column vector solution   */
        for(int k=row_-1; k>=0; k--)
        {
            x[k] = temp[k] / LD_inplace_.at(k,k);
            for(size_type j=k+1; j<row_; ++j)
            {
                x[k] -= LD_inplace_.at(j,k) * x[j];
            }
        }
        return true;
    }


    template<typename MatrixType>
    bool 
    LD<MatrixType>::solve(const MatrixType & rhs_mat, MatrixType & x_mat)
    {
        array_value_type temp(row_);

        #pragma omp parallel for
        for(size_type col_idx=0; col_idx<col_; ++col_idx)
         {



            /*   Solve for intermediate column vector solution  */
            for(size_type i=0; i<row_; ++i)
            {
                // temp[i] = rhs[i];
                temp[i] = rhs_mat(i,col_idx);
                for (size_type k = i; k >= 1; --k) 
                {
                    temp[i] -= LD_(i,i-k) * temp[i-k];
                }
            }

            /*   Solve for final column vector solution   */
            for(int k=row_-1; k>=0; k--)
            {
                x_mat.set(k,col_idx, temp[k] / LD_(k,k));
                for(size_type j=k+1; j<row_; ++j)
                {
                    x_mat.set(k,col_idx, x_mat(k,col_idx)-LD_(j,k) * x_mat(j,col_idx) );
                }
            }

         }

        return true;
    }

} // end namespace yakutat
#endif