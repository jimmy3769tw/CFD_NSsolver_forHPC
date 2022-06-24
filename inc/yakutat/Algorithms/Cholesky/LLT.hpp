#pragma once
#ifndef __LLT_H__
#define	__LLT_H__

#include <yakutat/matrix_definition.hpp>
#include <yakutat/mpi/matrix_definition.hpp>
#include <yakutat/backend/Matrix_traits.hpp>

namespace yakutat
{
    using namespace yakutat::backend;

    template<typename MatrixType>
    class LLT
        : private MatrixType
    {
        // Member types
        using self_type = LLT<MatrixType>;

        public:
            // Member types
            using size_type = typename MatrixTraits<MatrixType>::size_type;
            using value_type = typename MatrixTraits<MatrixType>::value_type;
            using array_size_type = typename MatrixTraits<MatrixType>::array_size_type;
            using array_value_type = typename MatrixTraits<MatrixType>::array_value_type;
            using matops_type = MatOps<MatrixType>;
            
           
            LLT();
            virtual ~LLT();

            bool factorize(const MatrixType & matrix);
            MatrixType & matrixL(void);
            MatrixType & matrixLT(void);
            bool solve(const array_value_type & rhs, array_value_type & x);
            bool solve(const MatrixType & rhs_mat, MatrixType & x_mat);

        private:
            MatrixType L_, LT_;
            size_type row_{0}, col_{0};
            matops_type matops_;

    };  

    //##################################################//
    //                Implementations                   //
    //##################################################//


    template<typename MatrixType>
	inline LLT<MatrixType>::LLT() {}

    template<typename MatrixType>
    inline LLT<MatrixType>::~LLT() {}

    template<typename MatrixType>
    inline bool LLT<MatrixType>::factorize(const MatrixType & matrix)
    {
        if (matrix.row() != matrix.col()) 
        {
			throw std::invalid_argument("LLT factorize: matrix.row() != matrix.col()");
		}
        value_type sum{0.0};
        row_ = matrix.row();
        col_ = matrix.col();
        L_ = matrix; // todo
        

        for(size_type j=0; j<row_; ++j)
        {
            sum = 0.0;
            for(size_type k=0; k<j; ++k)
            {
                sum += L_(j,k) * L_(j,k);
            } 
            L_.set(j,j, sqrt(matrix(j,j)- sum)  );

            
            for(size_type i=j+1; i<row_; ++i)
            {
                sum = 0.0;
                for(size_type k=0; k<j; ++k)
                {
                    sum += L_(i,k) * L_(j,k);
                } 
                
                L_.set(i,j, (1.0/L_(j,j)) * (matrix(i,j)-sum)  );
            }
            
        }
        return true;
    }

    template<typename MatrixType>
    inline MatrixType & LLT<MatrixType>::matrixLT(void)
    {
        LT_ = matops_.transpose(L_);
        return LT_;
    }

    template<typename MatrixType>
    inline MatrixType & LLT<MatrixType>::matrixL(void)
    {
        return L_;
    }

    template<typename MatrixType>
    bool LLT<MatrixType>::solve(const array_value_type & rhs, array_value_type & x)
    {
        array_value_type temp(rhs);

        for(size_type i=0; i< row_; ++i)
        {
            x[i] = rhs[i];
            for(size_type j=0; j< i; ++j)
            {
                x[i] -= L_(i,j) * x[j];
            }
            x[i] /= L_(i,i);
            temp[i] = x[i];
        }

        for(int i=int(row_-1); i>=0; --i)
        {
            x[i] = temp[i];
            for(int j=i+1; j< int(col_); ++j)
            {
                x[i] -= L_(j,i) * x[j];
            }
            x[i] /= L_(i,i);
        }
        return true;
    }


    template<typename MatrixType>
    bool 
    LLT<MatrixType>::solve(const MatrixType & rhs_mat, MatrixType & x_mat)
    {

        array_value_type temp(row_);

        for(size_type col_idx=0; col_idx<col_; ++col_idx)
        {

            for(size_type i=0; i< row_; ++i)
            {
                x_mat.set(i,col_idx, rhs_mat(i,col_idx));
                for(size_type j=0; j< i; ++j)
                {
                    x_mat.set(i,col_idx, x_mat(i,col_idx)-L_(i,j) * x_mat(j,col_idx) );
                }
                x_mat.set(i,col_idx, x_mat(i,col_idx)/L_(i,i));
                temp[i] = x_mat(i,col_idx);
            }

            for(int i=int(row_-1); i>=0; --i)
            {
                x_mat.set(i,col_idx, temp[i]);
                for(int j=i+1; j< int(col_); ++j)
                {
                    x_mat.set(i,col_idx, x_mat(i,col_idx)-L_(j,i) * x_mat(j,col_idx) );
                }
                x_mat.set(i,col_idx, x_mat(i,col_idx)/L_(i,i));
            }


        }
        return true;
    }

} // end namespace yakutat
#endif
