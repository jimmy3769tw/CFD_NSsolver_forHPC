#ifndef __MATRIX_OPERATIONS_H__
#define	__MATRIX_OPERATIONS_H__
#include <yakutat/matrix_definition.hpp>
#include <yakutat/mpi/matrix_definition.hpp>
#include <yakutat/backend/Matrix_traits.hpp>



namespace yakutat::backend
{
    template<typename MatrixType>
    class MatOps
        : private MatrixType
    {
        // Member types
        using self_type = MatOps<MatrixType>;

        public:
            // Member types
            using size_type = typename MatrixTraits<MatrixType>::size_type;
            using value_type = typename MatrixTraits<MatrixType>::value_type;
            using array_size_type = typename MatrixTraits<MatrixType>::array_size_type;
            using array_value_type = typename MatrixTraits<MatrixType>::array_value_type;


            MatOps();
            virtual ~MatOps();


            // --- Set row values to zero ---
            bool clearRowVal(MatrixType & matrix, size_type row_idx);
            // --- Set col values to zero ---
            bool clearColVal(MatrixType & matrix, size_type col_idx);
            // --- Get the diagonal value of sparse matrix ---
            value_type getDiaValue(MatrixType & matrix, size_type dia_idx);
            // --- Set the diagonal value of sparse matrix ---
            bool setDiaValue(MatrixType & matrix, size_type dia_idx, value_type val); 
            // --- Get diagonal matrix of sparse matrix  ---
            MatrixType diagonal(MatrixType & matrix);
            // --- Transpose the sparse matrix  ---
            MatrixType transpose(MatrixType & matrix);
            // --- Get the max element of sparse matrix  ---
            value_type maxCoeff(MatrixType & matrix);
            // --- Get the min element of sparse matrix  ---
            value_type minCoeff(MatrixType & matrix); 


        private:

    };

    //##################################################//
    //                Implementations                   //
    //##################################################//

    template<typename MatrixType>
	inline MatOps<MatrixType>::MatOps() {}


    template<typename MatrixType>
    inline MatOps<MatrixType>::~MatOps() {}

    template<typename MatrixType>
    inline bool 
    MatOps<MatrixType>::clearRowVal(MatrixType & matrix, size_type row_idx)
    {
        for(size_type j=0; j<matrix.col(); ++j)
            matrix.set(row_idx, j, 0);
        return true;
    }

    template<typename MatrixType>
    inline bool 
    MatOps<MatrixType>::clearColVal(MatrixType & matrix, size_type col_idx)
    {
        for(size_type i=0; i<matrix.row(); ++i)
            matrix.set(i,col_idx, 0);
        return true;
    }

    template<typename MatrixType>
    inline typename MatOps<MatrixType>::value_type 
    MatOps<MatrixType>::getDiaValue(MatrixType & matrix, size_type dia_idx)
    {
        return matrix(dia_idx, dia_idx);
    }

    template<typename MatrixType>
    inline bool 
    MatOps<MatrixType>::setDiaValue(MatrixType & matrix, size_type dia_idx, value_type val)
    {
        matrix.set(dia_idx, dia_idx, val);
        return true;
    }

    template<typename MatrixType>
    inline MatrixType
    MatOps<MatrixType>::diagonal(MatrixType & matrix)
    {
        size_type row, col;
        row = matrix.row();
        col = matrix.col();
        MatrixType result(row, col);

        #pragma omp simd
        for(size_type i=0; i<row; ++i)
            result.set(i, i, matrix(i,i));
        return result;
    }

    template<typename MatrixType>
    inline MatrixType 
    MatOps<MatrixType>::transpose(MatrixType & matrix)
    {
        size_type row, col;
        row = matrix.row();
        col = matrix.col();
        MatrixType result(col, row);
        
        for(size_type i=0; i<row; ++i)
        {
            #pragma omp simd
            for(size_type j=0; j<col; ++j)
            {
                if(matrix(i,j)!=0)
                {
                    result.set(j, i, matrix(i,j));
                }
            }
        }
        return result;
    }

    template<typename MatrixType>
    inline typename MatOps<MatrixType>::value_type 
    MatOps<MatrixType>::maxCoeff(MatrixType & matrix)
    {
        value_type result(matrix);
        size_type row{matrix.row()}, col{matrix.col()};
        result.zero();

        #pragma omp parallel for collapse(2)
        for(size_t i=0; i<row; ++i)
        {
            for(size_t j=0; j<col; ++j)
            {
                if(result < matrix(i,j))
                    result = matrix(i,j);
            }
        }
        return result;

    }

    template<typename MatrixType>
    inline typename MatOps<MatrixType>::value_type 
    MatOps<MatrixType>::minCoeff(MatrixType & matrix)
    {
        value_type result(matrix);
        size_type row{matrix.row()}, col{matrix.col()};
        result.zero();
        #pragma omp parallel for collapse(2)
        for(size_t i=0; i<row; ++i)
        {
            for(size_t j=0; j<col; ++j)
            {
                if(result > matrix(i,j))
                    result = matrix(i,j);
            }
        }
        return result;
    }
} // end namespace yakutat::backend

#endif