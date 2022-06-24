#pragma once
#ifndef __MPI_DYNAMICMATRIX_H__
#define	__MPI_DYNAMICMATRIX_H__

#include <vector>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <mpi.h>
#include <omp.h>
#include <math.h> 
#include <algorithm>

#include <yakutat/backend/mpi_tools.hpp>
namespace yakutat::mpi
{
    using namespace yakutat::backend;


    template<typename T>
    class DynamicMatrix
    {
        // Member types
        using self_type = DynamicMatrix<T>;
        using base_type = std::vector<T>;

        public:
            
            // Member types
            using size_type = typename base_type::size_type;
            using value_type = T;
            using array_size_type = std::vector<size_type>;
            using array_value_type = std::vector<value_type>;
            using data_type = base_type;

        // --- Constructor & Destructor --- matrix partition ---
        //--------------------------------------------- 
            DynamicMatrix();
            DynamicMatrix(size_type rows, size_type cols, value_type val);
            DynamicMatrix(size_type rows, size_type cols);
            DynamicMatrix(size_type n);
            DynamicMatrix(const self_type & matrix);
            virtual ~DynamicMatrix();
            bool partition(size_type myid, size_type nproc);
            bool partition(mpi_tools & world);
        //--------------------------------------------- 

        //##################################################//
        //               Basic matrix tools                 //
        //##################################################//
        //--------------------------------------------- 
        // --- Get matrix info ---
            size_type row(void) const noexcept;
            size_type col(void) const noexcept;
            value_type at(size_type row, size_type col) const; 

        // --- Add element ---
            self_type & set(size_type row, size_type col, value_type val); 
        // --- Accumulate element ---
            self_type & accumulate(size_type row, size_type col, value_type val); 
        // --- Clear contents ---
            bool clear(void);       
        // --- Resize matrix ---
            bool resize(size_type rows, size_type cols, value_type val); 
            bool resize(size_type rows, size_type cols); 
            bool resize(size_type n);
        // --- Output matrix data ---
            data_type & data() noexcept;
        // --- Fill matrix elements ---
            bool fill(value_type val);
        // --- Test whether matrix is empty ---
            bool empty(void);

        //--------------------------------------------- // matrix operators
            // --- mat1 + mat2 ---
            self_type add(const self_type & matrix) const; 
            // --- mat1 - mat2 ---
            self_type subtract(const self_type & matrix) const; 
            // --- mat * vector --- 
            array_value_type multiply(const array_value_type & x) const; 
            // --- mat1 * mat2 --- 
            self_type multiply(const self_type & matrix) const;  
        //--------------------------------------------- 

        //--------------------------------------------- // operations
            self_type & operator = (const self_type & matrix);
            self_type operator + (const self_type & matrix) const;
            self_type operator - (const self_type & matrix) const;
            array_value_type operator * (const array_value_type & x) const;
            self_type operator * (const self_type & matrix) const;
            value_type operator () (const size_type & row, const size_type & col) const;
        //---------------------------------------------

        //--------------------------------------------- // friend function
            template<typename X>
				friend std::ostream & operator << (std::ostream & os, const DynamicMatrix<X> & matrix);
        //---------------------------------------------

        //--------------------------------------------- // MPI
            self_type & MPI_allgather(size_type itag);
        //---------------------------------------------

        private:
            size_type row_{0};
            size_type col_{0}; 
            array_value_type values_;
            mpi_tools world_;

            //--------------------------------------------- // Basic function
            bool copy(const self_type & matrix);
            //---------------------------------------------//

    };

    //--------------------------------------------- // Constructor & Destructor --- matrix partition 

    template<typename T>
	inline DynamicMatrix<T>::DynamicMatrix() {}

    template<typename T>
	inline DynamicMatrix<T>::DynamicMatrix(size_type rows, size_type cols, value_type val)
    : row_(rows), col_(cols)
    {
        values_.resize(row_*col_, val);
    }

    template<typename T>
	inline DynamicMatrix<T>::DynamicMatrix(size_type rows, size_type cols)
    : row_(rows), col_(cols)
    {
        values_.resize(row_*col_);
    }

    template<typename T>
	inline DynamicMatrix<T>::DynamicMatrix(size_type n)
    : row_(n), col_(n)
    {
        values_.resize(row_*col_);
    }

    template<typename T>
	inline DynamicMatrix<T>::DynamicMatrix(const self_type & matrix)
    {
        copy(matrix);
    }

    template<typename T>
    inline bool DynamicMatrix<T>::partition(size_type myid, size_type nproc)
    {
        if(row_==0)
        {
            throw std::invalid_argument("DynamicMatrix partition: Your row of matrix is zero.");
            return false;
        }
        world_.MPI_division(row_, myid, nproc);
        return true;
    }

    template<typename T>
    inline bool DynamicMatrix<T>::partition(mpi_tools & world)
    {
        if(world.length==0)
        {
            throw std::invalid_argument("DynamicMatrix partition: Your world no division.");
            return false;
        }

        world_ = world;
        return true;
    }

    template<typename T>
    inline DynamicMatrix<T>::~DynamicMatrix() {}
    //--------------------------------------------- //

    //--------------------------------------------- //  Basic matrix tools
    template<typename T>
    inline DynamicMatrix<T> & DynamicMatrix<T>::set(size_type row, size_type col, value_type val)
    {
        if (row < 0 || col < 0 || row > row_-1 || col > col_-1) 
        {
			throw std::invalid_argument("DynamicMatrix set: Coordination out of range.");
		}

        values_[row * col_ + col] = val;
        return *this;
    }

    template<typename T>
    inline typename DynamicMatrix<T>::value_type DynamicMatrix<T>::at(size_type row, size_type col) const
    {
        if (row < 0 || col < 0 || row > row_-1 || col > col_-1) 
        {
			throw std::invalid_argument("DynamicMatrix at: Coordination out of range.");
		}

        return values_[row * col_ + col];
    }

    template<typename T>
    inline typename DynamicMatrix<T>::size_type DynamicMatrix<T>::row(void) const noexcept
    {
        return row_;
    }

    template<typename T>
    inline typename DynamicMatrix<T>::size_type DynamicMatrix<T>::col(void) const noexcept
    {
        return col_;
    }

    template<typename T>
    inline typename DynamicMatrix<T>::self_type & DynamicMatrix<T>::accumulate(size_type row, size_type col, value_type val)
    {
        if (row < 0 || col < 0 || row > row_-1 || col > col_-1) 
        {
			throw std::invalid_argument("DynamicMatrix accumulate: Coordination out of range.");
		}

        values_[row * col_ + col] += val;
        return *this;
    }

    template<typename T>
    inline bool DynamicMatrix<T>::clear(void)
    {
        std::fill(values_.begin(), values_.end(), 0);
        values_.clear();
        row_ = 0;
        col_ = 0;
        return true;
    }

    template<typename T>
    inline bool DynamicMatrix<T>::resize(size_type rows, size_type cols, value_type val)
    {
        values_.resize(rows * cols, val);
        row_ = rows;
        col_ = cols;
        return true;
    }

    template<typename T>
    inline bool DynamicMatrix<T>::resize(size_type rows, size_type cols)
    {
        values_.resize(rows * cols);
        row_ = rows;
        col_ = cols;
        return true;
    }

    template<typename T>
    inline bool DynamicMatrix<T>::resize(size_type n)
    {
        values_.resize(n *n);
        row_ = n;
        col_ = n;
        return true;
    }

    template<typename T>
    inline typename DynamicMatrix<T>::data_type & DynamicMatrix<T>::data() noexcept
    {
        return values_;
    }

    template<typename T>
    inline bool DynamicMatrix<T>::fill(value_type val)
    {
        std::fill(values_.begin(), values_.end(), val);
        return true;
    }

    template<typename T>
    inline bool DynamicMatrix<T>::empty(void)
    {
        return values_.empty();
    }
    //--------------------------------------------- //

    //--------------------------------------------- // matrix operators
    template<typename T>
    inline typename DynamicMatrix<T>::self_type DynamicMatrix<T>::add(const self_type & matrix) const
    {
        if(row_ != matrix.row_ || col_ != matrix.col_)
        {
            throw std::invalid_argument("DynamicMatrix Cannot add: szie of matrices don't match.");
        }
        DynamicMatrix<value_type> result(row_, col_);

        for(size_type i=world_.start; i<world_.end+1; ++i)
        {
            for(size_type j=0; j<col_; ++j)
            {
                result.values_[i*col_+j] = values_[i*col_+j] + matrix.values_[i*col_+j];
            }
        }
        return result;
    }

    template<typename T>
    inline typename DynamicMatrix<T>::self_type DynamicMatrix<T>::subtract(const self_type & matrix) const
    {
        if(row_ != matrix.row_ || col_ != matrix.col_)
        {
            throw std::invalid_argument("DynamicMatrix Cannot add: szie of matrices don't match.");
        }
        DynamicMatrix<value_type> result(row_, col_);

        for(size_type i=world_.start; i<world_.end+1; ++i)
        {
            for(size_type j=0; j<col_; ++j)
            {
                result.values_[i*col_+j] = values_[i*col_+j] - matrix.values_[i*col_+j];
            }
        }
        return result;
    }

    template<typename T>
    inline typename DynamicMatrix<T>::array_value_type DynamicMatrix<T>::multiply(const array_value_type & x) const
    {
        if( this->row_ != x.size() )
        {
			throw std::invalid_argument("DynamicMatrix col of matrix != vector.size()");
        }

        array_value_type result(row_, 0.0);
        size_type row_idx{0}, vec_idx{0};
        value_type sum{0.0};

        for(size_type i=world_.start; i<world_.end+1; ++i)
        {
            sum = 0.0;
            row_idx = i * col_;
            vec_idx = 0;
            #pragma omp simd reduction(+:sum)
            for(size_type j=row_idx; j<row_idx+col_; ++j)
            {
                if(values_[j]!=0 || x[vec_idx]!=0)
                sum += values_[j] * x[vec_idx];
                ++vec_idx;
            }
            result[i] = sum;
        }
        return result;
    }

    template<typename T>
    inline typename DynamicMatrix<T>::self_type DynamicMatrix<T>::multiply(const self_type & matrix) const
    {
        if( col_ != matrix.row_ )
        {
			throw std::invalid_argument("SparseMatrixELL Cannot multiply: col of left matrix != row of right matrix");
        }
        self_type result(row_, matrix.col_);
        value_type sum{0.0};

        for(size_type i=world_.start; i<world_.end+1; ++i)
        {
           for(size_type j=0; j<matrix.col_; ++j)
           {
                sum = 0;
                #pragma omp simd reduction(+:sum)
                for(size_type k=0; k<matrix.col_; ++k)
                {
                    if(values_[i * col_ + k]!=0 || matrix.values_[k * matrix.col_ + j]!=0)
                        sum += values_[i * col_ + k] * matrix.values_[k * matrix.col_ + j];
                }
                result.set(i, j, sum);
           }
        }
        return result;
    }

    //--------------------------------------------- //

    //--------------------------------------------- // operators
    template<typename T>
    inline typename DynamicMatrix<T>::self_type & DynamicMatrix<T>::operator = (const self_type & matrix)
    {   
        if(&matrix != this)
        {
            copy(matrix);
        }
        return *this;
    }

    template<typename T>
    inline typename DynamicMatrix<T>::value_type DynamicMatrix<T>::operator () (const size_type & row, const size_type & col) const
    {
        return at(row, col);
    }

    template<typename T>
	inline typename DynamicMatrix<T>::self_type DynamicMatrix<T>::operator + (const self_type & matrix) const
	{
		return add(matrix);
	}

    template<typename T>
    inline typename DynamicMatrix<T>::self_type DynamicMatrix<T>::operator - (const self_type & matrix) const
    {
        return subtract(matrix);
    }

    template<typename T>
    inline typename DynamicMatrix<T>::self_type DynamicMatrix<T>::operator * (const self_type & matrix) const
    {
        return multiply(matrix);
    }

    template<typename T>
    inline typename DynamicMatrix<T>::array_value_type DynamicMatrix<T>::operator * (const array_value_type & x) const
    {
        return multiply(x);
    }
    //--------------------------------------------- //


    //--------------------------------------------- // friend function
    template<typename T> 
	inline std::ostream & operator << (std::ostream & os, const DynamicMatrix<T> & matrix)
	{
        os << std::endl;
        
		for (size_t i = 0; i < matrix.row_; ++i) 
        {
			for (size_t j = 0; j < matrix.col_; ++j) 
            {
				if (j != 0) 
                {
					os << "   ";
				}

				os << std::setprecision(10) << std::setw(10) << matrix.at(i, j);
			}

			
            os << std::endl;
		}

		return os;
	}
    //---------------------------------------------//

    //--------------------------------------------- // MPI
    template<typename T>
    inline typename DynamicMatrix<T>::self_type & DynamicMatrix<T>::MPI_allgather(size_type itag)
    {
        MPI_Status istat[8];
        size_type master{0};
        size_type istart{0};
        size_type icount{0};  
        size_type row_idx{0};

        
        icount = world_.count; // row cut
        istart = world_.start;

        row_idx = istart * col_;
        icount *= col_;
        if(world_.myid > master)
        {
            MPI_Send( (void *)&values_[row_idx], icount, MPI_DOUBLE, master, itag, MPI_COMM_WORLD );
        }
        else if(world_.myid == master)
        {
            for(size_type i = 1 ; i < world_.nproc; ++i)
            {
                MPI_Recv( (void *)&values_[ world_.start_list[i]*col_ ], world_.count_list[i]*col_ , MPI_DOUBLE, i, itag, MPI_COMM_WORLD, istat );
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        icount = row_ * col_;
        MPI_Bcast(  (void *)&values_[0], icount, MPI_DOUBLE, master, MPI_COMM_WORLD   );
        return *this;
    }
    //---------------------------------------------


    //--------------------------------------------- // Basic function
    template<typename T>
    inline bool DynamicMatrix<T>::copy(const self_type & matrix)
    {
        row_ = matrix.row();
        col_ = matrix.col();
        values_ = matrix.values_;
        return true;
    }
    //---------------------------------------------//

} // end namespace yakutat::mpi

#endif
