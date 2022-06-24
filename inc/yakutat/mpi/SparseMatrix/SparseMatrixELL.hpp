#pragma once
#ifndef __MPI_SPARSEMATRIXELL_H__
#define	__MPI_SPARSEMATRIXELL_H__

#include <vector>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <mpi.h>
#include <omp.h>
#include <math.h> 

#include <yakutat/backend/mpi_tools.hpp>
#include <yakutat/backend/Matrix_traits.hpp>
#include <yakutat/mpi/DynamicMatrix.hpp>

namespace yakutat::mpi
{
    using namespace yakutat::backend;

    // T: vaule type; MemorySize: size of data length each row
    template<typename T, size_t MemorySize>
    class SparseMatrixELL
    {
        // Member types
        using self_type = SparseMatrixELL<T, MemorySize>;

            
        public:

            // Member types
            using MatrixType = DynamicMatrix<T>;
            using size_type = typename MatrixTraits<MatrixType>::size_type;
            using value_type = T;
            using array_size_type = typename MatrixTraits<MatrixType>::array_size_type;
            using array_value_type = std::vector<value_type>;
            using data_type = std::tuple<std::vector<int>&, array_value_type&>;
        
        // --- Constructor & Destructor --- matrix partition ---
        //--------------------------------------------- 
            SparseMatrixELL();
            SparseMatrixELL(size_type rows, size_type cols); 
            SparseMatrixELL(size_type n); 
            SparseMatrixELL(const self_type & matrix); 
            virtual ~SparseMatrixELL();
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
            size_type max_row_nnz(void) const noexcept;
            value_type at(size_type row, size_type col) const;  



        // --- Add element ---
            self_type & set(size_type row, size_type col, value_type val); 
        // --- Accumulate element ---
            self_type & accumulate(size_type row, size_type col, value_type val); 
        // --- Clear contents ---
            bool clear(void);       
        // --- fill zero ---
            bool zero(void);     
        // --- Resize matrix ---
            bool resize(size_type rows, size_type cols); 
            bool resize(size_type n); 
        // --- Output matrix data ---
            data_type data() noexcept;
        // --- Show matrix profile ---
            bool show_DataSturcture(void);

        //--------------------------------------------- // matrix operators
            // --- mat1 + mat2 ---
            self_type add(const self_type & matrix) const; 
            // --- mat1 - mat2 ---
            self_type subtract(const self_type & matrix) const; 
            // --- mat * vector ---  
            array_value_type  multiply(const array_value_type & x) const; 
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
            template<typename X, size_t Size>
				friend std::ostream & operator << (std::ostream & os, const SparseMatrixELL<X, Size> & matrix);
        //---------------------------------------------

        //--------------------------------------------- // MPI
            self_type & MPI_allgather(size_type itag);
        //---------------------------------------------



        private:
            size_type row_{0}, col_{0},  max_row_nnz_{0}; 
            
           // ell format
            std::vector<int> ell_col_indexes_;
            array_value_type ell_values_;

            //MPI variables
            mpi_tools world_;

    

        //--------------------------------------------- // Basic function
            bool insert(size_type row, size_type col, T val);
            bool remove(size_type row, size_type col);
            bool copy(const self_type & matrix);
        //---------------------------------------------//

    };


    //##################################################//
    //                Implementations                   //
    //##################################################//



    //--------------------------------------------- // Constructor & Destructor --- matrix partition 
    template<typename T, size_t MemorySize>
	inline SparseMatrixELL<T, MemorySize>::SparseMatrixELL() {}

    template<typename T, size_t MemorySize>
	inline SparseMatrixELL<T, MemorySize>::SparseMatrixELL(size_type rows, size_type cols)
    : row_(rows), col_(cols), max_row_nnz_(MemorySize)
	{
        ell_col_indexes_.resize(max_row_nnz_*rows, -1);
        ell_values_.resize(max_row_nnz_*rows, 0);
	}

    template<typename T, size_t MemorySize>
    inline SparseMatrixELL<T, MemorySize>::SparseMatrixELL(size_type n)
    : row_(n), col_(n), max_row_nnz_(MemorySize)
    {
        ell_col_indexes_.resize(max_row_nnz_*n, -1);
        ell_values_.resize(max_row_nnz_*n, 0);

    }

    template<typename T, size_t MemorySize>
    inline bool SparseMatrixELL<T, MemorySize>::partition(size_type myid, size_type nproc)
    {
        if(row_==0)
        {
            throw std::invalid_argument("SparseMatrixELL partition: Your row of matrix is zero.");
            return false;
        }
        
        world_.MPI_division(row_, myid, nproc);
        return true;
    }

    template<typename T, size_t MemorySize>
    inline bool SparseMatrixELL<T, MemorySize>::partition(mpi_tools & world)
    {
        if(world.length==0)
        {
            throw std::invalid_argument("SparseMatrixELL partition: Your world no division.");
            return false;
        }

        world_ = world;
        return true;
    }


    template<typename T, size_t MemorySize>
    inline SparseMatrixELL<T, MemorySize>::~SparseMatrixELL()  {}
    //--------------------------------------------- //

    
    
    //--------------------------------------------- //  Basic matrix tools
    template<typename T, size_t MemorySize>
	inline typename SparseMatrixELL<T, MemorySize>::self_type & SparseMatrixELL<T, MemorySize>::set(size_type row, size_type col, value_type val)
	{

        if (row < 0 || col < 0 || row > row_-1 || col > col_-1) 
        {
			throw std::invalid_argument("SparseMatrixELL set: Coordination out of range.");
		}

        if( !(val==0) )
        {
            if( at(row, col)==0 )
            {
                insert(row, col, val); 
            }
            else
            {
                size_type row_idx{0};
                row_idx = row * max_row_nnz_;
                
                for(size_type i=row_idx; i<row_idx+max_row_nnz_; ++i)
                {
                    if(ell_col_indexes_[i] == (int)col)
                    {
                        ell_values_[i] = val;
                        break;
                    }
                }
        
            }
        }
        else if( !(at(row, col)==0) && val==0 )
        {
            remove(row, col);
        }
		return *this;
	}



    template<typename T, size_t MemorySize>
	inline typename SparseMatrixELL<T, MemorySize>::value_type SparseMatrixELL<T, MemorySize>::at(size_type row, size_type col) const
	{
        if (row < 0 || col < 0 || row > row_-1 || col > col_-1) 
        {
			throw std::invalid_argument("SparseMatrixELL at: Coordination out of range.");
		}

        size_type row_idx = row * max_row_nnz_;
        value_type result{0.0};
    

        for(size_type i=row_idx; i<row_idx+max_row_nnz_; ++i)
        {
            if(ell_col_indexes_[i] == (int)col)
            {
                result = ell_values_[i];
                break;
            }
        }
        return result;
	}

    template<typename T, size_t MemorySize>
	inline typename SparseMatrixELL<T, MemorySize>::self_type & SparseMatrixELL<T, MemorySize>::accumulate(size_type row, size_type col, value_type val)
    {
        if (row < 0 || col < 0 || row > row_-1 || col > col_-1) 
        {
			throw std::invalid_argument("SparseMatrixELL accumulate: Coordination out of range.");
		}

        set(row, col, at(row, col) + val);
        return *this;

    }


    template<typename T, size_t MemorySize>
    inline bool SparseMatrixELL<T, MemorySize>::clear(void)
    {
        std::fill(ell_values_.begin(), ell_values_.end(),0);
        std::fill(ell_col_indexes_.begin(), ell_col_indexes_.end(),-1);

        ell_values_.clear();
        ell_col_indexes_.clear();
        
        row_ = 0;
        col_ = 0;
        return true;
    }

    template<typename T, size_t MemorySize>
    inline bool SparseMatrixELL<T, MemorySize>::zero(void)
    {
        std::fill(ell_values_.begin(), ell_values_.end(),0);
        std::fill(ell_col_indexes_.begin(), ell_col_indexes_.end(),-1);
        return true;
    }

    template<typename T, size_t MemorySize>
    inline bool SparseMatrixELL<T, MemorySize>::resize(size_type rows, size_type cols)
    {
        row_ = rows;
        col_ = cols;

        ell_col_indexes_.resize(max_row_nnz_*row_, -1);
        ell_values_.resize(max_row_nnz_*row_, 0);
        return true;
    }

    template<typename T, size_t MemorySize>
    inline bool SparseMatrixELL<T, MemorySize>::resize(size_type n)
    {
        row_ = n;
        col_ = n;

        ell_col_indexes_.resize(max_row_nnz_*row_, -1);
        ell_values_.resize(max_row_nnz_*row_, 0);
        return true;
    }

    
    template<typename T, size_t MemorySize>
	inline SparseMatrixELL<T, MemorySize>::SparseMatrixELL(const self_type & matrix)
	{
		copy(matrix);
	}
    template<typename T, size_t MemorySize>
    inline typename SparseMatrixELL<T, MemorySize>::data_type SparseMatrixELL<T, MemorySize>::data() noexcept
    {
        return {ell_col_indexes_, ell_values_};
    }

    
    template<typename T, size_t MemorySize>
    inline bool SparseMatrixELL<T, MemorySize>::show_DataSturcture(void)
    {

        size_type row_idx{0};
        std::cout << "ELL col_indexes" << std::endl;
        for(size_type i=0; i<row_; ++i)
        {
            row_idx = max_row_nnz_*i;
            for(size_type j=row_idx; j<row_idx+max_row_nnz_; ++j)
            {
                std::cout << std::setw(4) << ell_col_indexes_[j] << " " ;
            }
            std::cout << std::endl;
        }

        row_idx = 0;
        std::cout << "ELL value" << std::endl;
        for(size_type i=0; i<col_; ++i)
        {
            row_idx = max_row_nnz_*i;
            for(size_type j=row_idx; j<row_idx+max_row_nnz_; ++j)
            {
                std::cout << std::setw(4) << ell_values_[j] << " " ;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;


        return true;
    }
    


    template<typename T, size_t MemorySize>
    inline typename SparseMatrixELL<T, MemorySize>::size_type SparseMatrixELL<T, MemorySize>::row(void) const noexcept
    {
        return row_;
    } 

    template<typename T, size_t MemorySize>
    inline typename SparseMatrixELL<T, MemorySize>::size_type SparseMatrixELL<T, MemorySize>::col(void) const noexcept
    {
        return col_;
    }


    template<typename T, size_t MemorySize>
    inline typename SparseMatrixELL<T, MemorySize>::size_type SparseMatrixELL<T, MemorySize>::max_row_nnz(void) const noexcept
    {
        return max_row_nnz_;
    }

    template<typename T, size_t MemorySize>
    inline typename SparseMatrixELL<T, MemorySize>::self_type SparseMatrixELL<T, MemorySize>::add(const self_type & matrix) const
    {
        if(row_ != matrix.row_ || col_ != matrix.col_)
        {
            throw std::invalid_argument("SparseMatrixELL Cannot add: szie of matrices don't match.");
        }

        self_type result(row_, col_);
        

        for(size_type i=world_.start; i<world_.end+1; ++i)
        {
            for(size_type j=0; j<col_; ++j)
            {
               if(at(i,j)!=0 || matrix(i,j)!=0)
                    result.set(i, j, at(i,j) + matrix(i,j) );


            }
        }
        
        return result;
    }

    template<typename T, size_t MemorySize>
    inline typename SparseMatrixELL<T, MemorySize>::self_type SparseMatrixELL<T, MemorySize>::subtract(const self_type & matrix) const
    {
        if(row_ != matrix.row_ || col_ != matrix.col_)
        {
            throw std::invalid_argument("SparseMatrixELL Cannot subtract: szie of matrices don't match.");
        }

        self_type result(row_, col_);
        

        for(size_type i=world_.start; i<world_.end+1; ++i)
        {
            for(size_type j=0; j<col_; ++j)
            {
               if(at(i,j)!=0 || matrix(i,j)!=0)
                    result.set(i, j, at(i,j) - matrix(i,j) );
            }
        }
        
        return result;
    }

    template<typename T, size_t MemorySize> 
    inline typename SparseMatrixELL<T, MemorySize>::array_value_type SparseMatrixELL<T, MemorySize>::multiply(const array_value_type & x) const
    {
        if( this->row_ != x.size() )
        {
			throw std::invalid_argument("SparseMatrixELL col of matrix != vector.size()");
        }


        array_value_type result(row_, 0.0);
        value_type sum{0.0};
        size_type row_idx{0};

       
        
        for(size_type i=world_.start; i<world_.end+1; ++i)
        {
            sum = 0;
            row_idx = i * max_row_nnz_;

            #pragma omp simd reduction(+:sum)
            for(size_type j=row_idx; j<row_idx+max_row_nnz_; ++j)
            {
                if(ell_col_indexes_[j] >=0)
                    sum += ell_values_[j] * x[ ell_col_indexes_[j] ];
            }

            result[i] = sum;
        }
        return result;
    }

    template<typename T, size_t MemorySize>
    inline typename SparseMatrixELL<T, MemorySize>::self_type SparseMatrixELL<T, MemorySize>::multiply(const self_type & matrix) const
    {

        if( col_ != matrix.row_ )
        {
			throw std::invalid_argument("SparseMatrixELL Cannot multiply: col of left matrix != row of right matrix");
        }
        

        value_type sum{0.0};
        size_type row_idx{0};

        self_type result(row_, matrix.col_);
        
        for(size_type i=world_.start; i<world_.end+1; ++i)
        {
           for(size_type j=0; j<matrix.col_; ++j)
           {

                sum = 0.0;
                row_idx = i * max_row_nnz_;

                #pragma omp simd reduction(+:sum)
                for(size_type k=row_idx; k<row_idx+max_row_nnz_; ++k)
                {
                    if(ell_col_indexes_[k] >=0)
                        sum += ell_values_[k] * matrix(ell_col_indexes_[k],j);
                }
                result.set(i, j, sum);
           }
        }
        return result;
    }
    //--------------------------------------------- //



    //--------------------------------------------- // operations

    template<typename T, size_t MemorySize>
    inline typename SparseMatrixELL<T, MemorySize>::self_type & SparseMatrixELL<T, MemorySize>::operator = (const self_type & matrix)
    {   
        if(&matrix != this)
        {
            copy(matrix);
        }
        return *this;
    }

    template<typename T, size_t MemorySize>
	inline typename SparseMatrixELL<T, MemorySize>::self_type SparseMatrixELL<T, MemorySize>::operator + (const self_type & matrix) const
	{
		return add(matrix);
	}

    template<typename T, size_t MemorySize>
    inline typename SparseMatrixELL<T, MemorySize>::self_type SparseMatrixELL<T, MemorySize>::operator - (const self_type & matrix) const
    {
        return subtract(matrix);
    }

    template<typename T, size_t MemorySize>
    inline typename SparseMatrixELL<T, MemorySize>::array_value_type SparseMatrixELL<T, MemorySize>::operator * (const array_value_type & x) const
    {
        return multiply(x);
    }
    template<typename T, size_t MemorySize>
    inline typename SparseMatrixELL<T, MemorySize>::self_type SparseMatrixELL<T, MemorySize>::operator * (const self_type & matrix) const
    {
        return multiply(matrix);
    }

    template<typename T, size_t MemorySize>
    inline typename SparseMatrixELL<T, MemorySize>::value_type SparseMatrixELL<T, MemorySize>::operator () (const size_type & row, const size_type & col) const
    {
        return at(row, col);
    }
    //--------------------------------------------- //





    //--------------------------------------------- //  friend founction

    template<typename T, size_t MemorySize> 
	inline std::ostream & operator << (std::ostream & os, const SparseMatrixELL<T, MemorySize> & matrix)
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


    //--------------------------------------------- //




    //--------------------------------------------- // MPI
    template<typename T, size_t MemorySize>
    inline typename SparseMatrixELL<T, MemorySize>::self_type & SparseMatrixELL<T, MemorySize>::MPI_allgather(size_type itag)
    {
        MPI_Status istat[8];
        size_type master{0};
        size_type istart{0};
        size_type icount{0};  
        size_type row_idx{0};

        
        icount = world_.count; // row cut
        istart = world_.start;

        row_idx = istart * max_row_nnz_;
        icount *= max_row_nnz_;
        if(world_.myid > master)
        {
            MPI_Send( (void *)&ell_values_[row_idx], icount, MPI_DOUBLE, master, itag, MPI_COMM_WORLD );
            MPI_Send( (void *)&ell_col_indexes_[row_idx], icount, MPI_INT, master, itag, MPI_COMM_WORLD );
        }
        else if(world_.myid == master)
        {
            for(size_type i = 1 ; i < world_.nproc; ++i)
            {
                MPI_Recv( (void *)&ell_values_[ world_.start_list[i]*max_row_nnz_ ], world_.count_list[i]*max_row_nnz_ , MPI_DOUBLE, i, itag, MPI_COMM_WORLD, istat );
                MPI_Recv( (void *)&ell_col_indexes_[ world_.start_list[i]*max_row_nnz_ ], world_.count_list[i]*max_row_nnz_ , MPI_INT, i, itag, MPI_COMM_WORLD, istat );
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        icount = row_ * max_row_nnz_;
        MPI_Bcast(  (void *)&ell_values_[0], icount, MPI_DOUBLE, master, MPI_COMM_WORLD   );
        MPI_Bcast(  (void *)&ell_col_indexes_[0], icount, MPI_INT, master, MPI_COMM_WORLD   );

        


        return *this;

  
    }
    //---------------------------------------------

    //--------------------------------------------- // Basic function for sparse matrix

   

    template<typename T, size_t MemorySize>
    inline bool SparseMatrixELL<T, MemorySize>::insert(size_type row, size_type col, T val)
    {
        size_type row_idx = row * max_row_nnz_;
        for(size_type i=row_idx; i<row_idx+max_row_nnz_; ++i)
        {
            
            if(ell_col_indexes_[i] < 0)
            {
                ell_col_indexes_[i] = col;
                ell_values_[i] = val;
                return true;
            }
        }
        return false;
    }



    template<typename T, size_t MemorySize>
	inline bool SparseMatrixELL<T, MemorySize>::remove(size_type row, size_type col)
	{
        size_type row_idx = row * max_row_nnz_;
        for(size_type i=row_idx; i<row_idx+max_row_nnz_; ++i)
        {
            if(ell_col_indexes_[i] == (int)col)
            {
                ell_col_indexes_[i] = -1;
                ell_values_[i] = 0;
                return true;
            }
        }
        return false;
		
	}

    template<typename T, size_t MemorySize>
    inline bool SparseMatrixELL<T, MemorySize>::copy(const self_type & matrix)
    {
        row_ = matrix.row_; 
        col_ = matrix.col_; 
        max_row_nnz_ = MemorySize;

        ell_col_indexes_ = matrix.ell_col_indexes_;
        ell_values_ = matrix.ell_values_;
        return true;
    }

    
    //--------------------------------------------- //


} // end namespace yakutat::mpi


#endif
