#pragma once 

#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <numeric>
#include <utility>  // (since C++11) std::pair
#include <algorithm> // (until C++11)
#include <tuple>

// ! row major


#ifdef ELL0_MPI_ON
    #include "MAT_mpi.hpp"
#endif 

#include <omp.h>

    // *        | col0 | col1 | col2 | 
    // *  ______|____________________|
    // * | row0 |      |      |      | 
    // * | row1 |      |      |      |
    // * | row2 |      |      |      |
    // * | row3 |      |      |      |

namespace mat
{
    template<typename T>
    class ELL_matrix0
    {

        public:
        // ---------------------
        using CSR_type = 
            // <ptr, indices, values>
            std::tuple< std::vector<int>, std::vector<int>, std::vector<T> >; 
        // ---------------------


        //! ------ Constructor & Destructor ---------
        ELL_matrix0(){}
        ELL_matrix0(int rows, int cols, int MaxCol) { construct(rows, cols, MaxCol); }
        ELL_matrix0(int n, int MaxCol) { construct(n, n, MaxCol); }
        virtual ~ELL_matrix0(){ destruct(); };

        //* ------------------------------------------

        // ---------------------------
        void resize(int rows, int cols, int MaxCol){ construct(rows, cols, MaxCol); }
        void resize(int n, int MaxCol) { construct(n, n, MaxCol); }
        // ---------------------------

        // ---------------------------
        int row() const { return rows_; }
        int col() const { return cols_; } 
        int MaxRow()  { return max_idx_size_; }
        // ---------------------------

        std::vector<std::tuple<int, int, double> > get_mmt();


        // ! ----------------- set and get
        T at(int row, int col) const;
        ELL_matrix0<T> & set(int row, int col, T val);
        ELL_matrix0<T> & set(CSR_type c);




        void Show_ELL(void);
        // ! ------------------------------- Spmv
        std::vector<T> operator * (const std::vector<T> & x) const
        { return multiply_omp(x); }

        std::vector<T> multiply(const std::vector<T> & x) const;
        bool multiply(const std::vector<T> & x, std::vector<T> & r);

        std::vector<T> multiply_omp(const std::vector<T> & x) const;
        bool multiply_omp(const std::vector<T> & x, std::vector<T> & r);
        bool multiply_omp_Point(const std::vector<T> & x, std::vector<T> & r);
        


        #if defined(ELL0_MPI_ON)
        bool multiply_mpi( std::vector<T> & x, std::vector<T> & r);
        #endif

        // * -------------------------------


        #if defined (ELL0_MPI_ON)
            void mpi_init(MPI_Comm comm_world, int total)
                { mpi_.init(comm_world, total); }

            void mpi_init(MPI_Comm comm_world, int total, std::vector<int> beg, std::vector<int> end)
                { mpi_.init(comm_world, total, beg, end); }

            auto get_mpi()const{return mpi_;}
        #endif

        bool analyse();

        //--------------------------------------------- // friend function
        template<typename X>
            friend std::ostream & operator << (std::ostream & os, const ELL_matrix0<X> & matrix);
        //---------------------------------------------


        private:
            int rows_ , cols_, max_idx_size_, init_idx_;
            std::vector<int> Idx_;
            std::vector<double> Val_;
            std::vector<int> curr_idx_size;

            //---------------------------------------------
            #if defined (ELL0_MPI_ON)
                mat::mpiT mpi_;
            #endif
            //---------------------------------------------

        // * private function
        void construct(int rows, int cols, int Max_colIdx);


        // ----------- 
        void construct(int Maxoff, int Max_colIdx);
        void destruct(void){}
        // ----------- 
    };



    // ! ------------------------- construct and resize -------------------------

    template<typename T>
    inline void ELL_matrix0<T>::construct(int rows, int cols, int Max_colIdx)
    {
        if (rows < 0 || cols < 0) {
			throw std::invalid_argument("Matrix dimensions cannot be negative.");
		}

        rows_ = rows;
        cols_ = cols;
        max_idx_size_ = Max_colIdx;
        init_idx_ = cols-1;
        Idx_.resize(Max_colIdx*rows,  init_idx_);
        Val_.resize(Max_colIdx*rows,  0);
        curr_idx_size.resize(rows, 0);
    }


    // ------------------------------------------------------------------
    template<typename T>
    inline std::ostream & operator << (std::ostream & os, const ELL_matrix0<T> & matrix){
        // ----------------
        os << std::endl;

        // ----------------

        for(int i = 0, iter; i < matrix.rows_;++i){
            for(int j = 0; j < matrix.cols_;++j){
                const int pt = i*matrix.max_idx_size_;
                os << std::setw(4);
                for( iter=0; iter < matrix.curr_idx_size[i] ;++iter ){

                    if ( j == matrix.Idx_.at(pt+iter) ){
                        os << matrix.Val_.at(pt+iter);
                        break;
                    }
                }
                if (iter == matrix.curr_idx_size.at(i)  ){
                    os <<"x" ;
                }
            }
            os << std::endl;
        }
		return os;
    }
    // ------------------------------------------------------------------


    // * ------- set and get
    template<typename T>inline ELL_matrix0<T> & ELL_matrix0<T>::
    set(ELL_matrix0::CSR_type c){
        // --------------
        auto [prt, idx, val] = c;
        // --------------
        rows_ = prt.size()-1;
        // --------------

        max_idx_size_ = 0;

        // find max index
        for (int row = 0; row < rows_; ++row){
           max_idx_size_  = std::max( prt[row+1] - prt[row], max_idx_size_);
        }

        resize(rows_, max_idx_size_);

        // --------------
        for (int row = 0; row < rows_; ++row){
            for (int j = prt[row]; j < prt[row+1] ;++j){
                set(row, idx[j], val[j]);
            }
        }
        // --------------

        // ! We assume that your matrix is a n by n matrix for our convienious.
        cols_ = rows_;  

        return *this;
    }


    template<typename T>
    inline ELL_matrix0<T> & ELL_matrix0<T>::
    set(int row, int col, T val) {

        if (row < 0 || col < 0) {
			throw std::invalid_argument("Matrix dimensions cannot be negative.");
		}


        if (col >= cols_ ){
            cout << col << " >= " << cols_ << std::endl;
			throw std::invalid_argument("colIdx >= max_idx_size_ ");
        }


        if (row >= rows_)
			throw std::invalid_argument("row >= rows_");

        const int pt = row*max_idx_size_;

        if (val == T()){
            // bool NEW = true;
            for(int i = 0; i < curr_idx_size[row] ;++i ){

                if (col == Idx_[pt+i]){
                    // NEW = false; --> remove 
                    // -----------------------
                    for(int j = i ; j < curr_idx_size[row];++j){
                        Idx_[pt+j] = Idx_[pt+j+1];
                        Val_[pt+j] = Idx_[pt+j+1];
                    }
                    // -----------------------
                    Idx_[pt+curr_idx_size[row]] = init_idx_;
                    Val_[pt+curr_idx_size[row]] = T();
                    curr_idx_size[row]--;

                    return *this;
                }
            }
            return *this;
        }
        else{

            bool NEW = true;
            for(int i = 0; i < curr_idx_size[row] ;++i ){
                
                if ( col == Idx_[pt+i] ){
                    NEW = false;

                    Idx_[pt+i] = col;

                    Val_[pt+i] = val;
                    return *this;
                }
            }

            // -----------------------------------

            // -----------------------------------
            if (NEW){

                Idx_[pt+curr_idx_size[row]] = col;
                Val_[pt+curr_idx_size[row]] = val;  
                curr_idx_size[row]++;
            }
            // -----------------------------------
            

            // -----------------------------------
            if (curr_idx_size[row] > max_idx_size_)
            {
                throw std::invalid_argument(
                    "[ELL0::set] -> curr_idx_size[row] >max_idx_size_");
            }
            // -----------------------------------

            return *this;
        }
    }

    template<typename T>
    inline bool ELL_matrix0<T>::analyse(void) {

        // -----------------------------------
        std::vector<int> a(max_idx_size_+1, 0);
        for (auto v:curr_idx_size){ ++a[v]; }
        // -----------------------------------
        int i{0};
        // -----------------------------------
        for (auto v:a){
            cout << "i=" << i++ << ", " << v << "\n"; 
        }
        // -----------------------------------
        return true;
    }



    template<typename T>
    inline void ELL_matrix0<T>::Show_ELL(void){
        cout << "   |val";
        for (int i = 0;i < max_idx_size_ ;++i){ 
              cout<< std::setw(4) << "";
        }

        cout << "    |Idx\n";

        for (int i = 0;i < rows_ ;++i){

            // ----------------------------------------
            
            cout << std::setw(3) << i << "|";
            
            // ----------------------------------------
            for (int j = 0;j < max_idx_size_ ;++j){  
                cout<< std::setw(4) << Val_[i*max_idx_size_+j] << " " ;
            }
            // ----------------------------------------

            cout << "|";

            // ----------------------------------------
            for (int j = 0;j < max_idx_size_ ;++j){ 
                cout<< std::setw(4) << Idx_[i*max_idx_size_+j] << " " ;
            }
            // ----------------------------------------

            cout << std::endl;
        }
    }



    //  ! ---------------------------- Spmv ----------------------------

    template<typename T> inline bool ELL_matrix0<T>::
    multiply(const std::vector<T> & x, std::vector<T> & r) 
    {
        for(int i = 0; i < rows_;++i){
            const int pt = i*max_idx_size_;
            // for(int j = 0; j < curr_idx_size[i];++j){
            for(int j = 0; j < max_idx_size_;++j){
                r[i] += Val_[pt+j] * x[Idx_[pt+j]];
            }
        }
        return true;
    }

    template<typename T> inline std::vector<T> ELL_matrix0<T>::
    multiply(const std::vector<T> & x) const
    {
        std::vector<T> r(rows_, T());
        multiply(x, r);
        return r;
    }

    template<typename T> inline std::vector<T> ELL_matrix0<T>::
    multiply_omp(const std::vector<T> & x) const
    {
        std::vector<T> r(rows_, T());
        multiply_omp(x, r);
        return r;
    }



    template<typename T> inline bool ELL_matrix0<T>::
    multiply_omp(const std::vector<T> & x, std::vector<T> & r) 
    {
        #pragma omp parallel for
        for(int i = 0; i < rows_;++i){
            const int pt = i*max_idx_size_;
            r[i] = T();
            // for(int j = 0; j < curr_idx_size[i];++j){
            for(int j = 0; j < max_idx_size_;++j){
                r[i] += Val_[pt+j] * x[Idx_[pt+j]];
            }
        }
        return true;
    }


    template<typename T> inline bool ELL_matrix0<T>::
    multiply_omp_Point(const std::vector<T> & x, std::vector<T> & r) 
    {
        auto VAL = Val_.data();
        auto IDX = Idx_.data();
        auto X = x.data();
        auto R = r.data();

        #pragma omp parallel for default(none) \
        firstprivate(max_idx_size_, IDX, VAL,R, X)
        for(int i = 0; i < rows_;++i){
            const int pt = i*max_idx_size_;
            *(R+i) = T();
            // for(int j = 0; j < curr_idx_size[i];++j){
            for(int j = 0; j < max_idx_size_;++j){
                *(R+i) += *(VAL+pt+j) * *(X + *(IDX+pt+j) );
            }
        }
        return true;
    }



    #if defined (ELL0_MPI_ON)

    template<typename T> inline bool ELL_matrix0<T>::
    multiply_mpi( std::vector<T> & x, std::vector<T> & r) 
    {
        mpi_.allocate_vec(x);
        for(auto i = mpi_.beg(); i < mpi_.end();++i){
            const int pt = i*max_idx_size_;
            r[i] = T();

            for(int j = 0; j < max_idx_size_;++j){
                r[i] += Val_[pt+j] * x[Idx_[pt+j]];
            }
        }
        return true;
    }

    #endif
}
