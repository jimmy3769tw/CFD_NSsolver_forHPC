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





#ifdef ELL_MPI_ON
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
    using namespace std;
    template<typename T>
    class ELL_matrix
    {
        public:
        // ---------------------
        using CSR_type = 
            // <ptr, indices, values>
            std::tuple< std::vector<int>, std::vector<int>, std::vector<T> >; 
        // ---------------------


        //! ------ Constructor & Destructor ---------
        ELL_matrix(){}
        ELL_matrix(int rows, int cols, int max_idx_size) { construct(rows, cols, max_idx_size); }
        ELL_matrix(int n, int max_idx_size) { construct(n, n, max_idx_size); }

        virtual ~ELL_matrix(){ destruct(); };

        //* ------------------------------------------

        // ---------------------------
        void resize(int rows, int cols, int max_idx_size){ construct(rows, cols, max_idx_size); }
        void resize(int n, int max_idx_size) { construct(n, n, max_idx_size); }
        // ---------------------------

        // ---------------------------
        int row() const { return rows_; }
        int col() const { return cols_; } 
        int MaxRow()  { return max_idx_size_; }
        // ---------------------------

        std::vector<std::tuple<int, int, double> > get_mmt();


        // ! ----------------- set and get
        T& at(int row, int col);
        T& operator ()(const  int row, const  int col);

        ELL_matrix<T> & set(int row, int col, T val);
        ELL_matrix<T> & set(CSR_type c);

        void Show_ELL(void);
        // ! ------------------------------- Spmv
        std::vector<T> operator * ( std::vector<T> & x) 
        { return multiply_omp(x); }

        std::vector<T> multiply(std::vector<T> & x);
        bool multiply( std::vector<T> & x, std::vector<T> & r);

        std::vector<T> multiply_omp(std::vector<T> & x);
        bool multiply_omp( std::vector<T> & x, std::vector<T> & r);

        #if defined(ELL_MPI_ON)
        bool multiply_mpi( std::vector<T> & x, std::vector<T> & r);
        #endif
        // * -------------------------------

        // ! ------------------------------- SpmSp
        // ELL_matrix<T> operator * (const ELL_matrix<T> & matrix) const;
        // ELL_matrix<T> multiply(const ELL_matrix<T> & matrix) const;
        // * -------------------------------


        bool analyse();

        //--------------------------------------------- // friend function
        template<typename X>
            friend std::ostream & operator << (std::ostream & os, const ELL_matrix<X> & matrix);
        //---------------------------------------------


        //---------------------------------------------
        void sort_Idx(void);
        CSR_type get_CSR(void);
        //---------------------------------------------


        #if defined (ELL_MPI_ON)
            void mpi_init(MPI_Comm comm_world, int total)
                { mpi_.init(comm_world, total); }

            void mpi_init(MPI_Comm comm_world, int total, std::vector<int> beg, std::vector<int> end)
                { mpi_.init(comm_world, total, beg, end); }

            auto get_mpi()const{return mpi_;}
        #endif

        private:
            int rows_ , cols_, max_idx_size_, init_idx_;
            std::vector<std::pair<int,T>> Idx_val_;
            std::vector<int> curr_idx_size;

            //---------------------------------------------
            #if defined (ELL_MPI_ON)
                mat::mpiT mpi_;
            #endif
            //---------------------------------------------

        void construct(int rows, int cols, int max_idx_size);

        // ----------- 
        void construct(int Maxoff, int max_idx_size);

        void destruct(void){}
        // ----------- 
    };



    // ! ------------------------- construct and resize -------------------------


    template<typename T>
    inline void ELL_matrix<T>::
    construct(int rows, int cols, int max_idx_size)
    {
        if (rows < 0 || cols < 0) {
			throw std::invalid_argument("Matrix dimensions cannot be negative.");
		}

        rows_ = rows;
        cols_ = cols;
        max_idx_size_ = max_idx_size;
        init_idx_ = cols-1;
        Idx_val_.resize(max_idx_size*rows, std::make_pair(init_idx_, 0.0));
        curr_idx_size.resize(rows, 0);
    }

    // * ------------------------------------------------------------------------

    template<typename T>
    inline T& ELL_matrix<T>::operator()(int row, int col){
        const int pt = row*max_idx_size_;
        for(int j = 0; j < cols_;++j){
            auto & [IDX, VAL] = Idx_val_.at(pt);
            if (col == IDX)  return VAL;    
        }
    }

    // template<typename T>
    // inline ELL_matrix<T>::
    // set(int row, int col){
    //     if (rows < 0 || cols < 0) {
	// 		throw std::invalid_argument("Matrix dimensions cannot be negative.");
	// 	}
    //     return *this[row, col];
    // }

    // ------------------------------------------------------------------
    template<typename T>
    inline std::ostream & operator << (std::ostream & os, const ELL_matrix<T> & matrix){
        // ----------------
        os << std::endl;

        // ----------------

        for(int i = 0, iter; i < matrix.rows_;++i){
            for(int j = 0; j < matrix.cols_;++j){
                const int pt = i*matrix.max_idx_size_;
                os << std::setw(4);
                for( iter=0; iter < matrix.curr_idx_size[i] ;++iter ){
                    auto [IDX,VAL] = matrix.Idx_val_.at(pt+iter);
                    if ( j == IDX ){
                        os << VAL;
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
    template<typename T>inline ELL_matrix<T> & ELL_matrix<T>::
    set(ELL_matrix::CSR_type c){

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

        // --------------------
        resize(rows_, max_idx_size_);
        // --------------------

        for (int row = 0; row < rows_; ++row){
            for (int j = prt[row]; j < prt[row+1] ;++j){
                set(row, idx[j], val[j]);
            }
        }

        // --------------
        cols_ = rows_;
        // --------------

        return *this;
    }

    template<typename T>inline ELL_matrix<T> & ELL_matrix<T>::
    set(int row, int col, T val) {
        // ---------------------
        if (row < 0 || col < 0) {
			throw std::invalid_argument(
                "Matrix dimensions cannot be negative.");
		}

        if (col >= cols_ ){
			throw std::invalid_argument(
                "col >= cols_");
        }

        if (row >= rows_){
			throw std::invalid_argument(
                "row >= rows_");
        }
        // ---------------------


        // --------------------- if val == T()
        const int pt = row*max_idx_size_;

        if (val == T()){  

            // bool NEW = true;
            for(int i = 0; i < curr_idx_size[row] ;++i ){
                // -------------------
                int IDX;
                tie(IDX, ignore) = Idx_val_[pt+i];
                // -------------------
                if (col == IDX){
                    // NEW = false; --> shift(remove)

                    // -----------------------
                    for(int j = i ; j < curr_idx_size[row];++j){
                       Idx_val_[pt+j] = Idx_val_[pt+j+1]; 
                    }
                    // -----------------------

                    Idx_val_[pt+curr_idx_size[row]] = make_pair(init_idx_,T()); 
                    curr_idx_size[row]--;

                    return *this;
                }
            }
            return *this;
        }
        else{ // --------------------- val != T()
            bool NEW = true;
            for(int i = 0; i < curr_idx_size[row] ;++i ){

                auto [IDX, VAL] = Idx_val_[pt+i];
                if ( col == IDX ){

                    NEW = false;

                    Idx_val_[pt+i] = make_pair(col, val);

                    return *this;
                }
            }
        
            if (NEW){ Idx_val_[pt+(curr_idx_size[row]++)] = make_pair(col, val); }

            // ------------------------
            if (curr_idx_size[row] > max_idx_size_){
                throw invalid_argument(
                    "[ELL::set] -> curr_idx_size[row] >max_idx_size_");
            }
            // ------------------------

            return *this;
        }
    }

    template<typename T>inline bool ELL_matrix<T>::
    analyse() {

        std::vector<int> a(max_idx_size_+1, 0);
        for (auto v:curr_idx_size){
            ++a[v];
        }

        int i = 0;

        for (auto v:a){
            std::cout << "i=" << i++ << ", " << v << "\n"; 
        }
        return true;
    }


    template<typename T>inline void ELL_matrix<T>::
    sort_Idx(){

        for (int row = 0; row < rows_ ; ++row){
            const int pt = row*max_idx_size_;
            for (int i = 0; i < curr_idx_size[row] ; ++i){
                std::sort(Idx_val_.begin() + pt , Idx_val_.begin()+ pt +curr_idx_size[row]);
            }
        }
    }


    template<typename T>inline typename ELL_matrix<T>::CSR_type ELL_matrix<T>::
    get_CSR(){
        // --------------
        sort_Idx();
        // --------------

        // --------------
        std::vector<int> CSR_ptr(rows_+1,0);
        std::vector<int> CSR_idx;
        std::vector<T> CSR_data;
        // --------------

        int iter = 0;

        for (int row = 0; row < rows_ ; ++row){
            const int pt = row*max_idx_size_;
            for (int i = 0; i < curr_idx_size[row] ; ++i){
                auto [IDX,VAL] = Idx_val_[pt+i];
                CSR_idx.push_back(IDX);
                CSR_data.push_back(VAL);
                ++iter;
            }
            CSR_ptr.at(row+1) = iter;
        }

        return std::make_tuple(CSR_ptr, CSR_idx, CSR_data);
    }



    template<typename T>
    inline std::vector<std::tuple<int, int , double> > ELL_matrix<T>::
    get_mmt(){

        sort_Idx();

        std::vector<std::tuple<int, int , double> > mmt;
        // -----------------------
        for (int row = 0; row < rows_ ; ++row){
            // -----------------------
            const int pt = row*max_idx_size_;
            // -----------------------

            // -----------------------
            for (int i = 0; i < curr_idx_size[row] ; ++i){

                auto [IDX,VAL] = Idx_val_[pt+i];
                mmt.push_back(std::tie(row, IDX, VAL));
            }
            // -----------------------
        }
        // -----------------------

        return mmt;
    }



    template<typename T>inline void ELL_matrix<T>::
    Show_ELL(){
        std::cout << "   |val";
        for (int i = 0;i < max_idx_size_ ;++i){ 
              std::cout<< std::setw(4) << "";
        }

        std::cout << "    |Idx\n";

        for (int i = 0;i < rows_ ;++i){   // row
            std::cout << std::setw(3) << i << "|";
            for (int j = 0;j < max_idx_size_ ;++j){  //col
                T VAL;
                std::tie(std::ignore, VAL) = Idx_val_.at(i*max_idx_size_+j);
                std::cout<< std::setw(4) << VAL << " " ;
            }
            std::cout << "|";
            for (int j = 0;j < max_idx_size_ ;++j){  //col
                int IDX;
                std::tie(IDX,std::ignore) = Idx_val_.at(i*max_idx_size_+j);
                std::cout<< std::setw(4) << IDX << " " ;
            }

            std::cout << std::endl;
        }
    }



    //  ! ---------------------------- Spmv ----------------------------

    template<typename T> inline std::vector<T> ELL_matrix<T>::
    multiply(std::vector<T> & x)
    {
        std::vector<T> r(rows_, T());
        multiply(x, r);
        return r;
    }

    template<typename T> inline std::vector<T> ELL_matrix<T>::
    multiply_omp(std::vector<T> & x)
    {
        std::vector<T> r(rows_, T());
        multiply_omp(x, r);
        return r;
    }





    template<typename T> inline bool ELL_matrix<T>::
    multiply(std::vector<T> & x, std::vector<T> & r) 
    {

        for(int i = 0; i < rows_;++i){

            const int pt = i*max_idx_size_;
            // for(int j = 0; j < curr_idx_size[i];++j){ // for one of otpmization 
            for(int j = 0; j < max_idx_size_;++j){
                auto [IDX,VAL]  = Idx_val_[pt+j];
                r[i] += VAL * x[IDX];
            }
        }
        return true;
    }


    template<typename T> inline bool ELL_matrix<T>::
    multiply_omp(std::vector<T> & x, std::vector<T> & r) 
    {
        #pragma omp parallel for 
        for(int i = 0; i < rows_;++i){
            const int pt = i*max_idx_size_;
            r[i] = T();
            // for(int j = 0; j < curr_idx_size[i];++j){ // for one of otpmization 
            for(int j = 0; j < max_idx_size_;++j){

                auto [IDX,VAL]  = Idx_val_[pt+j];
                r[i] += VAL * x[IDX];
            }
        }
        return true;
    }




    #if defined (ELL_MPI_ON)

    template<typename T> inline bool ELL_matrix<T>::
    multiply_mpi( std::vector<T> & x, std::vector<T> & r) 
    {
        mpi_.allocate_vec(x);
        for(auto i = mpi_.beg(); i < mpi_.end();++i){
            const int pt = i*max_idx_size_;
            r[i] = T();

            for(int j = 0; j < max_idx_size_;++j){

                auto [IDX,VAL]  = Idx_val_[pt+j];
                r[i] += VAL * x[IDX];
            }
        }
        return true;
    }

    #endif
}
