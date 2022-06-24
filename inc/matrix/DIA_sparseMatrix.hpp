#pragma once 

#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>


namespace mat
{

    template<typename T>
    class DIA_matrix
    {
        public:
        // --- Constructor & Destructor ---
        DIA_matrix(int Maxoff, size_t Maxi);
        virtual ~DIA_matrix(void);
    
        // ---------- matrix infomation
        size_t sizeof_offset(void) const;  
        size_t sizeof_i(void) const;
        std::vector<T> multiply(const std::vector<T> & x) const;


        // ----------get and set 
        T at(int off,size_t i) const;
        T atCo(size_t row,size_t col) const;
        DIA_matrix<T> & set(int off, size_t i, T val);
        DIA_matrix<T> & setCo(size_t row, size_t col, T val);



        // ---------- operations
        DIA_matrix<T> & operator = (const DIA_matrix<T> &m);
        std::vector<T> operator * (const std::vector<T> & x) const;
        // ----------- 



        // ---------- friend funciton
        template<typename X>
        friend std::ostream & operator << (std::ostream & os, const DIA_matrix<T> & martix);


        // private: //=========================
    
        // ----------- Local variables
        int Maxoff_{0};
        size_t Maxi_{0}; // m:row, n:col
        size_t Sizeof_Maxoff_;
        std::vector<T> val_;
        std::vector<size_t> i_;
        // ----------- 

        // ----------- construct & destruct
        void construct(int Maxoff, size_t Maxi);
        void destruct(void);
        // ----------- 

    };


    // ==================== Cpp ==================

   // ------------ Constructor & Destructor
    template<typename T>
    inline DIA_matrix<T>::DIA_matrix(int Maxoff, size_t Maxi){
        this->construct(Maxoff, Maxi);
    }


    template<typename T>
    inline DIA_matrix<T>::~DIA_matrix(void){
        this->destruct();
    }



    template<typename T>
    inline void DIA_matrix<T>::construct(int Maxoff, size_t Maxi){
        if (Maxoff > Maxi ) {
			throw std::invalid_argument("Maxoff > Maxi.");
		}
        if (Maxi < 0){
			throw std::invalid_argument("Maxi < 0");
        }
        this -> Maxoff_ = Maxoff;
        this -> Sizeof_Maxoff_ = (Maxoff*2+1);
        this -> Maxi_ = Maxi;
        i_.resize(Maxi);
        val_.resize(this->Sizeof_Maxoff_*Maxi);
    }


    // ------------ get and set

    template<typename T>
    inline DIA_matrix<T> & DIA_matrix<T>::set(int off, size_t i, T val){
        if ( i < 0 || i >= this->Maxi_){
			throw std::invalid_argument("Coordination out of range(i < 0 || i >= this->Maxi_)).");
		}

        if (-off > this->Maxoff_ || off > this->Maxoff_ ){
			throw std::invalid_argument("out of range(off < -Maxoff_ ||  off > this->Maxoff_).");
        }

        const size_t idx = Sizeof_Maxoff_*i + (this->Maxoff_ + off);
        val_[idx] = val;
    }

    // | row 0, col 0 | row 0, col 1 |
    // | row 1, col 0| row 1 , col 1 |


    template<typename T>
    inline T DIA_matrix<T>::at(int off, size_t i) const{
        if ( i < 0 || i >= this->Maxi_){
			throw std::invalid_argument("Coordination out of range(i < 0 || i >= this->Maxi_)).");
		}

        if (-off > this->Maxoff_ || off > this->Maxoff_ ){
			throw std::invalid_argument("out of range(off < -Maxoff_ ||  off > this->Maxoff_).");
        }
        const size_t idx = this -> Sizeof_Maxoff_*i + (this->Maxoff_ + off);
        return val_[idx];
    }

    template<typename T>
    inline DIA_matrix<T> & DIA_matrix<T>::setCo(size_t row, size_t col, T val){
        if (row > this->Maxi_ || row < 0 || col < 0 || col > this ->Maxi_){
            throw std::invalid_argument("row > this->Maxi_; col > this ->Maxi_");
        }

        const int off = col-row; //j - i
        const int i = std::min(col,row);
        this->set(off,i, val);
    }

    template<typename T>
    inline T DIA_matrix<T>::atCo(size_t row, size_t col) const{
        if (row > this->Maxi_ || row < 0 || col < 0 || col > this ->Maxi_){
            throw std::invalid_argument("row > this->Maxi_; col > this ->Maxi_");
        }
        const int off = col-row; //j - i
        const int i = std::min(col,row);
        return this->at(off,i);
    }

    template<typename T>
    inline void DIA_matrix<T>::destruct(void)
    {
        // void
    }

    // ------------friend fuction
    template<typename T>
    inline std::ostream & operator << (std::ostream & os, const DIA_matrix<T> & matrix){
        os << std::endl;
        for(int i = 0; i < matrix.Maxi_;++i){
            for(int j = 0; j < matrix.Maxi_;++j){
                os << std::setw(4)  ;
                if (j <= i+matrix.Maxoff_ && j >= i-matrix.Maxoff_){
                    os  << matrix.at(j-i,std::min(j,i));
                }else{
                    os << 0. ;
                }
            }
            os << std::endl;
        }
    }

    template<typename T>
    inline std::vector<T> DIA_matrix<T>::operator * (const std::vector<T> & x) const
    {
        return this->multiply(x);
    }


    template<typename T> 
    inline std::vector<T> DIA_matrix<T>::multiply(const std::vector<T> & x) const
    {
        if (this -> Maxi_ != x.size() ){
            throw std::invalid_argument("col of matrix != vector.size()");
        }
        std::vector<T> result(this-> Maxi_, T()); //parallel vector ??

        for(size_t i = 0; i < this->Maxi_;++i){
            for(size_t j = 0; j < this->Maxi_;++j){
                result[i] += this->atCo(i,j) * x[j];
            }
        }
        return result;
    }

}


