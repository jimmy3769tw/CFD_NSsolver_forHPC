#pragma once
#include <cstdlib>  // for div
#include <string>
#include <vector>
#include <tuple>
#include <iostream>
#include <stdexcept>
#include "mpi.h"
using namespace std;

namespace mat{


    class mpiT{

        public:

        inline int beg(){return beg_;}
        inline int end(){return end_;}

        void init(MPI_Comm comm_world, int total_len, std::vector<int> beg_in, std::vector<int> end_in){

            com_ = comm_world;

            total_len_ = total_len;

            if (beg_in.size()!= end_in.size() ) 
            throw std::runtime_error(" didn't match");

            MPI_Comm_rank(comm_world, &rank_);

            MPI_Comm_size(comm_world, &size_ );

            init_local(beg_in, end_in);

            std::cout <<" \nbeg(), end()"<< beg() << ", " << end();

            init_global(beg_in, end_in);


            int sum{0};
            for (auto & x: table_len_ )
            {
                sum+=x; 
            }

            if (sum != total_len ) throw std::runtime_error("wrong the total length didn't match");
        }

        void init(MPI_Comm comm_world, int total_len){

            com_ = comm_world;

            total_len_ = total_len;

            MPI_Comm_rank(comm_world, &rank_);

            MPI_Comm_size(comm_world, &size_ );

            init_local();

            init_global();
        }

        void allocate_vec(std::vector<double> &x){

            MPI_Barrier(com_);

            for(size_t i = 0; i < size_ ; ++i){
                MPI_Bcast(  (void *)&x[table_beg_[i]], table_len_[i], MPI_DOUBLE, i, com_);

                MPI_Barrier(com_);
            }

            MPI_Barrier(com_);
        }   

        double allocate_scale_sum(double &x){
            double temp_g{x};
            MPI_Allreduce(&x, &temp_g, 1, MPI_DOUBLE, MPI_SUM, com_);
            return temp_g;
        }

        int allocate_scale_sum(int &x){
            int temp_g{x};
            MPI_Allreduce(&x, &temp_g, 1, MPI_INT, MPI_SUM, com_);
            return temp_g;
        }


        inline double 
        dot(const vector<double> & a, const vector<double> & b)
        {
            double temp_l{0.0};
            double temp_g{0.0};
            auto A = a.data()+beg_;
            auto B = b.data()+beg_;

            #pragma omp parallel for simd reduction(+:temp_l) default(none) \
            firstprivate(A, B, len_) schedule(static)
            for(size_t i=0; i<len_ ; ++i)
            {  temp_l += *(A+i) * *(B+i); }


            MPI_Allreduce(&temp_l, &temp_g, 1, MPI_DOUBLE, MPI_SUM, com_);
            return temp_g;
        }

        // --------------------------
        inline double 
        L2Norm(const std::vector<double> & a)
        { return sqrt( dot(a,a) ); }
        // --------------------------


        inline void init(vector<double> &a, const double &val) {
            auto A = a.data()+beg_;

            #pragma omp parallel for simd default(none) \
            firstprivate(A, val, len_) schedule(static)
            for(size_t i=0; i<len_ ; ++i)
            {  *(A+i) = val; }
        }


        inline void init(vector<double> &a, const vector<double> & b) {
            auto A = a.data()+beg_;
            auto B = b.data()+beg_;

            #pragma omp parallel for simd default(none) \
            firstprivate(A, len_, B) schedule(static)
            for(size_t i=0; i<len_ ; ++i)
            {  *(A+i) = *(B+i); }

        }
        // --------------------------
        inline void 
        zero(vector<double> & a){  init(a, 0.0); }
        // --------------------------

        // --------------------------
        inline void
        copy(vector<double> & a, vector <double> & r) { init( r, a); }
        // --------------------------


        // --------------------------
        inline double 
        getError(vector<double> &a, vector<double> &b){

            double temp_g{0.0}, temp_l{0.0};
            auto A = a.data()+beg_;
            auto B = b.data()+beg_;

            #pragma omp parallel for simd default(none) \
            firstprivate(A, len_, B) schedule(static) reduction(+:temp_l)
            for(size_t i=0; i<len_ ; ++i)
            {  temp_l += abs(*(A+i) - *(B+i)); }


            MPI_Allreduce(&temp_l, &temp_g, 1, MPI_DOUBLE, MPI_SUM, com_);

            return temp_g;
        }
        // --------------------------


        private:

        MPI_Comm com_;

        int rank_, size_, total_len_, beg_, end_, len_;

        std::vector<int> table_end_, table_beg_, table_len_;




        inline void init_local(){
            // -------------------
            auto [q, r] = std::div(total_len_, size_);

            tie(beg_, end_) = get_beg_end(rank_, q, r);
            len_ = end_ - beg_;
            // -------------------
        }



        inline void init_local(std::vector<int> beg, std::vector<int> end );
        // ---------------------------------------------
        inline void init_global(std::vector<int> beg, std::vector<int> end );
        inline void init_global();
        // ---------------------------------------------
        inline std::tuple<int, int> get_beg_end( int pid, const int q, const int r);
        // ---------------------------------------------
    };

        // -----------------------------------------------------------------
        inline void mpiT::init_local(std::vector<int> beg, std::vector<int> end ){
            // -------------------
            beg_ = beg.at(rank_);
            end_ = end.at(rank_);
            len_ = end_ - beg_;
            // -------------------
        }
        // -----------------------------------------------------------------

        inline std::tuple<int, int> mpiT::get_beg_end(
                int pid, const int q, const int r
        ){
            int beg = pid*q;

            int end = (pid+1)*q;
            
            if (pid < r){
                beg += pid;
                end += pid+1;
            }
            else {
                beg += r;
                end += r;
            }
            
            return std::make_tuple(beg, end);
        }




        inline void mpiT::init_global(std::vector<int> beg, std::vector<int> end ){
            // -------------------
            table_len_.resize(size_);

            table_beg_.resize(size_);

            table_end_.resize(size_);
            // -------------------
            table_beg_ = beg;
            table_end_ = end;

            for (int i = 0; i < size_ ;++i){
                table_len_[i] = table_end_[i] - table_beg_[i];
            }
            // -------------------
        }


        inline void mpiT::init_global(){

            auto [q, r] = std::div(total_len_, size_);

            // -------------------
            table_len_.resize(size_);

            table_beg_.resize(size_);

            table_end_.resize(size_);
            // -------------------

            // -------------------
            for (int i = 0; i < size_ ;++i){
                tie(table_beg_[i], table_end_[i]) = get_beg_end(i, q, r);
            }

            for (int i = 0; i < size_ ;++i){
                table_len_[i] = table_end_[i] - table_beg_[i];
            }
            // -------------------
        }




}

