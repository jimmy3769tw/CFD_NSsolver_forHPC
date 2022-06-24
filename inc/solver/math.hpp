
#pragma once 
#include <complex>
#include <vector>
#include <tuple>
#include <omp.h>

#include<iostream>
namespace math{

    using namespace std;

    inline void init(vector<double> &a, const double &val) {
        auto A = a.data();
        auto len = a.size();

        #pragma omp parallel for simd default(none) \
        firstprivate(A, val, len) schedule(static)
        for(size_t i=0; i<len ; ++i)
        {  *(A+i) = val; }
    }


    inline void init(vector<double> &a, const vector<double> & b) {
        auto len = a.size();
        auto A = a.data();
        auto B = b.data();

        #pragma omp parallel for simd default(none) \
        firstprivate(A, len, B) schedule(static)
        for(size_t i=0; i<len ; ++i)
        {  *(A+i) = *(B+i); }

    }


    inline double 
    dot(const vector<double> & a, const vector<double> & b)
    {
        double r{0.0};
        auto len = a.size();
        auto A = a.data();
        auto B = b.data();

        #pragma omp parallel for simd reduction(+:r) default(none) \
        firstprivate(A, B, len) schedule(static)
        for(size_t i=0; i<len ; ++i)
        {  r += *(A+i) * *(B+i); }

        return r;
    }
    // --------------------------

    inline double 
    dotNsimd(const vector<double> & a, const vector<double> & b)
    {
        double r{0.0};
        auto len = a.size();
        auto A = a.data();
        auto B = b.data();

        #pragma omp parallel for reduction(+:r) default(none) \
        firstprivate(A, B, len) schedule(static)
        for(size_t i=0; i<len ; ++i)
        {  r += *(A+i) * *(B+i); }

        return r;
    }
    // --------------------------


    // --------------------------
    inline double 
    dot_Npoint(const vector<double> & a, const vector<double> & b)
    {
        double r{0.0};

        #pragma omp parallel for simd reduction(+:r) \
        default(none) shared(a, b) schedule(static)
        for(size_t i=0; i<a.size() ; ++i)
        {  r += a[i] * b[i]; }

        return r;
    }
    // --------------------------


    // --------------------------
    inline double 
    dot_Npoint_Nsimd(const vector<double> & a, const vector<double> & b)
    {
        double r{0.0};

        #pragma omp parallel for reduction(+:r) \
        default(none) shared(a, b) schedule(static)
        for(size_t i=0; i<a.size() ; ++i)
        {  r += a[i] * b[i]; }

        return r;
    }
    // --------------------------


    inline double 
    dotseq(const vector<double> & a, const vector<double> & b)
    {
        double r{0.0};

        for(size_t i = 0; i < a.size() ; ++i)
        {  r += a[i] * b[i]; }

        return r;
    }
    // --------------------------

    // inline double 
    // dot(const vector<double> & a, const vector<double> & b)
    // {
    //     double r{0.0};
    //     return std::inner_product(a, b, );
    // }


    // --------------------------
    inline double 
    L2Norm(const std::vector<double> & a)
    { return sqrt( dot(a,a) ); }
    // --------------------------


    // --------------------------
    inline void 
    zero(vector<double> & a){  init(a, 0.0); }
    // --------------------------



    // --------------------------
    inline void
    copy(vector<double> & a, vector <double> & r)
    { init( r, a); }
    // --------------------------


    // --------------------------
    inline double 
    getError(vector<double> &a, vector<double> &b){

        double r{0.0};
        auto A = a.data();
        auto len = a.size();
        auto B = b.data();

        #pragma omp parallel for simd default(none) \
        firstprivate(A, len, B) schedule(static) reduction(+:r)
        for(size_t i=0; i<len ; ++i)
        {  r += abs(*(A+i) - *(B+i)); }

        return r;
    }
    // --------------------------
}