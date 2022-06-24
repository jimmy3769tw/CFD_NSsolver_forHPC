#pragma once

#include"controlPanel.hpp"

// Pressure_transform_X_result(t1, Mx, gridA);


void Pressure_transform_X_result(
    pressure &t1, 
    MxClass &Mx ,
    calDomain& Lo,
    grid & gridA
){
    #pragma omp parallel firstprivate(Lo)
    {
        #pragma omp for
        for (size_t i = Lo.i_begin; i < Lo.i_endof; ++i)
        for (size_t j = Lo.j_begin; j < Lo.j_endof; ++j)
        for (size_t k = Lo.k_begin; k < Lo.k_endof; ++k)
        {
            t1.p[gridA.icel(i,j,k)] = Mx.X_result[gridA.icelCal(i,j,k)];
        }
    }
}






void Pressure_transform_X_result_Dir(
    pressure &t1, 
    MxClass &Mx ,
    calDomain& Lo,
    grid & gA
){
    #pragma omp parallel firstprivate(Lo)
    {
        #pragma omp for
        for (size_t i = Lo.i_begin; i < Lo.i_endof; ++i)
        for (size_t j = Lo.j_begin; j < Lo.j_endof; ++j)
        for (size_t k = Lo.k_begin; k < Lo.k_endof; ++k)
        {
            t1.p[gA.icel(i,j,k)] = Mx.X_result[gA.icelDir(i,j,k)];
        }
    }
}







std::tuple<double , double> getMax(
    const std::vector<double> &x,
    const calDomain& Lo,
    grid & gA
){

    double maxVal = - 1.0e10;
    double minVal = 1.0e10;

    #pragma omp parallel firstprivate(Lo) reduction(max : maxVal) reduction (min : minVal)
    {
        #pragma omp for
        for (auto i = Lo.i_begin; i < Lo.i_endof; ++i)
        for (auto j = Lo.j_begin; j < Lo.j_endof; ++j)
        for (auto k = Lo.k_begin; k < Lo.k_endof; ++k)
        {
            if (maxVal < x[gA.icelCal(i,j,k)]) 
            { maxVal =  x[gA.icelCal(i,j,k)];}


            if (minVal > x[gA.icelCal(i,j,k)]) 
            { minVal =  x[gA.icelCal(i,j,k)];}
        }
    }

    return std::make_pair(maxVal, minVal);

}



#ifdef EIGEN_ON

void Pressure_transform_x_Eigen(
    pressure &t1, 
    MxClass &Mx ,
    calDomain& Lo,
    grid & gridA
){
    const auto [nx, ny , nz] = gridA.nxyz;

    for (size_t i = Lo.i_begin; i < Lo.i_endof; ++i)
    for (size_t j = Lo.j_begin; j < Lo.j_endof; ++j)
    for (size_t k = Lo.k_begin; k < Lo.k_endof; ++k)
    {
        t1.p[gridA.icel(i,j,k)] = Mx.x_Eigen[gridA.icelCal(i,j,k)];
    }
}

#endif