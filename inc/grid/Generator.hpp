#pragma once

#include"controlPanel.hpp"


#include <chrono>
#include <vector>


bool firtTouch(
    grid& gA,
    calDomain& Lo
);

void generateGride(
    simuClass& simu,
    shareMenory& ShareM,
    calDomain& Lo,
    grid& gA
){

    auto [nx, ny, nz, gC] = gA.nxyzgC;
    gA.resize();

    if (simu.Locality == 1){firtTouch(gA, Lo); }

    // * ---------------------  generator ---------------
    if(gA.Gridder == "non_uniform-Part"){}
    else if(gA.Gridder == "non_uniform-SML")
    { gA.init_nonuniform_SML();}
    else if(gA.Gridder == "uniform")
    { gA.init_uniform(); }
    else{ throw std::runtime_error("grid ?"); }

    // //*  2 preparation for IO
    gA.s_init();
    #pragma omp parallel for schedule(static)
    for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
    for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
    for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
    {
        gA.initNb(i,j,k);
    }


    // * 3 preparation 
    #pragma omp parallel for schedule(static)
    for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
    for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
    for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
    {
        gA.delta(i,j,k) = gA.Dx[i] * gA.Dy[j] * gA.Dz[k];
    }


}


void generateGride(
    simuClass& simu,
    calDomain& Lo,
    grid& gA
){

    auto [nx, ny, nz, gC] = gA.nxyzgC;
    gA.resize();

    if (simu.Locality == 1){firtTouch(gA, Lo); }

    // * ---------------------  generator ---------------
    if(gA.Gridder == "non_uniform-Part"){}
    else if(gA.Gridder == "non_uniform-SML")
    { gA.init_nonuniform_SML();}
    else if(gA.Gridder == "uniform")
    { gA.init_uniform(); }
    else{ throw std::runtime_error("grid ?"); }

    // //*  2 preparation for IO

    gA.s_init();

    #pragma omp parallel for schedule(static)
    for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
    for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
    for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
    {
        gA.initNb(i,j,k);
    }

}




bool firtTouch(
    grid& gA,
    calDomain& Lo
){
    
    std::vector<int> ifomp = {1,0,0}; 

    auto [nx, ny, nz, gC] = gA.nxyzgC;

    #pragma omp parallel for schedule(static) if (ifomp[0])
    for (int i = Lo.i_begin ; i < Lo.i_endof ; ++i ){
        gA.X.at(i) = 0.0;  
        gA.Xc.at(i-gC) = 0.0; 
    }

    #pragma omp parallel for schedule(static) if (ifomp[1])
    for (size_t j = Lo.j_begin ; j < Lo.j_endof ; ++j ){
        gA.Y.at(j) = 0.0;
        gA.Yc.at(j-gC) = 0.0;
    }

    #pragma omp parallel for schedule(static) if (ifomp[2])
    for (int k = Lo.k_begin ; k < Lo.k_endof ; ++k ){
        gA.Z.at(k) = 0.0 ; 
        gA.Zc.at(k-gC)  = 0.0;   
    }


    // ! first touch -----------------
    #pragma omp parallel for schedule(static) if (ifomp[0])
    for (int i = Lo.i_begin ; i < Lo.i_endof ; ++i )
        gA.Dx.at(i-gC) = 0.0;

    #pragma omp parallel for schedule(static) if (ifomp[1])
    for (size_t j = Lo.j_begin ; j < Lo.j_endof ; ++j )
        gA.Dy.at(j-gC) = 0.0;

    #pragma omp parallel for schedule(static) if (ifomp[2])
    for (int k = Lo.k_begin ; k < Lo.k_endof ; ++k )
        gA.Dz.at(k-gC) = 0.0;


    #pragma omp parallel for schedule(static) if (ifomp[0])
    for (int i = Lo.i_begin ; i < Lo.i_endof ; ++i )
        gA.Dxs.at(i-gC) = 0.0;

    #pragma omp parallel for schedule(static) if (ifomp[1])
    for (size_t j = Lo.j_begin ; j < Lo.j_endof ; ++j )
        gA.Dys.at(j-gC) = 0.0;

    #pragma omp parallel for schedule(static) if (ifomp[2])
    for (int k = Lo.k_begin ; k < Lo.k_endof ; ++k )
        gA.Dzs.at(k-gC) = 0.0;

    return true;

}