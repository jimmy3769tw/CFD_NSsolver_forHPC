#pragma once
#include "controlPanel.hpp"

void createBMatrix_seq (
    velocity& T1,
    MxClass& Mx,
    simuClass& simu,
    calDomain& Lo,
    grid& gA
    )
{

    Mx.matB.resize(gA.iceltotCal);
    double ddt = 1.0 / simu.dt;
    const auto [nx, ny, nz] = gA.nxyz;

    #if defined(JACOBI_PC)

    for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
    for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
    for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
    {
        const auto icel = gA.icel(i,j,k);
        const auto row = gA.icelCal(i,j,k);

        Mx.matB[row] = -gA.jp[row]*ddt*( 
                        ( T1.u[icel] - T1.u[gA.icel(i-1,j,k)] )/gA.Dx[i]+
                        ( T1.v[icel] - T1.v[gA.icel(i,j-1,k)] )/gA.Dy[j]+
                        ( T1.w[icel] - T1.w[gA.icel(i,j,k-1)] )/gA.Dz[k] );
    }

    #else

    for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
    for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
    for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
    {
        const auto icel = gA.icel(i,j,k);
        const auto row = gA.icelCal(i,j,k);

        Mx.matB[row] = -ddt*( 
                        ( T1.u[icel] - T1.u[gA.icel(i-1,j,k)] )/gA.Dx[i]+
                        ( T1.v[icel] - T1.v[gA.icel(i,j-1,k)] )/gA.Dy[j]+
                        ( T1.w[icel] - T1.w[gA.icel(i,j,k-1)] )/gA.Dz[k] );
    }
    #endif
}



void createBMatrix_omp (
    velocity& T1,
    MxClass& Mx,
    simuClass& simu,
    calDomain& Lo,
    grid& gA
    )
{

    Mx.matB.resize(gA.iceltotCal);
    #pragma omp parallel firstprivate(Lo)
    {

        double ddt = 1.0 / simu.dt;
        const auto [nx, ny, nz] = gA.nxyz;

        #if defined(JACOBI_PC)
        #pragma omp for  
        for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
        for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
        for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
        {
            const auto icel = gA.icel(i,j,k);
            const auto row = gA.icelCal(i,j,k);

            Mx.matB[row] = -gA.jp[row]*ddt*( 
                            ( T1.u[icel] - T1.u[gA.icel(i-1,j,k)] )/gA.Dx[i]+
                            ( T1.v[icel] - T1.v[gA.icel(i,j-1,k)] )/gA.Dy[j]+
                            ( T1.w[icel] - T1.w[gA.icel(i,j,k-1)] )/gA.Dz[k] );
        }

        #else

        #pragma omp for  
        for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
        for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
        for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
        {
            const auto icel = gA.icel(i,j,k);
            const auto row = gA.icelCal(i,j,k);

            Mx.matB[row] = -ddt*( 
                            ( T1.u[icel] - T1.u[gA.icel(i-1,j,k)] )/gA.Dx[i]+
                            ( T1.v[icel] - T1.v[gA.icel(i,j-1,k)] )/gA.Dy[j]+
                            ( T1.w[icel] - T1.w[gA.icel(i,j,k-1)] )/gA.Dz[k] );
        }
        #endif



    }
}



// Dir boundary condiaiton 
void createBMatrixDir(
    velocity& T1,
    pressure &t1, 
    MxClass& Mx,
    simuClass& simu,
    calDomain& Lo,
    grid& gA
    )
{
    auto ddt = 1.0 / simu.dt;

    const auto [nx, ny, nz] = gA.nxyz;

    Mx.matB.resize(gA.nznynxDir);

    // define B matrix
    for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
    for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
    for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
    {
        const auto icel = gA.icel(i,j,k);
        const auto row = gA.icelDir(i,j,k);

        Mx.matB.at(row) = -ddt*( 
                        ( T1.u[icel] - T1.u[gA.icel(i-1,j,k)] )/gA.Dx[i]+
                        ( T1.v[icel] - T1.v[gA.icel(i,j-1,k)] )/gA.Dy[j]+
                        ( T1.w[icel] - T1.w[gA.icel(i,j,k-1)] )/gA.Dz[k] );
    }


    int 
     xBeg = gA.gC-1, yBeg = gA.gC-1, zBeg = gA.gC-1,
     xEnd = gA.nxCal+2, yEnd = gA.nyCal+2, zEnd = gA.nzCal+2;
    // define B matrix
    for (auto i = Lo.i_begin-1 ; i < Lo.i_endof+1 ; ++i )
    for (auto j = Lo.j_begin-1 ; j < Lo.j_endof+1 ; ++j )
    for (auto k = Lo.k_begin-1 ; k < Lo.k_endof+1 ; ++k )
    {
        if (i==xBeg || i==xEnd || j==yBeg || j==yEnd || k==zBeg || k == zEnd ){
            Mx.matB.at( gA.icelDir(i,j,k) ) =  t1.p[gA.icel(i,j,k)];
        }
    }

}



#if defined(EIGEN_ON)

void createBMatrix_Eigen (
    velocity& T1,
    MxClass& Mx,
    simuClass& simu,
    calDomain& Lo,
    grid& gA
    )
{

    Mx.matB.resize(gA.iceltotCal);
    #pragma omp parallel firstprivate(Lo)
    {

        double ddt = 1.0 / simu.dt;
        const auto [nx, ny, nz] = gA.nxyz;

        #if defined(JACOBI_PC)
        #pragma omp for  
        for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
        for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
        for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
        {
            const auto icel = gA.icel(i,j,k);
            const auto row = gA.icelCal(i,j,k);

            Mx.matB_Eigen[row] = -gA.jp[row]*ddt*( 
                            ( T1.u[icel] - T1.u[gA.icel(i-1,j,k)] )/gA.Dx[i]+
                            ( T1.v[icel] - T1.v[gA.icel(i,j-1,k)] )/gA.Dy[j]+
                            ( T1.w[icel] - T1.w[gA.icel(i,j,k-1)] )/gA.Dz[k] );
        }

        #else

        #pragma omp for  
        for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
        for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
        for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
        {
            const auto icel = gA.icel(i,j,k);
            const auto row = gA.icelCal(i,j,k);

            Mx.matB_Eigen[row] = -ddt*( 
                            ( T1.u[icel] - T1.u[gA.icel(i-1,j,k)] )/gA.Dx[i]+
                            ( T1.v[icel] - T1.v[gA.icel(i,j-1,k)] )/gA.Dy[j]+
                            ( T1.w[icel] - T1.w[gA.icel(i,j,k-1)] )/gA.Dz[k] );
        }
        #endif



    }
}

#endif
