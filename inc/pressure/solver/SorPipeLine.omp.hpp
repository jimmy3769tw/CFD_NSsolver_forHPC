#pragma once
# include "controlPanel.hpp"


#include "omp.h"

bool SorPipeLine_omp(
    SORcoefficient& Sor,
    simuClass& simu,
    velocity& T1,
    pressure& t1,
    calDomain& Lo,
    grid& gA
)
{
    // -------------------------
    double omega = simu.p_sor_omega;

    int itmax = simu.p_sor_iter_max;

    const auto [nx, ny, nz, gC] = gA.nxyzgC;

    double reverse_dt = 1.0 / simu.dt;

    double mChangeMax;

    double rdt = 1.0 / simu.dt;

    if (simu.firstSOR) {

        Sor.cf.resize(gA.iceltotCal * 8);

        //! ========================= if ( A matrix) =========================
        #pragma omp parallel for firstprivate(Lo)
        for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
        for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
        for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
        {
            const int ii =  gA.icelCal(i,j,k);
            Sor.cf.at(ii*8  ) = gA.Dy[j] * gA.Dz[k] / gA.Dxs[i]  ;
            Sor.cf.at(ii*8+1) = gA.Dy[j] * gA.Dz[k] / gA.Dxs[i-1];

            Sor.cf.at(ii*8+2) = gA.Dx[i] * gA.Dz[k] / gA.Dys[j]  ;
            Sor.cf.at(ii*8+3) = gA.Dx[i] * gA.Dz[k] / gA.Dys[j-1];

            Sor.cf.at(ii*8+4) = gA.Dx[i] * gA.Dy[j] / gA.Dzs[k]  ;
            Sor.cf.at(ii*8+5) = gA.Dx[i] * gA.Dy[j] / gA.Dzs[k-1];

            Sor.cf.at(ii*8+7) = - 1.0 / ( gA.coef(i,j,k));
        }
        simu.firstSOR == false;
    }  


    double norm_rhs = 0.0;

    #pragma omp parallel for firstprivate(Lo, rdt)
    for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
    for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
    for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
    {
        const int icel = gA.icel(i,j,k);
        const int ii8 = gA.icelCal(i,j,k)*8;

        auto mChange = ( T1.u[icel] - T1.u[gA.icel(i-1,j,k)] ) * gA.Dy[j] * gA.Dz[k] 
                     + ( T1.v[icel] - T1.v[gA.icel(i,j-1,k)] ) * gA.Dx[i] * gA.Dz[k] 
                     + ( T1.w[icel] - T1.w[gA.icel(i,j,k-1)] ) * gA.Dx[i] * gA.Dy[j] ;

        Sor.cf[ii8+6] = mChange * rdt;
        norm_rhs += std::pow( mChange * rdt, 2);
    }

    norm_rhs = sqrt(norm_rhs);

    

    double interTime = omp_get_wtime();

    // * --------- setting Peridic and Neunann --------- 

    double residual;
    for (simu.iters = 0;  simu.iters < simu.p_sor_iter_max; ++simu.iters)
    {
        residual = 0;
        #pragma omp parallel for reduction(+:residual) firstprivate(Lo, norm_rhs)
        for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
        for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
        for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
        {
            const int icel = gA.icel(i,j,k);
            const int ii8 = gA.icelCal(i,j,k)*8;

            auto pNEW =(- t1.p[gA.icel(i+1,j,k)] * Sor.cf[ii8  ]
                        - t1.p[gA.icel(i-1,j,k)] * Sor.cf[ii8+1]
                        - t1.p[gA.icel(i,j+1,k)] * Sor.cf[ii8+2]
                        - t1.p[gA.icel(i,j-1,k)] * Sor.cf[ii8+3]
                        - t1.p[gA.icel(i,j,k+1)] * Sor.cf[ii8+4]
                        - t1.p[gA.icel(i,j,k-1)] * Sor.cf[ii8+5]
                + Sor.cf[ii8+6]) * Sor.cf[ii8+7];

            double pChange = std::abs(pNEW - t1.p[icel]);

            t1.p[icel] += omega * (pNEW - t1.p[icel]);
            residual += pChange/norm_rhs;
        }

        if (sqrt(residual) < simu.p_criteria){break;}

        // cout  << simu.iters << "\t residual : " << sqrt(residual)  << ", " <<  simu.p_criteria << endl;
    }

    simu.error = sqrt(residual);



    return true;

}



void SorPipeLine_seq_unpeeling_omp(
    SORcoefficient& Sor,
    simuClass& simu,
    velocity& T1,
    pressure& t1,
    calDomain& Lo,
    grid& gA
)
{


    // -------------------------
    double omega = simu.p_sor_omega;

    int itmax = simu.p_sor_iter_max;

    const auto [nx, ny, nz, gC] = gA.nxyzgC;

    double reverse_dt = 1.0 / simu.dt;

    double mChangeMax;

    double rdt = 1.0 / simu.dt;

    static bool first = true;

    if (first ){

        Sor.cf.resize(gA.iceltotCal * 8);

        //! ========================= if ( A matrix) =========================
        #pragma omp parallel for firstprivate(Lo)
        for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
        for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
        for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
        {
            const int ii =  gA.icelCal(i,j,k);
            Sor.cf.at(ii*8  ) = gA.Dy[j] * gA.Dz[k] / gA.Dxs[i]  ;
            Sor.cf.at(ii*8+1) = gA.Dy[j] * gA.Dz[k] / gA.Dxs[i-1];

            Sor.cf.at(ii*8+2) = gA.Dx[i] * gA.Dz[k] / gA.Dys[j]  ;
            Sor.cf.at(ii*8+3) = gA.Dx[i] * gA.Dz[k] / gA.Dys[j-1];

            Sor.cf.at(ii*8+4) = gA.Dx[i] * gA.Dy[j] / gA.Dzs[k]  ;
            Sor.cf.at(ii*8+5) = gA.Dx[i] * gA.Dy[j] / gA.Dzs[k-1];

            Sor.cf.at(ii*8+7) = - 1.0 / ( gA.coef(i,j,k));
        }
    }  

    first = false;

    #pragma omp parallel for firstprivate(Lo, rdt)
    for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
    for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
    for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
    {
        const int icel = gA.icel(i,j,k);
        const int ii8 = gA.icelCal(i,j,k)*8;

        auto mChange = ( T1.u[icel] - T1.u[gA.icel(i-1,j,k)] ) * gA.Dy[j] * gA.Dz[k] 
                     + ( T1.v[icel] - T1.v[gA.icel(i,j-1,k)] ) * gA.Dx[i] * gA.Dz[k] 
                     + ( T1.w[icel] - T1.w[gA.icel(i,j,k-1)] ) * gA.Dx[i] * gA.Dy[j] ;

        Sor.cf[ii8+6] = mChange * rdt;
    }

    simu.iters = 0;

    double interTime = omp_get_wtime();

    // * --------- setting Peridic and Neunann --------- 

    stopWatch time_ ;
    time_.start();
    // * main loop
    double residual;
    do {
        
        residual = 0.0;

        #pragma omp parallel for reduction(+:residual) firstprivate(Lo)
        for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
        for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
        for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
        {
            const auto icel = gA.icel(i,j,k);
            const auto ii8  = gA.icelCal(i,j,k)*8;

            auto pNEW =(- t1.p[gA.icel(i+1,j,k)] * Sor.cf[ii8  ]
                        - t1.p[gA.icel(i-1,j,k)] * Sor.cf[ii8+1]
                        - t1.p[gA.icel(i,j+1,k)] * Sor.cf[ii8+2]
                        - t1.p[gA.icel(i,j-1,k)] * Sor.cf[ii8+3] 
                        - t1.p[gA.icel(i,j,k+1)] * Sor.cf[ii8+4]
                        - t1.p[gA.icel(i,j,k-1)] * Sor.cf[ii8+5]
                + Sor.cf[ii8+6]) * Sor.cf[ii8+7];

            double pChange = std::abs(pNEW - t1.p[icel]);

            t1.p[icel] += omega * (pNEW - t1.p[icel]);

            residual += pChange;

        }

        ++simu.iters;

    }while(residual > simu.p_criteria && simu.iters < itmax);

    simu.error = residual;

    --simu.iters;

    time_.stop();
}