#pragma once
# include "SorPipeLine.top.hpp"

bool SorPipeLine_seq(
    SORcoefficient& Sor,
    simuClass& simu,
    velocity& T1,
    pressure& t1,
    calDomain& Lo,
    grid& gA
)
{
    double omega = simu.p_sor_omega;

    double pNEW;

    int itmax = simu.p_sor_iter_max;

    const auto [nx, ny, nz, gC] = gA.nxyzgC;


    double rdt = 1.0 / simu.dt;
    // cout << "Present - zeta: " <<  simu.p_criteria << ", itmax: "  << itmax << ", omega: " << omega << ", dt: " << simu.dt << endl;


    if (simu.loop == 1){
        //! ========================= if (main_count == 0) =========================
        Sor.cf.resize(gA.iceltotCal * 8);

        for (int i = Lo.i_begin ; i < Lo.i_endof ; ++i )
        for (size_t j = Lo.j_begin ; j < Lo.j_endof ; ++j )
        for (int k = Lo.k_begin ; k < Lo.k_endof ; ++k )
        {
            const int ii =  gA.icelCal(i,j,k);
            Sor.cf.at(ii*8  ) = gA.Dy[j] * gA.Dz[k] / gA.Dxs[i]  ;
            Sor.cf.at(ii*8+1) = gA.Dy[j] * gA.Dz[k] / gA.Dxs[i-1];
            Sor.cf.at(ii*8+2) = gA.Dx[i] * gA.Dz[k] / gA.Dys[j]  ;
            Sor.cf.at(ii*8+3) = gA.Dx[i] * gA.Dz[k] / gA.Dys[j-1];
            Sor.cf.at(ii*8+4) = gA.Dx[i] * gA.Dy[j] / gA.Dzs[k]  ;
            Sor.cf.at(ii*8+5) = gA.Dx[i] * gA.Dy[j] / gA.Dzs[k-1];
            Sor.cf.at(ii*8+7)   = - 1.0 / ( gA.coef(i,j,k));

        }


    }  

    for (int i = Lo.i_begin ; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin ; j < Lo.j_endof ; ++j )
    for (int k = Lo.k_begin ; k < Lo.k_endof ; ++k )
    {
        const int icel = gA.icel(i,j,k);
        const int ii8 = gA.icelCal(i,j,k)*8;

        auto mChange = ( T1.u[icel]- T1.u[gA.icel(i-1,j,k)] ) * gA.Dy[j] * gA.Dz[k] 
                + ( T1.v[icel]- T1.v[gA.icel(i,j-1,k)] ) * gA.Dx[i] * gA.Dz[k] 
                + ( T1.w[icel]- T1.w[gA.icel(i,j,k-1)] ) * gA.Dx[i] * gA.Dy[j] ; //FIXME"

        Sor.cf[ii8+6] = mChange * rdt;
    }
        // should slow ~~~
    simu.iters = 0;
    double interTime = omp_get_wtime();

    // * --------- setting Peridic and Neunann --------- 
    double pChange;

    // * main loop
    double residual;
    do {
        
        residual = 0.0;

        for (int i = Lo.i_begin ; i < Lo.i_endof ; ++i )
        for (size_t j = Lo.j_begin ; j < Lo.j_endof ; ++j )
        for (int k = Lo.k_begin ; k < Lo.k_endof ; ++k )
        {
            const int icel = gA.icel(i,j,k);
            const int ii8 = gA.icelCal(i,j,k)*8;
        
            pNEW = (- t1.p[gA.icel(i+1,j,k)] * Sor.cf[ii8  ]
                    - t1.p[gA.icel(i-1,j,k)] * Sor.cf[ii8+1]
                    - t1.p[gA.icel(i,j+1,k)] * Sor.cf[ii8+2]
                    - t1.p[gA.icel(i,j-1,k)] * Sor.cf[ii8+3] 
                    - t1.p[gA.icel(i,j,k+1)] * Sor.cf[ii8+4]
                    - t1.p[gA.icel(i,j,k-1)] * Sor.cf[ii8+5]

                    + Sor.cf[ii8+6]) * Sor.cf[ii8+7];

            pChange = std::abs(pNEW - t1.p[icel]);

            t1.p[icel] += omega * (pNEW - t1.p[icel]);


            #ifdef SOR_ERROR_MAXVAL
                if( pChange > residual ){
                    residual = pChange;
                }
            #endif

            #ifdef SOR_ERROR_SUM
                residual += pChange;
            #endif
        }

        #ifdef BC_P_PERIODIC_Z
            BcPressurePeriodicZ(simu,Lo,t1,gridA);
        #endif
    
        ++simu.iters;

    }while(residual > simu.p_criteria && simu.iters < itmax);

    // https://en.wikipedia.org/wiki/Errors_and_residuals
    simu.error = residual;

    --simu.iters;


    return true;

}




void SorPipeLine_seq_unpeeling(
    shareMenory& ShareM,
    SORcoefficient& Sor,
    simuClass& simu,
    velocity& T1,
    pressure& t1,
    calDomain& Lo,
    grid& gA
)
{
    double omega = 1.8;

    int itmax = 3000;

    double dt = simu.dt;

    double zeta = simu.p_criteria;

    const auto [nx, ny, nz, gC] = gA.nxyzgC;

    int ik = 0;


    double pChangeMax;
    
    do{
        ik ++;
        pChangeMax = 0.0;

        for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
        for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
        for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
        {
            auto icel = gA.icel(i,j,k);

            auto mChange = 
                      ( T1.u[icel] - T1.u[gA.icel(i-1,j,k)] ) * gA.Dy[j] * gA.Dz[k] \
                    + ( T1.v[icel] - T1.v[gA.icel(i,j-1,k)] ) * gA.Dx[i] * gA.Dz[k] \
                    + ( T1.w[icel] - T1.w[gA.icel(i,j,k-1)] ) * gA.Dx[i] * gA.Dy[j];

            auto pNEW = 
                   (- t1.p[gA.icel(i+1,j,k)] * gA.Dy[j] * gA.Dz[k] / gA.Dxs[i] 
                    - t1.p[gA.icel(i-1,j,k)] * gA.Dy[j] * gA.Dz[k] / gA.Dxs[i-1] 
                    - t1.p[gA.icel(i,j+1,k)] * gA.Dx[i] * gA.Dz[k] / gA.Dys[j] 
                    - t1.p[gA.icel(i,j-1,k)] * gA.Dx[i] * gA.Dz[k] / gA.Dys[j-1] 
                    - t1.p[gA.icel(i,j,k+1)] * gA.Dx[i] * gA.Dy[j] / gA.Dzs[k] 
                    - t1.p[gA.icel(i,j,k-1)] * gA.Dx[i] * gA.Dy[j] / gA.Dzs[k-1] 
                    + mChange / dt) / (- gA.coef(i,j,k)) ;

            auto pChange = std::abs(pNEW - t1.p[icel]);

            t1.p[icel] += omega * (pNEW - t1.p[icel]);
            
            if(pChange > pChangeMax)
                pChangeMax = pChange;

        }
    }while(pChangeMax>zeta && ik <itmax);
    --ik;
}
