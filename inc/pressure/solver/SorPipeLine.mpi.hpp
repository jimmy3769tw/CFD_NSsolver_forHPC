#pragma once
# include "SorPipeLine.top.hpp"
# include "mpi_tool/init_mpi.hpp"



void SorPipeLine_mpi_hybrid_omp(
    MPI_Comm comm_world,
    std::vector<int> &mpi_neighborhood,
    SORcoefficient& Sor,
    simuClass& simu,
    velocity& T1,
    pressure& t1,
    calDomain& Lo,
    grid& gA
){

    const int word_size = Lo.word_size;
    double omega = simu.p_sor_omega;
    double pNEW;
    double pChange;
    double mChange;
    int itmax = simu.p_sor_iter_max;
    const auto nx = gA.nx;
    const auto ny = gA.ny;
    const auto nz = gA.nz;
    const auto gC = gA.gC;

    double reverse_dt = 1.0 / simu.dt;
    // cout << "Present - zeta: " <<  simu.p_criteria << ", itmax: "  << itmax << ", omega: " << omega << ", dt: " << simu.dt << endl;

    double mChangeMax;
    double pChangeMax;

    if (simu.loop == 1){
        //! ========================= if (main_count == 0) =========================
        Sor.cf.resize(gA.iceltotCal * 8);

        for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
        for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
        for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
        {
            const int ii =  8*gA.icelCal(i,j,k);
            Sor.cf[ii  ] = gA.Dy[j] * gA.Dz[k] / gA.Dxs[i]  ;
            Sor.cf[ii+1] = gA.Dy[j] * gA.Dz[k] / gA.Dxs[i-1];
            Sor.cf[ii+2] = gA.Dx[i] * gA.Dz[k] / gA.Dys[j]  ;
            Sor.cf[ii+3] = gA.Dx[i] * gA.Dz[k] / gA.Dys[j-1];
            Sor.cf[ii+4] = gA.Dx[i] * gA.Dy[j] / gA.Dzs[k]  ;
            Sor.cf[ii+5] = gA.Dx[i] * gA.Dy[j] / gA.Dzs[k-1];

            Sor.cf[ii+7]   =  1.0 / (
                        - gA.Dy[j] * gA.Dz[k] / gA.Dxs[i] 
                        - gA.Dy[j] * gA.Dz[k] / gA.Dxs[i-1]
                        - gA.Dx[i] * gA.Dz[k] / gA.Dys[j] 
                        - gA.Dx[i] * gA.Dz[k] / gA.Dys[j-1]
                        - gA.Dx[i] * gA.Dy[j] / gA.Dzs[k] 
                        - gA.Dx[i] * gA.Dy[j] / gA.Dzs[k-1] );

        }
    }  


    for (int i = Lo.i_begin ; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin ; j < Lo.j_endof ; ++j )
    for (int k = Lo.k_begin ; k < Lo.k_endof ; ++k )
    {
        const int  icel = gA.icel(i,j,k);
        mChange = ( T1.u[icel]- T1.u[gA.NEIBcell[icel*12 +0]] ) * gA.Dy[j] * gA.Dz[k] 
                + ( T1.v[icel]- T1.v[gA.NEIBcell[icel*12 +2]] ) * gA.Dx[i] * gA.Dz[k] 
                + ( T1.w[icel]- T1.w[gA.NEIBcell[icel*12 +4]] ) * gA.Dx[i] * gA.Dy[j] ; //FIXME"
        Sor.cf[gA.icelCal(i,j,k)*8+6] = mChange * reverse_dt;
    }
        // should slow ~~~
    simu.iters = 0;
    double interTime = omp_get_wtime();

    // * --------- setting Peridic and Neunann --------- 

    stopWatch time_ ;
    time_.start();
    // * main loop
    double residual;
    auto recvbuf = residual;
    do {

        // ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~  MPI no-blocking send & recv ~~~~~~~~~~~~~~~~~~
        #ifdef MPI_DEBUG
        mpi_iSR_double_x_debugger(word_size, simu.TID, comm_world, mpi_neighborhood, t1.p, Lo, gA);
        #else
        mpi_iSR_double_x(1, comm_world, mpi_neighborhood, t1.p, Lo, gA );
        #endif
        // ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~  MPI no-blocking send & recv ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MPI_Barrier(comm_world);

        residual = 0.0;

        for (int i = Lo.i_begin ; i < Lo.i_endof ; ++i )
        for (size_t j = Lo.j_begin ; j < Lo.j_endof ; ++j )
        for (int k = Lo.k_begin ; k < Lo.k_endof ; ++k )
        {
            const int icel = gA.icel(i,j,k);
            const int icel_12 = icel*12;
            const int icel_08 = gA.icelCal(i,j,k)*8;
            pNEW = (- t1.p[gA.NEIBcell[icel_12+1]] * Sor.cf[icel_08  ]
                    - t1.p[gA.NEIBcell[icel_12  ]] * Sor.cf[icel_08+1]
                    - t1.p[gA.NEIBcell[icel_12+3]] * Sor.cf[icel_08+2]
                    - t1.p[gA.NEIBcell[icel_12+2]] * Sor.cf[icel_08+3] 
                    - t1.p[gA.NEIBcell[icel_12+5]] * Sor.cf[icel_08+4]
                    - t1.p[gA.NEIBcell[icel_12+4]] * Sor.cf[icel_08+5] 
                    + Sor.cf[icel_08+6]) * Sor.cf[icel_08+7];

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

        auto sendbuf = residual;

        #ifdef SOR_ERROR_MAXVAL

            MPI_Allreduce(&sendbuf, &recvbuf,1, MPI_DOUBLE, MPI_MAX, comm_world);

        #endif

        #ifdef SOR_ERROR_SUM
            MPI_Allreduce(&sendbuf, &recvbuf,1, MPI_DOUBLE, MPI_SUM, comm_world);
        #endif 


        #ifdef BC_P_PERIODIC_Z
            BcPressurePeriodicZ(simu,Lo,t1,gA);
        #endif
    
        ++simu.iters;
        // cout << ", " << recvbuf ;
    }while(recvbuf > simu.p_criteria && simu.iters < itmax);

    // https://en.wikipedia.org/wiki/Errors_and_residuals
    simu.error = recvbuf;

    --simu.iters;
}