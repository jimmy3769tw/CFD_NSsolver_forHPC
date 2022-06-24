#pragma once


#include"controlPanel.hpp"


void CheckSteadyStatebyL2norm(
    simuClass& simu,
    const velocity& T0,
    const velocity& T3,
    calDomain& Lo,
    grid& gridA
){
    auto L2norm = CheckL2Norm(simu, T0, T3, Lo,gridA);
    #pragma omp critical
    {
        cout  << std::flush << "TID=" <<simu.TID << ", L2norm=" << L2norm << endl;
    }
}


void CheckSteadyStatebyMaxVal_omp(
    simuClass& simu,
    shareMenory& ShareM,
    const velocity& T0,
    const velocity& T3,
    const calDomain& Lo,
    grid& gridA
)
{   
    double uDif_Max_ = 0;
    double vDif_Max_ = 0;
    double wDif_Max_ = 0;
    const size_t nx = gridA.nx;
    const size_t ny = gridA.ny;
    const size_t nz = gridA.nz;
    const size_t gC = gridA.gC;

    for (size_t i = Lo.i_begin; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin; j < Lo.j_endof ; ++j )
    for (size_t k = Lo.k_begin; k < Lo.k_endof ; ++k )
    {
        const int icel = i*nz*ny + j*nz + k ;
        const double  Dif = (T0.u[icel] - T3.u[icel]);
        if (Dif > uDif_Max_) {
            uDif_Max_= Dif;
        }
    }
    
    for (size_t i = Lo.i_begin; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin; j < Lo.j_endof ; ++j )
    for (size_t k = Lo.k_begin; k < Lo.k_endof ; ++k )
    {
        const int icel = i*nz*ny + j*nz + k ;
        const double  Dif = (T0.v[icel] - T3.v[icel]);
        if (Dif > vDif_Max_) {
            vDif_Max_= Dif;
        }         
    }

    for (size_t i = Lo.i_begin; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin; j < Lo.j_endof ; ++j )
    for (size_t k = Lo.k_begin; k < Lo.k_endof ; ++k )
    {
        const int icel = i*nz*ny + j*nz + k ;
        const double  Dif = (T0.w[icel] - T3.w[icel]);
        if (Dif > wDif_Max_) {
            wDif_Max_= Dif;
        }         
    }


    ShareM.uDif_Max = uDif_Max_;
    ShareM.vDif_Max = vDif_Max_;
    ShareM.wDif_Max = wDif_Max_;

    #pragma omp barrier
    if (uDif_Max_ > ShareM.Out)
    {
        {
            if (uDif_Max_ > ShareM.uDif_Max)
                ShareM.uDif_Max = uDif_Max_;
        }
    }

    if (vDif_Max_ > ShareM.Out)
    {        
        #pragma omp critical
        {
            if (vDif_Max_ > ShareM.vDif_Max)
                ShareM.vDif_Max = vDif_Max_;
        }
    }

    if (wDif_Max_ > ShareM.Out)
    {
        #pragma omp critical
        {
            if (wDif_Max_ > ShareM.wDif_Max)
                ShareM.wDif_Max = wDif_Max_;
        }
    }

    #pragma omp single
    {
        cout << "uDif_Max: " << ShareM.uDif_Max ;
        cout << ", vDif_Max: " << ShareM.vDif_Max ;
        cout << ", wDif_Max: " << ShareM.wDif_Max << endl;
    }

}

std::tuple<double, double, double> CheckSteadyStatebyMaxVal_seq(
    simuClass& simu,
    shareMenory& ShareM,
    const velocity& T0,
    const velocity& T3,
    calDomain& Lo,
    grid& gridA
)
{   
    double uDif_Max_ = 0.0;
    double vDif_Max_ = 0.0;
    double wDif_Max_ = 0.0;
    const auto [nx, ny, nz, gC] = gridA.nxyzgC;

    for (size_t i = Lo.i_begin; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin; j < Lo.j_endof ; ++j )
    for (size_t k = Lo.k_begin; k < Lo.k_endof ; ++k )
    {
        const int icel = i*nz*ny + j*nz + k ;
        const double  Dif = (T0.u[icel] - T3.u[icel]);
        if (Dif > uDif_Max_) {
            uDif_Max_= Dif;
        }
    }
    
    for (size_t i = Lo.i_begin; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin; j < Lo.j_endof ; ++j )
    for (size_t k = Lo.k_begin; k < Lo.k_endof ; ++k )
    {
        const int icel = i*nz*ny + j*nz + k ;
        const double  Dif = (T0.v[icel] - T3.v[icel]);
        if (Dif > vDif_Max_) {
            vDif_Max_= Dif;
        }         
    }

    for (size_t i = Lo.i_begin; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin; j < Lo.j_endof ; ++j )
    for (size_t k = Lo.k_begin; k < Lo.k_endof ; ++k )
    {
        const int icel = i*nz*ny + j*nz + k ;
        const double  Dif = (T0.w[icel] - T3.w[icel]);
        if (Dif > wDif_Max_) {
            wDif_Max_= Dif;
        }         
    }

    return std::make_tuple(uDif_Max_, vDif_Max_, wDif_Max_);
}



std::tuple<double, double, double> 
CheckSteadyStatebyMaxVal_omp(
    simuClass& simu,
    shareMenory& ShareM,
    velocity& T0,
    velocity& T3,
    calDomain& Lo,
    grid& gA
)
{   
    double uDif_Max_ = 0.0;
    double vDif_Max_ = 0.0;
    double wDif_Max_ = 0.0;
    auto [nx, ny, nz, gC] = gA.nxyzgC;

    #pragma omp parallel for reduction(max : uDif_Max_)
    for (size_t i = Lo.i_begin; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin; j < Lo.j_endof ; ++j )
    for (size_t k = Lo.k_begin; k < Lo.k_endof ; ++k )
    {
        const int icel = gA.icel(i,j,k);
        const double  Dif = (T0.u[icel] - T3.u[icel]);
        if (Dif > uDif_Max_) {
            uDif_Max_= Dif;
        }
    }
    
    #pragma omp parallel for reduction(max : vDif_Max_)
    for (size_t i = Lo.i_begin; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin; j < Lo.j_endof ; ++j )
    for (size_t k = Lo.k_begin; k < Lo.k_endof ; ++k )
    {
        const int icel = gA.icel(i,j,k);
        const double  Dif = (T0.v[icel] - T3.v[icel]);
        if (Dif > vDif_Max_) {
            vDif_Max_= Dif;
        }         
    }


    #pragma omp parallel for reduction(max : wDif_Max_)
    for (size_t i = Lo.i_begin; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin; j < Lo.j_endof ; ++j )
    for (size_t k = Lo.k_begin; k < Lo.k_endof ; ++k )
    {
        const int icel = gA.icel(i,j,k);
        const double  Dif = (T0.w[icel] - T3.w[icel]);
        if (Dif > wDif_Max_) {
            wDif_Max_= Dif;
        }         
    }


    return std::make_tuple(uDif_Max_, vDif_Max_, wDif_Max_);
}





#ifdef MPI_ON

std::tuple<double, double, double> 
CheckSteadyStatebyMaxVal_mpi(
    MPI_Comm &comm,
    simuClass& simu,
    shareMenory& ShareM,
    const velocity& T0,
    const velocity& T3,
    calDomain& Lo,
    grid& gridA
)
{

    double Max_[3] = {0.0};

    const auto [nx, ny, nz, gC] = gridA.nxyzgC;

    for (size_t i = Lo.i_begin; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin; j < Lo.j_endof ; ++j )
    for (size_t k = Lo.k_begin; k < Lo.k_endof ; ++k )
    {
        const int icel = i*nz*ny + j*nz + k ;
        const double  Dif = (T0.u[icel] - T3.u[icel]);
        if (Dif > Max_[0]) {
            Max_[0]= Dif;
        }
    }
    
    for (size_t i = Lo.i_begin; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin; j < Lo.j_endof ; ++j )
    for (size_t k = Lo.k_begin; k < Lo.k_endof ; ++k )
    {
        const int icel = i*nz*ny + j*nz + k ;
        const double  Dif = (T0.v[icel] - T3.v[icel]);
        if (Dif > Max_[1]) {
            Max_[1]= Dif;
        }         
    }

    for (size_t i = Lo.i_begin; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin; j < Lo.j_endof ; ++j )
    for (size_t k = Lo.k_begin; k < Lo.k_endof ; ++k )
    {
        const int icel = i*nz*ny + j*nz + k ;
        const double  Dif = (T0.w[icel] - T3.w[icel]);
        if (Dif > Max_[2]) {
            Max_[2]= Dif;
        }         
    }
    double  Max[3];

    MPI_Allreduce(Max_, Max, 3, MPI_DOUBLE, MPI_SUM, comm);
    return std::make_tuple(uDif_Max_, vDif_Max_, wDif_Max_);
}

#endif
