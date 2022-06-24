#pragma once
#include"controlPanel.hpp"

void update_UandF_seq(
    DfibArray& Dfib,
    const simuClass& simu,
    const pressure& t1,
    const velocity& oldVel,
    velocity& curtVel,
    const calDomain& Lo,
    grid& gA
)
{
    /*
    * update the velocity accordion to the Dfib method 
    * The U solid shoun't be equatal to zero if the solid can move.
    */ 

    double uS = 0.0, vS = 0.0, wS = 0.0;
    const auto [nx, ny, nz] = gA.nxyz;
    const auto dt = simu.dt;

    int d = 0, one = 1;
    auto iceltot = gA.iceltot;

    const int  xshift = 0*gA.iceltotCal;
    const int  yshift = 1*gA.iceltotCal;
    const int  zshift = 2*gA.iceltotCal;

    // !------------------------ In x direction
    if(Lo.i_endof == nx-gA.gC){d = one;}
    else {d = 0;}

    for (auto i = Lo.i_begin ; i < Lo.i_endof-d ; ++i )
    for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
    for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
    {
        // -----------------------
        auto  ic = gA.icel(i,j,k); 
        auto  icP = gA.icel(i+1,j,k);
        auto  T2u = oldVel.u[ic] - dt*( t1.p[icP]-t1.p[ic] ) / gA.Dxs[i];
        auto  eta = Dfib.eta[ic]+(Dfib.eta[icP]-Dfib.eta[ic])*gA.x_INTRPL(i);
        // -----------------------

        // -----------------------
        curtVel.u[ic]= Dfib.eta[ic] * uS + (1.0-eta) * T2u;
        Dfib.f[gA.icelCal(i,j,k)+xshift] = (curtVel.u[ic] - T2u) / dt;
        // -----------------------
    }


    // !------------------------ In y direction
    if(Lo.j_endof == ny-gA.gC){d = one;}
    else {d = 0;}

    for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
    for (auto j = Lo.j_begin ; j < Lo.j_endof-d ; ++j )
    for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
    {
        // -----------------------
        auto  ic = gA.icel(i,j,k); auto icP = gA.icel(i,j+1,k);
        auto  T2v = oldVel.v[ic] - dt*( t1.p[icP]-t1.p[ic] ) / gA.Dys[j];
        auto  eta = Dfib.eta[ic]+(Dfib.eta[icP]-Dfib.eta[ic])*gA.y_INTRPL(j);
        // -----------------------

        // -----------------------
        curtVel.v[ic] = Dfib.eta[ic] * vS + (1.0-eta) * T2v;
        Dfib.f[gA.icelCal(i,j,k)+yshift]  = (curtVel.v[ic] - T2v) / dt;
        // -----------------------
    }

    // !------------------------ In z direction
    if(Lo.k_endof == nz-gA.gC){d = one;}
    else {d = 0;}

    for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
    for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
    for (auto k = Lo.k_begin ; k < Lo.k_endof-d ; ++k )
    {
        // -----------------------
        auto  ic = gA.icel(i,j,k); auto icP = gA.icel(i,j,k+1);
        auto  T2w = oldVel.w[ic]-dt*( t1.p[icP]-t1.p[ic]) / gA.Dzs[k];
        auto  eta = Dfib.eta[ic]+(Dfib.eta[icP]-Dfib.eta[ic])*gA.z_INTRPL(k);
        // -----------------------

        // -----------------------
        curtVel.w[ic] = Dfib.eta[ic] * wS + (1.0-eta) * T2w;
        Dfib.f[gA.icelCal(i,j,k)+zshift] = (curtVel.w[ic]- T2w) / dt;
        // -----------------------
    }

}


void update_UandF_omp(
    DfibArray& Dfib,
    const simuClass& simu,
    const pressure& t1,
    const velocity& oldVel,
    velocity& curtVel,
    const calDomain& Lo,
    grid& gA
)
{




    #pragma omp parallel firstprivate(simu, Lo)
    {

        const int  xshift = 0*gA.iceltotCal;
        const int  yshift = 1*gA.iceltotCal;
        const int  zshift = 2*gA.iceltotCal;

        const auto [nx, ny, nz] = gA.nxyz;

        const auto iceltot = gA.iceltot;
        double uS = 0.0, vS = 0.0, wS = 0.0;
        auto dt = simu.dt;

        int one = 1;
        int d = 0;

        // * ---------------- In x direction ----------------
        if(Lo.i_endof == nx-gA.gC){d = one;}
        else {d = 0;}

        #pragma omp for simd nowait
        for (int i = Lo.i_begin ; i < Lo.i_endof-d ; ++i )
        for (size_t j = Lo.j_begin ; j < Lo.j_endof ; ++j )
        for (int k = Lo.k_begin ; k < Lo.k_endof ; ++k )
        {
            // -----------------------
            auto  ic = gA.icel(i,j,k); auto  icP = gA.icel(i+1,j,k);
            auto  T2u = oldVel.u[ic] - dt*( t1.p[icP]-t1.p[ic] ) / gA.Dxs[i];
            auto  eta = Dfib.eta[ic]+(Dfib.eta[icP]-Dfib.eta[ic])*gA.x_INTRPL(i);
            // -----------------------

            // -----------------------
            curtVel.u[ic]= Dfib.eta[ic] * uS + (1.0-eta) * T2u;
            Dfib.f[gA.icelCal(i,j,k)+xshift] = (curtVel.u[ic] - T2u) / dt;
            // -----------------------
        }


        // * ---------------- In y direction ----------------
        if(Lo.j_endof == ny-gA.gC){d = one;}
        else {d = 0;}

        #pragma omp for simd
        for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
        for (auto j = Lo.j_begin ; j < Lo.j_endof-d ; ++j )
        for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
        {
            // -----------------------
            auto  ic = gA.icel(i,j,k); auto icP = gA.icel(i,j+1,k);
            auto  T2v = oldVel.v[ic] - dt*( t1.p[icP]-t1.p[ic] ) / gA.Dys[j];
            auto  eta = Dfib.eta[ic]+(Dfib.eta[icP]-Dfib.eta[ic])*gA.y_INTRPL(j);
            // -----------------------

            // -----------------------
            curtVel.v[ic] = Dfib.eta[ic] * vS + (1.0-eta) * T2v;
            Dfib.f[gA.icelCal(i,j,k)+yshift]  = (curtVel.v[ic] - T2v) / dt;
            // -----------------------
        }

        // * ---------------- In z direction ----------------
        if(Lo.k_endof == nz-gA.gC){d = one;}
        else {d = 0;}

        #pragma omp for 
        for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
        for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
        for (auto k = Lo.k_begin ; k < Lo.k_endof-d ; ++k )
        {
            // -----------------------
            auto  ic = gA.icel(i,j,k); auto icP = gA.icel(i,j,k+1);
            auto  T2w = oldVel.w[ic]-dt*( t1.p[icP]-t1.p[ic]) / gA.Dzs[k];
            auto  eta = Dfib.eta[ic]+(Dfib.eta[icP]-Dfib.eta[ic])*gA.z_INTRPL(k);
            // -----------------------

            // -----------------------
            curtVel.w[ic] = Dfib.eta[ic] * wS + (1.0-eta) * T2w;
            Dfib.f[gA.icelCal(i,j,k)+zshift] = (curtVel.w[ic]- T2w) / dt;
            // -----------------------
        }

    }


}




void removeVelocity(
    DfibArray& Dfib,
    velocity& TX,
    const calDomain& Lo,
    grid& gA
)
{
    #pragma omp parallel firstprivate( Lo )
    {
        #pragma omp for simd nowait
        for (int i = Lo.i_begin ; i < Lo.i_endof ; ++i )
        for (size_t j = Lo.j_begin ; j < Lo.j_endof ; ++j )
        for (int k = Lo.k_begin ; k < Lo.k_endof ; ++k )
        {
            auto ic = gA.icel(i,j,k);
            TX.u[ic] *= 0.5*( Dfib.eta[ic]+Dfib.eta[gA.icel(i+1,j,k)] );
            TX.v[ic] *= 0.5*( Dfib.eta[ic]+Dfib.eta[gA.icel(i,j+1,k)] );
            TX.w[ic] *= 0.5*( Dfib.eta[ic]+Dfib.eta[gA.icel(i,j,k+1)] );
        }
    }
}
