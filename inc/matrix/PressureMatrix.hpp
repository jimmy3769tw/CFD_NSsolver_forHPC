#pragma once

#include"controlPanel.hpp"

auto createPressureMatrix(
    MxClass& Mx,
    simuClass& simu,
    calDomain& Lo,
    grid& gA
)
{
    auto [nx, ny, nz, gC] = gA.nxyzgC;

    Mx.matA_ell.resize(gA.iceltotCal,7);

    gA.jp.resize(gA.iceltotCal);

    // -------------------------------------------------
    auto xp = [&](auto i){return 1.0/gA.Dx[i]/gA.Dxs[i];};
    auto yp = [&](auto j){return 1.0/gA.Dy[j]/gA.Dys[j];};
    auto zp = [&](auto k){return 1.0/gA.Dz[k]/gA.Dzs[k];};
    // -------------------------------------------------

    // -------------------------------------------------
    auto xm = [&](auto i){return 1.0/gA.Dx[i]/gA.Dxs[i-1];};
    auto ym = [&](auto j){return 1.0/gA.Dy[j]/gA.Dys[j-1];};
    auto zm = [&](auto k){return 1.0/gA.Dz[k]/gA.Dzs[k-1];};
    // -------------------------------------------------

    size_t  xBeg = gC,
            yBeg = gC,
            zBeg = gC,
            xEnd = nx - gC -1 ,
            yEnd = ny - gC -1 ,
            zEnd = nz - gC -1 ;

    for (size_t i = xBeg ; i <= xEnd ; ++i )
    for (size_t j = yBeg ; j <= yEnd ; ++j )
    for (size_t k = zBeg ; k <= zEnd ; ++k )
    {
        double  aa = 0.0, bb = 0.0, cc = 0.0;
    // -------------------------------------------------
    #if defined(BC_0_P_NEUMANN) && defined(BC_1_P_NEUMANN)
        if(i==xBeg){ aa = -xm(i); }
        if(i==xEnd){ aa = -xp(i); }
    #endif
    // -------------------------------------------------

    // -------------------------------------------------
    #if defined(BC_2_P_NEUMANN) && defined(BC_3_P_NEUMANN) 
        if(j==yBeg) {bb = -ym(j);}
        if(j==yEnd) {bb = -yp(j);}
    #endif
    // -------------------------------------------------

    // -------------------------------------------------
    #if defined(BC_4_P_NEUMANN) && defined(BC_5_P_NEUMANN) 
        if(k==zBeg) { cc = -zm(k); }
        if(k==zEnd) { cc = -zp(k); }
    #endif
    // -------------------------------------------------
        const int row = gA.icelCal(i,j,k);
        auto val = xm(i)+xp(i)+ym(j)+yp(j)+zm(k)+zp(k)+ aa + bb + cc;

        #if defined(JACOBI_PC)
        gA.jp[row] = 1./val;
        #elif
        gA.jp[rpw] = 1.;
        #endif 

        Mx.matA_ell.set(row, row, val * gA.jp[row]);
    }

    auto c = gA.iceltotCal;

    for (size_t i = xBeg ; i <= xEnd ; ++i )
    for (size_t j = yBeg ; j <= yEnd ; ++j )
    for (size_t k = zBeg ; k <= zEnd ; ++k )
    {
        const int row = gA.icelCal(i,j,k);
    // -------------------------------------------------
        if(i!=xBeg) { Mx.matA_ell.set(row, gA.icelCal(i-1,j,k), -xm(i) * gA.jp[row] ); }
        if(i!=xEnd) { Mx.matA_ell.set(row, gA.icelCal(i+1,j,k), -xp(i) * gA.jp[row] ); }

        if(j!=yBeg) { Mx.matA_ell.set(row, gA.icelCal(i,j-1,k), -ym(j) * gA.jp[row] ); }
        if(j!=yEnd) { Mx.matA_ell.set(row, gA.icelCal(i,j+1,k), -yp(j) * gA.jp[row] ); }

        if(k!=zBeg) { Mx.matA_ell.set(row, gA.icelCal(i,j,k-1), -zm(k) * gA.jp[row] ); }
        if(k!=zEnd) { Mx.matA_ell.set(row, gA.icelCal(i,j,k+1), -zp(k) * gA.jp[row] ); }
    // -------------------------------------------------

    // -------------------------------------------------
    #if defined(BC_0_P_PERIODIC) &&  defined(BC_1_P_PERIODIC)
        if(i==xBeg) { Mx.matA_ell.set(row, gA.icelCal(xEnd,j,k), -xm(i) * gA.jp[row] ); }
        if(i==xEnd) { Mx.matA_ell.set(row, gA.icelCal(xBeg,j,k), -xp(i) * gA.jp[row] ); }
    #endif
    // -------------------------------------------------


    // -------------------------------------------------
    #if defined(BC_2_P_PERIODIC) &&  defined(BC_3_P_PERIODIC)
        if(j==yBeg) { Mx.matA_ell.set(row, gA.icelCal(i,yEnd,k), -ym(j) * gA.jp[row]); }
        if(j==yEnd) { Mx.matA_ell.set(row, gA.icelCal(i,yBeg,k), -yp(j) * gA.jp[row]); }
    #endif
    // -------------------------------------------------

    // -------------------------------------------------
    #if defined(BC_4_P_PERIODIC) && defined(BC_5_P_PERIODIC)
        if(k==zBeg) { Mx.matA_ell.set(row, gA.icelCal(i,j,zBeg), -zm(k) * gA.jp[row]); }
        if(k==zEnd) { Mx.matA_ell.set(row, gA.icelCal(j,j,zEnd), -zp(k) * gA.jp[row]); }
    #endif
    // -------------------------------------------------

    }

    return true;
}
