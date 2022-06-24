#pragma once

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::make_pair;

#include "sourceTerm.hpp"
#include "controlPanel.hpp"

/*

* A present main , p is positve, n is negative, m is minus
* In the quick shcemem and non-nuiform case, there are some diferent between A and B.

*/

// ! ======================== Quick ========================
inline std::pair<double, double> QuickUnitA(
    const double &Up,const double &Un,
    const int &pp,const int &p,const int &ic,const int &n,const int &nn,
    const std::vector<double> &Di,
    const std::vector<double> &Ds,
    const int &Idx,
    std::vector<double> &phi
){
    double phi_P, phi_N;

    if (Up > 0.0)
        phi_P  = 0.5*(phi[ic] + phi[p] )
                -0.125*Di[Idx+1]*Di[Idx+1]/Ds[Idx] * (
                (phi[p] - phi[ic] ) / Di[Idx+1] - (phi[ic] - phi[n] ) / Di[Idx]) ;
    else
        phi_P  = 0.5*(phi[ic] + phi[p] )
                - 0.125*Di[Idx+1]*Di[Idx+1]/Ds[Idx+1] * ( 
               (phi[pp] - phi[p] ) / Di[Idx+2] - (phi[p] - phi[ic] ) / Di[Idx+1] );

    if (Un > 0.0)
        phi_N  = 0.5*(phi[n] + phi[ic] ) 
                -0.125*Di[Idx]*Di[Idx]/Ds[Idx-1] * ( 
               (phi[ic] - phi[n] ) / Di[Idx] - (phi[n]   - phi[nn] ) / Di[Idx-1] ) ;
    else
        phi_N  = 0.5*(phi[n] + phi[ic] ) 
                -0.125*Di[Idx]*Di[Idx]/Ds[Idx]* (
                 ( phi[p] - phi[ic] )/Di[Idx+1] - ( phi[ic] - phi[n] )/Di[Idx]);

    return make_pair(phi_P,phi_N);
}



inline std::pair<double, double> QuickUnitB(
    const double &Up,const double &Un,
    const int &pp, const int &p, const int &ic, const int &n, const int &nn,
    const std::vector<double> &Di,
    const std::vector<double> &Ds,
    const int &Idx,
    std::vector<double> &phi
){
    double phi_P, phi_N;
    // ----------------------------------
    if (Up > 0.0)
        phi_P  = 0.5*(phi[ic] + phi[p])
        -0.125*Ds[Idx]*Ds[Idx]/Di[Idx] * ( 
            (phi[p] - phi[ic] ) / Ds[Idx] - (phi[ic] - phi[n] ) / Ds[Idx-1]) ;
    else
        phi_P  = 0.5*(phi[ic] + phi[p])
        - 0.125*Ds[Idx]*Ds[Idx]/Di[Idx+1] * ( 
            (phi[pp] - phi[p] ) / Ds[Idx+1] - (phi[p] - phi[ic] ) / Ds[Idx] );
    // ----------------------------------

    // ----------------------------------
    if (Un > 0.0)
        phi_N  = 0.5*(phi[n] + phi[ic]) 
        -0.125*Ds[Idx-1]*Ds[Idx-1]/Di[Idx-1] * ( 
            (phi[ic] - phi[n] ) / Ds[Idx-1] - (phi[n]  - phi[nn]) / Ds[Idx-2]) ;
    else
        phi_N  = 0.5*(phi[n] + phi[ic]) 
        -0.125*Ds[Idx-1]*Ds[Idx-1]/Di[Idx] * (
             (phi[p] - phi[ic] ) / Ds[Idx] - (phi[ic] - phi[n] ) / Ds[Idx-1]);
    // ----------------------------------
        return make_pair(phi_P,phi_N);
}

// * ======================== Quick ========================


// ! ======================== Lud ========================
inline std::pair<double, double> LudUnitA(
    const double &Up,const double &Un,
    const int &pp,const int &p,
    const int &ic,
    const int &n,const int &nn,
    const std::vector<double> &Di,
    const std::vector<double> &Ds,
    const int &Idx,
    std::vector<double> &phi
){
    double phi_P, phi_N;

    // ----------------------------------
    if (Up > 0.0) { phi_P = (phi[ic] + 0.5 * Ds[Idx] * (phi[ic] - phi[n]) / Di[Idx]); }
    else { phi_P = (phi[p] + 0.5 * Ds[Idx+1] * (phi[p] - phi[pp]) / Di[Idx+1]); }
    // ----------------------------------

    // ----------------------------------
    if (Un > 0.0) {phi_N = (phi[n] + 0.5 * Ds[Idx-1] * (phi[n] - phi[nn]) / Di[Idx-1]);}
    else { phi_N = (phi[ic] + 0.5 * Ds[Idx] *(phi[ic] - phi[p])/ Di[Idx+1]); }
    // ----------------------------------
        

    return make_pair(phi_P, phi_N);
}

inline std::pair<double, double> LudUnitB(
    const double &Up,const double &Un,
    const int &pp,const int &p,
    const int &ic,
    const int &n,const int &nn,
    const std::vector<double> &Di,
    const std::vector<double> &Ds,
    const int &Idx,
    std::vector<double> &phi
){
    double phi_P, phi_N;

    // ----------------------------------
    if (Up > 0.0) { phi_P = (phi[ic] + 0.5 * Ds[Idx] * (phi[ic] - phi[n]) / Di[Idx]); }
    else {phi_P = (phi[p] + 0.5 * Ds[Idx+1] * (phi[p] - phi[pp]) / Di[Idx+1]); }
    // ----------------------------------
        
    // ----------------------------------
    if (Un > 0.0) { phi_N = (phi[n] + 0.5 * Ds[Idx-1] * (phi[n] - phi[nn]) / Di[Idx-1]); }
    else {phi_N = (phi[ic] + 0.5 * Ds[Idx] *(phi[ic] - phi[p])/ Di[Idx+1]); }
    // ----------------------------------

    return make_pair(phi_P, phi_N);
}
// * ======================== Lud ========================


// ! ======================== Ud ========================
inline std::pair<double, double> UdUnit(
    const double &Up,const double &Un,
    const int &p,const int &ic,const int &n,
    std::vector<double> &Di,
    std::vector<double> &Ds,
    const int &Idx,
    std::vector<double> &phi
){
    double phi_P, phi_N;
    // ----------------------------------
    if (Up > 0.0) { phi_P = phi[ic];}
    else {phi_P = phi[p];}
    // ----------------------------------

    // ----------------------------------
    if (Un > 0.0) {phi_N = phi[n];}
    else {phi_N = phi[ic] ;}
    // ----------------------------------

    return make_pair(phi_P,phi_N);
}
// * ======================== Ud ========================


bool ConvectionDifussion(
    simuClass& simu,
    velocity& oldVel,
    velocity& curtVel,
    const calDomain& cal,
    grid& gA
)
{
    // ----------------------
    #if    defined (TEMPORAL_DISCRETIZATION_2_ORDER)\
        || defined (TEMPORAL_DISCRETIZATION_3_ORDER)

    auto C = discretizeTemporal(simu); 

    #endif
    // ----------------------

    // ----------------------
    auto ii = [&](auto i, auto j, auto k){return gA.icel(i,j,k);};
    const int nx = gA.nx, ny = gA.ny, nz = gA.nz;
    int d = 0, one = 1 ;
    // ----------------------
// ! ---------------------------------- x  ----------------------------------


    if(cal.i_endof == nx-gA.gC){ d = one;}
    else{ d = 0;}

    #pragma omp parallel for
    for (auto i = cal.i_begin; i < cal.i_endof-d ; ++i )
    for (auto j = cal.j_begin; j < cal.j_endof   ; ++j )
    for (auto k = cal.k_begin; k < cal.k_endof   ; ++k )
    { //  =================== u =================== //

        int ic = gA.icel(i, j, k);
        // ! -------------- convection term --------------
        const auto [xm, xp, ym, yp, zm, zp, xmm, xpp, ymm, ypp, zmm, zpp] = gA.getNb12(ic);
        const auto IT = gA.x_INTRPL(i);

        const double Up = 0.5 * (oldVel.u[xp]+oldVel.u[ic]); //main 
        const double Un = 0.5 * (oldVel.u[xm]+oldVel.u[ic]); //main

        const double Vp = oldVel.v[ic]+(oldVel.v[xp]-oldVel.v[ic])*IT;
        const double Vn = oldVel.v[ym]+(oldVel.v[ym+xp-ic]-oldVel.v[ym])*IT;

        const double Wp = oldVel.w[ic]+(oldVel.w[xp]-oldVel.w[ic])*IT;
        const double Wn = oldVel.w[zm]+(oldVel.w[zm+xp-ic]-oldVel.w[zm])*IT;

        // * ----------------------------------------------

    #if defined (CONVECTION_DIFUSSION_QUICK)

        auto [uPx, uNx] = QuickUnitA (Up, Un, xpp, xp, ic, xm, xmm, gA.Dx, gA.Dxs, i, oldVel.u);
        auto [uPy, uNy] = QuickUnitB (Vp, Vn, ypp, yp, ic, ym, ymm, gA.Dy, gA.Dys, j, oldVel.u);
        auto [uPz, uNz] = QuickUnitB (Wp, Wn, zpp, zp, ic, zm, zmm, gA.Dz, gA.Dzs, k, oldVel.u);

    #elif defined (CONVECTION_DIFUSSION_LUD)

        auto [uPx, uNx] = LudUnitA (Up, Un, xpp, xp, ic, xm, xmm, gA.Dx, gA.Dxs, i, oldVel.u);
        auto [uPy, uNy] = LudUnitB (Vp, Vn, ypp, yp, ic, ym, ymm, gA.Dy, gA.Dys, j, oldVel.u);
        auto [uPz, uNz] = LudUnitB (Wp, Wn, zpp, zp, ic, zm, zmm, gA.Dz, gA.Dzs, k, oldVel.u);

    #elif defined (CONVECTION_DIFUSSION_UD)

        auto [uPx, uNx] = UdUnit (Up, Un, xp, ic, xm, gA.Dx, gA.Dxs, i, oldVel.u);
        auto [uPy, uNy] = UdUnit (Vp, Vn, yp, ic, ym, gA.Dy, gA.Dys, j, oldVel.u);
        auto [uPz, uNz] = UdUnit (Wp, Wn, zp, ic, zm, gA.Dz, gA.Dzs, k, oldVel.u);

    #endif

        double convection = -simu.dt*(
                         (uPx*uPx-uNx*uNx) / gA.Dxs[i]
                        +(uPy*Vp-uNy*Vn) / gA.Dy[j]
                        +(uPz*Wp-uNz*Wn) / gA.Dz[k]);
        // * -------------- convection term --------------

        // ! -------------- difussion term --------------

        #ifdef TERBULENCE_SMAGORINSKY
            auto nu =  simu.nu + oldVel.Viseff[gA.icelCal(i,j,k)];
        #else
            auto nu =  simu.nu;
        #endif

         double difussion = nu*simu.dt*(
                    ( (oldVel.u[xp]-oldVel.u[ic]) / gA.Dx[i+1] 
                    - (oldVel.u[ic]-oldVel.u[xm]) / gA.Dx[i]   ) / gA.Dxs[i]+
                    ( (oldVel.u[yp]-oldVel.u[ic]) / gA.Dys[j]   
                    - (oldVel.u[ic]-oldVel.u[ym]) / gA.Dys[j-1] ) / gA.Dy[j]+
                    ( (oldVel.u[zp]-oldVel.u[ic]) / gA.Dzs[k]   
                    - (oldVel.u[ic]-oldVel.u[zm]) / gA.Dzs[k-1] ) / gA.Dz[k]);

      // * -------------- difussion term --------------


    #if defined(TEMPORAL_DISCRETIZATION_1_ORDER)
        curtVel.u[ic] = oldVel.u[ic] + convection + difussion;    
    #elif 
        souce(curtVel.u, oldVel.su, oldVel.u, convection+difussion, C, ic, gA.icelCal(i,j,k));
    #endif 

    }


// ! ---------------------------------- y  ----------------------------------

    if(cal.j_endof == ny-gA.gC){ d = one;}
    else{ d = 0;}


        #pragma omp parallel for
        for (auto i = cal.i_begin; i < cal.i_endof ; ++i )
        for (auto j = cal.j_begin; j < cal.j_endof-d ; ++j )
        for (auto k = cal.k_begin; k < cal.k_endof ; ++k )
        { //  =================== v =================== //
            int ic = gA.icel(i,j,k);
            // ! -------------- convection term --------------
            const auto [xm, xp, ym, yp, zm, zp, xmm, xpp, ymm, ypp, zmm, zpp] 
                          = gA.getNb12(ic);
    
            const auto IT = gA.y_INTRPL(j);

            const double Vp = 0.5 * (oldVel.v[yp]+oldVel.v[ic]); //main
            const double Vn = 0.5 * (oldVel.v[ym]+oldVel.v[ic]); //main

            const double Up = oldVel.u[ic]+(oldVel.u[yp]-oldVel.u[ic])*IT;
            const double Un = oldVel.u[xm]+(oldVel.u[xm+yp-ic]-oldVel.u[xm])*IT;

            const double Wp = oldVel.w[ic]+(oldVel.w[yp]-oldVel.w[ic])*IT;
            const double Wn = oldVel.w[zm]+(oldVel.w[zm+yp-ic]-oldVel.w[zm])*IT;


        #if defined (CONVECTION_DIFUSSION_QUICK)
        
            auto [vPx, vNx] = QuickUnitB (Up, Un, xpp, xp, ic, xm, xmm, gA.Dx, gA.Dxs, i, oldVel.v);
            auto [vPy, vNy] = QuickUnitA (Vp, Vn, ypp, yp, ic, ym, ymm, gA.Dy, gA.Dys, j, oldVel.v);
            auto [vPz, vNz] = QuickUnitB (Wp, Wn, zpp, zp, ic, zm, zmm, gA.Dz, gA.Dzs, k, oldVel.v);
        
        #elif defined (CONVECTION_DIFUSSION_LUD)

            auto [vPx, vNx] = LudUnitB (Up, Un, xpp, xp, ic, xm, xmm, gA.Dx, gA.Dxs, i, oldVel.v);
            auto [vPy, vNy] = LudUnitA (Vp, Vn, ypp, yp, ic, ym, ymm, gA.Dy, gA.Dys, j, oldVel.v);
            auto [vPz, vNz] = LudUnitB (Wp, Wn, zpp, zp, ic, zm, zmm, gA.Dz, gA.Dzs, k, oldVel.v);
        
        #elif defined (CONVECTION_DIFUSSION_UD)
        
            auto [vPx, vNx] = UdUnit (Up, Un, xp, ic, xm, gA.Dx, gA.Dxs, i, oldVel.v);
            auto [vPy, vNy] = UdUnit (Vp, Vn, yp, ic, ym, gA.Dy, gA.Dys, j, oldVel.v);
            auto [vPz, vNz] = UdUnit (Wp, Wn, zp, ic, zm, gA.Dz, gA.Dzs, k, oldVel.v);
        
        #endif

            
            double convection =   -simu.dt*(
                                     (Up*vPx-Un*vNx) / gA.Dx[i] 
                                    +(vPy*vPy-vNy*vNy) / gA.Dys[j] 
                                    +(Wp*vPz-Wn*vNz) / gA.Dz[k]); 
        // * -------------- convection term --------------


        // ! -------------- difussion term --------------

        #ifdef TERBULENCE_SMAGORINSKY
            const auto nu =  simu.nu + oldVel.Viseff[gA.icelCal(i, j, k)];
        #else
            const auto nu =  simu.nu;
        #endif
            const double difussion = nu*simu.dt*( 
                                    ( (oldVel.v[xp] - oldVel.v[ic]) / gA.Dxs[i]
                                    - (oldVel.v[ic] - oldVel.v[xm]) / gA.Dxs[i-1] ) / gA.Dx[i] +
                                    ( (oldVel.v[yp] - oldVel.v[ic]) / gA.Dy[j+1] 
                                    - (oldVel.v[ic] - oldVel.v[ym]) / gA.Dy[j]   ) / gA.Dys[j] +
                                    ( (oldVel.v[zp] - oldVel.v[ic]) / gA.Dzs[k] 
                                    - (oldVel.v[ic] - oldVel.v[zm]) / gA.Dzs[k-1] ) / gA.Dz[k] );
        // * -------------- difussion term --------------
        
    
    #if defined(TEMPORAL_DISCRETIZATION_1_ORDER)
        curtVel.v[ic] = oldVel.v[ic] + convection+difussion;
    #elif 
        souce(curtVel.v, oldVel.sv, oldVel.v, convection+difussion, C, ic, gA.icelCal(i,j,k));
    #endif 

    }



// ! ---------------------------------- z  ----------------------------------

    if(cal.k_endof == nz-gA.gC){ d = one;}
    else{ d = 0;}

    #pragma omp parallel for
    for (auto i = cal.i_begin; i < cal.i_endof   ; ++i )
    for (auto j = cal.j_begin; j < cal.j_endof   ; ++j )
    for (auto k = cal.k_begin; k < cal.k_endof-d ; ++k )
    { //  =================== w =================== //
        const int ic = gA.icel(i, j, k);
        // ! -------------- convection term --------------
        auto [xm, xp, ym, yp, zm, zp, xmm, xpp, ymm, ypp, zmm, zpp]
            = gA.getNb12(ic);

        const auto IT = gA.z_INTRPL(k);

        const double Wp = 0.5 * ( oldVel.w[zp]+oldVel.w[ic] );//main
        const double Wn = 0.5 * ( oldVel.w[zm]+oldVel.w[ic] );//main

        const double Up = oldVel.u[ic]+(oldVel.u[zp]-oldVel.u[ic])*IT;
        const double Un = oldVel.u[xm]+(oldVel.u[xm+zp-ic]-oldVel.u[xm])*IT; 

        const double Vp = oldVel.v[ic]+(oldVel.v[zp]-oldVel.v[ic])*IT;;
        const double Vn = oldVel.v[ym]+(oldVel.v[ym+zp-ic]-oldVel.v[ym])*IT;;
        // *----------------------------------------------

    #if defined (CONVECTION_DIFUSSION_QUICK)

        auto [wPx, wNx] = QuickUnitB (Up, Un, xpp, xp, ic, xm, xmm, gA.Dx, gA.Dxs, i, oldVel.w);
        auto [wPy, wNy] = QuickUnitB (Vp, Vn, ypp, yp, ic, ym, ymm, gA.Dy, gA.Dys, j, oldVel.w);
        auto [wPz, wNz] = QuickUnitA (Wp, Wn, zpp, zp, ic, zm, zmm, gA.Dz, gA.Dzs, k, oldVel.w);
        
    #elif defined (CONVECTION_DIFUSSION_LUD)

        auto [wPx, wNx] = LudUnitB (Up, Un, xpp, xp, ic, xm, xmm, gA.Dx, gA.Dxs, i, oldVel.w);
        auto [wPy, wNy] = LudUnitB (Vp, Vn, ypp, yp, ic, ym, ymm, gA.Dy, gA.Dys, j, oldVel.w);
        auto [wPz, wNz] = LudUnitA (Wp, Wn, zpp, zp, ic, zm, zmm, gA.Dz, gA.Dzs, k, oldVel.w);


    #elif defined (CONVECTION_DIFUSSION_UD)
   
        auto [wPx, wNx] = UdUnit (Up, Un, xp, ic, xm, gA.Dx, gA.Dxs, i, oldVel.w);
        auto [wPy, wNy] = UdUnit (Vp, Vn, yp, ic, ym, gA.Dy, gA.Dys, j, oldVel.w);
        auto [wPz, wNz] = UdUnit (Wp, Wn, zp, ic, zm, gA.Dz, gA.Dzs, k, oldVel.w);

    #endif

        const double convection =   -simu.dt*(
                                     (Up*wPx-Un*wNx) / gA.Dx[i]
                                    +(Vp*wPy-Vn*wNy) / gA.Dy[j]
                                    +(wPz*wPz-wNz*wNz) / gA.Dzs[k]);
        // * -------------- convection term --------------

        // * -------------- difussion term --------------
    #ifdef TERBULENCE_SMAGORINSKY
        auto nu =  simu.nu + oldVel.Viseff[gA.icelCal(i, j, k)];
    #else
        auto nu =  simu.nu;
    #endif

    const double difussion =    nu*simu.dt*(
                                ( ( oldVel.w[xp]-oldVel.w[ic] ) / gA.Dxs[i] 
                                - ( oldVel.w[ic]-oldVel.w[xm] ) / gA.Dxs[i-1] ) / gA.Dx[i]+
                                ( ( oldVel.w[yp]-oldVel.w[ic] ) / gA.Dys[j] 
                                - ( oldVel.w[ic]-oldVel.w[ym] ) / gA.Dys[j-1] ) / gA.Dy[j]+
                                ( ( oldVel.w[zp]-oldVel.w[ic] ) / gA.Dz[k+1] 
                                - ( oldVel.w[ic]-oldVel.w[zm] ) / gA.Dz[k]   ) / gA.Dzs[k]) ;
    // * -------------- difussion term --------------


    #if defined(TEMPORAL_DISCRETIZATION_1_ORDER)
        // souce(curtVel.w, oldVel.w, convection+difussion, ic);
        curtVel.w[ic] = oldVel.w[ic] + convection+difussion;
    #elif 
        souce(curtVel.w, oldVel.sw, oldVel.w, convection+difussion, C, ic, gA.icelCal(i,j,k));
    #endif 


    }

    #pragma omp barrier
    if (simu.Second_s)
        simu.Second_s = false;

    if (simu.First_s){
        simu.First_s = false;
        simu.Second_s = true;
    }
    
    return true;
}