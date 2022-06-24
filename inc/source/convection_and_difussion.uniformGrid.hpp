// #pragma once

// using std::cout;
// using std::endl;
// using std::string;
// using std::vector;
// #include "sourceTerm.hpp"
// #include "controlPanel.hpp"
// /*

// * A present main , p is positve, n is negative, m is minus
// * In the quick shcemem and non-nuiform case, there are some diferent between A and B.

// */

// // ! ======================== Quick ========================
// inline std::pair<double, double> Quick_unif(
//     const double &Up,const double &Un,
//     const int &pp,const int &p,const int &icel,const int &n,const int &nn,
//     std::vector<double> &phi
// ){
//     double phi_P, phi_N;
//         if (Up > 0.0)
//             phi_P  = 0.75*phi[icel] + 0.375*phi[p] - 0.125*phi[n];
//         else
//             phi_P  = 0.75*phi[p] + 0.375*phi[icel] - 0.125*phi[pp];

//         if (Un > 0.0)
//             phi_N  = 0.75*phi[n] + 0.375*phi[icel] - 0.125*phi[nn];
//         else
//             phi_N  = 0.75*phi[icel] + 0.375*phi[n] - 0.125*phi[p];

//         return std::make_pair(phi_P,phi_N);
// }


// // ! ======================== Lud ========================
// inline std::pair<double, double> LudUnit(
//     const double &Up,const double &Un,
//     const int &pp,const int &p,
//     const int &icel,
//     const int &n,const int &nn,
//     std::vector<double> &phi
// ){
//     double phi_P, phi_N;
//         if (Up > 0.0)
//             phi_P = (phi[icel] + 0.5 * (phi[icel] - phi[n]) );
//         else
//             phi_P = (phi[p]    + 0.5 * (phi[p] - phi[pp])  );

//         if (Un > 0.0)
//             phi_N = (phi[n] + 0.5 * (phi[n] - phi[nn]) );
//         else
//             phi_N = (phi[icel]  + 0.5  *(phi[icel] - phi[p]));

//         return std::make_pair(phi_P,phi_N);
// }


// // ! ======================== Ud ========================
// inline std::pair<double, double> UdUnit(
//     const double &Up,const double &Un,
//     const int &p,const int &icel,const int &n,
//     std::vector<double> &phi
// ){
//     double phi_P, phi_N;

//     if (Up > 0.0)
//         phi_P = phi[icel] ;
//     else
//         phi_P = phi[p];

//     if (Un > 0.0)
//         phi_N = phi[n];
//     else
//         phi_N = phi[icel] ;

//     return std::make_pair(phi_P,phi_N);
// }
// // ! ======================== Ud ========================


// bool ConvectionDifussion_UniformG(
//     simuClass& simu,
//     velocity& T0,
//     velocity& T1,
//     calDomain& Lo,
//     grid& gA
// )
// {

//     auto C = discretizeTemporal(simu);

//     double ddx = 1./gA.Dx[2];

//     double ddy = 1./gA.Dy[2];

//     double ddz = 1./gA.Dz[2];


//     auto ii = [&](auto i, auto j, auto k){return gA.icel(i,j,k);};

//     const size_t nx = gA.nx,
//                  ny = gA.ny,
//                  nz = gA.nz;

//     #pragma omp parallel for
//     for (size_t i = Lo.i_begin; i < Lo.i_endof-1 ; ++i )
//     for (size_t j = Lo.j_begin; j < Lo.j_endof ; ++j )
//     for (size_t k = Lo.k_begin; k < Lo.k_endof ; ++k )
//     { //  =================== u =================== //

//         int icel = gA.icel(i, j, k);
//         // * -------------- convection term --------------
//         const auto [xm, xp, ym, yp, zm, zp, xmm, xpp, ymm, ypp, zmm, zpp] = gA.getNb12(icel);

//         const double Upx = 0.5 * (T0.u[xp]+T0.u[icel]); //main 
//         const double Unx = 0.5 * (T0.u[xm]+T0.u[icel]); //main
//         const double Upy = 0.5 * (T0.v[xp]+T0.v[icel]);
//         const double Uny = 0.5 * (T0.v[ym]+T0.v[ym+xp-icel]);
//         const double Upz = 0.5 * (T0.w[xp]+T0.w[icel]);
//         const double Unz = 0.5 * (T0.w[zm]+T0.w[zm+xp-icel]);


//     #if defined (CONVECTION_DIFUSSION_QUICK)

//         auto [uPx, uNx] = Quick_unif (Upx, Unx, xpp, xp, icel, xm, xmm, T0.u);
//         auto [uPy, uNy] = Quick_unif (Upy, Uny, ypp, yp, icel, ym, ymm, T0.u);
//         auto [uPz, uNz] = Quick_unif (Upz, Unz, zpp, zp, icel, zm, zmm, T0.u);

//     #elif defined (CONVECTION_DIFUSSION_LUD)

//         auto [uPx, uNx] = LudUnit (Upx, Unx, xpp, xp, icel, xm, xmm, T0.u);
//         auto [uPy, uNy] = LudUnit (Upy, Uny, ypp, yp, icel, ym, ymm, T0.u);
//         auto [uPz, uNz] = LudUnit (Upz, Unz, zpp, zp, icel, zm, zmm, T0.u);

//     #elif defined (CONVECTION_DIFUSSION_UD)

//         auto [uPx, uNx] = UdUnit (Upx, Unx, xpp, xp, icel, xm, xmm, T0.u);
//         auto [uPy, uNy] = UdUnit (Upy, Uny, ypp, yp, icel, ym, ymm, T0.u);
//         auto [uPz, uNz] = UdUnit (Upz, Unz, zpp, zp, icel, zm, zmm, T0.u);

//     #endif

//         double convection = -simu.dt*(
//                          (uPx*Upx-uNx*Unx) * ddx
//                         +(uPy*Upy-uNy*Uny) * ddy
//                         +(uPz*Upz-uNz*Unz) * ddz);

//         // * -------------- convection term --------------

//         // * -------------- difussion term --------------
        

//         #ifdef TERBULENCE_SMAGORINSKY
//             const auto nu =  simu.nu + T0.Viseff[gA.icelCal(i, j, k)];
//         #else
//             const auto nu =  simu.nu;
//         #endif

//          double difussion =    nu*simu.dt*(
//                                     ( (T0.u[xp]- 2*T0.u[icel] + T0.u[xm]) )  * ddx * ddx +
//                                     ( (T0.u[yp]- 2*T0.u[icel] + T0.u[ym])  ) * ddy * ddy +
//                                     ( (T0.u[zp]- 2*T0.u[icel] + T0.u[zm])  ) * ddz * ddz );
//         // * -------------- difussion term --------------
//         #if defined(TEMPORAL_DISCRETIZATION_1_ORDER)
//             souce(T1.u, T0.u, convection+difussion, icel);
//         #elif 
//             souce(T1.u, T0.su, T0.u, convection+difussion, C, icel, gA.icelCal(i,j,k));
//         #endif 
//     }


//     // ! v_star calculation (y component) 
//         #pragma omp parallel for
//         for (size_t i = Lo.i_begin; i < Lo.i_endof ; ++i )
//         for (size_t j = Lo.j_begin; j < Lo.j_endof-1 ; ++j )
//         for (size_t k = Lo.k_begin; k < Lo.k_endof ; ++k )
//         { //  =================== v =================== //
//             int icel = gA.icel(i,j,k);
//             // * -------------- convection term --------------
//             const auto [xm, xp, ym, yp, zm, zp, xmm, xpp, ymm, ypp, zmm, zpp] = gA.getNb12(icel);

//             const double Upx = 0.5 * (T0.u[yp]+T0.u[icel]);
//             const double Unx = 0.5 * (T0.u[xm]+T0.u[xm+yp-icel]);
//             const double Upy = 0.5 * (T0.v[yp]+T0.v[icel]); //main
//             const double Uny = 0.5 * (T0.v[ym]+T0.v[icel]); //main
//             const double Upz = 0.5 * (T0.w[yp]+T0.w[icel]);
//             const double Unz = 0.5 * (T0.w[zm]+T0.w[zm+yp-icel]);


//         #if defined (CONVECTION_DIFUSSION_QUICK)
        
//             auto [vPx, vNx] = Quick_unif (Upx, Unx, xpp, xp, icel, xm, xmm, T0.v);
//             auto [vPy, vNy] = Quick_unif (Upy, Uny, ypp, yp, icel, ym, ymm, T0.v);
//             auto [vPz, vNz] = Quick_unif (Upz, Unz, zpp, zp, icel, zm, zmm, T0.v);
        
//         #elif defined (CONVECTION_DIFUSSION_LUD)

//             auto [vPx, vNx] = LudUnit (Upx, Unx, xpp, xp, icel, xm, xmm, T0.v);
//             auto [vPy, vNy] = LudUnit (Upy, Uny, ypp, yp, icel, ym, ymm, T0.v);
//             auto [vPz, vNz] = LudUnit (Upz, Unz, zpp, zp, icel, zm, zmm, T0.v);
        
//         #elif defined (CONVECTION_DIFUSSION_UD)
        
//             auto [vPx, vNx] = UdUnit (Upx, Unx, xp, icel, xm, T0.v);
//             auto [vPy, vNy] = UdUnit (Upy, Uny, yp, icel, ym, T0.v);
//             auto [vPz, vNz] = UdUnit (Upz, Unz, zp, icel, zm, T0.v);
        
//         #endif

            
//             const double convection =   -simu.dt*(
//                                          (Upx*vPx-Unx*vNx) * ddx
//                                         +(Upy*vPy-Uny*vNy) * ddy
//                                         +(Upz*vPz-Unz*vNz) * ddz); 
//         // * -------------- convection term --------------


//         // * -------------- difussion term --------------

//         #ifdef TERBULENCE_SMAGORINSKY
//             const auto nu =  simu.nu + T0.Viseff[gA.icelCal(i, j, k)];
//         #else
//             const auto nu =  simu.nu;
//         #endif
//             const double difussion = nu*simu.dt*( 
//                                     ( (T0.v[xp] - 2 * T0.v[icel] + T0.v[xm]) ) * ddx* ddx +
//                                     ( (T0.v[yp] - 2 * T0.v[icel] + T0.v[ym]) ) * ddy* ddy +
//                                     ( (T0.v[zp] - 2 * T0.v[icel] + T0.v[zm]) ) * ddz* ddz );
//         // * -------------- difussion term --------------
        
                    

//         #if defined(TEMPORAL_DISCRETIZATION_1_ORDER)
//             souce(T1.v, T0.v, convection+difussion, icel);
//         #elif 
//             souce(T1.v, T0.sv, T0.v, convection+difussion, C, icel, gA.icelCal(i,j,k));
//         #endif 
//     }

//  //w_star calculation (z component)
//     #pragma omp parallel for
//     for (size_t i = Lo.i_begin; i < Lo.i_endof ; ++i )
//     for (size_t j = Lo.j_begin; j < Lo.j_endof ; ++j )
//     for (size_t k = Lo.k_begin; k < Lo.k_endof-1 ; ++k )
//     { //  =================== w =================== //
//         const int icel = gA.icel(i, j, k);
//         // * -------------- convection term --------------
//         auto [xm, xp, ym, yp, zm, zp, xmm, xpp, ymm, ypp, zmm, zpp] = gA.getNb12(icel);
        
//         const double Upx = 0.5 * ( T0.u[zp]+T0.u[icel] );
//         const double Unx = 0.5 * ( T0.u[xm]+T0.u[xm+1] ); 
//         const double Upy = 0.5 * ( T0.v[zp]+T0.v[icel] );
//         const double Uny = 0.5 * ( T0.v[ym]+T0.v[ym+1] );
//         const double Upz = 0.5 * ( T0.w[zp]+T0.w[icel] );//main
//         const double Unz = 0.5 * ( T0.w[zm]+T0.w[icel] );//main

//     #if defined (CONVECTION_DIFUSSION_QUICK)

//         auto [wPx, wNx] = Quick_unif (Upx, Unx, xpp, xp, icel, xm, xmm, T0.w);
//         auto [wPy, wNy] = Quick_unif (Upy, Uny, ypp, yp, icel, ym, ymm, T0.w);
//         auto [wPz, wNz] = Quick_unif (Upz, Unz, zpp, zp, icel, zm, zmm, T0.w);
        
//     #elif defined (CONVECTION_DIFUSSION_LUD)

//         auto [wPx, wNx] = LudUnit (Upx, Unx, xpp, xp, icel, xm, xmm, T0.w);
//         auto [wPy, wNy] = LudUnit (Upy, Uny, ypp, yp, icel, ym, ymm, T0.w);
//         auto [wPz, wNz] = LudUnit (Upz, Unz, zpp, zp, icel, zm, zmm, T0.w);


//     #elif defined (CONVECTION_DIFUSSION_UD)
   
//         auto [wPx, wNx] = UdUnit (Upx, Unx, xp, icel, xm, T0.w);
//         auto [wPy, wNy] = UdUnit (Upy, Uny, yp, icel, ym, T0.w);
//         auto [wPz, wNz] = UdUnit (Upz, Unz, zp, icel, zm, T0.w);

//     #endif


//         const double convection =   -simu.dt*(
//                                      (Upx*wPx-Unx*wNx)  * ddx
//                                     +(Upy*wPy-Uny*wNy)  * ddy
//                                     +(Upz*wPz-Unz*wNz)  * ddz);
//         // * -------------- convection term --------------

//         // * -------------- difussion term --------------
//     #ifdef TERBULENCE_SMAGORINSKY
//         auto nu =  simu.nu + T0.Viseff[gA.icelCal(i, j, k)];
//     #else
//         auto nu =  simu.nu;
//     #endif

//     const double difussion =    nu*simu.dt*(
//                                 ( ( T0.w[xp]-T0.w[icel] ) / gA.Dxs[i] 
//                                 - ( T0.w[icel]-T0.w[xm] ) / gA.Dxs[i-1] ) * ddx+
//                                 ( ( T0.w[yp]-T0.w[icel] ) / gA.Dys[j] 
//                                 - ( T0.w[icel]-T0.w[ym] ) / gA.Dys[j-1] ) * ddy+
//                                 ( ( T0.w[zp]-T0.w[icel] ) / gA.Dz[k+1] 
//                                 - ( T0.w[icel]-T0.w[zm] ) / gA.Dz[k]   ) * ddz) ;
//     // * -------------- difussion term --------------


//     #if defined(TEMPORAL_DISCRETIZATION_1_ORDER)
//         souce(T1.w, T0.w, convection+difussion, icel);
//     #elif 
//         souce(T1.w, T0.sw, T0.w, convection+difussion, C, icel, gA.icelCal(i,j,k));
//     #endif 

//     }


//     #pragma omp barrier

//     if (simu.Second_s)
//         simu.Second_s = false;

//     if (simu.First_s){
//         simu.First_s = false;
//         simu.Second_s = true;
//     }

    
//     return true;
// }