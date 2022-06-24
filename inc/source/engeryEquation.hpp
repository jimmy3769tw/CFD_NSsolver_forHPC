#pragma once



#include"controlPanel.hpp"

// ! ======================== QuickUnitA ========================
inline std::pair<double, double> EngeryQuickUnit(
    const double &Up,const double &Un,
    const int &pp,const int &p,const int &icel,const int &n,const int &nn,
    std::vector<double> &Di,
    std::vector<double> &Ds,
    const size_t &Idx,
    std::vector<double> &phi
){
    double phi_P, phi_N;
        if(Up > 0)
            phi_P = 0.5*( phi[p] + phi[icel] ) 
                        -0.125*Ds[Idx]*Ds[Idx]/Di[Idx]*
                        (( (phi[p]-phi[icel])/Ds[Idx] 
                        - (phi[icel] - phi[n] )/Ds[Idx-1]));
        else
            phi_P = 0.5*( phi[p] + phi[icel] ) 
                        -0.125*Ds[Idx]*Ds[Idx]/Di[Idx]*
                        (( (phi[pp]-phi[p])/Ds[Idx+1] 
                        - (phi[p] - phi[icel] )/Ds[Idx]));

        if(Un > 0)
            phi_N = 0.5*( phi[icel] + phi[n] ) 
                        -0.125*Ds[Idx]*Ds[Idx]/Di[Idx]*
                        (( (phi[icel]-phi[n])/Ds[Idx-1] 
                        - (phi[n] - phi[nn] )/Ds[Idx-2]));
        else
            phi_N = 0.5*( phi[icel] + phi[n] ) 
                        -0.125*Ds[Idx]*Ds[Idx]/Di[Idx]*
                        (( (phi[p]-phi[icel])/Ds[Idx] 
                        - (phi[icel] - phi[n] )/Ds[Idx-1]));

        return std::make_pair(phi_P,phi_N);
}
// ! ======================== QuickUnitA ========================


// //*  EngeryEquQuickscheme(simu,T3,t1,Lo,gridA);
// void EngeryEquQuickscheme(
//     simuClass& simu,
//     velocity& T3,
//     pressure& t1,
//     calDomain& Lo,
//     grid& gridA
// )
// {


//     const size_t nx = gridA.nx, ny = gridA.ny, nz = gridA.nz;
//     const size_t gC = gridA.gC;


//     for (size_t i = Lo.i_begin; i < Lo.i_endof ; ++i )
//     for (size_t j = Lo.j_begin; j < Lo.j_endof ; ++j )
//     for (size_t k = Lo.k_begin; k < Lo.k_endof ; ++k )
//     {
//         const int icel = i*nz*ny + j*nz + k;
//         const int xn  = gridA.NEIBcell[icel*12+0]; /*-x negative */ const int xnn = gridA.NEIBcell[icel*12+6]; //-2x
//         const int xp  = gridA.NEIBcell[icel*12+1]; /*+x positive */ const int xpp = gridA.NEIBcell[icel*12+7]; //+2x 
//         const int yn  = gridA.NEIBcell[icel*12+2]; /*-y          */ const int ynn = gridA.NEIBcell[icel*12+8]; //-2y
//         const int yp  = gridA.NEIBcell[icel*12+3]; /*+y          */ const int ypp = gridA.NEIBcell[icel*12+9]; //+2y
//         const int zn  = gridA.NEIBcell[icel*12+4]; /*-z          */ const int znn = gridA.NEIBcell[icel*12+10];//-2z
//         const int zp  = gridA.NEIBcell[icel*12+5]; /*+z          */ const int zpp = gridA.NEIBcell[icel*12+11];//+2z

//         const double Tpx = (T3.u[xp]); 
//         const double Tnx = (T3.u[icel]); 
//         const double Tpy = (T3.v[yp]);
//         const double Tny = (T3.v[icel]);
//         const double Tpz = (T3.v[zp]);
//         const double Tnz = (T3.w[icel]);

//         auto  [TpX, TnX] = EngeryQuickUnit (Tpx, Tnx, xpp, xp, icel, xn, xnn, gridA.Dx, gridA.Dxs, i, t1.T);
//         auto  [TpY, TnY] = EngeryQuickUnit (Tpy, Tny, ypp, yp, icel, yn, ynn, gridA.Dy, gridA.Dys, j, t1.T);
//         auto  [TpZ, TnZ] = EngeryQuickUnit (Tpz, Tnz, zpp, zp, icel, zn, znn, gridA.Dz, gridA.Dzs, k, t1.T);

//         auto convection =   -simu.dt*( (TpX*TpX-TnX*TnX) / gridA.Dx[i]     // FIXME   Convection Term :: Quick Scheme
//                                       +(TpY*TpY-TnY*TnY) / gridA.Dy[j]  
//                                       +(TpZ*TpZ-TnZ*TnZ) / gridA.Dz[k]);

//         // * ----------------------- convetion -----------------------

//         // * ----------------------- difussion -----------------------
//         auto difussion = simu.dt * simu.alpha *(
//                         ( (t1.T[xp] - t1.T[icel])/gridA.Dxs[i] 
//                         - (t1.T[icel] -t1.T[xn])/gridA.Dxs[i-1] ) / gridA.Dx[i] +
//                         ( (t1.T[yp] - t1.T[icel])/gridA.Dxs[j] 
//                         - (t1.T[icel] -t1.T[yn])/gridA.Dxs[j-1] ) / gridA.Dy[i] +
//                         ( (t1.T[zp] - t1.T[icel])/gridA.Dxs[k] 
//                         - (t1.T[icel] -t1.T[zn])/gridA.Dxs[k-1] ) / gridA.Dz[i] ) ;
//         // * ----------------------- difussion -----------------------

//         t1.T[icel] = t1.T[icel]+ difussion + difussion;

//     }
// }

// Fe_phi_e = u_P * std::max(Fe, 0) - u_E *0 std::max(-Fe, 0);



// !BACKUP
// // -----------------------x
//         if(Tpx > 0)
//             TpX = 0.5*( t1.T[xp] + t1.T[icel] ) 
//                         -0.125*gridA.Dxs[i]*gridA.Dxs[i]/gridA.Dx[i]*
//                         (( (t1.T[xp]-t1.T[icel])/gridA.Dxs[i] 
//                         - (t1.T[icel] - t1.T[xn] )/gridA.Dxs[i-1]));
//         else
//             TpX = 0.5*( t1.T[xp] + t1.T[icel] ) 
//                         -0.125*gridA.Dxs[i]*gridA.Dxs[i]/gridA.Dx[i]*
//                         (( (t1.T[xpp]-t1.T[xp])/gridA.Dxs[i+1] 
//                         - (t1.T[xp] - t1.T[icel] )/gridA.Dxs[i]));


//         if(Tnx > 0)
//             TnX = 0.5*( t1.T[icel] + t1.T[xn] ) 
//                         -0.125*gridA.Dxs[i]*gridA.Dxs[i]/gridA.Dx[i]*
//                         (( (t1.T[icel]-t1.T[zn])/gridA.Dxs[i-1] 
//                         - (t1.T[xn] - t1.T[xnn] )/gridA.Dxs[i-2]));
//         else
//             TnX = 0.5*( t1.T[icel] + t1.T[nx] ) 
//                         -0.125*gridA.Dxs[i]*gridA.Dxs[i]/gridA.Dx[i]*
//                         (( (t1.T[xp]-t1.T[icel])/gridA.Dxs[i] 
//                         - (t1.T[icel] - t1.T[xn] )/gridA.Dxs[i-1]));

// // -----------------------y

//         if(Tpy > 0)
//             TpY = 0.5*( t1.T[yp] + t1.T[icel] )
//                         -0.125*gridA.Dys[j]*gridA.Dys[j]/gridA.Dy[j]*
//                         (( (t1.T[yp]-t1.T[icel])/gridA.Dys[j] 
//                         - (t1.T[icel] - t1.T[yn] )/gridA.Dys[j-1]));
//         else
//             TpY = 0.5*( t1.T[yp] + t1.T[icel] )
//                         -0.125*gridA.Dys[j]*gridA.Dys[j]/gridA.Dy[j]*
//                         (( (t1.T[ypp]-t1.T[yp])/gridA.Dys[j+1] 
//                         - (t1.T[yp] - t1.T[icel] )/gridA.Dys[j]));


//         if(Tny > 0)
//             TnY = 0.5*( t1.T[icel] + t1.T[yn] )
//                         -0.125*gridA.Dys[j]*gridA.Dys[j]/gridA.Dy[j]*
//                         (( (t1.T[icel]-t1.T[yn])/gridA.Dys[j-1] 
//                         - (t1.T[yn] - t1.T[ynn] )/gridA.Dys[j-2]));
//         else
//             TnY = 0.5*( t1.T[icel] + t1.T[yn] )
//                         -0.125*gridA.Dys[j]*gridA.Dys[j]/gridA.Dy[j]*
//                         (( (t1.T[yp]-t1.T[icel])/gridA.Dys[j] 
//                         - (t1.T[icel] - t1.T[yn] )/gridA.Dys[j-1]));

// // -----------------------z

//         if(Tpz > 0)
//             TpZ = 0.5*( t1.T[zp] + t1.T[icel] )
//                         -0.125*gridA.Dzs[k]*gridA.Dzs[k]/gridA.Dz[k]*
//                         (( (t1.T[zp]-t1.T[icel])/gridA.Dzs[k] 
//                         - (t1.T[icel] - t1.T[zn] )/gridA.Dzs[k-1]));
//         else
//             TpZ = 0.5*( t1.T[zp] + t1.T[icel] )
//                         -0.125*gridA.Dzs[k]*gridA.Dzs[k]/gridA.Dz[k]*
//                         (( (t1.T[zpp]-t1.T[zp])/gridA.Dzs[k+1] 
//                         - (t1.T[zp] - t1.T[icel] )/gridA.Dzs[k]));


//         if(Tnz > 0)
//             TnZ = 0.5*( t1.T[icel] + t1.T[zn] )
//                         -0.125*gridA.Dzs[k]*gridA.Dzs[k]/gridA.Dz[k]*
//                         (( (t1.T[icel]-t1.T[zn])/gridA.Dzs[k-1] 
//                         - (t1.T[zn] - t1.T[znn] )/gridA.Dzs[k-2]));
//         else
//             TnZ = 0.5*( t1.T[icel] + t1.T[zn] )
//                         -0.125*gridA.Dzs[k]*gridA.Dzs[k]/gridA.Dz[k]*
//                         (( (t1.T[zp]-t1.T[icel])/gridA.Dzs[k] 
//                         - (t1.T[icel] - t1.T[zn] )/gridA.Dzs[k-1]));