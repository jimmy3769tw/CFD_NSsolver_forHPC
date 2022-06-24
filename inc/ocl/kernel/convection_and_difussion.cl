#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomic:enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomic:enable


__kernel
void QuickUnitA(
    __local double * phiP,
    __local double * phiN,
    const double Up,
    const double Un,
    const int pp, const int p, const int icel, const int n, const int nn,
    __global double * Di,
    __global double * Ds,
    const uint Idx,
    __global double * phi
){
    double phi_P, phi_N;
        if (Up > 0.0)
            phi_P  = 0.5*(phi[icel] + phi[p] )
                    -0.125*Di[Idx+1]*Di[Idx+1]/Ds[Idx]*
                    ( (phi[p] - phi[icel] ) / Di[Idx+1] 
                    - (phi[icel] - phi[n] ) / Di[Idx]) ;
        else
            phi_P  = 0.5*(phi[icel] + phi[p] )
                    - 0.125*Di[Idx+1]*Di[Idx+1]/Ds[Idx+1]*
                    ( (phi[pp] - phi[p] )   / Di[Idx+2] 
                    - (phi[p] - phi[icel] ) / Di[Idx+1] );

        if (Un > 0.0)
            phi_N  = 0.5*(phi[n] + phi[icel] ) 
                    -0.125*Di[Idx]*Di[Idx]/Ds[Idx-1] *
                    ( (phi[icel] - phi[n] ) / Di[Idx] 
                    - (phi[n]   - phi[nn] ) / Di[Idx-1] ) ;
        else
            phi_N  = 0.5*(phi[n] + phi[icel] ) 
                    -0.125*Di[Idx]*Di[Idx]/Ds[Idx]*
                    ( (phi[p] - phi[icel] ) / Di[Idx+1] 
                    - (phi[icel] - phi[n] ) / Di[Idx])  ;
}



__kernel
void QuickUnitB(
    __local double * phiP,
    __local double * phiN,
    const double Up,
    const double Un,
    const int pp, const int p, const int icel, const int n, const int nn,
    __global double * Di,
    __global double * Ds,
    const uint Idx,
    __global double * phi
){
    double phi_P, phi_N;
        if (Up > 0.0)
            phi_P  = 0.5*(phi[icel] + phi[p] )
                    -0.125*Ds[Idx]*Ds[Idx]/Di[Idx]*
                    ( (phi[p] - phi[icel] ) / Ds[Idx] 
                    - (phi[icel] - phi[n] ) / Ds[Idx-1]) ;
        else
            phi_P  = 0.5*(phi[icel] + phi[p] )
                    - 0.125*Ds[Idx]*Ds[Idx]/Di[Idx+1]*
                    ( (phi[pp] - phi[p] )   / Ds[Idx+1] 
                    - (phi[p] - phi[icel] ) / Ds[Idx] );

        if (Un > 0.0)
            phi_N  = 0.5*(phi[n] + phi[icel] ) 
                    -0.125*Ds[Idx-1]*Ds[Idx-1]/Di[Idx-1] *
                    ( (phi[icel] - phi[n] ) / Ds[Idx-1] 
                    - (phi[n]   - phi[nn] ) / Ds[Idx-2] );
        else
            phi_N  = 0.5*(phi[n] + phi[icel] ) 
                    -0.125*Ds[Idx-1]*Ds[Idx-1]/Di[Idx] *
                    ( (phi[p] - phi[icel] ) / Ds[Idx] 
                    - (phi[icel] - phi[n] ) / Ds[Idx-1]);
}



    
__kernel
void CenterQuickScheme_u(
    __global const double *oldVel_u,
    __global const double *oldVel_v,
    __global const double *oldVel_w,
    __global double *curtVel_u,
    __global double *curtVel_v,
    __global double *curtVel_w,
    __global const double *Dx,
    __global const double *Dy,
    __global const double *Dz,
    __global const double *Dxs,
    __global const double *Dys,
    __global const double *Dzs,
    __global int *NEIBcell,
    const double dt,const double nu,
    const int nx,const int ny ,const int nz,const int gCells
)
{
    uint i = get_global_id(0);
    uint j = get_global_id(1);
    uint k = get_global_id(2);
    if ( i < gCells || i >= nx-gCells-1 ) return;
    if ( j < gCells || j >= ny-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;

        const int icelOffSet =  (i-2)*(nz-4)*(ny-4) + (j-2)*(nz-4) + (k-2);
        int icel = i*nz*ny + j*nz + k;
        // * -------------- convection term --------------
        const int xm  = NEIBcell[icel*12+0]; /*-x negative */ const int xmm = NEIBcell[icel*12+6]; //-2x
        const int xp  = NEIBcell[icel*12+1]; /*+x positive */ const int xpp = NEIBcell[icel*12+7]; //+2x 
        const int ym  = NEIBcell[icel*12+2]; /*-y          */ const int ymm = NEIBcell[icel*12+8]; //-2y
        const int yp  = NEIBcell[icel*12+3]; /*+y          */ const int ypp = NEIBcell[icel*12+9]; //+2y
        const int zm  = NEIBcell[icel*12+4]; /*-z          */ const int zmm = NEIBcell[icel*12+10];//-2z
        const int zp  = NEIBcell[icel*12+5]; /*+z          */ const int zpp = NEIBcell[icel*12+11];//+2z

        const double IT = (Dxs[i] * 0.5) / Dx[i];
        const double Up = 0.5 * (oldVel_u[xp]+oldVel_u[ic]); //main 
        const double Un = 0.5 * (oldVel_u[xm]+oldVel_u[ic]); //main

        const double Vp = oldVel_v[ic]+(oldVel_v[xp]-oldVel_v[ic])*IT;
        const double Vn = oldVel_v[ym]+(oldVel_v[ym+xp-ic]-oldVel_v[ym])*IT;

        const double Wp = oldVel_w[ic]+(oldVel_w[xp]-oldVel_w[ic])*IT;
        const double Wn = oldVel_w[zm]+(oldVel_w[zm+xp-ic]-oldVel_w[zm])*IT;


        __local double uPx, uNx , uPy, uNy , uPz, uNz ;

        QuickUnitA (&uPx, &uNx, Up, Un, xpp, xp, icel, xm, xmm, Dx, Dxs, i, oldVel_u);
        QuickUnitB (&uPy, &uNy, Vp, Vn, ypp, yp, icel, ym, ymm, Dy, Dys, j, oldVel_u);
        QuickUnitB (&uPz, &uNz, Wp, Wn, zpp, zp, icel, zm, zmm, Dz, Dzs, k, oldVel_u);


        double convection = -dt*(
                         (uPx*uPx-uNx*uNx) / Dxs[i]
                        +(uPy*Vp-uNy*Vn) / Dy[j]
                        +(uPz*Wp-uNz*Wn) / Dz[k]);

        // * -------------- convection term --------------

        // * -------------- difussion term --------------
        double difussion =    nu*dt*(
                                    ( (oldVel_u[xp]-oldVel_u[icel]) / Dx[i+1] 
                                    - (oldVel_u[icel]-oldVel_u[xm]) / Dx[i]   ) / Dxs[i]+
                                    ( (oldVel_u[yp]-oldVel_u[icel]) / Dys[j]   
                                    - (oldVel_u[icel]-oldVel_u[ym]) / Dys[j-1] ) / Dy[j]+
                                    ( (oldVel_u[zp]-oldVel_u[icel]) / Dzs[k]   
                                    - (oldVel_u[icel]-oldVel_u[zm]) / Dzs[k-1] ) / Dz[k]);
        // * -------------- difussion term --------------

        curtVel_u[icel] = oldVel_u[icel] + difussion + convection ;
}


__kernel
void CenterQuickScheme_v(
    __global const double *oldVel_u,
    __global const double *oldVel_v,
    __global const double *oldVel_w,
    __global double *curtVel_u,
    __global double *curtVel_v,
    __global double *curtVel_w,
    __global const double *Dx,
    __global const double *Dy,
    __global const double *Dz,
    __global const double *Dxs,
    __global const double *Dys,
    __global const double *Dzs,
    __global int *NEIBcell,
    const double dt,const double nu,
    const int nx,const int ny ,const int nz,const int gCells
)
{
    uint i = get_global_id(0);
    uint j = get_global_id(1);
    uint k = get_global_id(2);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( j < gCells || j >= ny-gCells-1 ) return;
    if ( k < gCells || k >= nz-gCells ) return;

        const int icelOffSet =  (i-2)*(nz-4)*(ny-4) + (j-2)*(nz-4) + (k-2);
        int icel = i*nz*ny + j*nz + k;
        // * -------------- convection term --------------
        const int xm  = NEIBcell[icel*12+0]; /*-x negative */ const int xmm = NEIBcell[icel*12+6]; //-2x
        const int xp  = NEIBcell[icel*12+1]; /*+x positive */ const int xpp = NEIBcell[icel*12+7]; //+2x 
        const int ym  = NEIBcell[icel*12+2]; /*-y          */ const int ymm = NEIBcell[icel*12+8]; //-2y
        const int yp  = NEIBcell[icel*12+3]; /*+y          */ const int ypp = NEIBcell[icel*12+9]; //+2y
        const int zm  = NEIBcell[icel*12+4]; /*-z          */ const int zmm = NEIBcell[icel*12+10];//-2z
        const int zp  = NEIBcell[icel*12+5]; /*+z          */ const int zpp = NEIBcell[icel*12+11];//+2z


        const double IT = (Dys[y ] * 0.5) / Dy[i];

        const double Up = 0.5 * (oldVel_u[xp]+oldVel_u[ic]); //main 
        const double Un = 0.5 * (oldVel_u[xm]+oldVel_u[ic]); //main
        const double Vp = oldVel_v[ic]+(oldVel_v[xp]-oldVel_v[ic])*IT;
        const double Vn = oldVel_v[ym]+(oldVel_v[ym+xp-ic]-oldVel_v[ym])*IT;
        const double Wp = oldVel_w[ic]+(oldVel_w[xp]-oldVel_w[ic])*IT;
        const double Wn = oldVel_w[zm]+(oldVel_w[zm+xp-ic]-oldVel_w[zm])*IT;

        __local double  vPx, vNx, vPy, vNy, vPz, vNz;
        QuickUnitB (&vPx, &vNx, Up, Un, xpp, xp, icel, xm, xmm, Dx, Dxs, i, oldVel_v);
        QuickUnitA (&vPy, &vNy, Vp, Vn, ypp, yp, icel, ym, ymm, Dy, Dys, j, oldVel_v);
        QuickUnitB (&vPz, &vNz, Wp, Wn, zpp, zp, icel, zm, zmm, Dz, Dzs, k, oldVel_v);

        double convection =   -dt*(
                                 (Up*vPx-Un*vNx) / Dx[i] 
                                +(vPy*vPy-vNy*vNy) / Dys[j] 
                                +(Wp*vPz-Wn*vNz) / Dz[k]); 

    // * -------------- convection term --------------

    // * -------------- difussion term --------------
        const double difussion = nu*dt*( 
                                ( (oldVel_v[xp] - oldVel_v[icel]) / Dxs[i]
                                - (oldVel_v[icel] - oldVel_v[xm]) / Dxs[i-1] ) / Dx[i] +
                                ( (oldVel_v[yp] - oldVel_v[icel]) / Dy[j+1] 
                                - (oldVel_v[icel] - oldVel_v[ym]) / Dy[j] ) / Dys[j] +
                                ( (oldVel_v[zp] - oldVel_v[icel]) / Dzs[k] 
                                - (oldVel_v[icel] - oldVel_v[zm]) / Dzs[k-1])/ Dz[k] );
    // * -------------- difussion term --------------

        const double temp = convection + difussion;

        curtVel_v[icel] = oldVel_v[icel] + temp;
}



__kernel
void CenterQuickScheme_w(
    __global const double *oldVel_u,
    __global const double *oldVel_v,
    __global const double *oldVel_w,
    __global double *curtVel_u,
    __global double *curtVel_v,
    __global double *curtVel_w,
    __global const double *Dx,
    __global const double *Dy,
    __global const double *Dz,
    __global const double *Dxs,
    __global const double *Dys,
    __global const double *Dzs,
    __global int *NEIBcell,
    const double dt,const double nu,
    const int nx,const int ny ,const int nz,const int gCells
)
{

    uint i = get_global_id(0);
    uint j = get_global_id(1);
    uint k = get_global_id(2);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( j < gCells || j >= ny-gCells ) return;
    if ( k < gCells || k >= nz-gCells-1 ) return;

        const int icelOffSet =  (i-2)*(nz-4)*(ny-4) + (j-2)*(nz-4) + (k-2);
        const int icel = i*nz*ny + j*nz + k;
        const int xp  = NEIBcell[icel*12+1]; /*+x positive */ const int xpp = NEIBcell[icel*12+7]; //+2x 
        const int xm  = NEIBcell[icel*12+0]; /*-x negative */ const int xmm = NEIBcell[icel*12+6]; //-2x
        const int ym  = NEIBcell[icel*12+2]; /*-y          */ const int ymm = NEIBcell[icel*12+8]; //-2y
        const int yp  = NEIBcell[icel*12+3]; /*+y          */ const int ypp = NEIBcell[icel*12+9]; //+2y
        const int zm  = NEIBcell[icel*12+4]; /*-z          */ const int zmm = NEIBcell[icel*12+10];//-2z
        const int zp  = NEIBcell[icel*12+5]; /*+z          */ const int zpp = NEIBcell[icel*12+11];//+2z

        // * -------------- convection term --------------
        const double IT = (Dzs[k] * 0.5) / Dz[k];

        const double Wp = 0.5 * ( oldVel_w[zp]+oldVel_w[ic] );//main
        const double Wn = 0.5 * ( oldVel_w[zm]+oldVel_w[ic] );//main
        const double Up = oldVel_u[ic]+(oldVel_u[zp]-oldVel_u[ic])*IT;
        const double Un = oldVel_u[xm]+(oldVel_u[xm+zp-ic]-oldVel_u[xm])*IT; 
        const double Vp = oldVel_v[ic]+(oldVel_v[zp]-oldVel_v[ic])*IT;;
        const double Vn = oldVel_v[ym]+(oldVel_v[ym+zp-ic]-oldVel_v[ym])*IT;;

        __local double  wPx, wNx, wPy, wNy, wPz, wNz;

        QuickUnitB (&wPx, &wNx, Up, Un, xpp, xp, icel, xm, xmm, Dx, Dxs, i, oldVel_w);
        QuickUnitB (&wPy, &wNy, Vp, Vn, ypp, yp, icel, ym, ymm, Dy, Dys, j, oldVel_w);
        QuickUnitA (&wPz, &wNz, Wp, Wn, zpp, zp, icel, zm, zmm, Dz, Dzs, k, oldVel_w);


        const double convection =   -dt*(
                                     (Up*wPx-Un*wNx) / Dx[i]
                                    +(Vp*wPy-Vn*wNy) / Dy[j]
                                    +(wPz*wPz-wNz*wNz) / Dzs[k]);

        // * -------------- convection term --------------

        // * -------------- difussion term --------------
        const double difussion = nu*dt*(
                                    ( ( oldVel_w[xp]-oldVel_w[icel] ) / Dxs[i] 
                                    - ( oldVel_w[icel]-oldVel_w[xm] ) / Dxs[i-1] ) / Dx[i]+
                                    ( ( oldVel_w[yp]-oldVel_w[icel] ) / Dys[j] 
                                    - ( oldVel_w[icel]-oldVel_w[ym] ) / Dys[j-1] ) / Dy[j]+
                                    ( ( oldVel_w[zp]-oldVel_w[icel] ) / Dz[k+1] 
                                    - ( oldVel_w[icel]-oldVel_w[zm] ) / Dz[k]   ) / Dzs[k]) ;
        // * -------------- difussion term --------------

        curtVel_w[icel] = oldVel_w[icel] + convection + difussion;
}
