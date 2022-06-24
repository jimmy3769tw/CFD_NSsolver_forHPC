#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomic:enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomic:enable

void update_UandF_x(
    __global const double *old_u,
    __global const double *old_v,
    __global const double *old_w,
    __global double *Currt_u,
    __global double *Currt_v,
    __global double *Currt_w,
    __global const double *Dxs,
    __global const double *Dys,
    __global const double *Dzs,
    __global const double *pressure,
    __global double *Dfib_f,
    __global const double *Dfib_eta,
    const double u_solid,
    const double v_solid,
    const double w_solid,
    const double dt,
    const double nu,
    const uint nx, const uint ny,const uint nz,const uint gC, const uint iceltot
)
{
    uint i = get_global_id(0);
    uint j = get_global_id(1);
    uint k = get_global_id(2);

    if ( i < gC || i >= nx-gC-1 ) return;
    if ( j < gC || j >= ny-gC ) return;
    if ( k < gC || k >= nz-gC ) return;

    uint icel = i*nz*ny + j*nz + k;

    // ! ---------------  x  ---------------
    double T2u = old_u[icel] - dt*( pressure[icel+(nz*ny)]-pressure[icel] ) / Dxs[i];

    Currt_u[icel] = Dfib_eta[icel] * u_solid 
        + (1.0-0.5*( Dfib_eta[icel] + Dfib_eta[icel+(nz*ny)])) * T2u;

    Dfib_f[icel] = (Currt_u[icel] - T2u) / dt;
}


__kernel
void update_UandF_y(
    __global const double *old_u,
    __global const double *old_v,
    __global const double *old_w,
    __global double *Currt_u,
    __global double *Currt_v,
    __global double *Currt_w,
    __global const double *Dxs,
    __global const double *Dys,
    __global const double *Dzs,
    __global const double *pressure,
    __global double *Dfib_f,
    __global const double *Dfib_eta,
    const double u_solid,
    const double v_solid,
    const double w_solid,
    const double dt,
    const double nu,
    const uint nx, const uint ny,const uint nz,const uint gC, const uint iceltot
)
{
    uint i = get_global_id(0);
    uint j = get_global_id(1);
    uint k = get_global_id(2);

    if ( i < gC || i >= nx-gC ) return;
    if ( j < gC || j >= ny-gC-1 ) return;
    if ( k < gC || k >= nz-gC ) return;

    uint icel = i*nz*ny + j*nz + k;

    double T2v = old_v[icel] - dt*( pressure[icel+nz]-pressure[icel] ) / Dys[j];

    Currt_v[icel] = Dfib_eta[icel] * v_solid + (1.0-0.5*(Dfib_eta[icel]+Dfib_eta[icel+(1*nz)])) * T2v;

    Dfib_f[icel+iceltot]  = (Currt_v[icel] - T2v) / dt;
}


__kernel
void update_UandF_z(
    __global const double *old_u,
    __global const double *old_v,
    __global const double *old_w,
    __global double *Currt_u,
    __global double *Currt_v,
    __global double *Currt_w,
    __global const double *Dxs,
    __global const double *Dys,
    __global const double *Dzs,
    __global const double *pressure,
    __global double *Dfib_f,
    __global const double *Dfib_eta,
    const double u_solid,
    const double v_solid,
    const double w_solid,
    const double dt,
    const double nu,
    const uint nx, const uint ny,const uint nz,const uint gC, const uint iceltot
)
{
    uint i = get_global_id(0);
    uint j = get_global_id(1);
    uint k = get_global_id(2);

    if ( i < gC || i >= nx-gC ) return;
    if ( j < gC || j >= ny-gC ) return;
    if ( k < gC || k >= nz-gC-1 ) return;

    uint icel = i*nz*ny + j*nz + k;

    double T2w = old_w[icel] - dt*( pressure[icel+1]-pressure[icel] ) / Dzs[k];

    Currt_w[icel] = Dfib_eta[icel] * w_solid + (1.0-0.5*(Dfib_eta[icel]+Dfib_eta[icel+1])) * T2w;

    Dfib_f[icel*2*iceltot]  = (Currt_w[icel]- T2w) / dt;
}



