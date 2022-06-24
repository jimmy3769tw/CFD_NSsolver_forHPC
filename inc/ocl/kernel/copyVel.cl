// #pragma OPENCL EXTENSION cl_khr_global_int32_base_atomic:enable
// #pragma OPENCL EXTENSION cl_khr_local_int32_base_atomic:enable

__kernel
void copyVel(
    __global double *Curr_u,
    __global double *Curr_v,
    __global double *Curr_w,
    __global const double *old_u,
    __global const double *old_v,
    __global const double *old_w,
    const uint nx,const uint ny,const uint nz,
    const uint gCells
)
{
    uint i = get_global_id(0);
    uint j = get_global_id(1);
    uint k = get_global_id(2);

    if ( i < 0 || i >= nx ) return;
    if ( j < 0 || j >= ny ) return;
    if ( k < 0 || k >= nz ) return;
    const int icel = i*nz*ny + j*nz + k; 

    T0_u[icel] = T3_u[icel] ;
    T0_v[icel] = T3_v[icel] ;
    T0_w[icel] = T3_w[icel] ;
}
