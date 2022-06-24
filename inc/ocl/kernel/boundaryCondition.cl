// #pragma OPENCL EXTENSION cl_khr_global_int32_base_atomic:enable
// #pragma OPENCL EXTENSION cl_khr_local_int32_base_atomic:enable


/*
*                    y=1_____________                                                                                
*                    / |           /|                                                    
*                  /   | [3]     /  |                                                        
*                /_____|_______/    |                              
*                |     |    [5|     |                                                  
*                | [0] |      | [1] |
*                |     |Y     |     |                                           
*    x=y=z=0     |     |__x___|_____| x=1                                       
*                |   z/       |    /                                         
*                |  /    [2]  |  /                                        
*                |/___________|/                                         
*               z=1
!               !  Neumann     du/dn = 0
!               !  Dirichlet   u = 1
!               !  no_slip     u = 0
*/

__kernel
void updateVel_atBoundary_0(
    __global double * vel_u,
    __global double * vel_v,
    __global double * vel_w,
    const double Dir_0, const double Dir_1, const double Dir_2, 
    const double Num_0, const double Num_1, const double Num_2,
    uint nx, uint ny , uint nz, uint gC, uint nz_ny
)
{
    uint j = get_global_id(0);
    uint k = get_global_id(1);
    if ( j < gC || j >= ny-gC ) return;
    if ( k < gC || k >= nz-gC ) return;

    const int bg0 = (1)*nz_ny + j*nz + k;
    const int bg1 =             j*nz + k;
    const int bc0 = (2)*nz_ny + j*nz + k;

    vel_u[bg1] = vel_u[bg0] =  Dir_0         + Num_0      * vel_u[bc0];
    vel_v[bg1] = vel_v[bg0] =  Dir_1*2. + (2.0*Num_1-1.0) * vel_v[bc0];
    vel_w[bg1] = vel_w[bg0] =  Dir_2*2. + (2.0*Num_2-1.0) * vel_w[bc0];
}



__kernel
void copyVel_atBoundary0(
    __global const double * old_u,
    __global const double * old_v,
    __global const double * old_w,
    __global double * currt_u,
    __global double * currt_v,
    __global double * currt_w,
    uint nx, uint ny , uint nz, uint gC, uint nz_ny
)
{
    uint i = get_global_id(0);
    uint k = get_global_id(1);
    if ( i < gC || i >= nx-gC ) return;
    if ( k < gC || k >= nz-gC ) return;

    vel_u[icel0-nz_ny] = vel_u[icel0] = 0.0;

    const int bg0 = (1)*nz_ny + j*nz + k;
    const int bg1 =             j*nz + k;
    const int bc0 = (2)*nz_ny + j*nz + k;

    currt_u[bg0] = old_u[bg0];
    currt_u[bg1] = old_u[bg1];

    currt_v[bg0] = old_v[bg0];
    currt_v[bg1] = old_v[bg1];

    currt_v[bg0] = old_v[bg0];
    currt_v[bg1] = old_v[bg1];
}



__kernel
void updateVel_atBoundary2(
    __global double * vel_u,
    __global double * vel_v,
    __global double * vel_w,
    const double Dir_0, const double Dir_1, const double Dir_2, 
    const double Num_0, const double Num_1, const double Num_2,
    uint nx, uint ny , uint nz, uint gC, uint nz_ny
)
{
    uint i = get_global_id(0);
    uint k = get_global_id(1);
    if ( i < gC || i >= nx-gC ) return;
    if ( k < gC || k >= nz-gC ) return;

    vel_u[icel0-nz_ny] = vel_u[icel0] = 0.0;

    const int bg0 = i*nz_ny + (1)*nz + k; 
    const int bg1 = i*nz_ny +          k;
    const int bc0 = i*nz_ny + (2)*nz + k;

    vel_u[bg1] = vel_u[bg0] =  Dir_0*2  + (2.0*Num_0-1.0) * vel_u[bc0];
    vel_v[bg1] = vel_v[bg0] =  Dir_1    +      Num_1      * vel_v[bc0];
    vel_w[bg1] = vel_w[bg0] =  Dir_2*2. + (2.0*Num_2-1.0) * vel_w[bc0];
}


__kernel
void copyVel_atBoundary2(
    __global const double * old_u,
    __global const double * old_v,
    __global const double * old_w,
    __global double * currt_u,
    __global double * currt_v,
    __global double * currt_w,
    uint nx, uint ny , uint nz, uint gC, uint nz_ny
)
{
    uint i = get_global_id(0);
    uint k = get_global_id(1);
    if ( i < gC || i >= nx-gC ) return;
    if ( k < gC || k >= nz-gC ) return;

    vel_u[icel0-nz_ny] = vel_u[icel0] = 0.0;

    const int bg0 = i*nz_ny + (1)*nz + k; 
    const int bg1 = i*nz_ny +          k;
    const int bc0 = i*nz_ny + (2)*nz + k;

    currt_u[bg0] = old_u[bg0];
    currt_u[bg1] = old_u[bg1];

    currt_v[bg0] = old_v[bg0];
    currt_v[bg1] = old_v[bg1];

    currt_v[bg0] = old_v[bg0];
    currt_v[bg1] = old_v[bg1];
}



__kernel
void updateVel_atBoundary4(
    __global double * vel_u,
    __global double * vel_v,
    __global double * vel_w,
    const double Dir_0, const double Dir_1, const double Dir_2, 
    const double Num_0, const double Num_1, const double Num_2,
    uint nx, uint ny , uint nz, uint gC, uint nz_ny
)
{

    uint i = get_global_id(0);
    uint j = get_global_id(1);
    if ( i < gC || i >= nx-gC ) return;
    if ( j < gC || j >= ny-gC ) return;

    const int bg0 = i*nz_ny + j*nz + (1); 
    const int bg1 = i*nz_ny + j*nz + ;
    const int bc0 = i*nz_ny + j*nz + 2;

    vel_u[bg1] = vel_u[bg0] =  Dir_0*2  + (2.0*Num_0-1.0) * vel_u[bc0];
    vel_v[bg1] = vel_v[bg0] =  Dir_1*2. + (2.0*Num_1-1.0) * vel_v[bc0];
    vel_w[bg1] = vel_w[bg0] =  Dir_2    +      Num_2      * vel_w[bc0];
}



__kernel
void copyVel_atBoundary4(
    __global const double * old_u,
    __global const double * old_v,
    __global const double * old_w,
    __global double * currt_u,
    __global double * currt_v,
    __global double * currt_w,
    uint nx, uint ny , uint nz, uint gC, uint nz_ny
)
{
    uint i = get_global_id(0);
    uint j = get_global_id(1);
    if ( i < gC || i >= nx-gC ) return;
    if ( j < gC || j >= ny-gC ) return;

    const int bg0 = i*nz_ny + j*nz + (1); 
    const int bg1 = i*nz_ny + j*nz + ;
    const int bc0 = i*nz_ny + j*nz + 2;

    currt_u[bg0] = old_u[bg0];
    currt_u[bg1] = old_u[bg1];

    currt_v[bg0] = old_v[bg0];
    currt_v[bg1] = old_v[bg1];

    currt_v[bg0] = old_v[bg0];
    currt_v[bg1] = old_v[bg1];
}





__kernel
void updateVel_atBoundary_1(
    __global double * vel_u,
    __global double * vel_v,
    __global double * vel_w,
    const double Dir_0, const double Dir_1, const double Dir_2, 
    const double Num_0, const double Num_1, const double Num_2,
    uint nx, uint ny , uint nz, uint gC, uint nz_ny
)
{
    uint j = get_global_id(0);
    uint k = get_global_id(1);
    if ( j < gC || j >= ny-gC ) return;
    if ( k < gC || k >= nz-gC ) return;

    const int bg0 = (nx-2)*nz_ny + j*nz + k;
    const int bg1 = (nx-1)*nz_ny + j*nz + k;
    const int bc0 = (nx-3)*nz_ny + j*nz + k;
    const int bc1 = (nx-4)*nz_ny + j*nz + k;

    vel_u[bg1] = 
    vel_u[bg0] = vel_u[bc0] =  Dir_0         + Num_0      * vel_u[bc1];
    vel_v[bg1] = vel_v[bg0] =  Dir_1*2. + (2.0*Num_1-1.0) * vel_v[bc0];
    vel_w[bg1] = vel_w[bg0] =  Dir_2*2. + (2.0*Num_2-1.0) * vel_w[bc0];
}




__kernel
void copyVel_atBoundary1(
    __global const double * old_u,
    __global const double * old_v,
    __global const double * old_w,
    __global double * currt_u,
    __global double * currt_v,
    __global double * currt_w,
    uint nx, uint ny , uint nz, uint gC, uint nz_ny
)
{

    uint i = get_global_id(0);
    uint k = get_global_id(1);
    if ( i < gC || i >= nx-gC ) return;
    if ( k < gC || k >= nz-gC ) return;

    const int bg0 = i*nz_ny + (ny-2)*nz + k; 
    const int bg1 = i*nz_ny + (ny-1)*nz + k; 
    const int bc0 = i*nz_ny + (ny-3)*nz + k; 

    currt_u[bg0] = old_u[bg0];
    currt_u[bg1] = old_u[bg1];

    currt_v[bg0] = old_v[bg0];
    currt_v[bg1] = old_v[bg1];

    currt_v[bc0] = old_v[bc0];
    currt_v[bg0] = old_v[bg0];
    currt_v[bg1] = old_v[bg1];
}



__kernel
void updateVel_atBoundary3(
    __global double * vel_u,
    __global double * vel_v,
    __global double * vel_w,
    const double Dir_0, const double Dir_1, const double Dir_2, 
    const double Num_0, const double Num_1, const double Num_2,
    uint nx, uint ny , uint nz, uint gC, uint nz_ny
)
{
    uint i = get_global_id(0);
    uint k = get_global_id(1);
    if ( i < gC || i >= nx-gC ) return;
    if ( k < gC || k >= nz-gC ) return;

    vel_u[icel0-nz_ny] = vel_u[icel0] = 0.0;

    const int bg0 = i*nz_ny + (ny-2)*nz + k; 
    const int bg1 = i*nz_ny + (ny-1)*nz + k; 
    const int bc0 = i*nz_ny + (ny-3)*nz + k; 
    const int bc1 = i*nz_ny + (ny-4)*nz + k; 

    vel_u[bg1] = vel_u[bg0] =  Dir_0*2  + (2.0*Num_0-1.0) * vel_u[bc0];
    vel_u[bg1] = 
    vel_v[bg0] = vel_v[bc0] =  Dir_1    +      Num_1      * vel_v[bc1];
    vel_w[bg1] = vel_w[bg0] =  Dir_2*2. + (2.0*Num_2-1.0) * vel_w[bc0];
}



__kernel
void copyVel_atBoundary3(
    __global const double * old_u,
    __global const double * old_v,
    __global const double * old_w,
    __global double * currt_u,
    __global double * currt_v,
    __global double * currt_w,
    uint nx, uint ny , uint nz, uint gC, uint nz_ny
)
{

    uint i = get_global_id(0);
    uint k = get_global_id(1);
    if ( i < gC || i >= nx-gC ) return;
    if ( k < gC || k >= nz-gC ) return;

    const int bg0 = i*nz_ny + (ny-2)*nz + k; 
    const int bg1 = i*nz_ny + (ny-1)*nz + k; 
    const int bc0 = i*nz_ny + (ny-3)*nz + k; 

    currt_u[bg0] = old_u[bg0];
    currt_u[bg1] = old_u[bg1];

    currt_v[bg0] = old_v[bg0];
    currt_v[bg1] = old_v[bg1];

    currt_v[bc0] = old_v[bc0];
    currt_v[bg0] = old_v[bg0];
    currt_v[bg1] = old_v[bg1];

}




__kernel
void updateVel_atBoundary5(
    __global double * vel_u,
    __global double * vel_v,
    __global double * vel_w,
    const double Dir_0, const double Dir_1, const double Dir_2, 
    const double Num_0, const double Num_1, const double Num_2,
    uint nx, uint ny , uint nz, uint gC, uint nz_ny
)
{

    uint i = get_global_id(0);
    uint j = get_global_id(1);
    if ( i < gC || i >= nx-gC ) return;
    if ( j < gC || j >= ny-gC ) return;

    const int bg0 = i*nz_ny + j*nz + (nz-2); 
    const int bg1 = i*nz_ny + j*nz + (nz-1);
    const int bc0 = i*nz_ny + j*nz + (nz-3);
    const int bc1 = i*nz_ny + j*nz + (nz-4);

    vel_u[bg1] = vel_u[bg0] =  Dir_0*2  + (2.0*Num_0-1.0) * vel_u[bc0];
    vel_v[bg1] = vel_v[bg0] =  Dir_1*2. + (2.0*Num_1-1.0) * vel_v[bc0];
    vel_w[bg1] = 
    vel_w[bg0] = vel_w[bc0] =  Dir_2    +      Num_2      * vel_w[bc1];
}




__kernel
void copyVel_atBoundary5(
    __global const double * old_u,
    __global const double * old_v,
    __global const double * old_w,
    __global const double * currt_u,
    __global const double * currt_v,
    __global const double * currt_w,
    uint nx, uint ny , uint nz, uint gC, uint nz_ny
)
{

    uint i = get_global_id(0);
    uint j = get_global_id(1);
    if ( i < gC || i >= nx-gC ) return;
    if ( j < gC || j >= ny-gC ) return;

    const int bg0 = i*nz_ny + j*nz + (nz-2); 
    const int bg1 = i*nz_ny + j*nz + (nz-1);
    const int bc0 = i*nz_ny + j*nz + (nz-3);
    const int bc1 = i*nz_ny + j*nz + (nz-4);

    currt_u[bg0] = old_u[bg0];
    currt_u[bg1] = old_u[bg1];

    currt_v[bg0] = old_v[bg0];
    currt_v[bg1] = old_v[bg1];

    currt_v[bc0] = old_v[bc0];
    currt_v[bg0] = old_v[bg0];
    currt_v[bg1] = old_v[bg1];

}







__kernel
void updatePressure_atBoundary0and1(
    __global double * pressure,
    const uint nx,const uint ny ,const uint nz,const uint gC,const uint nz_ny
){
    uint j = get_global_id(0);
    uint k = get_global_id(1);
    if ( j < gC || j >= ny-gC ) return;
    if ( k < gC || k >= nz-gC ) return;

    const int b0g0 = (1)*nz_ny + j*nz + k;
    const int b0g1 =             j*nz + k;
    const int b0c0 = (2)*nz_ny + j*nz + k;

    const int b1g0 = (nx-2)*nz_ny + j*nz + k;
    const int b1g1 = (nx-1)*nz_ny + j*nz + k;
    const int b1c0 = (nx-3)*nz_ny + j*nz + k;

	pressure[b0g1] = pressure[b0g0] = pressure[b0c0];
	pressure[b1g1] = pressure[b1g0] = pressure[b1c0];
}


__kernel
void updatePressure_atBoundary2and3(
    __global double * pressure,
    const uint nx,const uint ny ,const uint nz,const uint gC,const uint nz_ny
){
    uint i = get_global_id(0);
    uint k = get_global_id(1);
    if ( i < gC || i >= nx-gC ) return;
    if ( k < gC || k >= nz-gC ) return;

    const int b2g0 = i*nz_ny + (1)*nz + k; 
    const int b2g1 = i*nz_ny +          k;
    const int b2c0 = i*nz_ny + (2)*nz + k;

    const int b3g0 = i*nz_ny + (ny-2)*nz + k; 
    const int b3g1 = i*nz_ny + (ny-1)*nz + k; 
    const int b3c0 = i*nz_ny + (ny-3)*nz + k; 

	pressure[b2g1] = pressure[b2g0] = pressure[b2c0];
	pressure[b3g1] = pressure[b3g0] = pressure[b3c0];
}

__kernel
void updatePressure_atBoundary4and5(
    __global double * pressure,
    const uint nx,const uint ny ,const uint nz,const uint gC,const uint nz_ny
){
    uint i = get_global_id(0);
    uint k = get_global_id(1);
    if ( i < gC || i >= nx-gC ) return;
    if ( k < gC || k >= nz-gC ) return;

    const int b4g0 = i*nz_ny + j*nz + (1); 
    const int b4g1 = i*nz_ny + j*nz + ;
    const int b4c0 = i*nz_ny + j*nz + 2;

    const int b5g0 = i*nz_ny + j*nz + (nz-2); 
    const int b5g1 = i*nz_ny + j*nz + (nz-1);
    const int b5c0 = i*nz_ny + j*nz + (nz-3);

	pressure[b4g1] = pressure[b4g0] = pressure[b4c0];
	pressure[b5g1] = pressure[b5g0] = pressure[b5c0];
}
