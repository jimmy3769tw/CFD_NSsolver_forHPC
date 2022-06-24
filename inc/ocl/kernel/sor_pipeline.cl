`// #pragma OPENCL EXTENSION cl_khr_global_int32_base_atomic:enable
// #pragma OPENCL EXTENSION cl_khr_local_int32_base_atomic:enable

#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics:enable

// TODO :1.get coef 

__kernel
void Sor_pipeline_getCoef(
    __global double *SorCoef,
    __global const double *Dx,
    __global const double *Dy,
    __global const double *Dz,
    __global const double *Dxs,
    __global const double *Dys,
    __global const double *Dzs,
    uint nx, uint ny , uint nz, uint gCells
){

    uint i = get_global_id(0);
    uint j = get_global_id(1);
    uint k = get_global_id(2);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( j < gCells || j >= ny-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;

    int icelOffSet =  (i-2)*(nz-4)*(ny-4) + (j-2)*(nz-4) + (k-2);
    SorCoef[icelOffSet*8  ] = Dy[j] * Dz[k] / Dxs[i]  ;
    SorCoef[icelOffSet*8+1] = Dy[j] * Dz[k] / Dxs[i-1];
    SorCoef[icelOffSet*8+2] = Dx[i] * Dz[k] / Dys[j]  ;
    SorCoef[icelOffSet*8+3] = Dx[i] * Dz[k] / Dys[j-1];
    SorCoef[icelOffSet*8+4] = Dx[i] * Dy[j] / Dzs[k]  ;
    SorCoef[icelOffSet*8+5] = Dx[i] * Dy[j] / Dzs[k-1];
    SorCoef[icelOffSet*8+7] =  1.0L / ( - Dy[j] * Dz[k] / Dxs[i] 
                                        - Dy[j] * Dz[k] / Dxs[i-1]
                                        - Dx[i] * Dz[k] / Dys[j] 
                                        - Dx[i] * Dz[k] / Dys[j-1]
                                        - Dx[i] * Dy[j] / Dzs[k] 
                                        - Dx[i] * Dy[j] / Dzs[k-1] );
}


// FIXME::
__kernel
void Sor_pipeline_getCoef_m(
    __global double *SorCoef,
    __global const double *T1_u,  
    __global const double *T1_v,
    __global const double *T1_w,
    __global const double *Dx,
    __global const double *Dy,
    __global const double *Dz,
    __global const double *Dxs,
    __global const double *Dys,
    __global const double *Dzs,
    __global const int *NEIBcell,
    const double reverse_dt,
    const uint nx, const uint ny , const uint nz, const uint gCells
){

    uint i = get_global_id(0);
    uint j = get_global_id(1);
    uint k = get_global_id(2);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( j < gCells || j >= ny-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;

    const int icelOffSet =  (i-2)*(nz-4)*(ny-4) + (j-2)*(nz-4) + (k-2);
    const int  icel = i*nz*ny + j*nz + k;
    double mChange =  ( T1_u[icel]- T1_u[NEIBcell[icel*12 +0]] ) * Dy[j] * Dz[k] 
                    + ( T1_v[icel]- T1_v[NEIBcell[icel*12 +2]] ) * Dx[i] * Dz[k] 
                    + ( T1_w[icel]- T1_w[NEIBcell[icel*12 +4]] ) * Dx[i] * Dy[j] ;
    SorCoef[icelOffSet*8+6] = mChange * reverse_dt;

}



__kernel
void myAtomicAddG(__global float *addr, float val) {
	union {
		uint u32;
		float f32;
	} current, expected, next;

	do {
		current.f32 = *addr;
		next.f32 = current.f32 + val;
		expected.u32 = current.u32;
		current.u32 = atomic_cmpxchg( (__global uint*) addr, expected.u32, next.u32 );
	} while( current.u32 != expected.u32 );
}


__kernel
void myAtomMax(__global double *max,double val){
	union {
		ulong u64;
		double f64;
	} val_;
    val_.f64 = val;
    atom_max( (__global ulong*) max, val_.u64);
}


#define MAX(a,b) ((a) > (b) ? (a) : (b))


__kernel 
void Sor_pipeline_MAIN(
    __global double *presure,
    __global const double *SorCoef,
    __global double *LocalMax,
    __global uint *iter,
    __global const int *NEIBcell,
    const uint itmax,const  double zeta,const  double omega,
    const uint nx, const  uint ny , const uint nz, const uint gCells
){
    uint i = get_global_id(0);
    uint j = get_global_id(1);
    uint k = get_global_id(2); 



    do {
        LocalMax[0] = 0.0;
        barrier(CLK_GLOBAL_MEM_FENCE);

        if (i >= gCells && i < nx-gCells){
        if (j >= gCells && j < ny-gCells){
        if (k >= gCells && k < nz-gCells){
            LocalMax[0] = 0;
            const int icel = i*nz*ny + j*nz + k;
            const int icelOffSet = (i-2)*(nz-4)*(ny-4) + (j-2)*(nz-4) + (k-2);
            const int icel_12 = icel*12;
            const int icel_08 = icelOffSet*8;

            double  pNEW = (- presure[NEIBcell[icel_12+1]] * SorCoef[icel_08  ]
                            - presure[NEIBcell[icel_12  ]] * SorCoef[icel_08+1]
                            - presure[NEIBcell[icel_12+3]] * SorCoef[icel_08+2]
                            - presure[NEIBcell[icel_12+2]] * SorCoef[icel_08+3] 
                            - presure[NEIBcell[icel_12+5]] * SorCoef[icel_08+4]
                            - presure[NEIBcell[icel_12+4]] * SorCoef[icel_08+5] 
                            + SorCoef[icel_08+6]) * SorCoef[icel_08+7];
    
            double pChange = fabs(pNEW - presure[icel]);
            
            presure[icel] += omega * (pNEW - presure[icel]);

            // if ( pChange > zeta){
                // attomic (reduction)
                if( pChange > LocalMax[0]){
                    myAtomMax(LocalMax, pChange);
                    // LocalMax[0] = pChange;
                }
            // }

        }
        }
        }
        // if (i == 0 && j == 0 && k == 0){
        if (i == 3 && j == 3 && k == 3 ){
            iter[0] += 1;
        }
        // printf("%0.1E ,",*LocalMax);
    }while( LocalMax[0] > zeta && iter[0] < itmax);
    // }while (iter[0] < 30);

    if (i == 3 && j == 3 && k == 3 ){
        iter[0] -= 1;
    }
    barrier(CLK_GLOBAL_MEM_FENCE);
    if (i == 3 && j == 3 && k == 3 ){
        printf("iter:%u \n",iter[0]);
    }

}



__kernel 
void Sor_pipeline_MAIN_unpeeling(
    __global double *presure,
    __global double *LocalMax,
    __global uint *iter,
    __global const int *NEIBcell,
    __global const double *T1_u,  
    __global const double *T1_v,
    __global const double *T1_w,
    __global const double *Dx,
    __global const double *Dy,
    __global const double *Dz,
    __global const double *Dxs,
    __global const double *Dys,
    __global const double *Dzs,
    const uint itmax,const  double zeta,const  double omega,
    const uint nx, const  uint ny , const uint nz, const uint gCells
)
{
    uint i = get_global_id(0);
    uint j = get_global_id(1);
    uint k = get_global_id(2); 

    double dt = 1.0e-3;

    do{

        if (i == 3 && j == 3 && k==3)
            ++iter[0];

        barrier(CLK_GLOBAL_MEM_FENCE);
        LocalMax[0] = 0.0;
        barrier(CLK_GLOBAL_MEM_FENCE);
        if (i >= gCells && i < nx-gCells){
        if (j >= gCells && j < ny-gCells){
        if (k >= gCells && k < nz-gCells){
            int icel = i*nz*ny + j*nz + k;
            int xm  = NEIBcell[0 + 12*icel];
            int xp  = NEIBcell[1 + 12*icel];
            int ym  = NEIBcell[2 + 12*icel];
            int yp  = NEIBcell[3 + 12*icel];
            int zm  = NEIBcell[4 + 12*icel];
            int zp  = NEIBcell[5 + 12*icel];

            double mChange = ( T1_u[icel] - T1_u[xm] ) * Dy[j] * Dz[k] \
                    + ( T1_v[icel] - T1_v[ym] ) * Dx[i] * Dz[k] \
                    + ( T1_w[icel] - T1_w[zm] ) * Dx[i] * Dy[j];

            double pNEW = (- presure[xp] * Dy[j] * Dz[k] / Dxs[i] 
                    - presure[xm] * Dy[j] * Dz[k] / Dxs[i-1] 
                    - presure[yp] * Dx[i] * Dz[k] / Dys[j] 
                    - presure[ym] * Dx[i] * Dz[k] / Dys[j-1] 
                    - presure[zp] * Dx[i] * Dy[j] / Dzs[k] 
                    - presure[zm] * Dx[i] * Dy[j] / Dzs[k-1] 
                    + mChange / dt) /
                       (- Dy[j] * Dz[k] / Dxs[i] - Dy[j] * Dz[k] / Dxs[i-1] 
                        - Dx[i] * Dz[k] / Dys[j] - Dx[i] * Dz[k] / Dys[j-1] 
                        - Dx[i] * Dy[j] / Dzs[k] - Dx[i] * Dy[j] / Dzs[k-1] );

            double pChange = fabs(pNEW - presure[icel]);

            presure[icel] += omega * (pNEW - presure[icel]);
            
            if(pChange > LocalMax[0] )
                LocalMax[0]  = pChange;




        }
        }
        }
    }while(LocalMax[0] > zeta && iter[0] <itmax);

    if (i == 3 && j == 3 && k==3)
        --iter[0];


    if (i == 3 && j == 3 && k == 3 ){
        printf("iter:%u \n",iter[0]);
    }

}






            // printf("pc:%E , ",pChange);
            // if (icel < nx * ny * nz){
            //     printf("[%u,%u,%u]: %e", i, j , k, presure[icel]);
            // }

            // double pChange = max(pNEW - presure[icel], presure[icel] - pNEW); //abs
            // printf("SorCoef[icel_08+7]:%E \n",SorCoef[icel_08+7]);
            // printf("SorCoef[icel_08+6]:%E \n",SorCoef[icel_08+6]);
            // printf("SorCoef[icel_08+5]:%E \n",SorCoef[icel_08+5]);
            // printf("SorCoef[icel_08+4]:%E \n",SorCoef[icel_08+4]);
            // printf("SorCoef[icel_08+2]:%E \n",SorCoef[icel_08+2]);
            // printf("SorCoef[icel_08+1]:%E \n",SorCoef[icel_08+1]);
            // printf("SorCoef[icel_08+0]:%E \n",SorCoef[icel_08+0]);
            // printf("presure[NEIBcell[icel_12+1]]:%E \n", presure[NEIBcell[icel_12+1]] );


            // double PNEW__ = (- presure[NEIBcell[icel_12+1]] * SorCoef[icel_08  ]
            //                 - presure[NEIBcell[icel_12  ]] * SorCoef[icel_08+1]
            //                 - presure[NEIBcell[icel_12+3]] * SorCoef[icel_08+2]
            //                 - presure[NEIBcell[icel_12+2]] * SorCoef[icel_08+3] 
            //                 );
            // printf("pNEW:%E \n",PNEW__);