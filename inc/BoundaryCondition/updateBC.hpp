#pragma once

#include "controlPanel.hpp"
#include "Boundary_Condition.hpp"

bool BC_staggered_copy(
    const calDomain& Lo,
    const velocity& oldVel,
    velocity& currtVel,
    grid& gA
)
{

    #pragma omp parallel firstprivate(Lo)
    {
            const auto [nx, ny, nz, gC] = gA.nxyzgC;

        static const std::vector<int> 
                        b0_g = {1,0}, b0_c = {2,3,4},

                        b2_g = {1,0}, b2_c = {2,3,4},

                        b4_g = {1,0}, b4_c = {2,3,4},
                            
                        b1_g = {nx-2,nx-1}, b1_c = {nx-3,nx-4,nx-5},
                        
                        b3_g = {ny-2,ny-1}, b3_c = {ny-3,ny-4,ny-5},
                        
                        b5_g = {nz-2,nz-1}, b5_c = {nz-3,nz-4,nz-5};
        int d = 0;

        #pragma omp for simd collapse(2) nowait
        for(int j = Lo.j_begin-d ; j < Lo.j_endof+d; ++j)
        for(int k = Lo.k_begin-d ; k < Lo.k_endof+d; ++k)
        {// ====================== BC_0 and BC_1 ======================
            currtVel.u[gA.icel(b0_g[0],j,k)] = oldVel.u[gA.icel(b0_g[0],j,k)];
            currtVel.u[gA.icel(b0_g[1],j,k)] = oldVel.u[gA.icel(b0_g[1],j,k)];

            currtVel.v[gA.icel(b0_g[0],j,k)] = oldVel.v[gA.icel(b0_g[0],j,k)];
            currtVel.v[gA.icel(b0_g[1],j,k)] = oldVel.v[gA.icel(b0_g[1],j,k)];

            currtVel.w[gA.icel(b0_g[0],j,k)] = oldVel.w[gA.icel(b0_g[0],j,k)];
            currtVel.w[gA.icel(b0_g[1],j,k)] = oldVel.w[gA.icel(b0_g[1],j,k)];

            currtVel.u[gA.icel(b1_g[1],j,k)] = oldVel.u[gA.icel(b1_g[1],j,k)]; //omission?
            currtVel.u[gA.icel(b1_c[0],j,k)] = oldVel.u[gA.icel(b1_c[0],j,k)];
            currtVel.u[gA.icel(b1_g[0],j,k)] = oldVel.u[gA.icel(b1_g[0],j,k)];

            currtVel.v[gA.icel(b1_g[0],j,k)] = oldVel.v[gA.icel(b1_g[0],j,k)];
            currtVel.v[gA.icel(b1_g[1],j,k)] = oldVel.v[gA.icel(b1_g[1],j,k)];

            currtVel.w[gA.icel(b1_g[0],j,k)] = oldVel.w[gA.icel(b1_g[0],j,k)];
            currtVel.w[gA.icel(b1_g[1],j,k)] = oldVel.w[gA.icel(b1_g[1],j,k)];
        }

        #pragma omp for simd collapse(2) nowait
        for(int i = Lo.i_begin-d ; i < Lo.i_endof+d; ++i )
        for(int k = Lo.k_begin-d ; k < Lo.k_endof+d; ++k )
        {// ====================== BC_2 and BC_3 ======================
            currtVel.u[gA.icel(i,b2_g[0],k)] = oldVel.u[gA.icel(i,b2_g[0],k)];
            currtVel.u[gA.icel(i,b2_g[1],k)] = oldVel.u[gA.icel(i,b2_g[1],k)];

            currtVel.v[gA.icel(i,b2_g[0],k)] = oldVel.v[gA.icel(i,b2_g[0],k)];
            currtVel.v[gA.icel(i,b2_g[1],k)] = oldVel.v[gA.icel(i,b2_g[1],k)];

            currtVel.w[gA.icel(i,b2_g[0],k)] = oldVel.w[gA.icel(i,b2_g[0],k)];
            currtVel.w[gA.icel(i,b2_g[1],k)] = oldVel.w[gA.icel(i,b2_g[1],k)];


            currtVel.u[gA.icel(i,b3_g[0],k)] = oldVel.u[gA.icel(i,b3_g[0],k)];
            currtVel.u[gA.icel(i,b3_g[1],k)] = oldVel.u[gA.icel(i,b3_g[1],k)];

            currtVel.v[gA.icel(i,b3_c[0],k)] = oldVel.v[gA.icel(i,b3_c[0],k)];
            currtVel.v[gA.icel(i,b3_g[0],k)] = oldVel.v[gA.icel(i,b3_g[0],k)];
            currtVel.v[gA.icel(i,b3_g[1],k)] = oldVel.v[gA.icel(i,b3_g[1],k)]; //omission?

            currtVel.w[gA.icel(i,b3_g[0],k)] = oldVel.w[gA.icel(i,b3_g[0],k)];
            currtVel.w[gA.icel(i,b3_g[1],k)] = oldVel.w[gA.icel(i,b3_g[1],k)];
        }

        #pragma omp for simd collapse(2) nowait
        for(int i = Lo.i_begin-d ; i < Lo.i_endof+d; ++i )
        for(int j = Lo.j_begin-d ; j < Lo.j_endof+d; ++j )
        {// ====================== BC_4 and BC_5 ======================
            currtVel.u[gA.icel(i,j,b4_g[0])] = oldVel.u[gA.icel(i,j,b4_g[0])];
            currtVel.u[gA.icel(i,j,b4_g[1])] = oldVel.u[gA.icel(i,j,b4_g[1])];

            currtVel.v[gA.icel(i,j,b4_g[0])] = oldVel.v[gA.icel(i,j,b4_g[0])];
            currtVel.v[gA.icel(i,j,b4_g[1])] = oldVel.v[gA.icel(i,j,b4_g[1])];

            currtVel.w[gA.icel(i,j,b4_g[0])] = oldVel.w[gA.icel(i,j,b4_g[0])];
            currtVel.w[gA.icel(i,j,b4_g[1])] = oldVel.w[gA.icel(i,j,b4_g[1])];
            

            currtVel.u[gA.icel(i,j,b5_g[0])] = oldVel.u[gA.icel(i,j,b5_g[0])];
            currtVel.u[gA.icel(i,j,b5_g[1])] = oldVel.u[gA.icel(i,j,b5_g[1])];
            
            currtVel.v[gA.icel(i,j,b5_g[0])] = oldVel.v[gA.icel(i,j,b5_g[0])];
            currtVel.v[gA.icel(i,j,b5_g[1])] = oldVel.v[gA.icel(i,j,b5_g[1])];

            currtVel.w[gA.icel(i,j,b5_c[0])] = oldVel.w[gA.icel(i,j,b5_c[0])];
            currtVel.w[gA.icel(i,j,b5_g[0])] = oldVel.w[gA.icel(i,j,b5_g[0])];
            currtVel.w[gA.icel(i,j,b5_g[1])] = oldVel.w[gA.icel(i,j,b5_g[1])]; //omission?
        }

    }

    return true;
}



bool BC_updateSlid(
    const calDomain& Lo,
    velocity& vel,
    grid& gA
)
{
    auto ii = [&](auto &i, auto &j, auto &k){return gA.icel(i,j,k);};

    auto [Num, Dir] = getNumDIr();

    #pragma omp parallel firstprivate(Lo)
    {
        const auto [nx, ny, nz, gC] = gA.nxyzgC;

        std::vector<int> 
                        b0_g = {1,0}, b0_c = {2,3,4},

                        b2_g = {1,0}, b2_c = {2,3,4},

                        b4_g = {1,0}, b4_c = {2,3,4};
        int d = 0;

        #pragma omp for simd collapse(2) nowait
        for(int j = Lo.j_begin-d ; j < Lo.j_endof+d; ++j)
        for(int k = Lo.k_begin-d ; k < Lo.k_endof+d; ++k)
        {// ====================== BC_0 ======================
            const int b = 0;
            vel.u[gA.icel(b0_g[1],j,k)] = vel.u[gA.icel(b0_g[0],j,k)] =  Dir[b][0]    +      Num[b][0]      * vel.u[gA.icel(b0_c[0],j,k)];
            vel.v[gA.icel(b0_g[1],j,k)] = vel.v[gA.icel(b0_g[0],j,k)] =  Dir[b][1]*2. + (2.0*Num[b][1]-1.0) * vel.v[gA.icel(b0_c[0],j,k)];
            vel.w[gA.icel(b0_g[1],j,k)] = vel.w[gA.icel(b0_g[0],j,k)] =  Dir[b][2]*2. + (2.0*Num[b][2]-1.0) * vel.w[gA.icel(b0_c[0],j,k)];
        }


        #pragma omp for simd collapse(2) nowait
        for(int i = Lo.i_begin-d ; i < Lo.i_endof+d; ++i )
        for(int k = Lo.k_begin-d ; k < Lo.k_endof+d; ++k )
        {// ====================== BC_2 ======================
            const int b = 2;
            vel.u[gA.icel(i,b2_g[1],k)] = vel.u[gA.icel(i,b2_g[0],k)] =  Dir[b][0]*2. + (2.0*Num[b][0]-1.0) * vel.u[gA.icel(i,b2_c[0],k)];
            vel.v[gA.icel(i,b2_g[1],k)] = vel.v[gA.icel(i,b2_g[0],k)] =  Dir[b][1]    +      Num[b][1]      * vel.v[gA.icel(i,b2_c[0],k)];
            vel.w[gA.icel(i,b2_g[1],k)] = vel.w[gA.icel(i,b2_g[0],k)] =  Dir[b][2]*2. + (2.0*Num[b][2]-1.0) * vel.w[gA.icel(i,b2_c[0],k)];
        }


        #pragma omp for simd collapse(2) nowait
        for(int i = Lo.i_begin-d ; i < Lo.i_endof+d; ++i )
        for(int j = Lo.j_begin-d ; j < Lo.j_endof+d; ++j )
        {// ====================== BC_4 ======================
            const int b = 4;
            vel.u[gA.icel(i,j,b4_g[1])] = vel.u[gA.icel(i,j,b4_g[0])] =  Dir[b][0]*2. + (2.0*Num[b][0]-1.0) * vel.u[gA.icel(i,j,b4_c[0])];
            vel.v[gA.icel(i,j,b4_g[1])] = vel.v[gA.icel(i,j,b4_g[0])] =  Dir[b][1]*2. + (2.0*Num[b][1]-1.0) * vel.v[gA.icel(i,j,b4_c[0])];
            vel.w[gA.icel(i,j,b4_g[1])] = vel.w[gA.icel(i,j,b4_g[0])] =  Dir[b][2]    +      Num[b][2]      * vel.w[gA.icel(i,j,b4_c[0])];
        }

    }

    return true;
}