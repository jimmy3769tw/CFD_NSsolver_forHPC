#pragma once
#include "controlPanel.hpp"

auto getNumDIr(){

    vector<vector<double> > Num (6,vector<double>(3, 0.0) );
    vector<vector<double> > Dir (6,vector<double>(3, 0.0) );

    // #include "ChannelFlow.hpp"
    // #include "FLOWpassing.hpp"
    #include "BC_Cavity_A.hpp"

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
   return std::make_pair(Num, Dir);
}





bool BC_staggered_main(
    calDomain& Lo,
    velocity& vel,
    pressure& pre,
    grid& gA
)
{
    auto [Num, Dir] = getNumDIr();

#pragma omp parallel firstprivate(Lo) shared(gA)
{ //? ======================= OMP -start  ===================== ?

            const auto [nx, ny, nz, gC] = gA.nxyzgC;
/*

* |----|-----|----|-----|----|--    --|----|---|---|---|------|-----|------|
* 0----1-----2----3-----|----|--    --|----|---|---|---n-4---n-3---n-2----n-1
! g1---g0---c0----c1----|----|--    --|----|---|---|---c1----c0----g0-----g1

*/  
            static const std::vector<int> 
                        b0_g = {1,0}, b0_c = {2,3,4},

                        b2_g = {1,0}, b2_c = {2,3,4},

                        b4_g = {1,0}, b4_c = {2,3,4},
                            
                        b1_g = {nx-2,nx-1}, b1_c = {nx-3,nx-4,nx-5},
                        
                        b3_g = {ny-2,ny-1}, b3_c = {ny-3,ny-4,ny-5},
                        
                        b5_g = {nz-2,nz-1}, b5_c = {nz-3,nz-4,nz-5};

        int d = 0;

        #pragma omp for nowait
        for(auto j = Lo.j_begin-d ; j < Lo.j_endof+d; ++j)
        for(auto k = Lo.k_begin-d ; k < Lo.k_endof+d; ++k)
        {// ====================== BC_0 ======================
            const int b = 0;
            vel.u[gA.icel(b0_g[1],j,k)] = vel.u[gA.icel(b0_g[0],j,k)] =  Dir[b][0]    +      Num[b][0]      * vel.u[gA.icel(b0_c[0],j,k)];
            vel.v[gA.icel(b0_g[1],j,k)] = vel.v[gA.icel(b0_g[0],j,k)] =  Dir[b][1]*2. + (2.0*Num[b][1]-1.0) * vel.v[gA.icel(b0_c[0],j,k)];
            vel.w[gA.icel(b0_g[1],j,k)] = vel.w[gA.icel(b0_g[0],j,k)] =  Dir[b][2]*2. + (2.0*Num[b][2]-1.0) * vel.w[gA.icel(b0_c[0],j,k)];
        }


        #pragma omp for nowait
        for(auto j = Lo.j_begin-d ; j < Lo.j_endof+d; ++j)
        for(auto k = Lo.k_begin-d ; k < Lo.k_endof+d; ++k)
        {// ====================== BC_1 ======================
            const int b = 1;
            vel.u[gA.icel(b1_g[1],j,k)] = 
            vel.u[gA.icel(b1_g[0],j,k)] = vel.u[gA.icel(b1_c[0],j,k)] = Dir[b][0]     +      Num[b][0]      * vel.u[gA.icel(b1_c[1],j,k)];

            vel.v[gA.icel(b1_g[1],j,k)] = vel.v[gA.icel(b1_g[0],j,k)] = Dir[b][1]*2.  + (2.0*Num[b][1]-1.0) * vel.v[gA.icel(b1_c[0],j,k)];
            vel.w[gA.icel(b1_g[1],j,k)] = vel.w[gA.icel(b1_g[0],j,k)] = Dir[b][2]*2.  + (2.0*Num[b][2]-1.0) * vel.w[gA.icel(b1_c[0],j,k)];
        }


        #pragma omp for nowait
        for(auto i = Lo.i_begin-d ; i < Lo.i_endof+d; ++i )
        for(auto k = Lo.k_begin-d ; k < Lo.k_endof+d; ++k )
        {// ====================== BC_2 ======================
            const int b = 2;
            vel.u[gA.icel(i,b2_g[1],k)] = vel.u[gA.icel(i,b2_g[0],k)] =  Dir[b][0]*2. + (2.0*Num[b][0]-1.0) * vel.u[gA.icel(i,b2_c[0],k)];
            vel.v[gA.icel(i,b2_g[1],k)] = vel.v[gA.icel(i,b2_g[0],k)] =  Dir[b][1]    +      Num[b][1]      * vel.v[gA.icel(i,b2_c[0],k)];
            vel.w[gA.icel(i,b2_g[1],k)] = vel.w[gA.icel(i,b2_g[0],k)] =  Dir[b][2]*2. + (2.0*Num[b][2]-1.0) * vel.w[gA.icel(i,b2_c[0],k)];
        }

        #pragma omp for nowait
        for(auto i = Lo.i_begin-d ; i < Lo.i_endof+d; ++i )
        for(auto k = Lo.k_begin-d ; k < Lo.k_endof+d; ++k )
        {// ====================== BC_3 ======================
            const int b = 3;
            vel.v[gA.icel(i,b3_g[1],k)] = 
            vel.v[gA.icel(i,b3_g[0],k)] = vel.v[gA.icel(i,b3_c[0],k)] =  Dir[b][1]    +      Num[b][1]      * vel.v[gA.icel(i,b3_c[1],k)];

            vel.u[gA.icel(i,b3_g[1],k)] = vel.u[gA.icel(i,b3_g[0],k)] =  Dir[b][0]*2. + (2.0*Num[b][0]-1.0) * vel.u[gA.icel(i,b3_c[0],k)];
            vel.w[gA.icel(i,b3_g[1],k)] = vel.w[gA.icel(i,b3_g[0],k)] =  Dir[b][2]*2. + (2.0*Num[b][2]-1.0) * vel.w[gA.icel(i,b3_c[0],k)];
        }



        #pragma omp for nowait
        for(auto i = Lo.i_begin-d ; i < Lo.i_endof+d; ++i )
        for(auto j = Lo.j_begin-d ; j < Lo.j_endof+d; ++j )
        {// ====================== BC_4 ======================
            const int b = 4;
            vel.u[gA.icel(i,j,b4_g[1])] = vel.u[gA.icel(i,j,b4_g[0])] =  Dir[b][0]*2. + (2.0*Num[b][0]-1.0) * vel.u[gA.icel(i,j,b4_c[0])];
            vel.v[gA.icel(i,j,b4_g[1])] = vel.v[gA.icel(i,j,b4_g[0])] =  Dir[b][1]*2. + (2.0*Num[b][1]-1.0) * vel.v[gA.icel(i,j,b4_c[0])];
            vel.w[gA.icel(i,j,b4_g[1])] = vel.w[gA.icel(i,j,b4_g[0])] =  Dir[b][2]    +      Num[b][2]      * vel.w[gA.icel(i,j,b4_c[0])];
        }



        #pragma omp for nowait
        for(auto i = Lo.i_begin-d ; i < Lo.i_endof+d; ++i )
        for(auto j = Lo.j_begin-d ; j < Lo.j_endof+d; ++j )
        {// ====================== BC_5 ======================
            const int b = 5;
            vel.w[gA.icel(i,j,b5_g[1])] = 
            vel.w[gA.icel(i,j,b5_g[0])] = vel.w[gA.icel(i,j,b5_c[0])] =  Dir[b][2]    +      Num[b][2]      * vel.w[gA.icel(i,j,b5_c[1])];
        
            vel.u[gA.icel(i,j,b5_g[1])] = vel.u[gA.icel(i,j,b5_g[0])] =  Dir[b][0]*2. + (2.0*Num[b][0]-1.0) * vel.u[gA.icel(i,j,b5_c[0])];
            vel.v[gA.icel(i,j,b5_g[1])] = vel.v[gA.icel(i,j,b5_g[0])] =  Dir[b][1]*2. + (2.0*Num[b][1]-1.0) * vel.v[gA.icel(i,j,b5_c[0])];
        }


        // * =================== pressure =================

        d = 0;
        #pragma omp for nowait
        for(auto j = Lo.j_begin-d ; j < Lo.j_endof+d; ++j )
        for(auto k = Lo.k_begin-d ; k < Lo.k_endof+d; ++k )
        {
            pre.p[gA.icel(b0_g[1],j,k)] = pre.p[gA.icel(b0_g[0],j,k)] = pre.p[gA.icel(b0_c[0],j,k)];
            pre.p[gA.icel(b1_g[1],j,k)] = pre.p[gA.icel(b1_g[0],j,k)] = pre.p[gA.icel(b1_c[0],j,k)];
        }


        #pragma omp for nowait
        for(auto i = Lo.i_begin-d ; i < Lo.i_endof+d; ++i )
        for(auto k = Lo.k_begin-d ; k < Lo.k_endof+d; ++k )
        {
            pre.p[gA.icel(i,b2_g[1],k)] = pre.p[gA.icel(i,b2_g[0],k)] = pre.p[gA.icel(i,b2_c[0],k)];
            pre.p[gA.icel(i,b3_g[1],k)] = pre.p[gA.icel(i,b3_g[0],k)] = pre.p[gA.icel(i,b3_c[0],k)];
        }

        #pragma omp for 
        for(auto i = Lo.i_begin-d ; i < Lo.i_endof+d; ++i )
        for(auto j = Lo.j_begin-d ; j < Lo.j_endof+d; ++j )
        {
            pre.p[gA.icel(i,j,b4_g[1])] = pre.p[gA.icel(i,j,b4_g[0])] = pre.p[gA.icel(i,j,b4_c[0])];
            pre.p[gA.icel(i,j,b5_g[1])] = pre.p[gA.icel(i,j,b5_g[0])] = pre.p[gA.icel(i,j,b5_c[0])];
        }
    
    } // Omp fork join 
    return true;
}


