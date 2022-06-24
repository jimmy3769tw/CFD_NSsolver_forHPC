#pragma once
#include"controlPanel.hpp"
#include <cmath>

void DfibInit(
    DfibArray& Dfib,
    calDomain& Lo,
    grid& gridA
)
{

    auto [nx, ny, nz, gC] = gridA.nxyzgC;

    // !first touch

    for (int i = Lo.i_begin ; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin ; j < Lo.j_endof ; ++j )
    for (int k = Lo.k_begin ; k < Lo.k_endof ; ++k )
    {

        Dfib.eta[gridA.icel(i,j,k)] = 0.0;

        if (i == 2) {
            Dfib.eta[ gridA.icel(i-1,j,k) ] = 0.0;
            Dfib.eta[ gridA.icel(i-2,j,k) ] = 0.0;
        }

        if (i == nx -3 ) {
            Dfib.eta[ gridA.icel(i+1,j,k) ] = 0.0;
            Dfib.eta[ gridA.icel(i+2,j,k) ] = 0.0;
        }

        if (j == 2) {
            Dfib.eta[ gridA.icel(i,j-1,k) ]  = 0.0;
            Dfib.eta[ gridA.icel(i,j-2,k) ]  = 0.0;
        }

        if (j == nx -3 ) {
            Dfib.eta[  gridA.icel(i,j+1,k) ]= 0.0;
            Dfib.eta[  gridA.icel(i,j+2,k) ]= 0.0;
        }

        if (k == 2) {
            Dfib.eta[ gridA.icel(i,j,k-1) ] = 0.0;
            Dfib.eta[ gridA.icel(i,j,k-2) ] = 0.0;
        }

        if (k == nz -3 ) {
            Dfib.eta[ gridA.icel(i,j,k+1) ] = 0.0;;
            Dfib.eta[ gridA.icel(i,j,k+2) ] = 0.0;;
        }
    }
}

void Dfib_Sphere(
    DfibArray& Dfib,
    calDomain& Lo,
    grid& gridA
)
{

    
}



// * for X direction
void DFIB_CylinderX(
    DfibArray& Dfib,
    calDomain& Lo,
    grid& gA
)
{

    auto [nx, ny, nz ] = gA.nxyz;

    double  center_z = gA.lz/2.0, 
            center_y = gA.ly/2.0;

    double Radius = gA.lz/2.0; // * setting 

    int nSub = 100;

    vector <double> sdy (nSub), sdz (nSub);
    
    int one = 1;
    int zero = 0;

    for(size_t j = 0 ; j < ny-4; ++j)
    for(size_t k = 0 ; k < nz-4; ++k)
    {
        const int icel =  gA.icel(zero,j,k);
        const double Diagonal = sqrt(//gridA.Dx[0]*gridA.Dx[0] + 
                                        gA.Dy[j]*gA.Dy[j] + 
                                        gA.Dz[k]*gA.Dz[k]  ) * 0.5;

        const double Distance = sqrt( pow(((gA.Z[k] + gA.Z[k+1]) *0.5) - center_z, 2) 
                                    + pow(((gA.Y[j] + gA.Y[j+1]) *0.5) - center_y, 2) ); 

        if (abs (Distance - Radius) < Diagonal)
        {
            size_t xi = 0;
            const double dyg = gA.Dy[j] / double (nSub) ;
            const double dzg = gA.Dz[k] / double (nSub) ;

            for (size_t jj; jj < nSub ; jj++)
            for (size_t kk; kk < nSub ; kk++)
            { 
                sdy[jj] = gA.Z[j] + (jj) * dyg;
                sdz[jj] = gA.Z[k] + (kk) * dzg;
            }

            for (size_t jj; jj < nSub-1 ; jj++)
            for (size_t kk; kk < nSub-1 ; kk++)
            {
                auto Sdist = sqrt( sdy[jj] + sdy[jj+1]  );
                if (Sdist <= Radius){ ++xi;}
            }
            Dfib.eta[icel] = double(xi) * pow(nSub,-2);

        }
        else if(Distance > Radius){ Dfib.eta[icel] = 0.; }
        else{ Dfib.eta[icel] = 1.; }
    }  

    for(size_t i = 1 ; i < nx-1; ++i)
    for(size_t j = 0 ; j < ny-1; ++j)
    for(size_t k = 0 ; k < nz-1; ++k)
    {
        const int icel = gA.icel(i,j,k);
        Dfib.eta[icel]  = gA.icel(one,j,k);
    }
}


void DFIB_CylinderZ(
    DfibArray& Dfib,
    calDomain& Lo,
    grid& gA
)
{
    auto [nx, ny , nz, gC] = gA.nxyzgC;

    // * ------------------------- setting
    double  center_x    = 6.5, 
            center_y    = gA.ly/2.0,
            Radius      = 0.5;
    Dfib.cylinderDimension = Radius*2;
    Dfib.cylinderCenter.resize(2);
    Dfib.cylinderCenter[0] = center_x;
    Dfib.cylinderCenter[1] = center_y;
    int     nSubGrids   = 100;
    // * ------------------------- setting

    vector<double> sdy (nSubGrids+1);
    vector<double> sdx (nSubGrids+1);
    double Sdist;

    for(size_t i = gC ; i < nx - gC; ++i)
    for(size_t j = gC ; j < ny - gC; ++j)
    {
        const int icel = i*nz*ny + j*nz;
        const double Diagonal = sqrt(gA.Dx[i]*gA.Dx[i] + gA.Dy[j]*gA.Dy[j]) *0.5; //FIXME :: this is three way??
        const double Distance = sqrt(pow(((gA.X[i] + gA.X[i+1]) *0.5) -center_x,2) + pow(((gA.Y[j] + gA.Y[j+1]) *0.5)-center_y,2) );

        if (std::abs(Distance - Radius) < Diagonal) 
        {
            size_t xi = 0;
            const double dyg = gA.Dy[j] / double (nSubGrids) ;
            const double dxg = gA.Dx[i] / double (nSubGrids) ;
                // * | | | |
            for (size_t ii = 0; ii < nSubGrids+1 ; ii++)
            for (size_t jj = 0; jj < nSubGrids+1 ; jj++)
            { 
                sdy[jj] = gA.Y[j] + (jj) * dyg;
                sdx[ii] = gA.X[i] + (ii) * dxg;
            }

            for (size_t ii = 0; ii < nSubGrids ; ii++)
            for (size_t jj = 0; jj < nSubGrids ; jj++)
            {
                auto Sdist = sqrt(  pow((sdx[ii] + sdx[ii+1]) * 0.5 - center_x, 2) + 
                                    pow((sdy[jj] + sdy[jj+1]) * 0.5 - center_y, 2) );

                if (Sdist <= Radius){
                    ++xi;
                }
            }
            Dfib.eta[icel] = double(xi) * pow(nSubGrids,-2);

        }
        else if(Distance > Radius){
            Dfib.eta[icel] = 0.;
        }
        else{
            Dfib.eta[icel] = 1.;
        }
    }  

    for(size_t i = 0 ; i < nx-1; ++i)
    for(size_t j = 0 ; j < ny-1; ++j)
    for(size_t k = 1 ; k < nz-1; ++k)
        Dfib.eta[gA.icel(i,j,k)]  = Dfib.eta[ gA.icel(i,j,0)];
}

void PotentialFlow(
    simuClass& simu,
    DfibArray& Dfib,
    velocity& vel,
    pressure& pre,
    calDomain& Lo,
    grid& gA
)
{
    auto [nx, ny , nz, gC] = gA.nxyzgC;

    double U_infty = 1;
    double pi = acos(-1);

    auto D_c = Dfib.cylinderDimension;
    auto x_D = Dfib.cylinderCenter[0];
    auto y_D = Dfib.cylinderCenter[1];
    for(auto i = gA.gC ; i < gA.nx - gA.gC; ++i)
    for(auto j = gA.gC ; j < gA.ny - gA.gC; ++j)
    {
        auto tX = (gA.X[i] + gA.X[i+1]) *0.5 - x_D;
        auto tY = (gA.Y[j] + gA.Y[j+1]) *0.5 - y_D;
        auto tempX = pow((tX),2);
        auto tempY = pow((tY),2);
        // vel.u[gA.icel(i,j,gC)] += -D_c * (tempX + tempY) / (2 * pi) * (tempX - tempY) ;
        // vel.v[gA.icel(i,j,gC)] += -D_c * (tempX + tempY) / (2.* pi) * (2*(tX*tY));
        // pre.p[gA.icel(i,j,1)] = -0.5* simu.nu *(
        //     pow(vel.u[gA.icel(i,j,1)], 2.0) + pow(vel.v[gA.icel(i,j,1)], 2.0) );
    }

    for(auto i = gC+0 ; i < nx-2; ++i)
    for(auto j = gC+0 ; j < ny-2; ++j)
    for(auto k = gC+1 ; k < nz-2; ++k)
    {
        auto ick1 = gA.icel(i,j,gC);
        auto ic = gA.icel(i,j,k);
        vel.u[ic] = vel.u[ick1];
        pre.p[ic] = pre.p[ick1];
        vel.v[ic] = vel.v[ick1];
    }
}