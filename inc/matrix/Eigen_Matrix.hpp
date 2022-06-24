c#pragma once


#include "controlPanel.hpp"



void createPressureMatrix_Eigen (
    MxClass& Bicg,
    simuClass& simu,
    calDomain& Lo,
    grid& gridA
)
{
    Bicg.matA.resize(gridA.iceltotCal,7);
    const size_t    [nx, ny, nz, gC] = gridA.nxyzgC;

    for (int i = Lo.i_begin ; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin ; j < Lo.j_endof ; ++j )
    for (int k = Lo.k_begin ; k < Lo.k_endof ; ++k )
    {   
        double aa, bb, cc;
        const int icelOffSet =  (i-2)*(nz-4)*(ny-4) + (j-2)*(nz-4) + (k-2);
        if(i==gC)
            aa = -1.0/gridA.Dx.at(i)/gridA.Dxs.at(i-1);
         if(i==(nx-gC-1))
            aa = -1.0/gridA.Dx.at(i)/gridA.Dxs.at(i);
        else 
            aa = 0;
        // ------------------
        if(j==gC)     
            bb = -1.0L/gridA.Dy[j]/gridA.Dys[j-1];
        else if(j==(ny-gC-1))
            bb = -1.0L/gridA.Dy[j]/gridA.Dys[j];
        else
            bb = 0;
        // ------------------
        if(k==gC)
            cc = -1.0L/gridA.Dz[k]/gridA.Dzs[k-1];
        else if(k==(nz-gC-1))
            cc = -1.0L/gridA.Dz[k]/gridA.Dzs[k];
        else
            cc = 0.0L;

        const double coef = 
                    1./gridA.Dx[i]/gridA.Dxs[i] + 
                    1./gridA.Dx[i]/gridA.Dxs[i-1] + 
                    1./gridA.Dy[j]/gridA.Dys[j]+
                    1./gridA.Dy[j]/gridA.Dys[j-1] +
                    1./gridA.Dz[k]/gridA.Dzs[k] + 
                    1./gridA.Dz[k]/gridA.Dzs[k-1]
                    + aa + bb + cc;
        Bicg.matA_Eigen.insert(icelOffSet,icelOffSet) = coef;
    }

    // * ------------------- shift
    const int shift1 = 1;                     //x
    const int shift2 = nz-2*gC;               //y
    const int shift3 = (nz-2*gC)*(ny-2*gC);   //z

    for (int i = Lo.i_begin ; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin ; j < Lo.j_endof ; ++j )
    for (int k = Lo.k_begin ; k < Lo.k_endof ; ++k )
    {
        const int icelOffSet =  (i-2)*(nz-4)*(ny-4) + (j-2)*(nz-4) + (k-2);
        if(i!=gC)
            Bicg.matA_Eigen.insert(icelOffSet, icelOffSet-shift3) =
                                    -1/gridA.Dx[i]/gridA.Dxs[i-1] ;

        if(i!=(nx-gC-1))
            Bicg.matA_Eigen.insert(icelOffSet, icelOffSet+shift3) = 
                                    -1/gridA.Dx[i]/gridA.Dxs[i];

        if(j!=gC)
            Bicg.matA_Eigen.insert(icelOffSet, icelOffSet-shift2 ) =
                                    -1/gridA.Dy[j]/gridA.Dys[j-1] ;
        
        if(j!=(ny-gC-1))
            Bicg.matA_Eigen.insert(icelOffSet, icelOffSet+shift2 ) =
                                    -1/gridA.Dy[j]/gridA.Dys[j];

        if(k!=gC)
            Bicg.matA_Eigen.insert(icelOffSet, icelOffSet-shift1 ) =  
                                    -1/gridA.Dz[k]/gridA.Dzs[k-1];
        
        if(k!=(nz-gC-1))
            Bicg.matA_Eigen.insert(icelOffSet, icelOffSet+shift1 ) =
                                    -1/gridA.Dz[k]/gridA.Dzs[k];
    }

}


void createBMatrix_Eigen (
    velocity& T1,
    MxClass& Bicg,
    simuClass& simu,
    calDomain& Lo,
    grid& gA
    )
{
    double ddt = 1.0 / simu.dt;

    const size_t nx = gA.nx;
    const size_t ny = gA.ny;
    const size_t nz = gA.nz;

    Bicg.matB_Eigen.resize(gA.iceltotCal);

    for (int i = Lo.i_begin ; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin ; j < Lo.j_endof ; ++j )
    for (int k = Lo.k_begin ; k < Lo.k_endof ; ++k )
    {
        const int icel = i*nz*ny + j*nz + k;
        const int icelOffSet =  (i-2)*(nz-4)*(ny-4) + (j-2)*(nz-4) + (k-2);
        const int xn  = gA.NEIBcell[icel*12];
        const int yn  = gA.NEIBcell[icel*12+2];
        const int zn  = gA.NEIBcell[icel*12+4];

        Bicg.matB_Eigen[icelOffSet] = = -ddt*( 
                        ( T1.u[icel] - T1.u[gA.icel(i-1,j,k)] )/gA.Dx[i]+
                        ( T1.v[icel] - T1.v[gA.icel(i,j-1,k)] )/gA.Dy[j]+
                        ( T1.w[icel] - T1.w[gA.icel(i,j,k-1)] )/gA.Dz[k] );
    }

}

#endif