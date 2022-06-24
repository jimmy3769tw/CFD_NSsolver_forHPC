#pragma once
#include <vector>
#include <string>
#include <tuple>
#include <fstream>
#include <iostream>
#include <sys/types.h>
#include <dirent.h>
#include <cmath>
#include "controlPanel.hpp"
using namespace std;

// namespace::structureGrid
class grid{

        std::vector<double>  delta_;

        public:


        grid(){resize();}

        vector <string> DirectionName{"[x]", "[Y]", "[Z]"};

        // const string Gridder = "uniform" ;
        const std::string Gridder = "non_uniform-SML";

        // const double  lx = 5.0 ,      ly = 1.0 ,      lz = 1.0 ;
        // const double  lx = 1.0 ,      ly = 1.0 ,      lz = 1.0 ;
        // const double  lx = 10.0 ,      ly = 1.0 ,      lz = 1.0 ;

        // *   Accuracy  ------------
        // yousa 2020.03.18 scheme
        const double  lx = 17.0 ,      ly = 12.0 ,      lz = 3.0 ;

        // const double lx = 10, ly = 1, lz = 1;
        // const double lx = 1., ly = 1., lz = 1.;

        // *int
        const int       gC = 2;


        // yousa 2020.03.18 scheme
        const int       nxCal = 190,     // nx should > 40
                        nyCal = 210,
                        nzCal = 40;


        // const int       nxCal = 100,     // nx should > 40
        //                 nyCal = 20,
        //                 nzCal = 20;


        // const int       nxCal = 100,     // nx should > 40
        //                 nyCal = 30,
        //                 nzCal = 10;


        // const int       nxCal = 128,     // nx should > 40
        //                 nyCal = 60,
        //                 nzCal = 40;


        // Flow passing cylinder
        // const int       nxCal = 128,     // nx should > 40
        //                 nyCal = 61,
        //                 nzCal = 30;

        // const int       nxCal = 3,     // nx should > 40
        //                 nyCal = 3,
        //                 nzCal = 3;



        // const int       nxCal = 40*5,     // nx should > 40
        //                 nyCal = 40,
        //                 nzCal = 40;
        // const int       nxCal = 40,     // nx should > 40
        //                 nyCal = 40*4,
        //                 nzCal = 40*4;

        // const int       nxCal = 40,     // nx should > 40
        //                 nyCal = 40,
        //                 nzCal = 40;

        // ! -------------------------------------------- !
        // ------------------------
        const int iceltotCal =  (nxCal) * (nyCal) * (nzCal) ;
        // ------------------------

        // ------------------------
        const int       nx = nxCal + 2*gC, 
                        ny = nyCal + 2*gC,
                        nz = nzCal + 2*gC;

        const int iceltot = nx * ny * nz;
        // ------------------------

        std::vector<double>  Xc, Yc, Zc, Dx, Dy, Dz, Dxs, Dys, Dzs, X, Y, Z;



        vector<int> NEIBcell;

        vector<double> jp;// jacobi preconditioner

        // * ----------------------------------------------



        bool resize(){
                delta_.resize(iceltotCal);
                // ------------
                NEIBcell.resize(iceltot*12);
                // ------------

                // ------------
                Xc.resize(nxCal);
                Yc.resize(nyCal);
                Zc.resize(nzCal);
                // ------------

                // ------------
                Dx.resize(nx);
                Dy.resize(ny);
                Dz.resize(nz);
                // ------------

                // ------------
                Dxs.resize(nx);
                Dys.resize(ny);
                Dzs.resize(nz);
                // ------------

                // ------------
                X.resize(nx+1);
                Y.resize(ny+1);
                Z.resize(nz+1);
                // ------------
                return true;
        }

        // ---------------------------------
        const int       nzny = nz*ny, 
                        nxny = nx*ny, 
                        nynz = ny*nz;
        // ---------------------------------

        // ---------------------------------
        const std::tuple<int, int, int> nxyz{nx, ny, nz};

        const std::tuple<int, int, int, int> nxyzgC{nx, ny, nz, gC};

        const std::tuple<int, int, int> nxyzCal{nxCal, nyCal, nzCal};
        // ---------------------------------

        // ---------------------------------
        template<typename T1, typename T2, typename T3>
        inline size_t icel(T1 i, T2 j, T3 k){
                return i*nzny + j*nz +k;
        }
        // ---------------------------------


        // ---------------------------------
        int nzDir = (nzCal+2);
        int nznyDir = (nzCal+2)*(nyCal+2);
        int nznynxDir = (nzCal+2)*(nyCal+2)*(nxCal+2);
        // ---------------------------------

        // ---------------------------------
        template<typename T1, typename T2, typename T3>
        inline int icelDir(T1 i, T2 j, T3 k){
                return (i-gC)*nznyDir + (j-gC)*nzDir +(k-gC);
        } 
        // ---------------------------------


        // ---------------------------------
        int nznyCal = nzCal*nyCal;
        template<typename T1, typename T2, typename T3>
        inline size_t icelCal(T1 i, T2 j, T3 k){

                return (i-gC)*nznyCal + (j-gC)*nzCal + (k-gC);
        }
        // ---------------------------------


        // ---------------------------------
        template<typename T1, typename T2>
        inline size_t icel_z(T1 i, T2 j){
                return i*ny + j;
        }
        // ---------------------------------


        // ---------------------------------
        template<typename T1, typename T2, typename T3>
        inline size_t icelOffset(T1 i, T2 j, T3 k){
                return i*nznyCal + j*nzCal + k;
        }
        // ---------------------------------


        template<typename T1, typename T2, typename T3>
        double& delta(T1 i, T2 j, T3 k)
        { return delta_[icelCal(i, j, k)];}


        // ---------------------------------
        inline void showVel(std::vector<double>& which_U){

                cout << "\nicel ";
                for (size_t i = gC; i < nx-gC; ++i)
                {
                        cout << "\n\n" << "slide "<< i ;
                        for (size_t j = gC; j < ny-gC; ++j)
                        {
                                cout << "\n";
                                for (size_t k = gC; k < nz-gC; ++k)
                                {
                                        cout << which_U.at(icel(i,j,k)) << ", ";
                                }
                        }
                }
        }
        // ---------------------------------




        // ---------------------------------
        inline void showVelt(std::vector<double>& which_U){

                cout << "\nicel ";
                for (size_t i = 0; i < nx; ++i)
                {
                        cout << "\n\n" << "slide "<< i ;
                        for (size_t j = 0; j < ny; ++j)
                        {
                                cout << "\n";
                                for (size_t k = 0; k < nz; ++k)
                                {
                                        cout << which_U.at(icel(i,j,k)) << ", ";
                                }
                        }
                }
        }
        // ---------------------------------



        // ---------------------------------
        inline void showCalVector(std::vector<double>& which_U){

                cout << "\nicel ";
                for (size_t i = 0; i < nxCal; ++i)
                {
                        cout << "\n\n" << "slide "<< i ;
                        for (size_t j = 0; j < nyCal; ++j)
                        {
                                cout << "\n";
                                for (size_t k = 0; k < nzCal; ++k)
                                {
                                        cout << which_U.at(icelOffset(i,j,k)) << ", ";
                                }
                        }
                }
        }
        // ---------------------------------





        // ---------------------------------
        template<typename T1, typename T2, typename T3>
        inline double coef(T1 i, T2 j, T3 k){

                return
                  Dy[j] * Dz[k] / Dxs[i] 
                + Dy[j] * Dz[k] / Dxs[i-1]
                + Dx[i] * Dz[k] / Dys[j] 
                + Dx[i] * Dz[k] / Dys[j-1]
                + Dx[i] * Dy[j] / Dzs[k] 
                + Dx[i] * Dy[j] / Dzs[k-1] ;
        }
        // ---------------------------------


        // * init value ---------------------
        bool s_init(){

                for (size_t i = 0; i < nx-2*gC; ++i){
                        Xc.at(i) = ( X.at(i+gC) + X.at(i+gC+1) ) / 2.0;
                }

                for (size_t j = 0; j < ny-2*gC; ++j){
                        Yc.at(j) = ( Y.at(j+gC)+ Y.at(j+gC+1) ) / 2.0;
                }

                for (size_t k = 0; k < nz-2*gC; ++k){
                        Zc.at(k) = ( Z.at(k+gC)+ Z.at(k+gC+1) ) / 2.0;
                }
                return true;
        }


        std::string show(){

                static auto to_s = [&](auto a){ return std::to_string(a);};
                // ---------------------
                std::string sh_s;
                sh_s += ", [nx, ny, nz]=";
                sh_s += to_s(nxCal);
                sh_s += ", " ;
                sh_s += to_s(nyCal);
                sh_s += ", " ;
                sh_s += to_s(nzCal);
                // ---------------------

                return sh_s;
        }


        auto unit_SML(vector<double> &a, 
                const double lSml, const int nS, const int off_nS, 
                const double l ,const int n
        ){
                cout << "\n-------------unit_SML -------------\n";

        // * | goast Cell | Part[0](L) | Part[1].lx | Part[2].lx| Part[3].lx | Part[4].lx|goast Cell|
        // * | 0 1        | Part[0].nx | dySml                                          | nx-2, nx-1|
                auto nCal = n-2*gC;

                const double    dS = lSml / nS , 
                                d  = (l-lSml) / (nCal-nS);

                const double    SBgn = 0.5*(l - lSml),
                                SEnd = 0.5*(l + lSml);

                a.at(gC) = 0.0;

                if (nS != 0)
                {
                        // * Calculation the grid
                        if ( (n-nS)%2 != 0 ) throw std::invalid_argument("(n-nSml)%2 != 0");

                        int nlEnd[2]; 
                        int nsEnd[1];

                        nlEnd[0]  = (nCal-nS) / 2;
                        nsEnd[0]  = nlEnd[0] + nS;
                        nlEnd[1]  = nsEnd[0] + nlEnd[0] ;

                        nlEnd[0] += off_nS;
                        nsEnd[0] += off_nS;

                        cout << "nlEnd[0] " << nlEnd[0] << "\n";
                        cout << "nsEnd[0] " << nsEnd[0] << "\n";
                        cout << "nlEnd[1] " << nlEnd[1] << "\n";

                        if ( nlEnd[0] <= 0) throw std::invalid_argument("nL[0]");

                        
                        for(size_t i = 0; i < nlEnd[0]; i++)
                                a.at(i+gC+1) =  a.at(i+gC) + d;

                        for(size_t i = nlEnd[0]; i < nsEnd[0]; ++i)
                                a.at(i+gC+1) =  a.at(i+gC) + dS;

                        for(size_t i = nsEnd[0]; i < nlEnd[1]; i++)
                                a.at(i+gC+1) =  a.at(i+gC) + d;

                }
                else{

                        for(size_t i = gC; i < n-gC; i++)
                                a.at(i+1) =  a.at(i) + l/nCal;

                }
                
                init_gC_dx();

                return true;
        }

        bool init_nonuniform_SML(){

        // * | goast Cell | Part[0](L) | Part[1].lx | Part[2].lx| Part[3].lx | Part[4].lx|goast Cell|
        // * | 0 1        | Part[0].nx | dySml| nx-2, nx-1|

        // *Setting
                const double    lxS = 0.3 * lx,       lyS = 0.3 * ly,       lzS = 0.5 * lz;

                const int       nxS =  100,       nyS =  100 ,            nzS = 0 ;

                const int       off_nS_x = -4    ,       off_nS_y = 0,           off_nS_z = 0;

                unit_SML (X, lxS, nxS, off_nS_x, lx, nx);
                unit_SML (Y, lyS, nyS, off_nS_y, ly, ny);
                unit_SML (Z, lzS, nzS, off_nS_z, lz, nz);

                return true;
        }


        // * init value ---------------------
        bool init_uniform(){

                double dx = lx /( nx - 2.0*gC); 
                double dy = ly /( ny - 2.0*gC);
                double dz = lz /( nz - 2.0*gC);
                double Xi = 0.0, Yi = 0.0, Zi = 0.0;

                X.at(0) = Xi-2.0*dx;
                for(size_t i = 1; i < nx+1; ++i) {
                        X.at(i) = X.at(i-1) + dx;
                }

                Y[0] = Yi-2.0*dy;
                for(size_t j = 1; j < ny+1; ++j) {
                        Y.at(j) = Y.at(j-1) + dy;
                }

                Z[0] = Zi-2.0*dz;
                for(size_t k = 1; k < nz+1; ++k) {
                        Z.at(k) = Z.at(k-1) + dz;
                }

                fill(Dx.begin(), Dx.end(), dx);
                fill(Dy.begin(), Dy.end(), dy);
                fill(Dz.begin(), Dz.end(), dz);

                fill(Dxs.begin(), Dxs.end(), dx);
                fill(Dys.begin(), Dys.end(), dy);
                fill(Dzs.begin(), Dzs.end(), dz);
                
                return true;
        }

        inline auto x_INTRPL(int i){ return (Dxs[i] * 0.5) / Dx[i];}
        inline auto y_INTRPL(int j){ return (Dys[j] * 0.5) / Dy[j];}
        inline auto z_INTRPL(int k){ return (Dzs[k] * 0.5) / Dz[k];}

        // *  NEIBcell ----------------------------------
        // const auto [xm, xp, ym, yp, zn, zp, xmm, xpp, ymm, ypp, znn, zpp] = getNb12(icel);
        template<typename T1> 
        std::tuple<int, int, int , int , int, int,
                   int, int, int, int, int, int> getNb12(T1 &icel){
                const auto i  = icel *12;
                return make_tuple(
                        NEIBcell[i+0 ], /* - 1x + */ NEIBcell[i+1 ],
                        NEIBcell[i+2 ], /* - 1y + */ NEIBcell[i+3 ],   
                        NEIBcell[i+4 ], /* - 1z + */ NEIBcell[i+5 ], 
                        NEIBcell[i+6 ], /* - 2x + */ NEIBcell[i+7 ],
                        NEIBcell[i+8 ], /* - 2y + */ NEIBcell[i+9 ], 
                        NEIBcell[i+10], /* - 2z + */ NEIBcell[i+11] );
        }


        template<typename T1> 
        std::tuple<int, int, int , int , int, int> getNb6(T1 &icel){
        
                const auto ii  = icel *12;
                return make_tuple(
                        NEIBcell[ii+0], /* - x + */ NEIBcell[ii+1],
                        NEIBcell[ii+2], /* - y + */ NEIBcell[ii+3], 
                        NEIBcell[ii+4], /* - z + */ NEIBcell[ii+5] );
        }

        template<typename T1> 
        std::tuple<int, int, int , int , int, int> getNb6(T1 &i, T1 &j, T1 &k){

                const auto ii  = (i*nzny + j*nz + k)*12;
                return make_tuple(
                        NEIBcell[ii+0], /* - x + */ NEIBcell[ii+1],
                        NEIBcell[ii+2], /* - y + */ NEIBcell[ii+3], 
                        NEIBcell[ii+4], /* - z + */ NEIBcell[ii+5] );
        }

        template<typename T1, typename T2, typename T3>
        bool initNb(T1 &i, T2 &j, T3 &k){
                const auto iC  =  i*nzny + j*nz +k;
                const auto ii  = (i*nzny + j*nz +k)*12;

                NEIBcell.at(ii+0 ) = iC - nzny;     //-x 
                NEIBcell.at(ii+1 ) = iC + nzny;     //+x 
                NEIBcell.at(ii+2 ) = iC - nz;       //-y 
                NEIBcell.at(ii+3 ) = iC + nz;       //+y 
                NEIBcell.at(ii+4 ) = iC - 1;        //-z 
                NEIBcell.at(ii+5 ) = iC + 1;        //+z 

                NEIBcell.at(ii+6 ) = iC - 2*nzny;  //-2x
                NEIBcell.at(ii+7 ) = iC + 2*nzny;  //+2x
                NEIBcell.at(ii+8 ) = iC - 2*nz;    //-2y
                NEIBcell.at(ii+9 ) = iC + 2*nz;    //+2y
                NEIBcell.at(ii+10) = iC - 2;       //-2z
                NEIBcell.at(ii+11) = iC + 2;       //+2z
                return true;
        }




        bool init_gC_dx(){

                X[1] = 2*X[2] - X[3]; X[0] = 2*X[1] - X[2];
                Y[1] = 2*Y[2] - Y[3]; Y[0] = 2*Y[1] - Y[2];
                Z[1] = 2*Z[2] - Z[3]; Z[0] = 2*Z[1] - Z[2];

                X[nx-1] = 2*X[nx-2] - X[nx-3];
                Y[ny-1] = 2*Y[ny-2] - Y[ny-3];
                Z[nz-1] = 2*Z[nz-2] - Z[nz-3];

                X[nx] = 2*X[nx-1] - X[nx-2];
                Y[ny] = 2*Y[ny-1] - Y[ny-2];
                Z[nz] = 2*Z[nz-1] - Z[nz-2];

                for (auto i = gC ; i < nx-gC-1; ++i)
                        Dxs[i] = ( X[i+2] - X[i] ) / 2.0;

                for (auto j = gC; j < ny-gC-1; ++j)
                        Dys[j] = ( Y[j+2] - Y[j] ) / 2.0;

                for (auto k = gC; k < nz-gC-1; ++k)
                        Dzs[k] = ( Z[k+2] - Z[k] ) / 2.0;


                for (auto i = gC ; i < nx-gC-1; ++i)
                        Dx[i] = ( X[i+1] - X[i] );

                for (auto j = gC; j < ny-gC-1; ++j)
                        Dy[j] = ( Y[j+1] - Y[j] );

                for (auto k = gC; k < nz-gC-1; ++k)
                        Dz[k] = ( Z[k+1] - Z[k] );


                // Ghost boundary grid 
                Dx[0]     = Dx[1]     = Dx[2];
                Dy[0]     = Dy[1]     = Dy[2];
                Dz[0]     = Dz[1]     = Dz[2];

                Dxs[0]    = Dxs[1]    = Dxs[2];
                Dys[0]    = Dys[1]    = Dys[2];
                Dzs[0]    = Dzs[1]    = Dzs[2];


                Dx[nx-1]  = Dx[nx-2]  = Dx[nx-3]  = X[nx-2] - X[nx-3];
                Dy[ny-1]  = Dy[ny-2]  = Dy[ny-3]  = Y[ny-2] - Y[ny-3];
                Dz[nz-1]  = Dz[nz-2]  = Dz[nz-3]  = Z[nz-2] - Z[nz-3];

                Dxs[nx-1] = Dxs[nx-2] = Dxs[nx-3] = Dxs[nx-4];
                Dys[ny-1] = Dys[ny-2] = Dys[ny-3] = Dys[ny-4];
                Dzs[nz-1] = Dzs[nz-2] = Dzs[nz-3] = Dzs[nz-4];

                return true;
        }


auto createPressureMatrixCSR_forDir(){
        // TODO: Debug this scope!!
    auto xp = [&](auto i){return 1.0/Dx.at(i+gC)/Dxs.at(i+gC);};
    auto yp = [&](auto j){return 1.0/Dy.at(j+gC)/Dys.at(j+gC);};
    auto zp = [&](auto k){return 1.0/Dz.at(k+gC)/Dzs.at(k+gC);};

    auto xm = [&](auto i){return 1.0/Dx.at(i+gC)/Dxs.at(i+gC);};
    auto ym = [&](auto j){return 1.0/Dy.at(j+gC)/Dys.at(j+gC);};
    auto zm = [&](auto k){return 1.0/Dz.at(k+gC)/Dzs.at(k+gC);};

    auto n3 = (nyCal+gC-1) * (nyCal+gC-1) * (nyCal+gC-1);
    std::vector<int>  ptr, idx;
    std::vector<double>  val;

    ptr.clear(); ptr.reserve(n3 + 1); ptr.push_back(0);
    idx.clear(); idx.reserve(n3 * 7); // We use 7-point stencil, so the matrix
    val.clear(); val.reserve(n3 * 7); // will have at most n3 * 7 nonzero elements.
    // Interior point. Use 7-point finite difference stencil.


    int 
     xBeg = 0, yBeg = 0, zBeg = 0,
     xEnd = nxCal+1, yEnd = nyCal+1, zEnd = nzCal+1,
     sX = nznyDir, sY = nzDir, sZ = 1;

    for(int i = 0, iC = 0; i < nxCal+2; ++i) 
    {
        for(int j = 0 ; j < nyCal+2; ++j) 
        {
            for(int k = 0; k < nzCal+2; ++k, ++iC) 
            {
                if ( i==xBeg || i==xEnd  || j==yBeg || j==yEnd || k==zBeg || k == zEnd ){
                        idx.push_back(iC);
                        val.push_back(1.0);
                }else{
                        // ----------------------
                        idx.push_back(iC - sX);
                        val.push_back(-xm(i));

                        idx.push_back(iC - sY);
                        val.push_back(-ym(j));

                        idx.push_back(iC - sZ);
                        val.push_back(-zm(k));
                        // ----------------------

                        // ----------------------
                        idx.push_back(iC + sZ);
                        val.push_back(-zp(k));

                        idx.push_back(iC + sY);
                        val.push_back(-yp(j));

                        idx.push_back(iC + sX);
                        val.push_back(-xp(i));
                        // ----------------------
                }
                ptr.push_back(idx.size());
             }
        }
    }



    return std::make_tuple(ptr, idx, val);

}

auto createPressureMatrixCSR(){


    auto xp = [&](auto i){return 1.0/Dx.at(i+gC)/Dxs.at(i+gC);};
    auto yp = [&](auto j){return 1.0/Dy.at(j+gC)/Dys.at(j+gC);};
    auto zp = [&](auto k){return 1.0/Dz.at(k+gC)/Dzs.at(k+gC);};

    auto xm = [&](auto i){return 1.0/Dx.at(i+gC)/Dxs.at(i+gC-1);};
    auto ym = [&](auto j){return 1.0/Dy.at(j+gC)/Dys.at(j+gC-1);};
    auto zm = [&](auto k){return 1.0/Dz.at(k+gC)/Dzs.at(k+gC-1);};

    auto n3 = iceltotCal;
    std::vector<int> ptr, idx;
    std::vector<double>  val;
    jp.resize(n3);

    ptr.clear(); ptr.reserve(n3 + 1); ptr.push_back(0);
    idx.clear(); idx.reserve(n3 * 7); // We use 7-point stencil, so the matrix
    val.clear(); val.reserve(n3 * 7); // will have at most n3 * 7 nonzero elements.
    // Interior point. Use 7-point finite difference stencil.

    int 
     xBeg = 0, yBeg = 0, zBeg = 0,
     xEnd = nxCal-1, yEnd = nyCal-1, zEnd = nzCal-1;


    for(int i = 0, iCont = 0; i < nxCal; ++i) 
    {
        for(int j = 0 ; j < nyCal; ++j) 
        {
            for(int k = 0; k < nzCal; ++k, ++iCont) 
            {
                // Boundary point. Use Num condition.
                double  aa = 0.0, bb = 0.0, cc = 0.0;
                if(i==xBeg){aa = -xm(i);}
                if(i==xEnd){aa = -xp(i);}
                    
                if(j==yBeg){bb = -ym(j);}
                if(j==yEnd){bb = -yp(j);}

                if(k==zBeg){cc = -zm(k);}
                if(k==zEnd){cc = -zp(k);}
                auto ijk = xm(i) + ym(j) + zm(k) + 
                           xp(i) + yp(j) + zp(k) +
                            aa   + bb    + cc ;
                #if defined (JACOBI_PC)
                jp[iCont] = 1.0/ijk;        // Jacobi precondition on
                #else
                jp[iCont] = 1.0;            // Jacobi precondition off
                #endif
                // ----------------------
                if (i!=xBeg){
                idx.push_back(iCont - nznyCal);
                val.push_back(-xm(i)*jp[iCont]);
                }

                if (j!=yBeg){
                idx.push_back(iCont - nzCal);
                val.push_back(-ym(j)*jp[iCont]);
                }

                if (k!=zBeg){
                idx.push_back(iCont - 1);
                val.push_back(-zm(k)*jp[iCont]);
                }

                // ----------------------

                // ----------------------
                idx.push_back(iCont);
                val.push_back(ijk*jp[iCont]);
                // ----------------------

                // ----------------------

                if (k!=zEnd){
                idx.push_back(iCont + 1);
                val.push_back(-zp(k)*jp[iCont]);
                }

                if (j!=yEnd){
                idx.push_back(iCont + nzCal);
                val.push_back(-yp(j)*jp[iCont]);
                }

                if (i!=xEnd){
                idx.push_back(iCont + nznyCal);
                val.push_back(-xp(i)*jp[iCont]);
                }
                // ----------------------

                ptr.push_back(idx.size());
            }

        }
    }

    return std::make_tuple(ptr, idx, val);
}




void checkgC(){
        // check == X[nx - gC] == lx
        // -----------------------
        if (abs(X[nx - gC] - lx) > 1.0e-5){
                cout << nx - gC << "\t" << X[nx - gC] << "\t" << lx << std::endl;
                throw std::invalid_argument("X[nx - gC]");
        }

        if (abs(Y[ny - gC] - ly) > 1.0e-5){
                cout << ny - gC << "\t" << Y[ny - gC] << "\t" << ly << std::endl;
                throw std::invalid_argument("Y[ny - gC]");
        }

        if (abs(Z[nz - gC] - lz) > 1.0e-5){
                cout << nz - gC << "\t" << Z[nz - gC] << "\t" << lz << std::endl;
                throw std::invalid_argument("Z[nz - gC]");
        }
        // -----------------------

        // -----------------------
        for (auto i = 0; i < nx ; ++i ) {
                if (X[i] >= X[i+1])
                        throw std::invalid_argument("X[i] i:");
        }

        for (auto j = 0; j < ny ; ++j ) {
                if (Y[j] >= Y[j+1])
                        throw std::invalid_argument("Y[j]");
        }

        for (auto k = 0; k < nz ; ++k ) {
                if (Z[k] >= Z[k+1])
                        throw std::invalid_argument("Z[k]");
        }
}



bool io_csv(string filename){

        std::ofstream file;

        file.open (filename,std::ios::out); // | ios::app

        // * --------------------- x---------------------
        for (size_t i = 0; i < nx+1 ; ++i ) {

                file << "X[" << i << "] ," ; 
                
                file << X[i] << ", ";
                
                if (i < nx-gC+2) {
                
                        if (i < 2){file << "Xc[gC], not define, ";}
                        else if(i < nx-gC){ file << "Xc[" << i-2 << "] ," << Xc[i-2] << ", ";}
                        else { file << "Xc[gC], not define, "; }
                                
                }else{ file<< " , , ";}
                        
                if (i < nx) {
                        file << "Dx[" << i << "] ," << Dx[i] << ", "
                             << "Dxs[" << i << "] ," << Dxs[i] ;

                }else{ file<< " , , ";}

                file << std::endl ;

        }


                file << std::endl ;
        // * --------------------- y---------------------
        for (size_t j = 0; j < ny+1 ; ++j ) {

                file << "Y[" << j << "] ," ; 
                file << Y[j] << ", ";

                if (j <  ny-gC+2) {

                        if (j < 2){file << "Yc[gC], not define, ";}
                        else if (j < ny - gC){ file << "Yc[" << j-2 << "] ," <<  Yc[j-2] << ", "; }
                        else{ file << "Yc[gC], not define, "; }

                }
                else{ file << " , , ";}


                if (j < ny ) {

                        file    << "Dy[" << j << "] ,"
                                << Dy[j] << ", "
                                << "Dys[" << j << "] ," 
                                << Dys[j];
                }else{  file << " , , ";}

                file << std::endl ;

        }

                file << std::endl ;
        // *---------------------------- z --------------------
        for (size_t k = 0; k < nz+1 ; ++k ) {

                file << "Z[" << k << "] ," << Z[k] << ", ";

                if ( k < nz-gC+2 ) {

                        if (k < 2){ file << "Zc[gC], not define,"; }
                        else if (k < nz-gC){ file << "Zc[" << k-2 << "] ," << Zc[k-2] << ", "; }
                        else { file << "Zc[gC], not define,";}

                }else{ file<< " , , ";}



                if ( k < nz ) {
                        file    << "Dz[" << k << "] ," << Dz[k] << ", "
                                << "Dzs[" << k << "] ," << Dzs[k];
                }
                else{ file<< " , , "; }

                file << std::endl;
        }

        checkgC();
        return true;

}


};

