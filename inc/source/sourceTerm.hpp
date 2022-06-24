#pragma once

#include "controlPanel.hpp"
using std::cout;
using std::endl;
using std::string;
using std::vector;

#ifndef TEMPORAL_DISCRETIZATION_1_ORDER

inline vector<double> discretizeTemporal(simuClass& simu){

    // !----------------------------------------------
    #if defined (TEMPORAL_DISCRETIZATION_2_ORDER)

    if(simu.First_s) {T0.resize_s(gridA.iceltotCal);}


    double h1 = simu.dt_pre[0];
    double h2 = simu.dt_pre[1];

    vector<double> C(2);

    C[0] =  0.5*h1*(h1+2.0*h2)/h2;
    C[1] = -0.5*h1*h1/h2;
    return C;

    #elif defined (TEMPORAL_DISCRETIZATION_3_ORDER)
    // !----------------------------------------------
    auto p = [&](auto v, auto n){return std::pow(v, n);};

    if (simu.First_s){ T0.resize_s(gridA.iceltotCal*2); }

    simu.record_previous_dt();

    double h1 = simu.dt_pre[0];
    double h2 = simu.dt_pre[1];
    double h3 = simu.dt_pre[2];

    vector<double> C(3);

    C[0] =  2.0*p(h1,3.0) + 3.0*p(h1,2.0)*(2.0*h2+h3) + 6.0*h1*h2*(h2+h3)
            / (6.0*h2 *(h2+h3));

    C[1] = -( 2.0*p(h1,3.0) + 3.0*p(h1,2.0)*(h2 + h3) )
            / (6.0* h2 * h3);

    C[2] = (2.0*p(h1,3.0) + 3.0*p(h1,2.0)*h2)
            / (6.0 * h3 *(h2+h3));

    if (simu.Second_s){
        C[0] = 0.5*h1*(h1+2.0*h2)/h2;
        C[1] = -0.5*h1*h1/h2;
    }


    return C;

    // ----------------------------------------------
    #endif
}
    #endif


inline auto souce(  vector<double> &Unew, const vector<double> &U, 
                    double sTemp, int icel) { Unew[icel] = U[icel] + sTemp; }


inline auto souce(  vector<double> &Unew,
                    const  vector<double> &U, 
                    vector<double> &Us,
                    double sTemp, 
                    vector<double> &C, 
                    int icel, int icelCal){

    #if defined (TEMPORAL_DISCRETIZATION_2_ORDER)
    // !----------------------------------------------

        Us[icelCal] = sTemp;
        if (simu.loop == 1) { return U[icel] + sTemp; }
        else { Unew[icel] = U[icel] + C[0] * sTemp + C[1] * Us[icelCal];}


    #elif defined (TEMPORAL_DISCRETIZATION_3_ORDER)
    // !----------------------------------------------

        Us[icelCal*2 + 1] = Us[icelCal*2];
        Us[icelCal*2 ] = temp;

        if (simu.First_s)
        {
            Unew[icel] = U[icel] + sTemp;
        }
        else if (simu.Second_s)
        {
            Unew[icel] =  U[icel] 
                        + C[0] * sTemp 
                        + C[1] * Us[icelCal*2];
        }
        else
        {
            Unew[icel] =  U[icel] 
                        + C[0] * sTemp 
                        + C[1] * Us[icelCal*2]
                        + C[2] * Us[icelCal*2 + 1];
        }


    #endif

}
