#pragma once

#include"controlPanel.hpp"

inline auto virtualF_Int(
    DfibArray& Dfib,
    calDomain& Lo,
    grid& gA,
    const int whichDirection
)
{
    if (Dfib.ValSum.size() < 3)
        Dfib.ValSum.resize(3);

    const int ptr = gA.iceltotCal*whichDirection;
    double sum = 0.0;

    // #pragma omp parallel for firstprivate(Lo, ptr) reduction(+ : sum)
    // for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
    // for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
    // for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
    // {
    //    sum += Dfib.f[ ptr + gA.icelCal(i,j,k) ] * (gA.Dx[i]*gA.Dy[j]*gA.Dz[k]) ;
    // }

    #pragma omp parallel for firstprivate(Lo, ptr) reduction(+ : sum)
    for (auto i = Lo.i_begin ; i < Lo.i_endof ; ++i )
    for (auto j = Lo.j_begin ; j < Lo.j_endof ; ++j )
    for (auto k = Lo.k_begin ; k < Lo.k_endof ; ++k )
    {
       sum += Dfib.f[ ptr + gA.icelCal(i,j,k) ] * gA.delta(i,j,k);
    }

    Dfib.ValSum[whichDirection] = sum;

    return 0;
}



std::pair<double, double>
noDelay_IO_CD_CL(
    simuClass& simu,
    calDomain& Lo,
    DfibArray& Dfib,
    grid& gA
)
{
    virtualF_Int(Dfib, Lo, gA, 0);
    virtualF_Int(Dfib, Lo, gA, 1);

    auto area_cD    = gA.lz;
    auto area_cL    = gA.lz;

    auto cD = -2.0 * Dfib.ValSum[0] / area_cD;
    auto cL = -2.0 * Dfib.ValSum[1] / area_cL;


    if (simu.PID == 0)
    {

    cout << "[cD, cL] : " << cD << ", " << cL << endl; 

    std::ofstream file;

    std::string name = "Information/Time_cDcL";

    name += ".dat";

    file.open (name, std::ios::out|ios::app);

    std::string tab = " ";

    if (simu.loop == 1 ){

        std::vector<std::string> variables;
        variables.push_back("simulation time");
        variables.push_back("C<sub>D</sub>");
        variables.push_back("C<sub>L</sub>");

            file 
            << "TITLE     = \"\"\n"
            << "VARIABLES = \""
            << variables.at(0)
            << "\",\""
            << variables.at(1)
            << "\",\""
            << variables.at(2)
            << "\"\n"
            << "ZONE T=\""
            << simu.ZONE()
            << "\"";
    }

    file << "\n"<<  simu.time << tab << cD << tab  << cL ;

    file.close();
    }

    return make_pair(cD, cL);
}





std::pair<double, double>
Delay_IO_CD_CL(
    simuClass& simu,
    calDomain& Lo,
    DfibArray& Dfib,
    grid& gA
)
{

    std::string name = "Information/Time_cDcL";
    name += ".dat";


    virtualF_Int(Dfib, Lo, gA, 0);
    virtualF_Int(Dfib, Lo, gA, 1);

    auto area_cD = gA.lz;
    auto area_cL = gA.lz;

    auto cd = -2.0 * Dfib.ValSum[0] / area_cD;
    auto cl = -2.0 * Dfib.ValSum[1] / area_cL;

    simu.cd[ simu.loop % simu.delayIO ]  = cd;
    simu.cl[ simu.loop % simu.delayIO ]  = cl;

    bool IOfale = simu.loop % simu.delayIO ;

    if (simu.PID != 0 )     return make_pair(cd, cl);

    
    if (simu.loop == 1 )
    {
        std::ofstream file;
        file.open (name, std::ios::out|ios::app);
        std::vector<std::string> variables;
        variables.push_back("simulation time");
        variables.push_back("C<sub>D</sub>");
        variables.push_back("C<sub>L</sub>");

            file 
            << "TITLE     = \"\"\n"
            << "VARIABLES = \""
            << variables.at(0)
            << "\",\""
            << variables.at(1)
            << "\",\""
            << variables.at(2)
            << "\"\n"
            << "ZONE T=\""
            << simu.ZONE()
            << "\"";
    }


    if (!(IOfale))
    {
        std::ofstream file;
        file.open (name, std::ios::out|ios::app);
        std::string tab = " ";
        for (int i = 1; i < simu.delayIO+1 ; ++i)
        file << "\n"
         << simu.time
         << tab << simu.cd[simu.loop % simu.delayIO]
         << tab  << simu.cl[simu.loop % simu.delayIO] ;

        file.close();
    }


    return make_pair(cd, cl);
}