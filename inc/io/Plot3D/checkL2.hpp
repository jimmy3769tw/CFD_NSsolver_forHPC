#pragma once
#include"controlPanel.hpp"
#include "qfileRead.hpp"



void ReadPlot3DadnCheckL2(
    shareMenory& ShareM,
    int main_count,
    DfibArray& Dfib,
    simuClass& simu,
    pressure& t1,
    velocity& T3,
    calDomain& Lo,
    grid& gA
){


    // * read the file *.q
    ShareM.TheFileIsOpen = 1;

    atuo [nx, ny, nz , gC]  = gA.nxyzgC;



	std::string filename = "mx_in/P3D";

	filename += std::to_string(io);

	filename += ".q";


    double time = double(main_count) * simu.dt;

    vector<float> pD(gA.iceltotCal), ucD(gA.iceltotCal), vcD(gA.iceltotCal), wcD(gA.iceltotCal);


    auto [Nblock_read, temp_read, mach_read, alpha_read, reyn_read, 
            pPre, ucPre, vcPre, wcPre, EtaPre] = QfileRead(gA, filename);



    // * 2. processor the file 

    for (size_t i = Lo.i_begin; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin; j < Lo.j_endof ; ++j )
    for (size_t k = Lo.k_begin; k < Lo.k_endof ; ++k )
        pD[gA.icelCal(i,j,k)] = pPre[gA.icelCal(i,j,k)] 
                                - float(t1.p[gA.icel(i,j,k)]);



    for (size_t i = Lo.i_begin; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin; j < Lo.j_endof ; ++j )
    for (size_t k = Lo.k_begin; k < Lo.k_endof ; ++k )
        ucD[gA.icelCal(i,j,k)] = ucPre[gA.icelCal(i,j,k)] 
                                -float(0.5 * (T3.u[gA.icel(i,j-1,k)] + T3.u[gA.icel(i,j,k)]));


    for (size_t i = Lo.i_begin; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin; j < Lo.j_endof ; ++j )
    for (size_t k = Lo.k_begin; k < Lo.k_endof ; ++k )
        vcD[gA.icelCal(i,j,k)] =  vcPrevious[gA.icelCal(i,j,k)] 
                                 - float(0.5 * (T3.v[gA.icel(i,j,k)] + T3.v[gA.icel(i,j,k)]));

    for (size_t i = Lo.i_begin; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin; j < Lo.j_endof ; ++j )
    for (size_t k = Lo.k_begin; k < Lo.k_endof ; ++k )
        wcD[gA.icelCal(i,j,k)] = wcPrevious[gA.icelCal(i,j,k)] 
                                 - float(0.5 * (T3.w[gA.icel(i,j,k-1)] + T3.w[gA.icel(i,j,k)]));



    // * 3. Calculate the L2 norm
    std::vector<float>  L2_temp(4, 0.0);

    #pragma omp parallel for reduction(+: L2_temp) 
    for (size_t i = Lo.i_begin; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin; j < Lo.j_endof ; ++j )
    for (size_t k = Lo.k_begin; k < Lo.k_endof ; ++k )
    {

        L2_temp[0] += std::pow(pD[gA.icelCal(i,j,k)],2);
        L2_temp[1] += std::pow(ucD[gA.icelCal(i,j,k)],2);
        L2_temp[2] += std::pow(vcD[gA.icelCal(i,j,k)],2);
        L2_temp[3] += std::pow(wcD[gA.icelCal(i,j,k)],2);
    }


    cout<< ", L2norm (uc, vc, wc):" ;

    for(auto x :L2){

        cout << x / gA.iceltotCal << ", ";
    }

// * 4. plot the differernt


    std::ofstream file;

    std::string filename = "mx_Diffrent/P3D";

    filename += std::to_string(simu.loop);

    filename += ".q";	

    file.open(filename,std::ofstream::binary);


    if (file.is_open()){

        //*  mach, alpha, reyn, time //
        float temp = simu.time, mach = 0, alpha = 0, reyn = simu.Re;

        file.write((char *)(&Nblock), sizeof(int));  
        file.write((char *)(&nxCal),  sizeof(int));
        file.write((char *)(&nyCal),  sizeof(int));
        file.write((char *)(&nzCal),  sizeof(int));

        file.write((char *)(&mach),  sizeof(float));
        file.write((char *)(&alpha), sizeof(float));
        file.write((char *)(&reyn),  sizeof(float));
        file.write((char *)(&time),  sizeof(float));

        for (size_t k = 0; k < nzCal; ++k) 
        for (size_t j = 0; j < nyCal; ++j) 
        for (size_t i = 0; i < nxCal; ++i)
            file.write((char *)(&pD[ii(i,j,k)]), sizeof(float));


        for (size_t k = 0; k < nzCal; ++k) 
        for (size_t j = 0; j < nyCal; ++j) 
        for (size_t i = 0; i < nxCal; ++i)
            file.write((char *)(&ucD[ii(i,j,k)]), sizeof(float));


        for (size_t k = 0; k < nzCal; ++k) 
        for (size_t j = 0; j < nyCal; ++j) 
        for (size_t i = 0; i < nxCal; ++i)
            file.write((char *)(&vcD[ii(i,j,k)]), sizeof(float));


        for (size_t k = 0; k < nzCal; ++k) 
        for (size_t j = 0; j < nyCal; ++j) 
        for (size_t i = 0; i < nxCal; ++i)
            file.write((char *)(&wcD[ii(i,j,k)]), sizeof(float));


        for (size_t k = 0; k < nzCal; ++k) 
        for (size_t j = 0; j < nyCal; ++j) 
        for (size_t i = 0; i < nxCal; ++i)
            file.write((char *)(&EtaPre[ii(i,j,k)]), sizeof(float));
    

        file.close();

    }

    return true;

}