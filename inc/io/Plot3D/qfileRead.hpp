#pragma once

#include "controlPanel.hpp"


auto QfileRead(
    grid& gA,
    std::string filename
)
{
    //*  prepare the data to read
    float mach, alpha, reyn, time;
    int Nblock, nxCal, nyCal, nzCal;


    uint nodata = 5;

    std::vector<std::vector<float > >
            data(nodata, std::vector<float>(gA.iceltotCal, 0.0));

    std::vector<std::vector<double > >
            data_new(nodata, std::vector<double>(gA.iceltot, 0.0));


    // * read data;

    std::ifstream file;
	file.open(filename,std::ifstream::binary);


    auto ii = [&](auto i, auto j, auto k ){
            return i*nzCal*nyCal + j*nzCal +k;};

    if (file.is_open()){

        file.read((char *)(&Nblock), sizeof(int));  
        file.read((char *)(&nxCal), sizeof(int));
        file.read((char *)(&nyCal), sizeof(int));
        file.read((char *)(&nzCal), sizeof(int));

        file.read((char *)(&mach), sizeof(float));
        file.read((char *)(&alpha), sizeof(float));
        file.read((char *)(&reyn), sizeof(float));
        file.read((char *)(&time), sizeof(float));

        for (auto & x:data){

            for (size_t k = 0; k < nzCal; ++k)
            for (size_t j = 0; j < nyCal; ++j)
            for (size_t i = 0; i < nxCal; ++i) 
                file.read((char *)(&x[ii(i,j,k)]), sizeof(float));

        }

        file.close();
    }
    else{

         file.close();
         std::cout << "I can't Read The File." << std::endl;
    }



    for (size_t k = gA.gC; k < nzCal + gA.gC; ++k) 
    for (size_t j = gA.gC; j < nyCal + gA.gC; ++j) 
    for (size_t i = gA.gC; i < nxCal + gA.gC; ++i) 
    {
        data_new[0][gA.icel(i,j,k)] = data[0][gA.icelCal(i,j,k)];
        data_new[4][gA.icel(i,j,k)] = data[4][gA.icelCal(i,j,k)];
    }

    for (size_t k = gA.gC; k < nzCal + gA.gC; ++k) 
    for (size_t j = gA.gC; j < nyCal + gA.gC; ++j) 
    for (size_t i = gA.gC; i < nxCal + gA.gC; ++i) 
    {

        if (nxCal + gA.gC -1 == i)
            data_new[1][gA.icel(i,j,k)] = data[1][ gA.icelCal(i,j,k)] ; 
        else
            data_new[1][gA.icel(i,j,k)] = (data[1][ gA.icelCal(i,j,k)] + data[1][ gA.icelCal(i+1,j,k)]) / 2.0; 


        if (nyCal + gA.gC -1 == j)
            data_new[2][gA.icel(i,j,k)] = data[2][ gA.icelCal(i,j,k)];
        else
            data_new[2][gA.icel(i,j,k)] = (data[2][ gA.icelCal(i,j,k)] + data[2][ gA.icelCal(i,j+1,k)]) / 2.0;


        if (nzCal + gA.gC -1 == k)
            data_new[3][gA.icel(i,j,k)] = data[3][ gA.icelCal(i,j,k)] ;
        else
            data_new[3][gA.icel(i,j,k)] = (data[3][ gA.icelCal(i,j,k)] + data[3][ gA.icelCal(i,j,k+1)]) / 2.0;

    }
    

    return std::make_tuple(Nblock, mach, alpha, reyn, time,
            data_new[0], data_new[1], data_new[2], data_new[3], data_new[4]);
}

