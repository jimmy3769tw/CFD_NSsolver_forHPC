#pragma once

#include "controlPanel.hpp"

#include "mpi.h"
//   auto [mpi_word_size, mpi_word_rank, mpi_coord, mpi_neighborhood] = 
//      mpi_init(reorder, XYZ, argc, argv);

std::tuple< MPI_Comm, int, int, std::vector<int>, std::vector<int> >
mpi_init(
    bool &reorder,
    std::vector<int> &XYZ,
    int argc, 
    char **argv
)
{

    MPI_Init(&argc, &argv);

    int mpi_word_rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_word_rank);

    int mpi_word_size;

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_word_size);

    // int Check = 1;

    // for (auto dims:XYZ){
    //   Check *= dims;
    // }
   
    // if (Check!=mpi_word_size) {
    //     MPI_Barrier(comm_world);
    //     MPI_Finalize();
    //    throw std::invalid_argument("MPI size doesn't correct.");
    // }

    // * virtual process topology functions 
    // * number of dimensions of cartesian grid
    const int ndims = XYZ.size();

    MPI_Comm comm_world;

    int periods[ndims] = {0}; 

    MPI_Cart_create(MPI_COMM_WORLD, ndims, &XYZ[0], periods, reorder, &comm_world);

    std::vector<int> mpi_coord(ndims);
    
    MPI_Cart_coords(comm_world, mpi_word_rank, mpi_coord.size(), &mpi_coord[0]);

    std::vector<int> mpi_neighborhood(2*ndims);

    auto disp = 1;
    for (size_t i = 0; i < ndims ; ++i){
      MPI_Cart_shift(comm_world, i, disp, 
        &mpi_neighborhood[i*2], &mpi_neighborhood[i*2+1]);
    }

    cout << "[ank, coord ] " << mpi_word_rank <<", {";
    for(size_t i = 0; i < mpi_coord.size(); i++){
      cout << mpi_coord.at(i) ;
      if (i != mpi_coord.size()-1)
        cout << ", ";
    }

    for (size_t i = 0 ; i < ndims ;++i){

      std::cout << "}>> Neberhood | " 
        << mpi_neighborhood[i*2]<< ", "
        << mpi_neighborhood[i*2+1] << endl;
    }

    // ! ## INFO
    char name[1024];
    int length=1024, minor;
    MPI_Get_processor_name(name, &length);

    int major;
    MPI_Get_version(&major, &minor);
    
    // MPI_Request req[10];

    if (mpi_word_rank == 0){
        cout << "\nMPI Version " << major << "." << minor << endl;
        cout << "This Project is from " << name << endl;
    }


    return std::tie(comm_world, mpi_word_size, mpi_word_rank, mpi_coord, mpi_neighborhood);
}

