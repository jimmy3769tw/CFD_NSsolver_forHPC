
// #define MPI_NONBLOCKING_SR
#define MPI_BLOCKING_SR


// mx_comm_world, t1.p, localDomain, gridA)

void mpi_iSR_double_x(
  int nolayers,
  MPI_Comm &comm_world,
  std::vector<int> &mpi_neighborhood,
  std::vector<double> &v,
  calDomain& localDomain,
  grid& gridA
){

  if(nolayers < 0){
    throw std::invalid_argument("nolayers < 0 is not suporting!!");
  }

  MPI_Barrier(comm_world);

  MPI_Request requests[4];

  int nynz = gridA.nz * gridA.ny;

  int itag[2];
  
  int nolayers_nynz = nolayers * nynz;
  
  int shift_recv[2];

  int shift_send[2];

  itag[0] = 110;

  itag[1] = 220;
  // make a tuple 
  
  shift_send[0] = (localDomain.i_begin) * nynz;
  shift_recv[0] = (localDomain.i_endof) * nynz;// !`+`

    // const int icel_endof = (localDomain.i_endof_table.at(i_rank)-1)*nz*ny + (ny-1)*nz + (nz-1);
    // const int count = 1- icel_begin + icel_endof;

  shift_send[1] = (localDomain.i_endof - nolayers) * nynz;
  shift_recv[1] = (localDomain.i_begin - nolayers) * nynz;// !`-`

  #ifdef MPI_NONBLOCKING_SR
    // !`+`
    MPI_Isend(&v[shift_send[0]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(0), itag[0], comm_world, requests+0);
    MPI_Irecv(&v[shift_recv[0]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(1), itag[0], comm_world, requests+1);
    
    // !`-`
    MPI_Isend(&v[shift_send[1]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(1), itag[1], comm_world, requests+2);
    MPI_Irecv(&v[shift_recv[1]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(0), itag[1], comm_world, requests+3);
   
    MPI_Status status[4];

    MPI_Waitall(4, requests, status);
    
  #endif

  #ifdef MPI_BLOCKING_SR
    MPI_Status status[2];
    // !`+`
    MPI_Send(&v[shift_send[0]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(0), itag[0], comm_world);
    MPI_Recv(&v[shift_recv[0]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(1), itag[0], comm_world, status);
    
    // !`-`
    MPI_Send(&v[shift_send[1]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(1), itag[1], comm_world);
    MPI_Recv(&v[shift_recv[1]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(0), itag[1], comm_world, status+1);
  #endif

  MPI_Barrier(comm_world);
}



void mpi_iSR_double_x_half(
  int nolayers,
  MPI_Comm &comm_world,
  std::vector<int> &mpi_neighborhood,
  std::vector<double> &v,
  calDomain& localDomain,
  grid& gridA
){

  MPI_Barrier(comm_world);

  int itag[2];

  const auto abs_nolayers = std::abs(nolayers);
  
  const int nolayers_nynz = abs_nolayers* gridA.nz *gridA.ny;
  
  const int nynz = gridA.nz * gridA.ny;

  int shift_recv[2];

  int shift_send[2];

  itag[0] = 100;
  itag[1] = 200;

  shift_send[0] = (localDomain.i_begin) * nynz;
  shift_recv[0] = (localDomain.i_endof) * nynz;// !`+`
  shift_send[1] = (localDomain.i_endof-abs_nolayers) * nynz;
  shift_recv[1] = (localDomain.i_begin-abs_nolayers) * nynz;// !`-`

  if (nolayers > 0){

  #ifdef MPI_NONBLOCKING_SR
    MPI_Request requests[2];
    // * --------------------- `+`*
    MPI_Isend(&v[shift_send[0]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(0), itag[0], comm_world, requests+0);
    MPI_Irecv(&v[shift_recv[0]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(1), itag[0], comm_world, requests+1);
    }
    else{
    // * --------------------- `-`*
    MPI_Request requests[2];
    MPI_Isend(&v[shift_send[1]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(1), itag[1], comm_world, requests+0);
    MPI_Irecv(&v[shift_recv[1]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(0), itag[1], comm_world, requests+1);
  #endif



  #ifdef MPI_BLOCKING_SR
    MPI_Status status[1];
    // * --------------------- `+`*
    MPI_Send(&v[shift_send[0]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(0), itag[0], comm_world);
    MPI_Recv(&v[shift_recv[0]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(1), itag[0], comm_world, status);
    }
    else{
    // * --------------------- `-`*
    MPI_Status status[1];
    MPI_Send(&v[shift_send[1]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(1), itag[1], comm_world);
    MPI_Recv(&v[shift_recv[1]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(0), itag[1], comm_world, status);
  #endif


  }
  #ifdef MPI_NONBLOCKING_SR
  MPI_Status status[2];

  MPI_Waitall(2, requests, status);
  #endif
  
  MPI_Barrier(comm_world);
}



void mpi_iSR_double_x_Collect_to_Master(
  MPI_Comm &comm_world,
  const int rank, const int rank_size,
  std::vector<double> &v,
  calDomain& localDomain,
  grid& gridA
){

  // ---------------------------------
  MPI_Barrier(comm_world); 
  int itag = 100;   
  size_t nynz = gridA.nz *gridA.ny;
  const int master = 0;
  // ---------------------------------




  // #ifdef MPI_NONBLOCKING_SR
  // MPI_Request requests[rank_size];
  // #endif

  MPI_Status status[rank_size];

  for (size_t i_rank = 1; i_rank < rank_size ; ++ i_rank){

    itag += 10;
    
    int shift = nynz * (localDomain.i_begin_table.at(i_rank)) ;
    int count = nynz * (localDomain.i_length_table.at(i_rank));

    if (i_rank == rank_size-1) count += nynz * gridA.gC;

    // #ifdef MPI_NONBLOCKING_SR
      // if (rank != master){
      //   if (rank == i_rank) MPI_Isend(&v[shift], count, MPI_DOUBLE, master, itag, comm_world, requests+i_rank*2+0);
      // }
      // else{
      //   MPI_Irecv(&v[shift], count, MPI_DOUBLE, i_rank, itag, comm_world, requests+i_rank*2+1);
      // }
    // #endif

    // #ifdef MPI_BLOCKING_SR
      if (rank != master){

        if (rank == i_rank) MPI_Send((void *)&v[shift], count, MPI_DOUBLE, master, itag, comm_world);

      }
      else{

        MPI_Recv((void *)&v[shift], count, MPI_DOUBLE, i_rank, itag, comm_world, status+i_rank);

      }
    // #endif

  }

  MPI_Barrier(comm_world);

  // #ifdef MPI_NONBLOCKING_SR
    // MPI_Waitall(rank_size*2, requests, status);
  // #endif

}




void mpi_Bcast(
  MPI_Comm &comm_world,
  grid& gridA,
  calDomain& localDomain,
  std::vector<double> &v
){

  int nx, ny, nz, shift{0};

  // ---------------------
  if (v.size() == gridA.iceltot) { 
    std::tie(nx, ny, nz) = gridA.nxyz;
    shift = 0;
  }

  else if (v.size() == gridA.iceltotCal) { 
    std::tie(nx, ny, nz) = gridA.nxyzCal; 
    shift = gridA.gC;
  }
  else{ 
    cout << "\nv.size()" << v.size() ;
    throw std::invalid_argument("checking V!!!"); 
  }
  // ---------------------
  MPI_Barrier(comm_world);

  int nynz = ny*nz;

  for(size_t i = 0; i < localDomain.i_begin_table.size() ; ++i){
      MPI_Bcast(  (void *)&v[(localDomain.i_begin_table[i]-shift)* nynz], 
                    localDomain.i_length_table[i] * nynz, 
                    MPI_DOUBLE, i, comm_world);

      MPI_Barrier(comm_world);
  }

  MPI_Barrier(comm_world);
}




void mpi_debug_file(
  calDomain& localDomain,
  simuClass& simu,
  int vec_size,
  int mpi_word_size,
  MPI_Comm &mx_comm_world,
  std::vector<int> &mpi_neighborhood,
  grid& gridA
)
{
  std::vector<double> deVec1(vec_size);

  std::fill(deVec1.begin(), deVec1.end(),simu.PID);

  for (size_t i = 0; i < gridA.nx ; ++i)
  for (size_t j = 0; j < gridA.ny ; ++j)
  for (size_t k = 0; k < gridA.nz ; ++k)
  {
    const int icel = i*gridA.nz*gridA.ny + j*gridA.nz + k;
    deVec1.at(icel) *= 10000;
    deVec1.at(icel) += icel;
  } 

  auto deVec2 = deVec1;
  mpi_Bcast(mx_comm_world,gridA, localDomain,  deVec1);

  // mpi_iSR_double_x(2, mx_comm_world, mpi_neighborhood, deVec2, localDomain, gridA );

  std::cout <<  std::flush ;
  MPI_Barrier(mx_comm_world);


	ofstream file;
  std::string filename = "Information/MPI_file";
	filename += std::to_string(simu.PID);
	filename += ".dat";	

	file.open(filename);

	for (size_t i = 0; i < gridA.nx; ++i)
	for (size_t j = 0; j < gridA.ny; ++j) 
	for (size_t k = 0; k < gridA.nz; ++k) 
	{
    const int icel = i*gridA.nz*gridA.ny + j*gridA.nz + k;
		file <<"(" <<  icel << ")\ti\t" << i << "\tj\t" << j << "\tk\t" << k 
				<< "\t" <<  deVec1.at(icel) << "\t" << deVec2.at(icel)<< "\t" << deVec1.at(icel) -  deVec2.at(icel) <<endl;

		// file.write((char *)(&T3.u.at(icel)), sizeof(double));
	}
	file.close();

      // for (size_t i = 0; i < gridA.nx ; ++i)
      // for (size_t j = 0; j < gridA.ny ; ++j)
      // for (size_t k = 0; k < gridA.nz ; ++k)
      // {
      //   const int icel = i*gridA.nz*gridA.ny + j*gridA.nz + k;
      //   // if (abs(deVec1.at(icel) - deVec2.at(icel)) > 1.0e-4){
      //   cout << std::flush;
      //   MPI_Barrier(mx_comm_world);
      //   if (simu.PID == 0 ){
      //   cout << endl;
      //   cout << std::flush  << i << ", " 
      //                       << j << ", " 
      //                       << k << "::";
      //   }

      //   for(size_t ii = 0; ii < mpi_word_size ; ii++ ){
      //     cout << std::flush ;
      //     MPI_Barrier(mx_comm_world);
      //     if (simu.PID == ii ){
      //       cout <<  " | " << deVec1.at(icel) << " | " << deVec2.at(icel);
      //     } 
      //   }
      // }

}





