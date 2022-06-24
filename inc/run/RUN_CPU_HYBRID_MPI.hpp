# include "controlPanel.hpp"

# include "run/general.hpp"

# include <mpi.h>

# include "run/RUN_CPU.hpp"
# include "mpi_tool/init_mpi.hpp"
# include "dfib/virtualForceIntergrator_mpi.hpp"
# include "mpi_tool/pointTOpoint/blocking_p2p.hpp"
# include "mpi_tool/pointTOpoint/nonBlocking_p2p.hpp"
# include "pressure/solver/SorPipeLine.mpi.hpp"
# include "solver/mpi/bicgstabRe2_mpi.hpp"

// #define MPI_DEBUG

bool Run_Hy_MPI_OpenMP(
    grid &gA,
    clockstruct &timer,
    simuClass &simu,
    int argc, char **argv
){

  timer.beginNew.start();

  //  * struct ---------------------------
  DfibArray Dfib;

  pressure t1;

  velocity T0 , T1, T3;

  SORcoefficient Sor;

  shareMenory ShareM;

  MxClass Mx;


  //  * struct ---------------------------
  DfibArray Dfib_loc;

  pressure t1_loc;

  velocity T0_loc , T1_loc, T3_loc;

  SORcoefficient Sor_loc;

  MxClass Mx_loc;
  //  * struct ---------------------------

  //  * struct ---------------------------

  printf("pid: %d\n", getpid());

  // ~~~~~~~~~~~~~~~~~~~~~~~~ openMP ~~~~~~~~~~~~~~~~~~~~~~~~
  #ifdef _OPENMP
    std::cout << "OpenMP: " << _OPENMP << std::endl;
  #else
    std::cout << "Your compiler does not support OpenMP." << std::endl;
  #endif
  printf("pid: %d\n", getpid());

  int ompThreads = omp_get_max_threads() ;  // Using Inv OMP_NUM_THREADS
  omp_set_num_threads(ompThreads);
  // ~~~~~~~~~~~~~~~~~~~~~~~~ openMP ~~~~~~~~~~~~~~~~~~~~~~~~

  // ! ============================  divid Domain ============================


  calDomain globalDomain;
  calDomain localDomain;

  std::vector<int> grid_size{gA.nx, gA.ny, gA.nz};

  std::vector<int> dims{1, 1, 1};
  std::vector<int> dims_loc{3};

  globalDomain.initTable(grid_size, dims);

  globalDomain.initLocal(0, 0, 0);

  bool reorder = true;

  auto [mx_comm_world, 
        mpi_word_size, 
        mpi_word_rank, 
        mpi_coord, 
        mpi_neighborhood]= mpi_init(reorder, dims_loc, argc, argv);


  for (size_t i = 0 ; i < (4 - dims_loc.size()) ; ++i){
    dims_loc.push_back(1);
    mpi_coord.push_back(0);
  }

  localDomain.initTable(grid_size, dims_loc);
  localDomain.initLocal(mpi_coord.at(0), mpi_coord.at(1), mpi_coord.at(2));

  // * ========================================================================

  // check calDomain
  simu.PID = mpi_word_rank;


  if (simu.PID == 0)
  {
    if(opendir("Information") == NULL)
    { if (system("mkdir Information") != 0){ return 1; } }

    if(opendir("mx_out") == NULL)
    { if (system("mkdir mx_out") != 0){ return 1; } }

    if(opendir("Information/Chronograph") == NULL)
    { if (system("mkdir Information/Chronograph") != 0){ return 1; } }
  }


  for(size_t i = 0; i < mpi_word_size ; ++i){
    MPI_Barrier(mx_comm_world);
    if (mpi_word_rank == i){
      cout  << "[Rank: {begin,endof}] = "<< mpi_word_rank 
            << ": [{ "  << localDomain.i_begin
            << ", "     << localDomain.i_endof

            << "},{ "   << localDomain.j_begin
            << ", "     << localDomain.j_endof

            << "},{ "   << localDomain.k_begin
            << ", "     << localDomain.k_endof

            << "}]"     << endl;
    }
  }




  MPI_Barrier(mx_comm_world);
  
  if (simu.PID == 0) cout << "\n<i_begin_table>,<i_endof_table>";
  
  for(size_t j = 0; j < mpi_word_size ; ++j){
    cout << std::flush;
    MPI_Barrier(mx_comm_world);

    if (mpi_word_rank == j){
      cout << "\n===== " << j << " =====" << endl;

      for(size_t i = 0; i < mpi_word_size ; ++i){
      cout << "{" << localDomain.i_begin_table.at(i) << ", ";
      cout << localDomain.i_endof_table.at(i) << "}, ";
    }

    }
  }


  cout << std::flush;

  MPI_Barrier(mx_comm_world);

  if (simu.PID == 0) cout << "\n<i_length_table>";

  for(size_t j = 0; j < mpi_word_size ; ++j){
    MPI_Barrier(mx_comm_world);
    if (mpi_word_rank == j){
      cout << "\n===== " << j << " =====" << endl;

      for(size_t i = 0; i < mpi_word_size ; ++i){
        cout << localDomain.i_length_table.at(i) << ", ";
      }
    }
  }

  // for (auto &x:T0.U(0)){
  //   x = simu.PID;
  // }


  // for(size_t j = 0; j < mpi_word_size ; ++j)
  // {
  //   MPI_Barrier(mx_comm_world);
  //   if (mpi_word_rank == j){
  //     cout << "\n===== " << j << " =====" << endl;
  //      T0.showVelt(0, gA);
  //   }
  // }

  // mpi_Bcast(mx_comm_world, gA, localDomain, T0.U(0));
  // mpi_iSR_double_x(2, mx_comm_world, mpi_neighborhood, T0.U(0), localDomain, gA );
  // mpi_iSR_double_x_Collect_to_Master(mx_comm_world, simu.PID, mpi_word_size, T0.U(0), localDomain, gA );

  // for(size_t j = 0; j < mpi_word_size ; ++j)
  // {
  //   MPI_Barrier(mx_comm_world);
  //   if (mpi_word_rank == j){
  //     cout << "\n===== " << j << " =====" << endl;
  //      T0.showVelt(0, gA);
  //   }
  // }


  // return 0;

  // check calDomain

  // *!##  Init the variables (First policy or data Locality) and Generate grids -------------




  // ! ## Init the variables (First policy or data Locality) and Generate grids 

  resize_variable(gA, t1, T0, T1, T3, Dfib); 
  resize_variable(gA, t1_loc, T0_loc, T1_loc, T3_loc, Dfib_loc); 
  {
    double ui = 0.0, vi = 0.0, wi = 0.0;

    // T0.iniU_omp(ui, vi, wi);
    // T1.iniU_omp(ui, vi, wi);
    // T3.iniU_omp(ui, vi, wi);
    // t1.init_p(0.0);

    T0_loc.iniU(ui, vi, wi);
    T1_loc.iniU(ui, vi, wi);
    T3_loc.iniU(ui, vi, wi);
    t1_loc.init_p(0.0);
  }

  generateGride(simu, ShareM, globalDomain, gA);

  if (simu.PID == 0)
  {
    gA.io_csv("Information/gA.csv");
    OutputPlot3D_Xfile(simu, gA);
  } 

  // ! ## Read the exist data .-------------
  // {
  // auto [Nblock_read, mach_read, alpha_read, reyn_read, time_read, 
  //         pPre, ucPre, vcPre, wcPre, EtaPre] = QfileRead(gA, "mx_in/in.q");

  // t1.p = pPre;
  // T1.u = T3.u = T0.u = ucPre;
  // T1.v = T3.v = T0.v = vcPre;
  // T1.w = T3.w = T0.w = wcPre;
  // PotentialFlow(simu, Dfib, T0, t1, globalDomain, gA);

  // simu.set_Re(reyn_read);

  // simu.set_time(time_read);
  // }
  // * --------------------------------------

  // ! ## The first BC -------------
  // T1.u = T3.u = T0.u ;
  // T1.v = T3.v = T0.v ;
  // T1.w = T3.w = T0.w ;


  T1_loc.u = T3_loc.u = T0_loc.u ;
  T1_loc.v = T3_loc.v = T0_loc.v ;
  T1_loc.w = T3_loc.w = T0_loc.w ;


  // BC_staggered_main( globalDomain, T0, t1, gA );
  // BC_staggered_copy( globalDomain, T0, T1, gA );

  BC_staggered_main( localDomain, T0_loc, t1_loc, gA );
  BC_staggered_copy( localDomain, T0_loc, T1_loc, gA );

  // !## prepare for poisson equation --------------------------
    timer.pressure.start(); 


  #if  defined (P_SOLVER_ELL) \
    || defined (P_SOLVER_CSR) \
    || defined (P_SOLVER_SPE) \
    || defined (P_SOLVER_AMGCL_BUILTIN) \

    // Mx.X_result.resize(gA.iceltotCal, 0.0);
    Mx_loc.X_result.resize(gA.iceltotCal, 0.0);
  #endif



  // !### ---------------------------------------- ELL (BICG)
  #if defined (P_SOLVER_ELL)

  // ---------------------------
  // createPressureMatrix(Mx, simu, localDomain, gA);
  // ---------------------------

  // ---------------------------
  // solver::bicgstabRe2<ELL_type> pSolver(Mx.matA_ell);
  // ---------------------------


  // ---------------------------
  createPressureMatrix(Mx_loc, simu, globalDomain, gA);
  // ---------------------------
  // ---------------------------
  auto [beg, end] = localDomain.get_begEnd(gA.nyCal * gA.nzCal);
  Mx_loc.matA_ell.mpi_init(mx_comm_world, gA.iceltotCal, beg, end);
  // ---------------------------

  // ---------------------------
  solver::bicgstabRe2_mpi<ELL_type> pSolver_loc(Mx_loc.matA_ell);
  // ---------------------------
  // !### ---------------------------------------- CSR
  #elif defined (P_SOLVER_CSR)

  //  ----------------------------------------------------------------
  // Mx.matA_csr.set( gA.createPressureMatrixCSR() );
  //  ----------------------------------------------------------------

  //  ----------------------------------------------------------------
  // solver::bicgstabRe2<CSR_type> pSolver(Mx.matA_csr);
  //  ----------------------------------------------------------------


  //  ----------------------------------------------------------------
  Mx_loc.matA_csr.set( gA.createPressureMatrixCSR() );
  //  ----------------------------------------------------------------
  //  ----------------------------------------------------------------
  {
  auto [beg, end] = localDomain.get_begEnd(gA.nyCal * gA.nzCal);
  Mx_loc.matA_csr.mpi_init(mx_comm_world, gA.iceltotCal, beg,  end);
  }
  //  ----------------------------------------------------------------
  //  ----------------------------------------------------------------
  solver::bicgstabRe2_mpi<CSR_type> pSolver_loc(Mx_loc.matA_csr);
  //  ----------------------------------------------------------------


  // !### ---------------------------------------- SPE
  #elif defined(P_SOLVER_SPE)

  // --------------------------------------
  // Mx.matA_spe.resize(gA.nxCal, gA.nyCal, gA.nzCal);
  // --------------------------------------

  // --------------------------------------
  // Mx.matA_spe.setupPressure( gA.gC, gA.Dx, gA.Dy, gA.Dz,gA.Dxs, gA.Dys, gA.Dzs);
  // --------------------------------------

  // #if defined(JACOBI_PC)
  // gA.jp =  Mx.matA_spe.set_Jp_pc(); 
  // #endif


  // --------------------------------------
  // solver::bicgstabRe2<SPE_type> pSolver(Mx.matA_spe);
  // --------------------------------------


  // --------------------------------------
  Mx_loc.matA_spe.resize(gA.nxCal, gA.nyCal, gA.nzCal);
  // --------------------------------------

  // --------------------------------------
  {
  auto [beg, end] = localDomain.get_begEnd(gA.nyCal * gA.nzCal);
  Mx_loc.matA_spe.mpi_init(mx_comm_world, gA.iceltotCal, beg, end);
  }
  // --------------------------------------

  // --------------------------------------
  #if defined(JACOBI_PC)
  gA.jp =  Mx_loc.matA_spe.set_Jp_pc(); 
  #endif
  // --------------------------------------

  //  ----------------------------------------------------------------
  solver::bicgstabRe2_mpi<SPE_type> pSolver_loc(Mx_loc.matA_spe);
  //  ----------------------------------------------------------------
  #endif 


  #ifndef JACOBI_PC
  // gA.jp.clear();
  #endif


  #if  defined (P_SOLVER_ELL) \
    || defined (P_SOLVER_CSR) \
    || defined (P_SOLVER_SPE)

    // pSolver.setTolerance(simu.p_criteria);
    pSolver_loc.setTolerance(simu.p_criteria);

  #endif  // P_SOLVER_ELL  P_SOLVER_CSR  P_SOLVER_SPE



  timer.pressure.stop();
  // * ----------------------------------------------------

  // std::fill(Dfib_loc.eta.begin(), Dfib_loc.eta.end(), 0.0);
  // if (simu.DfibMethod == "OFF"){}
  // else if (simu.DfibMethod == "DFIB_Cylinder-X")
  //   DFIB_CylinderX(Dfib_loc, globalDomain, gA);
  // else if (simu.DfibMethod == "DFIB_Cylinder-Z")
  //   DFIB_CylinderZ(Dfib_loc, globalDomain, gA);
  // else
  //   throw std::invalid_argument("DFIB method??");


  // * ------------------------------------------------------

    // ---------------------------------
    OutputPLOT3D_Qfile(Dfib, simu, t1, T0, gA); 
    // ---------------------------------

    // ---------------------------------
    #ifdef TERBULENCE_SMAGORINSKY
    // T0.Viseff.resize(gA.iceltotCal, 0.0);
    T0_loc.Viseff.resize(gA.iceltotCal, 0.0);
    T1_loc.Viseff.resize(gA.iceltotCal);
    #endif
    // ---------------------------------


  auto perLoopWTime = timer.beginNew.elapsedTime();
  // ! ==================  main loop begin ==================
  for (simu.loop = 1; simu.get_finishloop(); simu.finishloop())
  {

    // ! ~~~~~~~~~~~~~~~~~~~~~~~~~~  send & recv ~~~~~~~~~~~~~~~~~~~~~~~~~~
    for(int whichDim = 0; whichDim < 3 ; whichDim++){
      mpi_iSR_double_x(2, mx_comm_world, mpi_neighborhood, T0_loc.U(whichDim), localDomain, gA );
    }
    // * ~~~~~~~~~~~~~~~~~~~~~~~~~~  send & recv ~~~~~~~~~~~~~~~~~~~~~~~~~~

    // !## 1. get the T1/u* --------------------
    timer.convectionDifussion.start();
    // ---------------------------------
    #if defined (TERBULENCE_SMAGORINSKY)
    // SmagorinskyModel(ShareM, simu, T0, t1, globalDomain, gA);
    SmagorinskyModel(ShareM, simu, T0_loc, t1_loc, localDomain, gA);
    #endif
    // ---------------------------------

    // ---------------------------------
    // ConvectionDifussion(simu, T0, T1, globalDomain, gA);
    ConvectionDifussion(simu, T0_loc, T1_loc, localDomain, gA);
    // ---------------------------------

    // --------------------------------- the error could be happen if you call BC_updateSlid().
    // BC_updateSlid(globalDomain, T1, gA);
    // BC_updateSlid(localDomain, T1_loc, gA);
    // ---------------------------------
    timer.convectionDifussion.stop();
    //*---------------------------------------------------



    // ! ~~~~~~~~~~~~~~~~~~~~~~~~~~  send & recv ~~~~~~~~~~~~~~~~~~~~~~~~~~
    for(int whichDim = 0; whichDim < 3 ; whichDim++){
      mpi_iSR_double_x(1, mx_comm_world, mpi_neighborhood, T1_loc.U(whichDim), localDomain, gA );
    }
    // * ~~~~~~~~~~~~~~~~~~~~~~~~~~  send & recv ~~~~~~~~~~~~~~~~~~~~~~~~~~
    // -----------------------
    // get_Cfl(T1, globalDomain, gA, simu);
    // -----------------------

    // -----------------------
    // get_Cfl(T1_loc, localDomain, gA, simu);
    // -----------------------


    // -----------------------
    // simu.timestepper();  // todo : debug
    // -----------------------


    // !## 2. Get the pressure poisson equation's result.
    timer.pressure.start(); 


    // !### ---------------------------------------- SOR method 

    #if defined (P_SOLVER_SOR)
    if (simu.loop > pLoopIni) {
    // ---------------------------------
    // SorPipeLine_omp(Sor, simu, T1, t1, globalDomain, gA);
    SorPipeLine_mpi_hybrid_omp(mx_comm_world, mpi_neighborhood, Sor, simu, T1_loc, t1_loc, localDomain, gA);
    // ---------------------------------
    }

    // !### ---------------------------------------- CSR/SPE/ELL (BICG)(CG) / AMGCL
    #elif defined (P_SOLVER_BICG_CSR) \
      ||  defined (P_SOLVER_BICG_SPE) \
      ||  defined (P_SOLVER_BICG_ELL) \
      ||  defined (P_SOLVER_AMGCL_BUILTIN)

    int pLoopIni  = 0;

    if (simu.loop > pLoopIni)
    {
      // ---------------------------------
      // createBMatrix_omp( T1, Mx, simu, globalDomain, gA);
      // ---------------------------------

      // ---------------------------------
      // std::tie(simu.iters, simu.error) = pSolver(Mx.matB, Mx.X_result);
      // ---------------------------------

      // ---------------------------------
      // Pressure_transform_X_result(t1, Mx, globalDomain, gA);
      // ---------------------------------

      // ---------------------------------
      createBMatrix_omp(T1_loc, Mx_loc, simu, localDomain, gA);
      mpi_Bcast(mx_comm_world, gA, localDomain, Mx_loc.matB);
      // ---------------------------------

      // ---------------------------------
      std::tie(simu.iters, simu.error) = pSolver_loc(Mx_loc.matB, Mx_loc.X_result);
      // ---------------------------------

      // ---------------------------------
      Pressure_transform_X_result(t1_loc, Mx_loc, localDomain, gA);
      // ---------------------------------
    }

    // ------------------------------------------------------------------
    auto [min_p , max_p] = getMax(Mx_loc.X_result, localDomain, gA);
    {
      auto min_pTEMP = min_p;
      auto max_pTEMP = max_p;
      MPI_Allreduce(&min_pTEMP, &min_p, 1, MPI_DOUBLE, MPI_MIN, mx_comm_world);
      MPI_Allreduce(&max_pTEMP, &max_p, 1, MPI_DOUBLE, MPI_MAX, mx_comm_world);
    }
    if (simu.PID == 0) std::cout  << "[Min_p, Max_p]:"<< min_p  << ", "<< max_p << std::endl;
    // ------------------------------------------------------------------

    #endif 
    
    if (simu.PID == 0) simu.printInfo();
    timer.pressure.stop(); 
    // * -----------------------------------------------------------


    // ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ send & recv ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      mpi_iSR_double_x(1, mx_comm_world, mpi_neighborhood, t1_loc.p, localDomain, gA );
    // * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ send & recv ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // !## 3. Implement the DFIB method, the velocity(T1) be update to velocity(T3).
    timer.updateT1toT3.start();
    // ----------------------
    // update_UandF_omp(Dfib, simu, t1, T1, T3, globalDomain, gA); 
    // ----------------------

    // ----------------------
    update_UandF_omp(Dfib, simu, t1_loc, T1_loc, T3_loc, localDomain, gA); 
    // ----------------------
    // BC_updateSlid(globalDomain, T3, gA); //Ahmad didn't update it !!
    timer.updateT1toT3.stop();
    // * -----------------------------------------------------------


    // ! Check the staedy state (L2 norm) -------------------------------
    timer.checkL2norm.start();
    // CheckSteadyStatebyMaxVal_omp(simu, ShareM, T0, T3, globalDomain, gA);
    // CheckSteadyStatebyMaxVal_mpi(mx_comm_world, simu, ShareM, T0, T3, localDomain, gA);
    // CheckSteadyStatebyL2norm(simu, T0_loc, T3_loc, localDomain, gA);
    timer.checkL2norm.stop();
    // *  ----------------------------------------------------------------

    // ! ## 4. copy the vel.
    timer.updateT3toT0.start();
    // T0.copy_omp(T3);
    T0_loc.copy_omp(T3_loc);
    timer.updateT3toT0.stop();
    // * --------------------


    // !## 5. Update vel and pressure at ghost cell.
    timer.BC.start();
    // BC_staggered_main( globalDomain, T0, t1, gA );
    // BC_staggered_copy( globalDomain, T0, T1, gA );
    BC_staggered_main( localDomain, T0_loc, t1_loc, gA);
    BC_staggered_copy( localDomain, T0_loc, T1_loc, gA );

    timer.BC.stop();
    // * -----------------------------------------

    // !## 6 Get both Cd and Cl.  -------------
    if (simu.DfibMethod != "OFF")
    {
      timer.getCdCl.start();
      // getCD_CL(simu, globalDomain, Dfib, gA);
      getCD_CL_mpi(mx_comm_world, simu, localDomain, Dfib_loc, gA);
      timer.getCdCl.stop();
    }
    // * --------------------------------



    // !## 7 write plot3D formart .-------------
    timer.beginNew.stop();
    if (simu.get_writefile()){

      // ! ~~~~~~~~~~~~~~~~~~~~~~~~~~  send & recv ~~~~~~~~~~~~~~~~~~~~~~~~~~
      for(int whichDim = 0; whichDim < 3 ; whichDim++){
        mpi_iSR_double_x_Collect_to_Master(mx_comm_world, simu.PID, mpi_word_size, T3_loc.U(whichDim), localDomain, gA );
      }
      mpi_iSR_double_x_Collect_to_Master(mx_comm_world, simu.PID, mpi_word_size, t1_loc.p, localDomain, gA );

      MPI_Barrier(mx_comm_world);

      // * ~~~~~~~~~~~~~~~~~~~~~~~~~~  send & recv ~~~~~~~~~~~~~~~~~~~~~~~~~~

      if (simu.PID == 0){
        // OutputPLOT3D_Qfile(Dfib, simu, t1, T3, gA);
        OutputPLOT3D_Qfile(Dfib, simu, t1_loc, T3_loc, gA);
      }

    }
    // * -----------------------------------------


    // !## 8. IO ------------------------------------------------
    MPI_Barrier(mx_comm_world);
    if (simu.PID == 0) 
    cout 
    <<     "[simu]{time, dt} ="  << simu.get_time()          
    <<     ", " << simu.dt       << endl
    <<     "[wall time] < "      << timer.beginNew.elapsedTime() - perLoopWTime 
    <<     " OF "                << timer.beginNew.elapsedTime()
    <<     " >" << endl

    // ! =============== NEXT time step ===============
    <<     "===================================================" << endl
    <<     "LOOP :"  << simu.loop + 1
    <<     gA.show() << ", file :" << simu.get_file() << endl;
    // * ----------------------------------------------------------
    perLoopWTime = timer.beginNew.elapsedTime();
    timer.beginNew.start();

    recordingTime(timer, simu, ShareM);
  }
  timer.beginNew.stop();
  // ! ==================  main loop begin ==================


  MPI_Barrier(mx_comm_world);
  
  MPI_Finalize();  
  
  return true;
}

