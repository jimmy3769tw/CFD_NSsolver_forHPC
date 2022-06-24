#pragma once

#include "controlPanel.hpp"

#include "run/general.hpp"

#include <algorithm>

bool runSeqOmp_predictionMethod(
    grid &gA,
    clockstruct &timer,
    simuClass &simu,
    int argc, char **argv
){

  #if  defined (PC_SEQ) || defined (PC_OMP)
  timer.beginNew.start();

  if(opendir("Information") == NULL)
  { if (system("mkdir Information") != 0){ return 1; } }

	if(opendir("mx_out") == NULL)
	{ if (system("mkdir mx_out") != 0){ return 1; } }

	if(opendir("Information/Chronograph") == NULL)
	{ if (system("mkdir Information/Chronograph") != 0){ return 1; } }


  //  * struct ---------------------------
  DfibArray Dfib;

  pressure t1;

  velocity T0 , T1, T3;

  SORcoefficient Sor;

  shareMenory ShareM;

  MxClass Mx;
  //  * struct ---------------------------

  printf("pid: %d\n", getpid());

#if defined (PC_OMP)
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
#endif

  // ! ============================  divid Domain ============================

  calDomain globalDomain;

  std::vector<int> grid_size{gA.nx, gA.ny, gA.nz};

  std::vector<int> dims{1, 1, 1};

  globalDomain.initTable(grid_size, dims);

  globalDomain.initLocal(0, 0, 0);

  // * ========================================================================

  // ! ## Init the variables (First policy or data Locality) and Generate grids 

  resize_variable(gA, t1, T0, T1, T3, Dfib); // ! resize shared memory

  double ui = 1.0, vi = 0.0, wi = 0.0;

  #if defined (PC_SEQ)
    T0.iniU(ui, vi, wi);
    T1.iniU(ui, vi, wi);
    T3.iniU(ui, vi, wi);
  #elif defined (PC_OMP)
    T0.iniU_omp(ui, vi, wi);
    T1.iniU_omp(ui, vi, wi);
    T3.iniU_omp(ui, vi, wi);
  #endif

  t1.init_p(0.0);

  generateGride(simu, ShareM, globalDomain, gA);

  gA.io_csv("Information/gA.csv");

  OutputPlot3D_Xfile(simu, gA);


  // ! ## Read the exist data. -------------
  // auto [Nblock_read, mach_read, alpha_read, reyn_read, time_read, 
  //         pPre, ucPre, vcPre, wcPre, EtaPre] = QfileRead(gA, "mx_in/in.q");

  // t1.p = pPre;
  // T1.u = T3.u = T0.u = ucPre;
  // T1.v = T3.v = T0.v = vcPre;
  // T1.w = T3.w = T0.w = wcPre;
  // PotentialFlow(simu, Dfib, T0, t1, globalDomain, gA);

  // simu.set_Re(reyn_read);

  // simu.set_time(time_read);
  // * --------------------------------------

  // ! ## The first BC -------------
  T1.u = T3.u = T0.u ;
  T1.v = T3.v = T0.v ;
  T1.w = T3.w = T0.w ;

  BC_staggered_main( globalDomain, T0, t1, gA );
  BC_staggered_copy( globalDomain, T0, T1, gA );

  // !## prepare for poisson equation --------------------------
    timer.pressure.start();


  #if  defined (P_SOLVER_ELL) \
    || defined (P_SOLVER_CSR) \
    || defined (P_SOLVER_SPE) \
    || defined (P_SOLVER_AMGCL_BUILTIN) \

    Mx.X_result.resize(gA.iceltotCal, 0.0);

  #endif



  // !### ---------------------------------------- ELL (BICG)
  #if defined (P_SOLVER_ELL)
  // ---------------------------

    timer.set_MaA.start();
    createPressureMatrix(Mx, simu, globalDomain, gA);
    timer.set_MaA.stop();
    std::cout << "ELL_setUptime : "<< timer.set_MaA.elapsedTime() << std::endl;
     Mx.matA_ell.analyse();
  // ---------------------------

  // ---------------------------
    solver::bicgstabRe2<ELL_type> pSolver(Mx.matA_ell);
  // ---------------------------

  // ---------------------------
  // !### ---------------------------------------- CSR
  #elif defined (P_SOLVER_CSR)

  //  ----------------------------------------------------------------
  Mx.matA_csr.set( gA.createPressureMatrixCSR() );
  //  ----------------------------------------------------------------

  //  ----------------------------------------------------------------
  solver::bicgstabRe2<CSR_type> pSolver(Mx.matA_csr);
  //  ----------------------------------------------------------------

  // !### ---------------------------------------- SPE
  #elif defined(P_SOLVER_SPE)

  // --------------------------------------
  Mx.matA_spe.resize(gA.nxCal, gA.nyCal, gA.nzCal);
  // --------------------------------------

  // --------------------------------------
  Mx.matA_spe.setupPressure( gA.gC, gA.Dx, gA.Dy, gA.Dz,gA.Dxs, gA.Dys, gA.Dzs);

  #if defined(JACOBI_PC)
  gA.jp =  Mx.matA_spe.set_Jp_pc(); 
  #endif // JACOBI_PC
  // --------------------------------------

  // --------------------------------------
  solver::bicgstabRe2<SPE_type> pSolver(Mx.matA_spe);
  // --------------------------------------

  // !### ---------------------------------------- EIGEN
  #elif defined (EIGEN_ON)

    // -------------------------------------------
    timer.set_MaA.start();
    Mx.x_Eigen.resize(gA.iceltotCal);
    Mx.matA_Eigen.resize(gA.iceltotCal, gA.iceltotCal);
    // -------------------------------------------

    {
        auto [prt, idx, val] = gA.createPressureMatrixCSR() ;
        std::vector<Eigen::Triplet<double>> coefficients;
        coefficients.reserve(gA.iceltotCal);

        for (int row = 0; row < prt.size()-1; ++row){
            for (int j = prt[row]; j < prt[row+1] ;++j){
            coefficients.push_back(Eigen::Triplet<double>( row, idx[j], val[j]));
            }
        }
      Mx.matA_Eigen.setFromTriplets(coefficients.begin(), coefficients.end());

      Mx.solver_Eigen.compute(Mx.matA_Eigen);
      Mx.matB_Eigen.resize(gA.iceltotCal);
      Mx.solver_Eigen.setTolerance(simu.p_criteria);
    }


    timer.set_MaA.stop();
    std::cout << "[EIGEN] setUp time : "<< timer.set_MaA.elapsedTime() << std::endl;
    // -------------------------------------------
  // !### ---------------------------------------- AMGCL(buildin)
  #elif defined (P_SOLVER_AMGCL_BUILTIN)

    // -------------------------------------------
    auto [ptr, idx, values] = gA.createPressureMatrixCSR();

    auto amgcl_mat = std::tie(gA.iceltotCal, ptr, idx, values);
    // -------------------------------------------
    SolverBuiltin::params prm;

    prm.solver.tol = simu.p_criteria;

    prm.solver.maxiter = simu.p_iterMax;

    prm.solver.ns_search = true;

    prm.solver.verbose = true;

    SolverBuiltin pSolver(amgcl_mat, prm);
    // ---------------------------------
    // SolverBuiltin pSolver(amgcl_mat);
    // -------------------------------------------
    // amgcl::io::mm_write("test_io_vec.mm", rhs.data(), n, 1);

    amgcl::io::mm_write("Information/This_Matrix.mtx", amgcl_mat);
    // -------------------------------------------

  #endif  // P_SOLVER

  #ifndef JACOBI_PC
  // gA.jp.clear();
  #endif

  #if  defined (P_SOLVER_ELL) \
    || defined (P_SOLVER_CSR) \
    || defined (P_SOLVER_SPE)

    pSolver.setTolerance(simu.p_criteria);

  #endif


  timer.pressure.stop();
  // * ----------------------------------------------------

  // ! ##  DfibGetEta(Dfib, globalDomain, gA) -------------
  std::fill(Dfib.eta.begin(), Dfib.eta.end(), 0.0);

  if (simu.DfibMethod == "OFF"){}
  else if (simu.DfibMethod == "DFIB_Cylinder-X")
    DFIB_CylinderX(Dfib, globalDomain, gA);
  else if (simu.DfibMethod == "DFIB_Cylinder-Z")
    DFIB_CylinderZ(Dfib, globalDomain, gA);
  else
    throw std::invalid_argument("DFIB method??");

  // * ------------------------------------------------------

    // ---------------------------------
    OutputPLOT3D_Qfile(Dfib, simu, t1, T0, gA); ///Init
    // ---------------------------------

    // ---------------------------------
    #ifdef TERBULENCE_SMAGORINSKY
    T0.Viseff.resize(gA.iceltotCal, 0.0);
    #endif
    // ---------------------------------
    auto source_termF = Dfib.f;
  // !  ------------- remove velocity (eta == 0) -------------
  // removeVelocity(Dfib, T0, globalDomain, gA);
  // removeVelocity(Dfib, T1, globalDomain, gA);
  // removeVelocity(Dfib, T3, globalDomain, gA);
    // * -------------------------------------------------------------
  auto perLoopWTime = timer.beginNew.elapsedTime();
  // ! ==================  main loop begin ==================
  for (simu.loop = 1; simu.get_finishloop(); simu.finishloop())
  {
    // ==============================================================================
    std::fill(source_termF.begin(), source_termF.end(), double());
    for (int inner_loop = 0 ; inner_loop < 3; ++inner_loop)
    {

      // !## PC method 
      // !## 1. get the T1/u* (PredictionMethod) --------------------
      timer.convectionDifussion.start();
      // ---------------------------------
      #if defined (TERBULENCE_SMAGORINSKY)
      SmagorinskyModel(ShareM, simu, T0, t1, globalDomain, gA);
      #endif
      // ---------------------------------

      // ---------------------------------
      ConvectionDifussion(simu, T0, T1, globalDomain, gA);
      // BC_updateSlid(globalDomain, T1, gA);  // If we call the BC_updateSlid() the error could be happen.
      // ---------------------------------


      for(int whichDim = 0; whichDim < 3 ; whichDim++){
        for (auto i = globalDomain.i_begin ; i < globalDomain.i_endof ; ++i )
        for (auto j = globalDomain.j_begin ; j < globalDomain.j_endof ; ++j )
        for (auto k = globalDomain.k_begin ; k < globalDomain.k_endof ; ++k )
        {
          T1.U(whichDim)[gA.icelCal(i,j,k)] += source_termF[gA.icelCal(i,j,k) + gA.iceltotCal*whichDim];
        }
      }

      timer.convectionDifussion.stop();
      //*---------------------------------------------------


      // !## 2. Get the pressure, form solving poisson equation 
      timer.pressure.start(); 

      // !### ---------------------------------------- SOR method 
      int pLoopIni  = 1;
      if (simu.loop > pLoopIni)
      {
      #if defined (P_SOLVER_SOR)

      // ---------------------------------
      #if defined (PC_OMP)
      SorPipeLine_omp(Sor, simu, T1, t1, globalDomain, gA);
      #elif defined (PC_SEQ)
      SorPipeLine_seq( Sor, simu, T1, t1, globalDomain, gA);
      #endif
      // ---------------------------------

      // !### ---------------------------------------- CSR/SPE/ELL (BICG) 
      #elif defined (P_SOLVER_BICG_CSR) \
        ||  defined (P_SOLVER_BICG_SPE) \
        ||  defined (P_SOLVER_BICG_ELL) \
        ||  defined (P_SOLVER_AMGCL_BUILTIN)

      // ---------------------------------
      #if defined (PC_SEQ)
      createBMatrix_seq(T1, Mx, simu, globalDomain, gA);
      #elif defined (PC_OMP)
      createBMatrix_omp(T1, Mx, simu, globalDomain, gA);
      #endif   // PC_SEQ
      // ---------------------------------

      // ---------------------------------
        std::tie(simu.iters, simu.error) = pSolver(Mx.matB, Mx.X_result);
      // ---------------------------------

      // ---------------------------------
        Pressure_transform_X_result(t1, Mx, globalDomain, gA);
      // ---------------------------------

      // auto [min_p , max_p] = getMax(Mx.X_result, globalDomain, gA);
      // std::cout  << "[Min_p, Max_p]:"<< min_p  << ", "<< max_p;

      // ---------------------------------------- EIGNE 
      #elif defined (P_SOLVER_EIGEN_CSR)

        createBMatrix_Eigen(T1, Mx, simu, globalDomain, gA);
      
        Mx.x_Eigen = Mx.solver_Eigen.solve(Mx.matB_Eigen);

        simu.iters = Mx.solver_Eigen.iterations();
        simu.error = Mx.solver_Eigen.error()  ;
      
      // ---------------------------------
        Pressure_transform_x_Eigen(t1, Mx, globalDomain, gA);
      // ---------------------------------

      #endif
      }



        // !## 3. Implement the DFIB method, the velocity(T1) be update to velocity(T3).
        timer.updateT1toT3.start();
        // ----------------------
        #if defined (PC_SEQ)
        update_UandF_seq(Dfib, simu, t1, T1, T3, globalDomain, gA);
        #elif defined (PC_OMP)
        update_UandF_omp(Dfib, simu, t1, T1, T3, globalDomain, gA); 
        #endif
        // ----------------------

        for(int whichDim = 0; whichDim < 3 ; whichDim++){
          for (auto i = globalDomain.i_begin ; i < globalDomain.i_endof ; ++i )
          for (auto j = globalDomain.j_begin ; j < globalDomain.j_endof ; ++j )
          for (auto k = globalDomain.k_begin ; k < globalDomain.k_endof ; ++k )
          {
            source_termF[gA.icelCal(i,j,k) + gA.iceltotCal*whichDim]
             += simu.dt * Dfib.f[gA.icelCal(i,j,k) + gA.iceltotCal*whichDim];
          }
        }

        {
            double sum = 0;

        #pragma omp parallel for reduction(+ : sum)
        for (auto i = globalDomain.i_begin ; i < globalDomain.i_endof ; ++i )
        for (auto j = globalDomain.j_begin ; j < globalDomain.j_endof ; ++j )
        for (auto k = globalDomain.k_begin ; k < globalDomain.k_endof ; ++k )
        {
            sum += Dfib.f[ gA.iceltotCal*0 + gA.icelCal(i,j,k) ] * gA.delta(i,j,k);
        }

        std::cout << "inner_loop\t" << inner_loop << "\tcD " << -2.0 * sum / gA.lz;

        sum  = 0;
        #pragma omp parallel for reduction(+ : sum)
        for (auto i = globalDomain.i_begin ; i < globalDomain.i_endof ; ++i )
        for (auto j = globalDomain.j_begin ; j < globalDomain.j_endof ; ++j )
        for (auto k = globalDomain.k_begin ; k < globalDomain.k_endof ; ++k )
        {
            sum += Dfib.f[ gA.iceltotCal*1 + gA.icelCal(i,j,k) ] * gA.delta(i,j,k);
        }

        sum  = 0;
        #pragma omp parallel for reduction(+ : sum)
        for (auto i = globalDomain.i_begin ; i < globalDomain.i_endof ; ++i )
        for (auto j = globalDomain.j_begin ; j < globalDomain.j_endof ; ++j )
        for (auto k = globalDomain.k_begin ; k < globalDomain.k_endof ; ++k )
        {
            sum += Dfib.f[ gA.iceltotCal*1 + gA.icelCal(i,j,k) ] * gA.delta(i,j,k);
        }


        std::cout << "\tcL " << -2.0 * sum / gA.lz  << "\t iters\t" << simu.iters << std::endl ;

        }

        timer.updateT1toT3.stop();
        // * -----------------------------------------------------------
      }  
      // end of prediction
      // ==============================================================================

    // !　## (PC) 
    // !## 1. get the T1/u* --------------------
    timer.convectionDifussion.start();
    // ---------------------------------
    #if defined (TERBULENCE_SMAGORINSKY)
    SmagorinskyModel(ShareM, simu, T0, t1, globalDomain, gA);
    #endif
    // ---------------------------------

    // ---------------------------------
    ConvectionDifussion(simu, T0, T1, globalDomain, gA);
    // BC_updateSlid(globalDomain, T1, gA);  // If we call the BC_updateSlid() the error could be happen.
    // ---------------------------------
    timer.convectionDifussion.stop();

    for(int whichDim = 0; whichDim < 3 ; whichDim++){
    for (auto i = globalDomain.i_begin ; i < globalDomain.i_endof ; ++i )
    for (auto j = globalDomain.j_begin ; j < globalDomain.j_endof ; ++j )
    for (auto k = globalDomain.k_begin ; k < globalDomain.k_endof ; ++k )
    {
        T1.U(whichDim)[gA.icelCal(i,j,k)] += source_termF[gA.icelCal(i,j,k) + gA.iceltotCal*whichDim];
    }
    }

    //*---------------------------------------------------

    // -----------------------
    // get_Cfl(T1, globalDomain, gA, simu);
    // simu.timestepper();  // todo  
    // -----------------------

    // !## 2. Get the pressure, form solving poisson equation 
    timer.pressure.start(); 

    // !### ---------------------------------------- SOR method 
    int pLoopIni  = 1;
    if (simu.loop > pLoopIni)
    {
    #if defined (P_SOLVER_SOR)


    // ---------------------------------
    #if defined (PC_OMP)
    SorPipeLine_omp(Sor, simu, T1, t1, globalDomain, gA);
    #elif defined (PC_SEQ)
    SorPipeLine_seq( Sor, simu, T1, t1, globalDomain, gA);
    #endif
    // ---------------------------------

    // !### ---------------------------------------- CSR/SPE/ELL (BICG) 
    #elif defined (P_SOLVER_BICG_CSR) \
      ||  defined (P_SOLVER_BICG_SPE) \
      ||  defined (P_SOLVER_BICG_ELL) \
      ||  defined (P_SOLVER_AMGCL_BUILTIN)



    // ---------------------------------
    #if defined (PC_SEQ)
    createBMatrix_seq(T1, Mx, simu, globalDomain, gA);
    #elif defined (PC_OMP)
    createBMatrix_omp(T1, Mx, simu, globalDomain, gA);
    #endif   // PC_SEQ
    // ---------------------------------

    // ---------------------------------
      std::tie(simu.iters, simu.error) = pSolver(Mx.matB, Mx.X_result);
    // ---------------------------------

    // ---------------------------------
      Pressure_transform_X_result(t1, Mx, globalDomain, gA);
    // ---------------------------------

    // auto [min_p , max_p] = getMax(Mx.X_result, globalDomain, gA);
    // std::cout  << "[Min_p, Max_p]:"<< min_p  << ", "<< max_p;

    // ---------------------------------------- EIGNE 
    #elif defined (P_SOLVER_EIGEN_CSR)

      createBMatrix_Eigen(T1, Mx, simu, globalDomain, gA);
    
      Mx.x_Eigen = Mx.solver_Eigen.solve(Mx.matB_Eigen);

      simu.iters = Mx.solver_Eigen.iterations();
      simu.error = Mx.solver_Eigen.error()  ;
    
    // ---------------------------------
      Pressure_transform_x_Eigen(t1, Mx, globalDomain, gA);
    // ---------------------------------


    #endif
    }




    // {

    //   double Min_p{1.0e10}, Max_p{-1.0e10};
    //   for (size_t i = globalDomain.i_begin; i < globalDomain.i_endof ; ++i )
    //   for (size_t j = globalDomain.j_begin; j < globalDomain.j_endof ; ++j )
    //   for (size_t k = globalDomain.k_begin; k < globalDomain.k_endof ; ++k )
    //   {
    //       double pressure_ = t1.p[gA.icel(i,j,k)];
    //       Max_p = std::max(pressure_, Max_p);
    //       Min_p = std::min(pressure_, Min_p);
    //   }

    //   std::cout  << "[p]:"<< Min_p  << ", "<< Max_p;
    
    // }
  
    simu.printInfo();
    // * -----------------------------------------------------------
    timer.pressure.stop(); 

    // !## 3. Implement the DFIB method, the velocity(T1) be update to velocity(T3).
    timer.updateT1toT3.start();
    // ----------------------
    #if defined (PC_SEQ)
    update_UandF_seq(Dfib, simu, t1, T1, T3, globalDomain, gA);
    #elif defined (PC_OMP)
    update_UandF_omp(Dfib, simu, t1, T1, T3, globalDomain, gA); 
    #endif
    // ----------------------
    // BC_updateSlid(globalDomain, T3, gA); //Ahmad didn't update it !!
    timer.updateT1toT3.stop();
    // * -----------------------------------------------------------


    // ! Check the staedy state (L2 norm) -------------------------------
    timer.checkL2norm.start();
    #if defined (PC_SEQ)
    auto [uDif, vDif, wDif] = CheckSteadyStatebyMaxVal_seq(simu, ShareM, T0, T3, globalDomain, gA); 
    #elif defined (PC_OMP)
    auto [uDif, vDif, wDif] = CheckSteadyStatebyMaxVal_omp(simu, ShareM, T0, T3, globalDomain, gA);
    #endif
    timer.checkL2norm.stop();
    // *  ----------------------------------------------------------------

    // ! ## 4. copy the vel.
    timer.updateT3toT0.start();
    T0.copy_omp(T3);
    timer.updateT3toT0.stop();
    // * --------------------


    // !## 5. Update vel and pressure at ghost cell.
    timer.BC.start();
    
    BC_staggered_main( globalDomain, T0, t1, gA );
    BC_staggered_copy( globalDomain, T0, T1, gA );

    timer.BC.stop();
    // * -----------------------------------------

    // !## 6 Get both Cd and Cl.  -------------
    timer.getCdCl.start();
    timer.beginNew.stop();
    if (simu.DfibMethod != "OFF")
    { Delay_IO_CD_CL(simu, globalDomain, Dfib, gA);}
    timer.getCdCl.stop();
    // * --------------------------------

    // !## 7 write plot3D formart. (IO)-------------
    if (simu.get_writefile()){
      OutputPLOT3D_Qfile(Dfib, simu, t1, T3, gA);
    }
    // * -----------------------------------------


    // !## 8. IO ------------------------------------------------
    timer.preIOclock.start();

    if (simu.loop %10)

    cout << (simu.loop + 1) << "\t" <<  simu.get_time()  << "\t" << timer.beginNew.elapsedTime() - perLoopWTime << "\n";

    timer.beginNew.start(); perLoopWTime = timer.beginNew.elapsedTime();
    recordingTime(timer, simu, ShareM); 
    timer.preIOclock.stop();
  }
  timer.beginNew.stop();
  // ! ==================  main loop begin ==================


  #endif
  return true;
}
