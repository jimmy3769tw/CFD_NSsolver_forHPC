#pragma once 

#include "controlPanel.hpp"
#include "ocl/oclClass.hpp"
#include "controlPanel.hpp"

#include "run/RUN_CPU.hpp"

#include "run/RUN_GPU.hpp"

// #include "run/general.hpp"

bool RunOnGpu(
    grid &gA,
    clockstruct &timer,
    simuClass &simu,
    int argc, 
    char **argv
){
}

//   timer.beginNew.start();

//   if(opendir("../Information") == NULL)
//   { if (system("mkdir ../Information") != 0){ return 1; } }


//   //  * struct ---------------------------
//   DfibArray Dfib;

//   pressure t1;

//   velocity T0 , T1, T3;

//   SORcoefficient Sor;

//   shareMenory ShareM;

//   MxClass Mx;

//   //  * struct ---------------------------

//   printf("pid: %d\n", getpid());

// #if defined (PC_OMP)
//   // ~~~~~~~~~~~~~~~~~~~~~~~~ openMP ~~~~~~~~~~~~~~~~~~~~~~~~
//   #ifdef _OPENMP
//     std::cout << "OpenMP: " << _OPENMP << std::endl;
//   #else
//     std::cout << "Your compiler does not support OpenMP." << std::endl;
//   #endif
//   printf("pid: %d\n", getpid());

//   int ompThreads = omp_get_max_threads() ;  // Using Inv OMP_NUM_THREADS
//   omp_set_num_threads(ompThreads);
//   // ~~~~~~~~~~~~~~~~~~~~~~~~ openMP ~~~~~~~~~~~~~~~~~~~~~~~~
// #endif

//   // ! ============================  divid Domain ============================

//   calDomain Lo;

//   std::vector<int> grid_size{gA.nx, gA.ny, gA.nz};

//   std::vector<int> dims{1, 1, 1};

//   Lo.initTable(grid_size, dims);

//   Lo.initLocal(0, 0, 0);

//   // * ========================================================================

//   // ! ## Init the variables (First policy or data Locality) and Generate grids 

//   resize_variable(gA, simu, t1, T0, T1, T3, Dfib); // ! resize shared memory

//   double ui = 1.0, vi = 0.0, wi = 0.0;

//   #if defined (PC_SEQ)
//     T0.iniU(ui, vi, wi);
//     T1.iniU(ui, vi, wi);
//     T3.iniU(ui, vi, wi);
//   #elif defined (PC_OMP)
//     T0.iniU_omp(ui, vi, wi);
//     T1.iniU_omp(ui, vi, wi);
//     T3.iniU_omp(ui, vi, wi);
//   #endif

//   t1.init_p(0.0);

//   generateGride(simu, ShareM, Lo, gA);

//   gA.io_csv("Information/gA.csv");

//   OutputPlot3D_Xfile(simu, gA);


//   // ! ## Read the exist data .-------------
//   // auto [Nblock_read, mach_read, alpha_read, reyn_read, time_read, 
//   //         pPre, ucPre, vcPre, wcPre, EtaPre] = QfileRead(gA, "mx_in/in.q");

//   // t1.p = pPre;
//   // T1.u = T3.u = T0.u = ucPre;
//   // T1.v = T3.v = T0.v = vcPre;
//   // T1.w = T3.w = T0.w = wcPre;
//   // PotentialFlow(simu, Dfib, T0, t1, Lo, gA);

//   // simu.set_Re(reyn_read);

//   // simu.set_time(time_read);
//   // * --------------------------------------

//   // ! ## The first BC -------------
//   T1.u = T3.u = T0.u ;
//   T1.v = T3.v = T0.v ;
//   T1.w = T3.w = T0.w ;

//   BC_staggered_main( Lo, T0, t1, gA );
//   BC_staggered_copy( Lo, T0, T1, gA );

//   // !## prepare for poisson equation --------------------------
//     timer.p.start();


//   #if  defined (P_SOLVER_ELL) \
//     || defined (P_SOLVER_CSR) \
//     || defined (P_SOLVER_SPE) \
//     || defined (P_SOLVER_AMGCL_BUILTIN) \

//     Mx.X_result.resize(gA.iceltotCal, 0.0);

//   #endif




//   // !### ---------------------------------------- ELL (BICG)
//   #if defined (P_SOLVER_ELL)
//   // ---------------------------
//     timer.set_MaA.start();
//     createPressureMatrix(Mx, simu, Lo, gA);
//     timer.set_MaA.stop();
//     std::cout << "ELL_setUptime : "<< timer.set_MaA.elapsedTime() << std::endl;
//      Mx.matA_ell.analyse();
//   // ---------------------------

//   // ---------------------------
//     solver::bicgstab<ELL_type> pSolver(Mx.matA_ell);
//   // ---------------------------

//   // ---------------------------
//   // !### ---------------------------------------- CSR
//   #elif defined (P_SOLVER_CSR)

//   //  ----------------------------------------------------------------
//   Mx.matA_csr.set( gA.createPressureMatrixCSR() );
//   //  ----------------------------------------------------------------

//   //  ----------------------------------------------------------------
//   solver::bicgstabRe2<CSR_type> pSolver(Mx.matA_csr);
//   //  ----------------------------------------------------------------

//   // !### ---------------------------------------- SPE
//   #elif defined(P_SOLVER_SPE)

//   // --------------------------------------
//   Mx.matA_spe.resize(gA.nxCal, gA.nyCal, gA.nzCal);
//   // --------------------------------------

//   // --------------------------------------
//   Mx.matA_spe.setupPressure( gA.gC, gA.Dx, gA.Dy, gA.Dz,gA.Dxs, gA.Dys, gA.Dzs);
//   // gA.jp =  Mx.matA_spe.set_Jp_pc(); 
//   // --------------------------------------

//   // --------------------------------------
//   solver::bicgstabRe2<SPE_type> pSolver(Mx.matA_spe);
//   // --------------------------------------

//   // !### ---------------------------------------- EIGEN
//   #elif defined (EIGEN_ON)


//     // -------------------------------------------
//     timer.set_MaA.start();
//     Mx.x_Eigen.resize(gA.iceltotCal);
//     Mx.matA_Eigen.resize(gA.iceltotCal, gA.iceltotCal);
//     // -------------------------------------------

//     createPressureMatrix_Eigen(Mx, simu, Lo, gA);

//     Mx.matA_Eigen.makeCompressed();

//     timer.set_MaA.stop();

//     std::cout << "[EIGEN] setUp time : "<< timer.set_MaA.elapsedTime() << std::endl;
//     // -------------------------------------------

//   // !### ---------------------------------------- AMGCL(buildin)
//   #elif defined (P_SOLVER_AMGCL_BUILTIN)

//     // -------------------------------------------
//     auto [ptr, idx, values] = gA.createPressureMatrixCSR();

//     auto amgcl_mat = std::tie(gA.iceltotCal, ptr, idx, values);
//     // -------------------------------------------
//     SolverBuiltin::params prm;

//     prm.solver.tol = simu.p_criteria;

//     prm.solver.maxiter = simu.p_iterMax;

//     prm.solver.ns_search = true;

//     prm.solver.verbose = true;
//     SolverBuiltin pSolver(amgcl_mat, prm);
//     // ---------------------------------
//     // SolverBuiltin pSolver(amgcl_mat);
//     // -------------------------------------------
//     // amgcl::io::mm_write("test_io_vec.mm", rhs.data(), n, 1);

//     amgcl::io::mm_write("Information/This_Matrix.mtx", amgcl_mat);
//     // -------------------------------------------

//   #endif



//   #if  defined (P_SOLVER_ELL) \
//     || defined (P_SOLVER_CSR) \
//     || defined (P_SOLVER_SPE)

//     pSolver.setTolerance(simu.p_criteria);

//   #endif


//   timer.p.stop();
//   // * ----------------------------------------------------

//   // ! ##  DfibGetEta(Dfib, Lo, gA) -------------
//   std::fill(Dfib.eta.begin(), Dfib.eta.end(), 0.0);

//   if (simu.DfibMethod == "OFF"){}
//   else if (simu.DfibMethod == "DFIB_Cylinder-X")
//     DFIB_CylinderX(Dfib, Lo, gA);
//   else if (simu.DfibMethod == "DFIB_Cylinder-Z")
//     DFIB_CylinderZ(Dfib, Lo, gA);
//   else
//     throw std::invalid_argument("DFIB method??");

//   // * ------------------------------------------------------

//     // ---------------------------------
//     OutputPLOT3D_Qfile(Dfib, simu, t1, T0, Lo, gA); ///Init
//     // ---------------------------------

//     // ---------------------------------
//     #ifdef TERBULENCE_SMAGORINSKY
//     T0.Viseff.resize(gA.iceltotCal, 0.0);
//     #endif
//     // ---------------------------------



//   // ***** perf sch
//   stopWatch Ptimer[20];



//     // =========================== GPU ===========================
//     // =========================== GPU ===========================
//     // =========================== GPU ===========================


//     // !## load the OCL program 

//     // size_t G = 1024;
//     // OCL.init("NVIDIA", CL_DEVICE_TYPE_GPU, 1*G); // Mesa, pocl

//     // std::vector<std::string> Source_File = {"inc/ocl/kernal/updateUandF.cl"
//     //                                         // ,"src/4_5_1_UpdateT1toT3.cl"
//     //                                         };

//     // OCL.SetKernelProgram_from_SourceFile_m(Source_File);
//     // // !## Set the Ocl buffer 

//     // // paramether
//     // ocl.LocalMax = cl::Buffer(ocl.ctx, CL_MEM_READ_WRITE, 1 * sizeof(cl_double));

//     // ocl.iter = cl::Buffer(ocl.ctx, CL_MEM_READ_WRITE, 1 * sizeof(cl_uint));

//     // OCL.LocalMax = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE, 1 * sizeof(cl_double));

//     // OCL.iter = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE, 1 * sizeof(cl_uint));

//     // OCL.T0_u = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
//     //                 T0.u.size() * sizeof(cl_double), T0.u.data());

//     // OCL.T0_v = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
//     //                 T0.v.size() * sizeof(cl_double), T0.v.data());

//     // OCL.T0_w = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
//     //                 T0.w.size() * sizeof(cl_double), T0.w.data());

//     // OCL.pressure = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
//     //                 t1.p.size() * sizeof(cl_double), t1.p.data());

//     // OCL.Dfib_eta = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
//     //                 Dfib.eta.size() * sizeof(cl_double), Dfib.eta.data());

//     // OCL.Dfib_f = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE, gridA.iceltot * sizeof(cl_double));


//     // OCL.Dx = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
//     //     gridA.Dx.size() * sizeof(cl_double), gridA.Dx.data());

//     // OCL.Dy= cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
//     //     gridA.Dy.size() * sizeof(cl_double), gridA.Dy.data());

//     // OCL.Dz = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
//     //     gridA.Dz.size()  * sizeof(cl_double), gridA.Dz.data());

//     // OCL.Dxs = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
//     //     gridA.Dys.size() * sizeof(cl_double), gridA.Dxs.data());

//     // OCL.Dys= cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
//     //     gridA.Dys.size() * sizeof(cl_double), gridA.Dys.data());

//     // OCL.Dzs = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
//     //     gridA.Dzs.size()  * sizeof(cl_double), gridA.Dzs.data());



//     // !## prepare for linear equations sovler.
//     // SorPopeLine_OCL_pre(OCL,ShareM, Sor, simu, T1, t1, Lo, gridA);




//     // =========================== GPU ===========================
//     // =========================== GPU ===========================
//     // =========================== GPU ===========================


//   // !  ------------- remove velocity (eta == 0) -------------
//   // removeVelocity(Dfib, T0, Lo, gA);
//   // removeVelocity(Dfib, T1, Lo, gA);
//   // removeVelocity(Dfib, T3, Lo, gA);
//     // * -------------------------------------------------------------
//   auto perLoopWTime = timer.beginNew.elapsedTime();
//   // ! ==================  main loop begin ==================
//   for (simu.loop = 1; simu.get_finishloop(); simu.finishloop())
//   {

//     // ---------------------------------
//     #if defined (TERBULENCE_SMAGORINSKY)
//         SmagorinskyModel(ShareM, simu, T0, t1, Lo, gA);
//     #endif
//     // ---------------------------------

//     // !## 1. get the T1/u* --------------------
//     /*  T1(Tsstart) */ timer.convectionDifussion.start();
//       ConvectionDifussion(simu, T0, T1, Lo, gA);
//       // BC_updateSlid(Lo, T1, gA);
//     /*  T1(Tsstart) */ timer.convectionDifussion.stop();
//     //*---------------------------------------------------

//     // -----------------------
//     get_Cfl(T1, Lo, gA, simu);
//     // simu.timestepper();
//     // -----------------------



//     // !## 2. get the pressure ------------ .2.
//      /* get the pressure */ timer.p.start(); 


//     // !### ---------------------------------------- SOR
//     #if defined (P_SOLVER_SOR)

//     if (simu.loop > 1)
//     {
//     // ---------------------------------
//     #if defined (PC_OMP)
//     SorPipeLine_omp(ShareM, Sor, simu, T1, t1, Lo, gA);
//     #elif defined (PC_SEQ)
//     SorPipeLine_seq(ShareM, Sor, simu, T1, t1, Lo, gA);
//     #endif
//     // ---------------------------------
//     }

//     // !### ---------------------------------------- CSR/SPE/ELL (BICG)(CG) / AMGCL
//     #elif defined (P_SOLVER_BICG_CSR) \
//       ||  defined (P_SOLVER_BICG_SPE) \
//       ||  defined (P_SOLVER_BICG_ELL) \
//       ||  defined (P_SOLVER_AMGCL_BUILTIN)

//     int pLoopIni  = 1;

//     if (simu.loop > pLoopIni)
//     {

//     // ---------------------------------
//     #if defined (PC_SEQ)
//     createBMatrix_seq(T1, Mx, simu, Lo, gA);
//     #elif defined (PC_OMP)
//     createBMatrix_omp(T1, Mx, simu, Lo, gA);
//     #endif
//     // ---------------------------------

//     // ---------------------------------
//       std::tie(simu.iters, simu.error) = pSolver(Mx.matB, Mx.X_result);
//     // ---------------------------------

//     // ---------------------------------
//       Pressure_transform_X_result(t1, Mx, Lo, gA);
//     // ---------------------------------
//     }

//     #endif 
//     // !### ---------------------------------------- endif
    

//     /* get the pressure */ timer.p.stop(); simu.printInfo();

//     auto [min_p , max_p] = getMax(Mx.X_result, Lo, gA);
//     std::cout  << "[Min_p, Max_p]:"<< min_p  << ", "<< max_p << std::endl;

//     // * .2. ------------- get the pressure ------------ .2.

//     // !## 3. update T1 to T3 (DFIB is included) -------------
//                             timer.updateT1toT3.start();
//     // ----------------------
//     #if defined (PC_SEQ)
//     update_UandF_seq(Dfib, simu, t1, T1, T3, Lo, gA);
//     #elif defined (PC_OMP)
//     update_UandF_omp(Dfib, simu, t1, T1, T3, Lo, gA); 
//     #endif
//     // ----------------------
//     // BC_updateSlid(Lo, T3, gA); //Ahmad didn't up it !!
//                             timer.updateT1toT3.stop();
//     // * 3 -------------------------------------------------------------


//     // !  ------------------ check staedy state (L2 norm) ------------------
//                              timer.checkL2norm.start();
//     #if defined (PC_SEQ)
//     CheckSteadyStatebyMaxVal_seq(simu, ShareM, T0, T3, Lo, gA); 
//     #elif defined (PC_OMP)
//     CheckSteadyStatebyMaxVal_omp(simu, ShareM, T0, T3, Lo, gA);
//     #endif
//                              timer.checkL2norm.stop();
//     // *  -------------------------------------------------------------------

//     // ! ## 4. cpoy ------------------
//               timer.updateT3toT0.start();
//     T0.copy_omp(T3);
//               timer.updateT3toT0.stop();
//     // * ----------------------------------------


//     // !## 5  Update the ghost cells  -------------
//     timer.BC.start();
    
//     BC_staggered_main( Lo, T0, t1, gA );
//     BC_staggered_copy( Lo, T0, T1, gA );

//     timer.BC.stop();
//     // * -----------------------------------------

//     // !## 6 get Cd and Cl  -------------
//     if (simu.DfibMethod != "OFF")
//     { getCD_CL(simu, Lo, Dfib, gA);}
//     // * --------------------------------

//     // !## 7 write plot3D formart .-------------
//     // if (simu.get_writefile()){
//     //   OutputPLOT3D_Qfile(Dfib, simu, t1, T3, Lo, gA);
//     // }
//     // * -----------------------------------------


//     // !## 8. IO loop------------------------------------------------
//     timer.beginNew.stop();
//     cout 
//     <<     "[simu]{time, dt} ="  << simu.get_time()          
//     <<     ", " << simu.dt       << endl
//     <<     "[wall time] < "      << timer.beginNew.elapsedTime() - perLoopWTime 
//     <<     " OF "                << timer.beginNew.elapsedTime()
//     <<     " >" << endl

//     // ! =============== NEXT time step ===============
//     <<     "===================================================" << endl
//     <<     "LOOP :"  << simu.loop + 1
//     <<     gA.show() << ", file :" << simu.get_file() << endl;
//     // * ----------------------------------------------------------


//     timer.beginNew.start(); perLoopWTime = timer.beginNew.elapsedTime();

//     recordingTime(timer, simu, ShareM);
//   }
//   timer.beginNew.stop();
  // ! ==================  main loop begin ==================
//   return true;
// }



    // CenterQuickScheme_OCL(OCL, So, simu, T0, T1, Lo, gridA);

    // SorPipeLine_OCL(OCL, ShareM, Sor, simu, T1, t1, Lo, gridA);

    // UpdateT1toT3_OCL(OCL, Dfib, simu, t1, T1, T3, Lo, gridA);

    // UpdateT3toT0_OCL(OCL, simu, T3, T0, Lo, gridA);
    
    // BoundaryCondtion_OCL(OCL,simu, Lo, T0, t1, gridA);


    // if (simu.get_writefile())
    // {
        // OCL.queue.enqueueReadBuffer(OCL.T3_u, CL_FALSE, 0,gridA.iceltot*sizeof(cl_double),T3.u.data());
        // OCL.queue.enqueueReadBuffer(OCL.T3_v, CL_FALSE, 0,gridA.iceltot*sizeof(cl_double),T3.v.data());
        // OCL.queue.enqueueReadBuffer(OCL.T3_w, CL_TRUE, 0,gridA.iceltot*sizeof(cl_double),T3.w.data());

        // OutputPLOT3D_Qfile(Dfib, simu, t1, T3, Lo, gridA);

      // ReadPlot3DadnCheckL2(ShareM,main_count,Dfib,simu,t1,T3,Lo,gridA);
      // getGhiaProfile(simu, T3, Lo, gridA);
    // }
