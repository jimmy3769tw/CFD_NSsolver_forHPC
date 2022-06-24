#include "controlPanel.hpp"

#if defined (PC_HYBRID_MPI_OMP)

#include "run/RUN_CPU_HYBRID_MPI.hpp"

#else

#include "run/RUN_CPU.hpp"
#include "run/RUN_CPU_PC_METHOD.hpp"

#endif 



using std::string, 
      std::vector, 
      std::abs, 
      std::pow, 
      std::sqrt, 
      std::cout, 
      std::endl; 

int main(int argc, char **argv){


    // * CFD_MX_struct
    clockstruct timer;

    simuClass simu;

    grid gA;

    // *   parameter  ------------
    simu.set_time_max(40);

    simu.set_dt(0.002); 

    simu.set_Re(40.0);

    simu.init_fileStep(1);
    // *   parameter  ------------


    // *   Accuracy  ------------
    std::cout << "\n----------------------------------------------------" 
              << "\n| Gridder :                       | "             << gA.Gridder
              << "\n| [Lx:Ly:Lz]" << gA.lx << " : " << gA.ly << " : " << gA.lz
              << "\n----------------------------------------------------" ;


    // *  >>>>>>>>>>>>>----- SETTING (RUNTIME)
    // DFIB_Cylinder-X, DFIB_Cylinder-Z, "DFIB_Cylinder-Y", "OFF"
    simu.DfibMethod = "DFIB_Cylinder-Z"; 
    simu.Dfib_boolT = true;

    simu.Locality = 1;

    // *  >>>>>>>>>>>>>----- SETTING (RUNTIME)


    // ! RUN ============
    int mpi_word_rank = 0;

    if (mpi_word_rank == 0){
        std::cout << "\n# MainLoop :" << simu.loop_max 
                  << "\n----------------------------------------------------\n";
    }
    
    // ! TEXT ============
    #if defined (PC_SEQ)

    cout << simu.Sequential_TEXT << endl;

    #elif defined (PC_OMP)


    #elif defined (PC_HYBRID_MPI_OMP)


    #elif defined (PC_OCL)

    #endif


    #if defined (PC_SEQ) || defined (PC_OMP)
    
    // runSeqOmp(gA, timer, simu, argc, argv);
    runSeqOmp_predictionMethod(gA, timer, simu, argc, argv);
    #elif defined (PC_HYBRID_MPI_OMP)

    Run_Hy_MPI_OpenMP(gA, timer, simu, argc, argv);

    #elif defined (PC_OCL)

    #endif

    // ! RUN RUN RUN ============

    std::cout << "Finish !\n";

    return 0;

}