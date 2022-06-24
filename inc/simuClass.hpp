#pragma once 
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>

class simuClass{


	public:
	simuClass(){
		cd.resize(delayIO);
		cl.resize(delayIO);
	}

	std::vector<double> cd;
	std::vector<double> cl;
	std::vector<std::string> io_string;

	int TID = 0;

	int PID = 0;

	int n = 1;

	int PS_bool_OMP_Manual = 0; 	// Manual loop parallelization

	int PS_bool_OMP_FOR = 0;

	int Locality;

	int delayIO = 100;


	// *--------------------------
	double Re;
	double nu;

	bool set_Re(double re){

		Re = re;
    	nu = 1.0 / re;
		return true;
	}

	// *--------------------------


	// ! ----------  convection and difussion term ----------

	double alpha, error{0.};

	bool First_s = true;

	bool Second_s = false;

	// * ----------  convection and difussion term ----------



	size_t accumulate = 0;

	// !-----------  Solver Pressure ----------- 
	size_t p_iterMax = 5000;

	int iters{0};

	bool firstSOR{true};

	double cfl_max = 0.0, P_max, P_min;

	std::vector<double> U_max, U_min;

	double p_criteria = 1.e-3; 
	// amgcl def 1e-8
	// ahmad def 1e-3

	double p_sor_omega = 1.8;

	size_t p_sor_iter_max = 5000;
	// *-----------  Solver Pressure ----------- 


	// ! Time Variance Authority  ------------------

	double dt,  dt_init, dt_next, dt_file;

	std::vector<double> dt_pre;

	size_t currentFile = 0;

	size_t loop_max, loop;

	double 	time_max = 0.0, time = 0.0, time_file = 0.0, time_file_next = 0.0;

	auto set_time(double t){
		time = t;
		time_file = t;
	}


	auto get_simuT(){ return time;}

	auto get_time(){ return time;}

	auto set_time_max(double t){ time_max= t;}

	bool set_dt(double t){

		if (time_max == 0)	std::cout << "time_max " << std::endl;

		dt_next = dt_init = dt = t;

		dt_pre.resize(4,t);

		loop_max = time_max / dt;

		return true;
	}

	// !------- io ------------

	void init_fileStep(double t){ dt_file = t; }

	double get_file(){ return ( time / dt_file);	}

	bool get_writefile(){
		time_file_next = time_file + dt_file;

		if (time > time_file_next ){

			currentFile++;

			time_file += dt_file;

			return true;

		}else{return false;}
	}



	bool get_finishloop(){
		if (loop < loop_max || time < time_max )
			return true;
		else
			return false;
	}
	// *------- io ------------


	bool finishloop(){

		size_t j = dt_pre.size()-1;

		for (size_t i=0; i < dt_pre.size()-1 ; ++i, --j){
			dt_pre[j] = dt_pre[j-1];
		}

		dt_pre[0] = dt;

		time += dt;

		loop++;

		dt = dt_next;

		return true;
	}


	inline void printInfo(){
		std::cout 
		<< "[" << ZONE()
		<< "] iter:" << iters 
		<< ", error " << error 
		<< std::endl;
	}


	// * data file 



	void timestepper()
	{

		dt_next = 1e10;
		for (auto &x:U_max){
			dt_next = std::min(dt_next, cflFactor/(Re*x*x));
		}

		for (auto &x:U_min){
			dt_next = std::min(dt_next, gridFoFactor*0.5*Re*x*x);
		}

		if (dt_next < dtMin) { dt_next = dtMin; }

	}


	
	std::string ZONE(){
		std::string A;

	#if defined (P_SOLVER_SOR)

		A += "SOR";

	#elif defined (P_SOLVER_BICG_CSR)

		A += "CSR (BICG)";

	#elif defined (P_SOLVER_AMGCL_BUILTIN)


	#elif defined (P_SOLVER_BICG_SPE)

		A += "SPE (BICG)";

	#elif defined (P_SOLVER_BICG_ELL)

		A += "ELL (BICG)";

	#elif defined (P_SOLVER_AMGCL_EIGEN)

		A += "EIGEN (AMGCL)";

	#elif defined (P_SOLVER_EIGEN_CSR)

		A += "EIGEN";


	#endif
		return A;
	}

	std::string TurbulenceModeling;

	std::string DfibMethod = "OFF";	
	bool Dfib_boolT = false;
	std::string Periodic;
	
	std::string Neumann;
	
	// std::string convection_and_difussion;
	

#ifdef TEMPORAL_DISCRETIZATION_1_ORDER

	std::string timeDiscretization = "explicit_Euler";

#endif


	std::string Sequential_TEXT =  R"(
┌──────────────────────────────────────────────────────────────────────┐
│                                                                      │
│                                             ▄      ▀           ▀▀█   │
│  ▄▄▄    ▄▄▄    ▄▄▄▄  ▄   ▄   ▄▄▄   ▄ ▄▄   ▄▄█▄▄  ▄▄▄     ▄▄▄     █   │
│ █   ▀  █▀  █  █▀ ▀█  █   █  █▀  █  █▀  █    █      █    ▀   █    █   │
│  ▀▀▀▄  █▀▀▀▀  █   █  █   █  █▀▀▀▀  █   █    █      █    ▄▀▀▀█    █   │
│ ▀▄▄▄▀  ▀█▄▄▀  ▀█▄██  ▀▄▄▀█  ▀█▄▄▀  █   █    ▀▄▄  ▄▄█▄▄  ▀▄▄▀█    ▀▄▄ │
│                   █                                                  │
│                   ▀                                                  │
└──────────────────────────────────────────────────────────────────────┘    
	)";





	std::string CFDMX_TEXT = R"(
┌──────────────────────────────────────────┐
│                                          │
│   [0;1;34;94m▄▄▄[0m  [0;34m▄▄▄▄▄▄[0m [0;34m▄▄▄▄[0m          [0;37m▄[0m    [0;37m▄[0m [0;37m▄[0m    [0;1;30;90m▄[0m│
│ [0;34m▄▀[0m   [0;34m▀[0m [0;34m█[0m      [0;34m█[0m   [0;37m▀▄[0m        [0;37m██[0m  [0;1;30;90m██[0m  [0;1;30;90m█[0m  [0;1;30;90m█[0m │
│ [0;34m█[0m      [0;37m█▄▄▄▄▄[0m [0;37m█[0m    [0;37m█[0m        [0;1;30;90m█[0m [0;1;30;90m██[0m [0;1;30;90m█[0m   [0;1;30;90m██[0m  │
│ [0;37m█[0m      [0;37m█[0m      [0;37m█[0m    [0;1;30;90m█[0m  [0;1;30;90m▀▀▀[0m   [0;1;30;90m█[0m [0;1;30;90m▀[0;1;34;94m▀[0m [0;1;34;94m█[0m  [0;1;34;94m▄▀▀▄[0m │
│  [0;37m▀▄▄▄▀[0m [0;1;30;90m█[0m      [0;1;30;90m█▄▄▄▀[0m         [0;1;34;94m█[0m    [0;1;34;94m█[0m [0;1;34;94m▄▀[0m  [0;34m▀▄[0m│
│                                          │
│                                          │
└──────────────────────────────────────────┘
	)";



	private:
		int time_Discretization = 0;

		double                              dtMin           {1e-12},
		/* Timstep adjustment: CFL       */ cflFactor       {0.05},
		/* Timstep adjustment: Grid Fo   */ gridFoFactor    {0.15};

};

