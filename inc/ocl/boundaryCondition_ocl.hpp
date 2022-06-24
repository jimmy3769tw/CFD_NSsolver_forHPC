#include <ocl/oclClass.hpp>


void BC_staggered_main_vel_ocl(
	cl::Buffer &currt_u,
	cl::Buffer &currt_v,
	cl::Buffer &currt_w,
    OCLstruct &OCL,
    simuClass& simu,
    calDomain& Lo,
    velocity& T0,
    pressure& t1,
    grid& gridA
)
{
    using cl_double = double;
    using cl_uint = uint;
    using cl_int = int32_t;

	static std::vector<std::string>functionName(6);

	functionName[0] = "updateVel_atBoundary0";
	functionName[1] = "updateVel_atBoundary1";
	functionName[2] = "updateVel_atBoundary2";
	functionName[3] = "updateVel_atBoundary3";
	functionName[4] = "updateVel_atBoundary4";
	functionName[5] = "updateVel_atBoundary5";

	const std::string sourcFileName = "inc/ocl/kernal/BoundaryCondition.cl";

	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				const cl_double,const cl_double ,const  cl_double,
				const cl_double,const cl_double ,const  cl_double,
				const cl_uint, const cl_uint, const cl_uint, 
				const cl_uint, const cl_uint >Kernel_BC_0(OCL.prg_m[sourcFileName],functionName[0]);

	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				const cl_double,const cl_double ,const  cl_double,
				const cl_double,const cl_double ,const  cl_double,
				const cl_uint, const cl_uint, const cl_uint, 
				const cl_uint, const cl_uint >Kernel_BC_1(OCL.prg_m[sourcFileName],functionName[1]);


	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				const cl_double,const cl_double ,const  cl_double,
				const cl_double,const cl_double ,const  cl_double,
				const cl_uint, const cl_uint, const cl_uint, 
				const cl_uint, const cl_uint >Kernel_BC_2(OCL.prg_m[sourcFileName],functionName[2]);

	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				const cl_double,const cl_double ,const  cl_double,
				const cl_double,const cl_double ,const  cl_double,
				const cl_uint, const cl_uint, const cl_uint, 
				const cl_uint, const cl_uint >Kernel_BC_3(OCL.prg_m[sourcFileName],functionName[3]);


	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				const cl_double,const cl_double ,const  cl_double,
				const cl_double,const cl_double ,const  cl_double,
				const cl_uint, const cl_uint, const cl_uint, 
				const cl_uint, const cl_uint >Kernel_BC_4(OCL.prg_m[sourcFileName],functionName[4]);

	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				const cl_double,const cl_double ,const  cl_double,
				const cl_double,const cl_double ,const  cl_double,
				const cl_uint, const cl_uint, const cl_uint, 
				const cl_uint, const cl_uint >Kernel_BC_5(OCL.prg_m[sourcFileName],functionName[5]);



	auto config2D_x = cl::EnqueueArgs(OCL.queue,{ (gridA.ny+OCL.D2[0]-1)/OCL.D2[0]*OCL.D2[0],
                                                  (gridA.nz+OCL.D2[1]-1)/OCL.D2[1]*OCL.D2[1]}, 
												   OCL.D2[0],OCL.D2[1]); 

	auto config2D_y = cl::EnqueueArgs(OCL.queue,{ (gridA.nx+OCL.D2[0]-1)/OCL.D2[0]*OCL.D2[0],
                                                  (gridA.nz+OCL.D2[1]-1)/OCL.D2[1]*OCL.D2[1]},
												   OCL.D2[0],OCL.D2[1]); 

	auto config2D_z = cl::EnqueueArgs(OCL.queue,{ (gridA.nx+OCL.D2[0]-1)/OCL.D2[0]*OCL.D2[0],
                                                  (gridA.ny+OCL.D2[1]-1)/OCL.D2[1]*OCL.D2[1]},
												   OCL.D2[0],OCL.D2[1]);

    // * ----------------- BoundaryCondtion
    static auto [Num, Dir] = getNumDIr();

	Kernel_BC_0(config2D_x, currt_u, currt_v, currt_w, Num[0], Num[1], Num[2], Dir[0], Dir[1], Dir[2], 
						gridA.nx, gridA.ny, gridA.nz, gridA.gC, gridA.nz_ny);

	Kernel_BC_1(config2D_x, currt_u, currt_v, currt_w, Num[0], Num[1], Num[2], Dir[0], Dir[1], Dir[2], 
						gridA.nx, gridA.ny, gridA.nz, gridA.gC, gridA.nz_ny);

	Kernel_BC_2(config2D_y, currt_u, currt_v, currt_w, Num[0], Num[1], Num[2], Dir[0], Dir[1], Dir[2], 
						gridA.nx, gridA.ny, gridA.nz, gridA.gC, gridA.nz_ny);

	Kernel_BC_3(config2D_y, currt_u, currt_v, currt_w, Num[0], Num[1], Num[2], Dir[0], Dir[1], Dir[2], 
						gridA.nx, gridA.ny, gridA.nz, gridA.gC, gridA.nz_ny);

	Kernel_BC_4(config2D_z, currt_u, currt_v, currt_w, Num[0], Num[1], Num[2], Dir[0], Dir[1], Dir[2], 
						gridA.nx, gridA.ny, gridA.nz, gridA.gC, gridA.nz_ny);

	Kernel_BC_5(config2D_z, currt_u, currt_v, currt_w, Num[0], Num[1], Num[2], Dir[0], Dir[1], Dir[2], 
						gridA.nx, gridA.ny, gridA.nz, gridA.gC, gridA.nz_ny);
}




void BC_staggered_main_vel_ocl(
	cl::Buffer &currt_u,
	cl::Buffer &currt_v,
	cl::Buffer &currt_w,
    OCLstruct &OCL,
    simuClass& simu,
    calDomain& Lo,
    velocity& T0,
    pressure& t1,
    grid& gridA
)
{
    using cl_double = double;
    using cl_uint = uint;
    using cl_int = int32_t;

	static std::vector<std::string>functionName(6);

	functionName[0] = "updateVel_atBoundary0";
	functionName[1] = "updateVel_atBoundary1";
	functionName[2] = "updateVel_atBoundary2";
	functionName[3] = "updateVel_atBoundary3";
	functionName[4] = "updateVel_atBoundary4";
	functionName[5] = "updateVel_atBoundary5";

	const std::string sourcFileName = "inc/ocl/kernal/BoundaryCondition.cl";

	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				const cl_double,const cl_double ,const  cl_double,
				const cl_double,const cl_double ,const  cl_double,
				const cl_uint, const cl_uint, const cl_uint, 
				const cl_uint, const cl_uint >Kernel_BC_0(OCL.prg_m[sourcFileName],functionName[0]);

	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				const cl_double,const cl_double ,const  cl_double,
				const cl_double,const cl_double ,const  cl_double,
				const cl_uint, const cl_uint, const cl_uint, 
				const cl_uint, const cl_uint >Kernel_BC_1(OCL.prg_m[sourcFileName],functionName[1]);


	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				const cl_double,const cl_double ,const  cl_double,
				const cl_double,const cl_double ,const  cl_double,
				const cl_uint, const cl_uint, const cl_uint, 
				const cl_uint, const cl_uint >Kernel_BC_2(OCL.prg_m[sourcFileName],functionName[2]);

	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				const cl_double,const cl_double ,const  cl_double,
				const cl_double,const cl_double ,const  cl_double,
				const cl_uint, const cl_uint, const cl_uint, 
				const cl_uint, const cl_uint >Kernel_BC_3(OCL.prg_m[sourcFileName],functionName[3]);


	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				const cl_double,const cl_double ,const  cl_double,
				const cl_double,const cl_double ,const  cl_double,
				const cl_uint, const cl_uint, const cl_uint, 
				const cl_uint, const cl_uint >Kernel_BC_4(OCL.prg_m[sourcFileName],functionName[4]);

	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				const cl_double,const cl_double ,const  cl_double,
				const cl_double,const cl_double ,const  cl_double,
				const cl_uint, const cl_uint, const cl_uint, 
				const cl_uint, const cl_uint >Kernel_BC_5(OCL.prg_m[sourcFileName],functionName[5]);



	auto config2D_x = cl::EnqueueArgs(OCL.queue,{ (gridA.ny+OCL.D2[0]-1)/OCL.D2[0]*OCL.D2[0],
                                                  (gridA.nz+OCL.D2[1]-1)/OCL.D2[1]*OCL.D2[1]}, 
												   OCL.D2[0],OCL.D2[1]); 

	auto config2D_y = cl::EnqueueArgs(OCL.queue,{ (gridA.nx+OCL.D2[0]-1)/OCL.D2[0]*OCL.D2[0],
                                                  (gridA.nz+OCL.D2[1]-1)/OCL.D2[1]*OCL.D2[1]},
												   OCL.D2[0],OCL.D2[1]); 

	auto config2D_z = cl::EnqueueArgs(OCL.queue,{ (gridA.nx+OCL.D2[0]-1)/OCL.D2[0]*OCL.D2[0],
                                                  (gridA.ny+OCL.D2[1]-1)/OCL.D2[1]*OCL.D2[1]},
												   OCL.D2[0],OCL.D2[1]);

    // * ----------------- BoundaryCondtion
    static auto [Num, Dir] = getNumDIr();

	Kernel_BC_0(config2D_x, currt_u, currt_v, currt_w, Num[0], Num[1], Num[2], Dir[0], Dir[1], Dir[2], 
						gridA.nx, gridA.ny, gridA.nz, gridA.gC, gridA.nz_ny);

	Kernel_BC_1(config2D_x, currt_u, currt_v, currt_w, Num[0], Num[1], Num[2], Dir[0], Dir[1], Dir[2], 
						gridA.nx, gridA.ny, gridA.nz, gridA.gC, gridA.nz_ny);

	Kernel_BC_2(config2D_y, currt_u, currt_v, currt_w, Num[0], Num[1], Num[2], Dir[0], Dir[1], Dir[2], 
						gridA.nx, gridA.ny, gridA.nz, gridA.gC, gridA.nz_ny);

	Kernel_BC_3(config2D_y, currt_u, currt_v, currt_w, Num[0], Num[1], Num[2], Dir[0], Dir[1], Dir[2], 
						gridA.nx, gridA.ny, gridA.nz, gridA.gC, gridA.nz_ny);

	Kernel_BC_4(config2D_z, currt_u, currt_v, currt_w, Num[0], Num[1], Num[2], Dir[0], Dir[1], Dir[2], 
						gridA.nx, gridA.ny, gridA.nz, gridA.gC, gridA.nz_ny);

	Kernel_BC_5(config2D_z, currt_u, currt_v, currt_w, Num[0], Num[1], Num[2], Dir[0], Dir[1], Dir[2], 
						gridA.nx, gridA.ny, gridA.nz, gridA.gC, gridA.nz_ny);
}

void BC_staggered_main_vel_ocl(
    OCLstruct &OCL,
    simuClass& simu,
    calDomain& Lo,
    velocity& T0,
    pressure& t1,
    grid& gridA
)
{

    using cl_double = double;
    using cl_uint = uint;
    using cl_int = int32_t;

	static const std::string sourcFileName = "inc/ocl/kernal/BoundaryCondition.cl";

	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				const cl_uint, const cl_uint, const cl_uint, 
				const cl_uint, const cl_uint >Kernel_updatePressure_atBoundary0and1(OCL.prg_m[sourcFileName],"updatePressure_atBoundary0and1");

	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				const cl_uint, const cl_uint, const cl_uint, 
				const cl_uint, const cl_uint >Kernel_updatePressure_atBoundary2and3(OCL.prg_m[sourcFileName],"updatePressure_atBoundary2and3");

	static cl::KernelFunctor<cl::Buffer&, 
				const cl_uint, const cl_uint, const cl_uint, 
				const cl_uint, const cl_uint >Kernel_updatePressure_atBoundary4and5(OCL.prg_m[sourcFileName],"updatePressure_atBoundary4and5");


	auto config2D_x = cl::EnqueueArgs(OCL.queue,{ (gridA.ny+OCL.D2[0]-1)/OCL.D2[0]*OCL.D2[0],
                                                  (gridA.nz+OCL.D2[1]-1)/OCL.D2[1]*OCL.D2[1]}, 
												   OCL.D2[0],OCL.D2[1]); 

	auto config2D_y = cl::EnqueueArgs(OCL.queue,{ (gridA.nx+OCL.D2[0]-1)/OCL.D2[0]*OCL.D2[0],
                                                  (gridA.nz+OCL.D2[1]-1)/OCL.D2[1]*OCL.D2[1]},
												   OCL.D2[0],OCL.D2[1]); 

	auto config2D_z = cl::EnqueueArgs(OCL.queue,{ (gridA.nx+OCL.D2[0]-1)/OCL.D2[0]*OCL.D2[0],
                                                  (gridA.ny+OCL.D2[1]-1)/OCL.D2[1]*OCL.D2[1]},
												   OCL.D2[0],OCL.D2[1]);

	Kernel_updatePressure_atBoundary0and1(config2D_x, OCL.pressure, Num[0], Num[1], Num[2], Dir[0], Dir[1], Dir[2], 
						gridA.nx, gridA.ny, gridA.nz, gridA.gC, gridA.nz_ny);

	Kernel_updatePressure_atBoundary2and3(config2D_y, OCL.pressure, Num[0], Num[1], Num[2], Dir[0], Dir[1], Dir[2], 
						gridA.nx, gridA.ny, gridA.nz, gridA.gC, gridA.nz_ny);

	Kernel_updatePressure_atBoundary4and5(config2D_z, OCL.pressure, Num[0], Num[1], Num[2], Dir[0], Dir[1], Dir[2], 
						gridA.nx, gridA.ny, gridA.nz, gridA.gC, gridA.nz_ny);

}
