#include <Boundary_Coundition.hpp>

void copy_vel_ocl(
    OCLstruct &OCL,
    simuClass& simu,
    velocity& T3,
    velocity& T0,
    calDomain& Lo,
    grid& gridA
){

	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&, 
                            const cl::Buffer&,const cl::Buffer&,const cl::Buffer&,
                            const int,const int, const  int, const  int >
                            Kernel_UpdateT3toT0(OCL.prg_m["inc/kernal/copyVel.cl"], "copyVel");

	auto config3D = cl::EnqueueArgs(OCL.queue,{ (gridA.nx+OCL.D3[0]-1)/OCL.D3[0]*OCL.D3[0],
												(gridA.ny+OCL.D3[1]-1)/OCL.D3[1]*OCL.D3[1],
												(gridA.nz+OCL.D3[2]-1)/OCL.D3[2]*OCL.D3[2]},
												{OCL.D3[0], OCL.D3[1], OCL.D3[2]}); 

	Kernel_UpdateT3toT0(config3D,OCL.T0_u, OCL.T0_v, OCL.T0_w,
                                 OCL.T3_u, OCL.T3_v, OCL.T3_w,
                                 gridA.nx, gridA.ny, gridA.nz, gridA.gC);

}



void update_UandF_ocl(
    OCLstruct& OCL, 
    DfibArray& Dfib,
    simuClass& simu,
    pressure& t1,
    velocity& T1,
    velocity& T3,
    calDomain& Lo,
    grid& gridA
){
    double u_solid = 0.0; 
    double v_solid = 0.0; 
    double w_solid = 0.0; 
    // !---------------------------------

    using cl_double = double;
    using cl_uint = uint;
    using cl_int =int32_t;
    string fileName = "inc/ocl/kernal/updateUandF.cl";
    // *---------------------------------

    // ! ---------------------------------
	static cl::KernelFunctor<const cl::Buffer&,const cl::Buffer&,const cl::Buffer&,
                             cl::Buffer&, cl::Buffer&, cl::Buffer&, 
                             const cl::Buffer&, const cl::Buffer&, const cl::Buffer&,
                             const cl::Buffer&, cl::Buffer&, const cl::Buffer&,
                             const cl_double,const cl_double, const cl_double,
                             const cl_double, const cl_double, 
                             const cl_uint, const cl_uint, const cl_uint, const cl_uint,
                             const cl_uint
                             >
                             kernel_update_UandF_x(OCL.prg_m[fileName], "update_UandF_x");


	static cl::KernelFunctor<const cl::Buffer&,const cl::Buffer&,const cl::Buffer&,
                             cl::Buffer&, cl::Buffer&, cl::Buffer&, 
                             const cl::Buffer&, const cl::Buffer&, const cl::Buffer&,
                             const cl::Buffer&, cl::Buffer&, const cl::Buffer&,
                             const cl_double,const cl_double, const cl_double,
                             const cl_double, const cl_double, 
                             const cl_uint, const cl_uint, const cl_uint, const cl_uint,
                             const cl_uint
                             >
                             Kernel_update_UandF_y(OCL.prg_m[fileName], "update_UandF_y");




	static cl::KernelFunctor<const cl::Buffer&,const cl::Buffer&,const cl::Buffer&,
                             cl::Buffer&, cl::Buffer&, cl::Buffer&, 
                             const cl::Buffer&, const cl::Buffer&, const cl::Buffer&,
                             const cl::Buffer&, cl::Buffer&, const cl::Buffer&,
                             const cl_double,const cl_double, const cl_double,
                             const cl_double, const cl_double, 
                             const cl_uint, const cl_uint, const cl_uint, const cl_uint,
                             const cl_uint
                             >
                             Kernel_update_UandF_z(OCL.prg_m[fileName], "update_UandF_z");
    // *---------------------------------



	auto config3D = cl::EnqueueArgs(OCL.queue,{ (gridA.nx+OCL.D3[0]-1)/OCL.D3[0]*OCL.D3[0],
												(gridA.ny+OCL.D3[1]-1)/OCL.D3[1]*OCL.D3[1],
												(gridA.nz+OCL.D3[2]-1)/OCL.D3[2]*OCL.D3[2]},
												{OCL.D3[0], OCL.D3[1], OCL.D3[2]}); 


	kernel_update_UandF_x(config3D,OCL.T1_u,OCL.T1_v, OCL.T1_w,
                                 OCL.T3_u,OCL.T3_v, OCL.T3_w,
                                 OCL.Dxs, OCL.Dys, OCL.Dzs,
                                 OCL.pressure, OCL.Dfib_f, OCL.Dfib_eta,
                                 u_solid, v_solid, w_solid,
                                 simu.dt, simu.nu,
                                 gridA.nx, gridA.ny, gridA.nz,gridA.gC,
                                 gridA.iceltot);


	Kernel_update_UandF_y(config3D,OCL.T1_u,OCL.T1_v, OCL.T1_w,
                                 OCL.T3_u,OCL.T3_v, OCL.T3_w,
                                 OCL.Dxs, OCL.Dys, OCL.Dzs,
                                 OCL.pressure, OCL.Dfib_f, OCL.Dfib_eta,
                                 u_solid, v_solid, w_solid,
                                 simu.dt, simu.nu,
                                 gridA.nx, gridA.ny, gridA.nz,gridA.gC,
                                 gridA.iceltot);


	Kernel_update_UandF_z(config3D,OCL.T1_u,OCL.T1_v, OCL.T1_w,
                                 OCL.T3_u,OCL.T3_v, OCL.T3_w,
                                 OCL.Dxs, OCL.Dys, OCL.Dzs,
                                 OCL.pressure, OCL.Dfib_f, OCL.Dfib_eta,
                                 u_solid, v_solid, w_solid,
                                 simu.dt, simu.nu,
                                 gridA.nx, gridA.ny, gridA.nz,gridA.gC,
                                 gridA.iceltot);

}