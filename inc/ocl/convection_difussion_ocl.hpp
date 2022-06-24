
#include "4_1_1_CenterQuickscheme_OCL.hpp"
#ifdef OCL_ON

void CenterQuickScheme_OCL(
    OCLstruct &OCL, 
    sourceTerm& So,
    simuClass& simu,
    velocity& T0,
    velocity& T1,
    calDomain& Lo,
    grid& gridA
)
{


	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&, cl::Buffer&,cl::Buffer&,cl::Buffer&,
                            const cl::Buffer&, const cl::Buffer&, const cl::Buffer&,const cl::Buffer&, const cl::Buffer&, const cl::Buffer&,
							const cl::Buffer&, double, double, const int, const int, const int, int>
                            Kernel_CenterQuickScheme_u(OCL.prg_m["src/4_1_1_CenterQuickscheme.cl"], "CenterQuickScheme_u");

	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&, cl::Buffer&,cl::Buffer&,cl::Buffer&,
                            const cl::Buffer&,const cl::Buffer&,const cl::Buffer&,const cl::Buffer&,const cl::Buffer&,const cl::Buffer&,
							const cl::Buffer&, double, double, const int, const int, const int, int>
                            Kernel_CenterQuickScheme_v(OCL.prg_m["src/4_1_1_CenterQuickscheme.cl"], "CenterQuickScheme_v");

	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&, cl::Buffer&,cl::Buffer&,cl::Buffer&,
                            const cl::Buffer&,const cl::Buffer&,const cl::Buffer&,const cl::Buffer&,const cl::Buffer&,const cl::Buffer&,
							const cl::Buffer&, double, double, const int, const int, const int, int>
                            Kernel_CenterQuickScheme_w(OCL.prg_m["src/4_1_1_CenterQuickscheme.cl"], "CenterQuickScheme_w");

	auto config3D = cl::EnqueueArgs(OCL.queue,{ (gridA.nx+OCL.D3[0]-1)/OCL.D3[0]*OCL.D3[0],
												(gridA.ny+OCL.D3[1]-1)/OCL.D3[1]*OCL.D3[1],
												(gridA.nz+OCL.D3[2]-1)/OCL.D3[2]*OCL.D3[2]},
												{OCL.D3[0], OCL.D3[1], OCL.D3[2]}); 

// FIXME:

	Kernel_CenterQuickScheme_u(config3D,  
                                    OCL.T0_u,OCL.T0_v, OCL.T0_w,
                                    OCL.T1_u,OCL.T1_v, OCL.T1_w,
                                    OCL.Dx, OCL.Dy, OCL.Dz,
                                    OCL.Dxs, OCL.Dys, OCL.Dzs,
                                    OCL.NEIBcell,
                                    simu.dt,simu.nu, gridA.nx,gridA.ny,gridA.nz, gridA.gC);


	Kernel_CenterQuickScheme_v(config3D, 
                                   OCL.T0_u,OCL.T0_v, OCL.T0_w,
                                   OCL.T1_u,OCL.T1_v, OCL.T1_w,
                                   OCL.Dx, OCL.Dy, OCL.Dz,
                                   OCL.Dxs, OCL.Dys, OCL.Dzs,
                                   OCL.NEIBcell,
                                   simu.dt, simu.nu, gridA.nx, gridA.ny, gridA.nz, gridA.gC);

    Kernel_CenterQuickScheme_w(config3D,  
                                   OCL.T0_u,OCL.T0_v, OCL.T0_w,
                                   OCL.T1_u,OCL.T1_v, OCL.T1_w,
                                   OCL.Dx, OCL.Dy, OCL.Dz,
                                   OCL.Dxs, OCL.Dys, OCL.Dzs,OCL.NEIBcell,
                                   simu.dt,simu.nu, gridA.nx, gridA.ny, gridA.nz,gridA.gC);
    
}



#endif
