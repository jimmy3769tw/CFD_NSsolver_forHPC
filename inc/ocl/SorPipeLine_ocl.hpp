

void prepare_for_SorPopeLine_OCL(
    OCLstruct &OCL,
    shareMenory& ShareM,
    SORcoefficient& Sor,
    simuClass& simu,
    velocity& T1,
    pressure& t1,
    calDomain& Lo,
    grid& gridA
)
{


	// * cl::Buffer
    OCL.SorCoef = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE, 8*gridA.iceltotCal*sizeof(cl_double));


	static cl::KernelFunctor<cl::Buffer&,
                            const cl::Buffer&,const cl::Buffer&,const cl::Buffer&,
                            const cl::Buffer&,const cl::Buffer&,const cl::Buffer&,
                            int, int, int, int >
                            Kernel_Sor_pipeline_getCoef(OCL.prg_m["inc/kernal/SorPipeLine.cl"], "Sor_pipeline_getCoef");


	auto config3D = cl::EnqueueArgs(OCL.queue,{ (gridA.nx+OCL.D3[0]-1)/OCL.D3[0]*OCL.D3[0],
												(gridA.ny+OCL.D3[1]-1)/OCL.D3[1]*OCL.D3[1],
												(gridA.nz+OCL.D3[2]-1)/OCL.D3[2]*OCL.D3[2]},
												{OCL.D3[0], OCL.D3[1], OCL.D3[2]}); 

	Kernel_Sor_pipeline_getCoef(config3D,OCL.SorCoef,
                                 OCL.Dx, OCL.Dy, OCL.Dz,
                                 OCL.Dxs, OCL.Dys, OCL.Dzs,
                                 gridA.nx, gridA.ny, gridA.nz, gridA.gC);

}



void SorPipeLine_OCL(
    OCLstruct &OCL,
    shareMenory& ShareM,
    SORcoefficient& Sor,
    simuClass& simu,
    velocity& T1,
    pressure& t1,
    calDomain& Lo,
    grid& gridA
)
{


    cout << "passing" << endl;
    // TODO add max value




	// OCL.queue.enqueueFillBuffer(OCL.LocalMax, 0, 0, 1*sizeof(cl_double));
	OCL.queue.enqueueFillBuffer(OCL.iter, 0, 0,1*sizeof(cl_uint));


	auto config3D = cl::EnqueueArgs(OCL.queue,{ (gridA.nx+OCL.D3[0]-1)/OCL.D3[0]*OCL.D3[0],
												(gridA.ny+OCL.D3[1]-1)/OCL.D3[1]*OCL.D3[1],
												(gridA.nz+OCL.D3[2]-1)/OCL.D3[2]*OCL.D3[2]},
												{OCL.D3[0], OCL.D3[1], OCL.D3[2]}); 
                                                
    double reverse_dt = 1.0 / simu.dt;

	static cl::KernelFunctor<cl::Buffer&,
                            const cl::Buffer&, const cl::Buffer&, const cl::Buffer&, 
                            const cl::Buffer&, const cl::Buffer&, const cl::Buffer&, 
                            const cl::Buffer&, const cl::Buffer&, const cl::Buffer&, const cl::Buffer&, 
							const double ,const int, const  int, const int, const int >
                            Kernel_Sor_pipeline_getCoef_m(OCL.prg_m["src/4_2_1_SorPipeLine.cl"], "Sor_pipeline_getCoef_m");


	static cl::KernelFunctor<cl::Buffer&, const cl::Buffer&,
                            cl::Buffer&,cl::Buffer&, const cl::Buffer&,
							 const int,const double, const double, 
                             const int, const  int, const  int, const int >
                            Kernel_Sor_pipeline_MAIN(OCL.prg_m["src/4_2_1_SorPipeLine.cl"], "Sor_pipeline_MAIN");



	static cl::KernelFunctor<cl::Buffer&, cl::Buffer&,
                            cl::Buffer&, const cl::Buffer&,
                            const cl::Buffer&, const cl::Buffer&, const cl::Buffer&, 
                            const cl::Buffer&, const cl::Buffer&, const cl::Buffer&, 
                            const cl::Buffer&, const cl::Buffer&, const cl::Buffer&, 
                            const int,const double, const double, 
                            const int, const  int, const  int, const int >
                            Kernel_Sor_pipeline_MAIN_unpeeling(OCL.prg_m["src/4_2_1_SorPipeLine.cl"], "Sor_pipeline_MAIN_unpeeling");



    double omega = 1.8;
    // Kernel_Sor_pipeline_MAIN_unpeeling(config3D, OCL.pressure, OCL.LocalMax,
    //                                     OCL.iter, OCL.NEIBcell,
    //                                     OCL.T1_u, OCL.T1_v, OCL.T1_w,
    //                                     OCL.Dx, OCL.Dy, OCL.Dz,
    //                                     OCL.Dxs, OCL.Dys, OCL.Dzs,
    //                                     simu.iterMax, simu.p_criteria , omega,
    //                                     gridA.nx, gridA.ny, gridA.nz, gridA.gC);


    Kernel_Sor_pipeline_getCoef_m(config3D, OCL.SorCoef,
                                 OCL.T1_u, OCL.T1_v, OCL.T1_w,
                                 OCL.Dx, OCL.Dy, OCL.Dz,
                                 OCL.Dxs, OCL.Dys, OCL.Dzs,OCL.NEIBcell,reverse_dt,
                                 gridA.nx, gridA.ny, gridA.nz, gridA.gC);

    Kernel_Sor_pipeline_MAIN(config3D,OCL.pressure, OCL.SorCoef,
                        OCL.LocalMax, OCL.iter, OCL.NEIBcell,
                        simu.iterMax, simu.p_criteria , omega,
                        gridA.nx, gridA.ny, gridA.nz, gridA.gC);

    // std::vector<double> PP(gridA.iceltot,2);
    // OCL.queue.enqueueReadBuffer(OCL.pressure, CL_TRUE, 0,t1.p.size()*  sizeof(cl_double),PP.data());

    // for (auto v : PP ){
    //     cout << v << "| ";
    // }



}