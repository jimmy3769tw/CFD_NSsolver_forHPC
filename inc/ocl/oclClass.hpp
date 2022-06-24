#pragma once

#include "controlPanel.hpp"


// // config cl2.hpp through Macros
// #define CL_HPP_CL_1_2_DEFAULT_BUILD
// #define CL_HPP_TARGET_OPENCL_VERSION 120
// #define CL_HPP_MINIMUM_OPENCL_VERSION 120
// #define CL_HPP_ENABLE_EXCEPTIONS
// #include <CL/cl2.hpp>



// config cl2.hpp through Macros
#define CL_HPP_CL_1_2_DEFAULT_BUILD
#define CL_HPP_TARGET_OPENCL_VERSION 200
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>





class OCLstruct{

	public:

    // template<typename T2, typename T1>
	// void init(T1 PlatformName, cl_device_type type,T2 globalMemoryMB ){
	// 	// * A1 Platform
	// 	platform = getPlatform(PlatformName);
	// 	// * A2 Device
	// 	device = getDevice(platform, type, globalMemoryMB);
	// 	// * A3 ctx
	// 	ctx =  cl::Context(device);

	// 	queue = cl::CommandQueue(ctx,device);
	// 	queue11 = queue();
	// }



    // auto SetKernelProgram_from_SourceFile_m(std::vector<std::string> Source_File){
	// 	int size = Source_File.size();


	// 	for (uint i = 0 ; i < size ; i++){

	// 		prg_m[Source_File[i]] = cl::Program(ctx,LoadProgram(Source_File[i]));
	// 		cout <<"\n load " << i;

	// 		try {
	// 			prg_m[Source_File[i]].build();

	// 		} catch (cl::Error &e) {
	// 			std::cout << "filename" << Source_File[i];

	// 			std::cerr << "\n" << prg_m[Source_File[i]].getBuildInfo<CL_PROGRAM_BUILD_LOG>(device);

	// 			throw cl::Error(CL_INVALID_PROGRAM,"Failed to build kernel");
	// 		}
	// 	}
	// }


	// std::vector<cl::size_type>D1 = {256};
	// std::vector<cl::size_type>D2 = {64,64};
	// std::vector<cl::size_type>D3 = {8,8,8};

	// * prg
	// std::map<std::string, cl::Program> prg_m;



	// cl::Program ;
	// * queue 
	// cl::CommandQueue queue;
	// cl_command_queue queue11;

	// cl::Context ctx;



	// cl::Buffer<cl::Program>
	// * cl::Buffer 
	// cl::Buffer SorCoef;
	// cl::Buffer T0_u, T0_v,T0_w;
	// cl::Buffer T1_u, T1_v,T1_w;
	// cl::Buffer T3_u, T3_v,T3_w;
	// cl::Buffer Dfib_eta, Dfib_f, pressure;
	// cl::Buffer Dx, Dy, Dz;
	// cl::Buffer Dxs, Dys, Dzs, NEIBcell;
	// cl::Buffer temperature;
	// cl::Buffer LocalMax;
	// cl::Buffer iter;

	// cl::Platform platform;
	// cl::Device device;


// #ifdef BoostCompute_ON
//     compute::vector<cl_float> bo_x;
//     compute::vector<cl_float> bo_y;
//     compute::vector<cl_float> bo_z;
//     compute::device bo_device_gpu;
//     compute::context bo_ctx;
//     compute::command_queue bo_que;
//     compute::program bo_program;
// #endif




	// void set_phiBufer(
	// 	const std::vector<double> &phi,
	// 	cl::Buffer &Buffer,
	// )
	// {
	// 	// Currt_pBuf = cl::Buffer(ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
	// 	// 				phi.size() * sizeof(cl_double), phi.data());
	// }



	private:

	// // A utility function for getting a specific platform based on vendor's name
	// auto getPlatform(const std::string& vendorNameFilter) {
	//     std::vector<cl::Platform> platforms;
	//     cl::Platform::get(&platforms);
	//     for(const auto& p: platforms) {
	//         if(p.getInfo<CL_PLATFORM_VENDOR>().find(vendorNameFilter) != std::string::npos) {
	//             return p;
	//         }
	//     }
	//     throw cl::Error(CL_INVALID_PLATFORM, "No platform has given vendorName");
	// }


	// // A utility function for getting a device based on the amount of global memory.
	// auto getDevice(cl::Platform& platform, cl_device_type type, size_t globalMemoryMB) {
	//     std::vector<cl::Device> devices;
	//     platform.getDevices(type, &devices);
	//     globalMemoryMB *= 1024 * 1024; // from MB to bytes
	//     for(const auto& d: devices) {
	//         if( d.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() >= globalMemoryMB ) return d;
	//     }
	//     throw cl::Error(CL_INVALID_DEVICE, "No device has needed global memory size");
	// }



    // template<typename T>
	// inline std::string LoadProgram(T filename)
	// {
	// 	std::vector<char> temp;
	// 	fstream file;
	// 	// char buffer[200];
	// 	char ch;
	// 	file.open(filename,ios::in);
	// 	if(!file){
	// 		throw invalid_argument("The file could't open.");
	// 	}
	// 	else
	// 	{
	// 		while (file.get(ch))
	// 		{
	// 			temp.push_back(ch);
	// 		}
	// 		file.close();
	// 	}
	// 	std::string s(temp.begin(), temp.end()) ;
	// 	// cout << s << endl;  //for check
	// 	return s;
	// }

};




// ! for Boost Compute ----------------
// #ifdef BoostCompute_ON

// #define CL_TARGET_OPENCL_VERSION 120


// #include <boost/compute/container/vector.hpp

// 	// #include <boost/compute.hpp>

// 		// #include <boost/compute/system.hpp>

// 	// #include <boost/compute/system.hpp>

// #include <boost/compute/algorithm/copy.hpp>

// #include <boost/compute/utility/dim.hpp>

// namespace compute = boost::compute;

// #endif
