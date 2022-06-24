#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>
#include <grid/structureGrid.hpp>
class velocity{

	public:
		// -------------------------------
		std::vector <double> u, v, w,
							 su, sv, sw,
							 Viseff, c ;
		// -------------------------------

		template< typename T>
		std::vector<double>& U(T whichDim){
			if (whichDim < 0 || whichDim > 2) throw std::runtime_error("out of dimensions!");
			if (whichDim == 0){return u;}
			else if (whichDim == 1){return v;}

			return w;
		}

		template< typename T>
		std::vector<double>& sU(T whichDim){
			if (whichDim < 0 || whichDim > 2) throw std::runtime_error("out of dimensions!");
			if (whichDim == 0){return su;}
			else if (whichDim == 1){return sv;}
			return sw;
		}

		bool resize_s(size_t n){
			// -----------
			su.resize(n);
			sv.resize(n);
			sw.resize(n);			
			// -----------
			return true;
		}


		bool resize(size_t n){

			for (size_t i = 0; i < 3; i++){
				U(i).resize(n);
				sU(i).resize(1);
			}

			return true;
		}

		bool iniU(double ui, double vi, double wi){
			std::fill(u.begin(), u.end(), ui);
			std::fill(v.begin(), v.end(), vi);
			std::fill(w.begin(), w.end(), wi);
			return true;
		}


		inline void fill_(std::vector<double> &v, double val){

			#pragma omp for schedule(static)
			for(size_t i = 0; i < v.size() ; ++i){
				v[i] = val;
			}
		}

		bool iniU_omp(double ui, double vi, double wi){

				fill_(u,ui);
				fill_(v,vi);
				fill_(w,wi);
			return true;
		}


		void copy_seq(const velocity &b){
			u = b.u;
			v = b.v;
			w = b.w;
		}

		void copy_omp(const velocity &b){

			#pragma omp parallel 
			{

				#pragma omp for simd schedule(static)
				for (size_t i = 0 ; i < u.size() ; ++i)
					u[i] = b.u[i];

				#pragma omp for simd  schedule(static)
				for (size_t i = 0 ; i < v.size() ; ++i)
					v[i] = b.v[i];

				#pragma omp for simd  schedule(static)
				for (size_t i = 0 ; i < w.size() ; ++i)
					w[i] = b.w[i];

			}
		}


		void coutDifferent(velocity &b){

			std::vector<double>accumulate(3,0.);
			for (size_t j = 0 ; j < 3 ; ++j){
				for (size_t i = 0 ; i < U(j).size() ; ++i){
					accumulate[j] += std::abs(U(j)[i] - b.U(j)[i]);
				}
			}
			std::cout << "\nDifferent" ;
			for (auto &x:accumulate){
				std::cout << x << ", "; 
			}
		}


		template< typename T>
		void showVel(T whichDim, grid thegird){
			thegird.showVel(U(whichDim));
		}


		template< typename T>
		void showVelt(T whichDim, grid thegird){
			thegird.showVelt(U(whichDim));
		}

	private:


};




class pressure{
	public:
		std::vector<double> p;
		std::vector<double> T;

		bool init_p(double pi){

			std::fill(p.begin(), p.end(), pi);
			return true;
		}

		bool iniP(double pi){
			std::fill(p.begin(), p.end(), pi);
			return true;
		}

	private:

};

struct DfibArray{
    std::vector<double> f, 
						eta, 
						Val_sumz,
						Val_sumz_sumy, 
						ValSum;

    double cylinderDimension;

	std::vector<double> cylinderCenter;
};


struct shareMenory{

	double pChangeMax;

	double Out;

	// * Check Max val 
	double uDif_Max;

	double vDif_Max;
	
	double wDif_Max;
};
