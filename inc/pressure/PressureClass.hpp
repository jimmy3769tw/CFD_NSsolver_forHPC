#pragma once
#include "controlPanel.hpp"

// ---------------------------------------
#include "matrix/ELL_sparseMatrix.hpp"
#include "matrix/ELL_sparseMatrix0.hpp"
// ---------------------------------------

// ---------------------------------------
#include "matrix/CSR_sparseMatrix.hpp"
#include "matrix/CSR_sparseMatrix0.hpp"
// ---------------------------------------

// ---------------------------------------
#include "matrix/SPE_sparseMatrix.hpp"
#include "matrix/SPE_sparseMatrix1.hpp"
#include "matrix/SPE_sparseMatrix0.hpp"
// ---------------------------------------

using ELL_type = typename mat::ELL_matrix<double>;
using CSR_type = typename mat::CSR_matrix<double>;
using SPE_type = typename mat::SPE_matrix1<double>;

struct MxClass
{
    ELL_type matA_ell;
    CSR_type matA_csr;
    SPE_type matA_spe;

    std::vector <double> matB;

	std::vector <double> X_result;

	//* Eigen ================================
#ifdef EIGEN_ON
    Eigen::SparseMatrix<double, Eigen::RowMajor> matA_Eigen;

	Eigen::VectorXd x_Eigen;

    Eigen::VectorXd matB_Eigen;

	Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor> >  solver_Eigen;
#endif
	//* Eigen ================================
};


struct SORcoefficient{
    std::vector<double> cf;
};