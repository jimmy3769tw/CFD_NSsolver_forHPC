#pragma once

#include "controlPanel.hpp"

#ifdef AMGCL_ON

#include <amgcl/amg.hpp>

#include <amgcl/io/mm.hpp>

#include <amgcl/backend/builtin.hpp>  // for build OpenMP
#include <amgcl/adapter/crs_tuple.hpp> // for build crs_tuple


// ! make_solver
// 1.1
#include <amgcl/make_solver.hpp>
// 1.2
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/coarsening/plain_aggregates.hpp>
#include <amgcl/coarsening/rigid_body_modes.hpp>

#include <amgcl/profiler.hpp>

#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/relaxation/chebyshev.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/bicgstabl.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/solver/lgmres.hpp> // 好像可以

  typedef amgcl::make_solver<
             amgcl::amg< 
/*PBackend */ amgcl::backend::builtin<double>,
              amgcl::coarsening::smoothed_aggregation,
              amgcl::relaxation::ilu0>,
/*SBackend */amgcl::solver::bicgstab<amgcl::backend::builtin<double> >
  > SolverBuiltin;




#ifdef EIGEN_ON
#include <amgcl/backend/eigen.hpp>

// #include <amgcl/coarsening/smoothed_aggregation.hpp>
// #include <amgcl/relaxation/spai0.hpp>
// #include <amgcl/make_solver.hpp>
// #include <amgcl/solver/bicgstab.hpp>

  typedef amgcl::make_solver<
      amgcl::amg<
          amgcl::backend::eigen<double>,                      
          amgcl::coarsening::smoothed_aggregation,
          amgcl::relaxation::spai0
          >,
      amgcl::solver::bicgstab<amgcl::backend::eigen<double> > //SBackend
      > SolverEigen;



#endif

#endif