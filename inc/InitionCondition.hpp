#pragma once
#include "controlPanel.hpp"


void Locality(
    simuClass& simu,
    calDomain& Lo,
    velocity& T0,
    velocity& T1,
    velocity& T3,
    pressure& t1,
    grid& gridA
)
{
   const auto [nx, ny, nz, gC] = gridA.nxyzgC; 

   for (size_t i = Lo.i_begin  ; i < Lo.i_endof ; ++i )
   for (size_t j = Lo.j_begin  ; j < Lo.j_endof ; ++j )
   for (size_t k = Lo.k_begin  ; k < Lo.k_endof ; ++k )
   {
      const int icel = i*nz*ny + j*nz + k;
      T0.u[icel] = 0.;         T0.v[icel] = 0.;         T0.w[icel] = 0.;
      T1.u[icel] = 0.;         T1.v[icel] = 0.;         T1.w[icel] = 0.;
      T3.u[icel] = 0.;         T3.v[icel] = 0.;         T3.w[icel] = 0.;
      t1.p[icel] = 0.;

      if (i == 2) {
         T0.u[icel - nz*ny] = 0.;
         T0.u[icel - 2*nz*ny] = 0.;
      }
      if (i == nx -3 ) {
         T0.u[icel + nz*ny] = 0.;
         T0.u[icel + 2*nz*ny] = 0.;
      }
      if (i == 2) {
         T0.u[icel-nz] = 0.;
         T0.u[icel-2*nz] = 0.;
      }
      if (i == nx -3 ) {
         T0.u[icel + nz] = 0.;
         T0.u[icel + 2*nz] = 0.;
      }
      if (k == 2) {
         T0.u[icel - 1] = 0.;
         T0.u[icel - 2] = 0.;
      }
      if (k == nz -3 ) {
         T0.u[ icel + 1] = 0.;
         T0.u[ icel + 2] = 0.;
      }
   }
}

