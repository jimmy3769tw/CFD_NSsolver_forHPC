#include "controlPanel.hpp"
#pragma once 
#ifdef OCL_ON

void SorPopeLine_OCL_pre(
    OCLstruct &OCL,
    shareMenory& ShareM,
    SORcoefficient& Sor,
    simuClass& simu,
    velocity& T1,
    pressure& t1,
    calDomain& Lo,
    grid& gridA
);


void SorPipeLine_OCL(
    OCLstruct &OCL,
    shareMenory& ShareM,
    SORcoefficient& Sor,
    simuClass& simu,
    velocity& T1,
    pressure& t1,
    calDomain& Lo,
    grid& gridA
);



#endif
