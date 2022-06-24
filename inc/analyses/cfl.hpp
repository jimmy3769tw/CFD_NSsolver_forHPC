#pragma once
// simuvariables
bool get_Cfl(
    velocity& vel,
    calDomain& Lo,
    grid& gA,
    simuClass& simu
){

    std::vector<double> c(gA.iceltotCal,0.0);

    simu.U_max.resize(3, -1e10);
    simu.U_min.resize(3,  1e10);

    #pragma omp parallel firstprivate(Lo)
    {
        int d = 0;

        auto [nx, ny, nz , gC] = gA.nxyzgC;

        if(Lo.i_endof == nx-gA.gC){d = 1;}
        else {d = 0;}

        #pragma omp for nowait
        for (auto i = Lo.i_begin ; i < Lo.i_endof-d ; ++i )
        for (auto j = Lo.j_begin ; j < Lo.j_endof   ; ++j )
        for (auto k = Lo.k_begin ; k < Lo.k_endof   ; ++k )
        {
            auto ic = gA.icel(i,j,k);
            c[gA.icelCal(i,j,k)] += vel.u[ic] / gA.Dxs[i];
            simu.U_max[0] = std::max(simu.U_max[0], vel.u[ic] );
            simu.U_min[0] = std::min(simu.U_min[0], vel.u[ic] );
        }


        if(Lo.j_endof == ny-gA.gC) {d = 1;}
        else {d = 0;}
            

        #pragma omp for nowait
        for (auto i = Lo.i_begin ; i < Lo.i_endof   ; ++i )
        for (auto j = Lo.j_begin ; j < Lo.j_endof-d ; ++j )
        for (auto k = Lo.k_begin ; k < Lo.k_endof   ; ++k )
        {
            auto ic = gA.icel(i,j,k);
            c[gA.icelCal(i,j,k)] += vel.v[ic] / gA.Dys[j];
            simu.U_max[1] = std::max(simu.U_max[1], vel.v[ic] );
            simu.U_min[1] = std::min(simu.U_min[1], vel.v[ic] );
        }

        if(Lo.k_endof == nz-gA.gC){d = 1;}
        else {d = 0;}

        #pragma omp for nowait
        for (auto i = Lo.i_begin ; i < Lo.i_endof   ; ++i )
        for (auto j = Lo.j_begin ; j < Lo.j_endof   ; ++j )
        for (auto k = Lo.k_begin ; k < Lo.k_endof-d ; ++k )
        {
            auto ic = gA.icel(i,j,k);
            c[gA.icelCal(i,j,k)] += vel.w[ic] / gA.Dzs[k];
            simu.U_max[2] = std::max(simu.U_max[2], vel.w[ic] );
            simu.U_min[2] = std::min(simu.U_min[2], vel.w[ic] );
        }

    }

    
    auto max_c = std::max_element(c.begin(), c.end());
    double cfl = *max_c * simu.dt;
    if (simu.PID == 0)
    cout << std::endl   <<  "[0;1;34;94m" 
                        << "[CFL] :"<< cfl <<  ", Max" <<  simu.cfl_max 
                        << "[0m" << std::endl;

    simu.cfl_max = std::max(simu.cfl_max, cfl);



    return true;
}

