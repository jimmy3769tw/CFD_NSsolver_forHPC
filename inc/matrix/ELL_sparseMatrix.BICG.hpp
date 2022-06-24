    //  ! ----------------------------------- un pre bicg -----------------------------------
    // A -- this mat (#rows, #cols)
    // B -- impute vector (n)
    using std::vector;
    template<typename T>
    inline std::pair<int, double> ELL_matrix<T>::npc_bicgstab(
        std::vector <double> &B,
        std::vector <double> &x
    ){


        auto m = this -> rows_;
        auto n = this -> cols_;

        if ( m != n){
            throw std::invalid_argument(" This mat don't suport BICG");
        }

        if (m != x.size() ){
            std::cout <<  m << ", " << x.size() << std::endl;
            throw std::invalid_argument(" m != x.size()_!!");
        }

        if (n != B.size() ){
            throw std::invalid_argument(" x.n != B.size()!!");
        }

        static vector<double> r(m);
        std::fill(r.begin(), r.end(), 0.0);

        std::fill(x.begin(),x.end(),0.0);

        // ! npc BiCGSTAB

        // //* 1.r0 = b − Ax0 // I.C. r @ i = 0
        auto p = multiply(x);
        // p is present the temp ;
        std::transform(B.begin(), B.end(), p.begin(), r.begin(), std::minus<double>()); 

        // //*2 r̂0 = r0 // r_hat[i] = const
        auto r_hat = r;

        // //*3 ρ0 = α = ω0 = 1
        double rho = 1.0, alpha = 1.0, omeaga = 1.0;

        // * 4. v0 = p0 = 0  //I.C.
        std::fill(p.begin(), p.end(), 0.0);

        static vector<double> v(m);
        // !-----------------
        static vector<double> s(m);
        static vector<double> t(m);
        static vector<double> h(m);
        static vector<double> rm(m);
        static vector<double> xm(m);
        // !-----------------
        std::fill(v.begin(), v.end(), 0.0);
        std::fill(s.begin(), s.end(), 0.0);
        std::fill(h.begin(), h.end(), 0.0);
        std::fill(t.begin(), t.end(), 0.0);
        std::fill(xm.begin(), xm.end(), 0.0);
        std::fill(rm.begin(), rm.end(), 0.0);

        double accurate = 1.0e-10;
        

        auto norm0 = std::inner_product(B.begin(), B.end(),B.begin() ,0.0);
        const auto tolerance = std::pow((this->accurate_bicg_ * norm0), 2);
        int iter_max = this->iter_max_bicg_;

       
        // //! 5. For i = 1, 2, 3, 
        int iter= 1; double error;
        for (double beta = 0.0; iter < iter_max ;++iter){
            xm.assign(x.begin(), x.end());
            // * 1.ρi = (r̂0, ri−1)
            double rhom = rho;
            rho = std::inner_product(r_hat.begin(), r_hat.end(),r.begin() ,0.0);
            // * ---------------------

            // * 2.β = (ρi/ρi−1)(α/ωi−1)
            beta = (rho / rhom) * (alpha / omeaga);
            
            // * ---------------------


            // * 3.pi = r_{i−1} + β(p_{i−1} − ω_{i−1} * v_{i−1})
            for (int i = 0 ; i < m ;++i ){
                p[i] = r[i] + beta * (p[i] - omeaga * v[i]);
            }
            // * ---------------------

            // * 4.vi = A p_i
            v = this->multiply(p);
            // * ---------------------


            // * 5.α = ρi/(r̂0, vi)
            alpha = rho / std::inner_product(r_hat.begin(), r_hat.end(),v.begin(),0.0);
            // * ---------------------

            // * 6. h = x_{i−1} + α pi
            for (int i = 0 ; i < m ;++i ){
                h[i] =  xm[i]+ alpha * p[i];
            }   
            // * ---------------------

            // // ! 7.If h is accurate enough, then set xi = h and quit
            // error = std::inner_product( h.begin(), h.end(), 
            //                             x.begin(), 0.0);
            // if (error < tolerance){
            //     x.swap(h);
            //     break;
            // }


            // * 8.s = ri−1 − αvi
            for (int i = 0 ; i < m ;++i ){
                s[i] = r[i] - alpha * v[i];
            }



            // * 9. t = As
            t = this->multiply(s);
            // * ---------------------

            // * 10. ωi = (t, s)/(t, t)

            omeaga = std::inner_product(t.begin(), t.end(),s.begin(),0.0)
                                                /
                     std::inner_product(t.begin(), t.end(),t.begin(),0.0);
            // * ---------------------

            // * 11 xi = h + ωi*s
            for (int i = 0 ; i < m ;++i ){
                x[i] = h[i] + omeaga * s[i];
            }
            // * ---------------------
            // ! 12. If xi is accurate enough, then quit
            error = std::inner_product( xm.begin(), xm.end(), 
                                        x.begin(), 0.0, 
                                        std::plus<>(),
                                        [](auto &a, auto &b){return std::pow((a-b),2);});

            if (error < tolerance){
                break;
            }
            // * ---------------------


            // * 13 ri = s − ω_i t
            for (int i = 0 ; i < m ;++i ){
                r[i] = s[i] - omeaga * t[i];
            }
            // * ---------------------
        }

        return std::make_pair( iter, std::sqrt(error) / norm0);
    }

// ---------------------------------- omp -------------------------

    template<typename T>
    inline std::pair<int, double> ELL_matrix<T>::npc_bicgstab_omp(
        std::vector <double> &B,
        std::vector <double> &x
    ){

        auto m = this -> rows_;
        auto n = this -> cols_;
        

        if ( m != n){
            throw std::invalid_argument(" This mat don't suport BICG");
        }

        if (m != x.size() ){
            throw std::invalid_argument(" x.size() != this->.row_!");
        }

        if (n != B.size() ){
            throw std::invalid_argument(" x.size() != this->.row_!");
        }

        // int TID = omp_get_thread_num();
        // int Tsize = omp_get_max_threads();
        // auto [q,r] = std::div(m, Tsize);


        static vector<double> r(m);
        static vector<double> p(m);

        #pragma omp parallel for schedule(static) if(ompBool)
        for (int i = 0 ; i < m ;++i )
            r[i] = 0.0;

        #pragma omp parallel for schedule(static) if(ompBool)
        for (int i = 0 ; i < m ;++i )
            p[i] = 0.0;

        #pragma omp parallel for schedule(static) if(ompBool)
        for (int i = 0 ; i < m ;++i )
            x[i] = 0.0;



        // ! npc BiCGSTAB

        // //* 1.r0 = b − Ax0 // I.C. r @ i = 0
        multiply_omp(x,p);

        #pragma omp parallel for schedule(static) if(ompBool)
        for (int i = 0 ; i < m ;++i )
            r[i] = B[i] - p[i];


        // //*2 r̂0 = r0 // r_hat[i] = const
        auto r_hat = r;

        // //*3 ρ0 = α = ω0 = 1
        double rho = 1.0, alpha = 1.0, omeaga = 1.0;

        // * 4. v0 = p0 = 0  //I.C.
        static vector <double> v(m);

        // !-----------------
        static vector<double> s(m);
        static vector<double> t(m);
        static vector<double> h(m);
        static vector<double> rm(m);
        static vector<double> xm(m);
        // !-----------------



        #pragma omp parallel for schedule(static) if(ompBool)
        for (int i = 0 ; i < m ;++i ){
            p[i]    = 0.0;
            v[i]    = 0.0;
            t[i]    = 0.0;
            s[i]    = 0.0;
            h[i]    = 0.0;
            rm[i]   = 0.0;
            xm[i]   = 0.0;
        }
        

        auto norm0 = inner_product_omp_for(B, B);
        auto tolerance = std::pow(this->accurate_bicg_ * norm0, 2);
        int iter_max = this->iter_max_bicg_;

        double error;

        // //! 5. For i = 1, 2, 3, 
        int iter = 1 ;
        for ( double beta = 0.0; iter < iter_max ;++iter){

            #pragma omp parallel for schedule(static) if(ompBool)
            for (int i = 0 ; i < m ;++i ){
                xm[i]= x[i];
            }

            // * 1.ρi = (r̂0, ri−1)
            double rhom = rho;
            rho = inner_product_omp_for(r_hat, r);
            // * ---------------------

            // * 2.β = (ρi/ρi−1)(α/ωi−1)
            beta = (rho / rhom) * (alpha / omeaga);
            // * ---------------------

            // * 3.pi = r_{i−1} + β(p_{i−1} − ω_{i−1} * v_{i−1})
            #pragma omp parallel for schedule(static) if(ompBool)
            for (int i = 0 ; i < m ;++i ){
                p[i] = r[i] + beta * (p[i] - omeaga * v[i]);
            }
            // * ---------------------
            // * 4.vi = A p_i
            
            multiply_omp(p, v);
            // * ---------------------


            // * 5.α = ρi/(r̂0, vi)
            double temp = inner_product_omp_for(r_hat, v);
            alpha = rho / temp;
            // * ---------------------

            // * 6. h = x_{i−1} + α pi
            #pragma omp parallel for schedule(static) if(ompBool)
            for (int i = 0 ; i < m ;++i ){
                h[i] =  xm[i]+ alpha * p[i];
            }   
            // * ---------------------


            // ! 7.If h is accurate enough, then set xi = h and quit
            error = inner_product_omp_for(h, x);
            if (error < tolerance){
                x.swap(h);
                break;
            }
        

            // * 8.s = ri−1 − αvi


            #pragma omp parallel for schedule(static) if(ompBool)
            for (int i = 0 ; i < m ;++i ){
                s[i] = r[i] - alpha * v[i];
            }


            // * ---------------------

            // * 9. t = As
            double time_1 = omp_get_wtime();

            multiply_omp(s, t);

            timer[1] +=  omp_get_wtime() - time_1; 
            // * ---------------------

            // * 10. ωi = (t, s)/(t, t)

            omeaga =    inner_product_omp_for(t,s)
                                    /
                        inner_product_omp_for(t,t);
            // * ---------------------

            // * 11 xi = h + ωi*s
            #pragma omp parallel for schedule(static) if(ompBool)
            for (int i = 0 ; i < m ;++i ){
                x[i] = h[i] + omeaga * s[i];
            }
            // * ---------------------

            // ! 12. If xi is accurate enough, then quit
            error = this->accumulate_omp(xm, x);

            if (error < tolerance){
                break;
            }
            // * ---------------------


            // * 13 ri = s − ω_i t
            #pragma omp parallel for schedule(static) if(ompBool)
            for (int i = 0 ; i < m ;++i ){
                r[i] = s[i] - omeaga * t[i];
            }


            // * ---------------------
        }
        std::cout << "inner product" << timer[0] << std::endl;
        std::cout << "multiply_omp " << timer[1] << std::endl;
        return std::make_pair( iter, std::sqrt(error) / norm0);
    }

    #define ELL_MPI_ON
    #ifdef ELL_MPI_ON
    template<typename T>
    inline std::pair<int, double> ELL_matrix<T>::npc_bicgstab_mpi(
        std::vector <double> &B,
        std::vector <double> &x
    ){
        auto m = this -> rows_;
        auto n = this -> cols_;

        if ( m != n){
            throw std::invalid_argument(" This mat don't suport BICG");
        }

        if (m != x.size() ){
            throw std::invalid_argument(" x.size() != this->.row_!");
        }

        if (n != B.size() ){
            throw std::invalid_argument(" x.size() != this->.row_!");
        }

        static vector<double> r(m);
        for (int i = 0 ; i < m ;++i )
            r[i] = 0.0;
    
        std::fill(x.begin(),x.end(),0.0);

        // ! npc BiCGSTAB

        // //* 1.r0 = b − Ax0 // I.C. r @ i = 0
        static vector <double> p(m);
        p = multiply_mpi(x);
        // p is present the temp ;
        std::transform(B.begin(), B.end(), p.begin(), r.begin(), std::minus<double>()); 

        // //*2 r̂0 = r0 // r_hat[i] = const
        auto r_hat = r;

        // //*3 ρ0 = α = ω0 = 1
        double rho = 1.0, alpha = 1.0, omeaga = 1.0;

        // * 4. v0 = p0 = 0  //I.C.
    
        static vector <double> v(m);
        // !-----------------
        static vector<double> s(m);
        static vector<double> t(m);
        static vector<double> h(m);
        static vector<double> rm(m);
        static vector<double> xm(m);
        // !-----------------
        // !-----------------
        for (int i = 0 ; i < m ;++i ){
            p[i]    = 0.0;
            v[i]    = 0.0;
            t[i]    = 0.0;
            s[i]    = 0.0;
            h[i]    = 0.0;
            rm[i]   = 0.0;
            xm[i]   = 0.0;
        }
        
        double accurate = 1.0e-10;
        
        auto norm0 = this->inner_product_mpi(B, B);

        auto tolerance = std::pow((this->accurate_bicg_ * norm0), 2);
        int iter_max = this->iter_max_bicg_;

       
        // //! 5. For i = 1, 2, 3, 
        int iter= 1; double error;
        for (double beta = 0.0; iter < iter_max ;++iter){
            xm.assign(x.begin(), x.end());
            // * 1.ρi = (r̂0, ri−1)
            double rhom = rho;
            rho = this->inner_product_mpi(r_hat, r);
            // * ---------------------

            // * 2.β = (ρi/ρi−1)(α/ωi−1)
            beta = (rho / rhom) * (alpha / omeaga);
            // * ---------------------


            // * 3.pi = r_{i−1} + β(p_{i−1} − ω_{i−1} * v_{i−1})
            for (int i = this->m_begin_ ; i < this->m_endof_ ; ++i ){
                p[i] = r[i] + beta * (p[i] - omeaga * v[i]);
            }
            // * ---------------------

            // * 4.vi = A p_i
            v = this->multiply_mpi(p);
            // * ---------------------


            // * 5.α = ρi/(r̂0, vi)
            alpha = rho / this->inner_product_mpi(r_hat, v);
            // * ---------------------

            // * 6. h = x_{i−1} + α pi
            for (int i = this->m_begin_ ; i < this->m_endof_ ; ++i ){
                h[i] =  xm[i] + alpha * p[i];
            } 
            // * ---------------------

            // ! 7.If h is accurate enough, then set xi = h and quit
            // error = std::inner_product( h.begin(), h.end(), 
            //                             x.begin(), 0.0);
            // if (error < tolerance){
            //     x.swap(h);
            //     break;
            // }


            // * 8.s = ri−1 − αvi
            for (int i = this->m_begin_ ; i < this->m_endof_ ; ++i ){
                s[i] = r[i] - alpha * v[i];
            }



            // * 9. t = As
            t = this->multiply_mpi(s);
            // * ---------------------

            // * 10. ωi = (t, s)/(t, t)

            omeaga = this->inner_product_mpi(t,s)
                                /
                     this->inner_product_mpi(t,t);
            // * ---------------------

            // * 11 xi = h + ωi*s
            for (int i = this->m_begin_ ; i < this->m_endof_ ; ++i ){
                x[i] = h[i] + omeaga * s[i];
            }
            // * ---------------------

            // ! 12. If xi is accurate enough, then quit
            error = this->accumulate_pow_mpi(xm, x);


            if (error < tolerance){
                break;
            }
            // * ---------------------


            // * 13 ri = s − ω_i t
            for (int i = this->m_begin_ ; i < this->m_endof_ ; ++i ){
                r[i] = s[i] - omeaga * t[i];
            }
            // * ---------------------
        }
        this -> allocate_vector_mpi(x);

        return std::make_pair( iter, std::sqrt(error) / norm0);
}
#endif
//  ! ----------------------------------- end of un pre bicg -----------------------------------
