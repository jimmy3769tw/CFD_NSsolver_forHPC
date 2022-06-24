#ifndef __MPI_TOOLS_H__
#define	__MPI_TOOLS_H__

#include <mpi.h>
#include <vector>

namespace yakutat::backend
{
    
    class mpi_tools
    {
        using self_type = mpi_tools;

        public:
            using size_type = size_t; 
            using value_type = double;
            using array_size_type = std::vector<size_type>;
            using array_value_type = std::vector<value_type>;

            mpi_tools();
            virtual ~mpi_tools();
            bool MPI_division(size_type size, size_type id, size_type proc);
            bool MPI_vector_collect(array_value_type & a, size_type itag);
        
            self_type & operator = (self_type & world);



            //MPI variables
            size_type myid{0}, nproc{0};
            size_type length;
            size_type start{0}, end{0}, count{0};
            array_size_type start_list, end_list, count_list;

        private:
            MPI_Status istat_[8];
            size_type master_{0};
            size_type istart_{0};
            size_type icount_{0};  

            bool copy(self_type & world);

    };

    inline mpi_tools::mpi_tools() {}
    inline mpi_tools::~mpi_tools() {}




    inline bool mpi_tools::MPI_division(size_type size, size_type id, size_type proc)
    {
        length = size;
        myid = id;
        nproc = proc;
      
        size_type Xdv = length / nproc;				    
        size_type Xr = length - Xdv * nproc;
        

        start_list.reserve(nproc);
        end_list.reserve(nproc);
        count_list.reserve(nproc);

        for(size_type i=0; i<nproc; ++i)
        {
            if(i < Xr)
            {
                start_list.push_back( i * (Xdv + 1) + 0 );			
                end_list.push_back( start_list[i]+ Xdv );
            }
            else
            {
                start_list.push_back( i * Xdv + Xr + 0 );		

                end_list.push_back( start_list[i] + Xdv - 1 );		

            }

            count_list.push_back( end_list[i] - start_list[i] + 1 );	

        }
        start = start_list[myid];
        end = end_list[myid];
        count = count_list[myid];
        

        return true;
    }


    inline bool mpi_tools::MPI_vector_collect(array_value_type & a, size_type itag)
    {
        icount_ = count;
        istart_ = start;

        if( myid > master_ )
        {
            MPI_Send(   (void *)&a[istart_], icount_, MPI_DOUBLE, master_, itag, MPI_COMM_WORLD  );
        }
        else if(myid == master_)
        {
            for(size_type i = 1 ; i < nproc ; ++i)
            {
                MPI_Recv(   (void *)&a[ start_list[i] ], icount_, MPI_DOUBLE, i, itag, MPI_COMM_WORLD, istat_  );
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        icount_ = length;
        MPI_Bcast(  (void *)&a[0], icount_, MPI_DOUBLE, master_, MPI_COMM_WORLD   );

        return true;
    }

    inline mpi_tools::self_type & mpi_tools::operator = (mpi_tools & world)
    {
        copy(world);
        return *this;
    }


    inline bool mpi_tools::copy(mpi_tools & world)
    {
        myid = world.myid;
        nproc = world.nproc;
        length = world.length;
        start = world.start;
        end = world.end;
        count = world.count;
        start_list.assign(world.start_list.begin(), world.start_list.end());
        end_list.assign(world.end_list.begin(), world.end_list.end());
        count_list.assign(world.count_list.begin(), world.count_list.end());


        return true;
    }

} // end namespace yakutat::backend

#endif
