#include "../method/Util.h"
#include "../mrc/mrcstack.h"
#include "Process.h"
#include "opts.h"

int main(int argc, char *argv[]){
    SysInfo info;
    MPI_Init(&argc, &argv); 
    MPI_Comm_rank(MPI_COMM_WORLD, &(info.id));
    MPI_Comm_size(MPI_COMM_WORLD, &(info.procs));
    MPI_Get_processor_name(info.processor_name, &(info.namelen));
    
    options opts;
    InitOpts(&opts);
    if (GetOpts(argc, argv, &opts) <= 0)
    {
        EX_TRACE("***WRONG INPUT.\n");
        return -1;
    }

    if (info.id == 0)
    {
        PrintOpts(opts);
    }

    cudaSetDevice(opts.GPU);
    Process process;
    process.DoIt(opts, info);

    // CalcErr calcerr;
    // calcerr.DoIt(opts, txt);

    // GLWT::lwt haar;
    // haar.test();
    // haar.test_time();
    // CalcReprojPnp a;
    // a.test();
    // test(3,2,2);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize(); 

    return 0;
}