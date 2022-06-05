#!/usr/bin/env python3
import sys
from bs4 import BeautifulSoup

collectives = {
    "Bcast",
    "Barrier",
    "Allreduce",
    "Allgather",
    "Allgatherv",
    "Alltoall",
    "Alltoallv",
    "Reduce",
    "Gather",
    "Gatherv",
    "Scan",
    "Exscan",
    "Scatter",
    "Scatterv",
    "Reduce_scatter",
    "Iallreduce",
    "Ibcast",
    "Ialltoall",
    "Iscatter",
    "Ibarrier",
}

prim_names={
    "Send",
    "Recv",
    "Isend",
    "Irecv",
    "Sendrecv",
    "Bcast",
    "Barrier",
    "Allreduce",
    "Allgather",
    "Allgatherv",
    "Alltoall",
    "Alltoallv",
    "Reduce",
    "Gather",
    "Gatherv",
    "Scan",
    "Exscan",
    "Scatter",
    "Scatterv",
    "Reduce_scatter",
    "Waitall",
    "Wait",
    "Waitany",
    "Test",
    "Iallreduce",
    "Ibcast",
    "Ialltoall",
    "Iscatter",
    "Ibarrier",
    "Testany"
}
if len(sys.argv) < 2 :
    print( "Insert file to read" )
    exit(1)

else:
    with open(sys.argv[1], 'rb') as f:
        Bs_data = BeautifulSoup(f.read(),"xml")
    # mpi_rank = Bs_data.find_all('func')
    # print (len(prim_names))
    for primitive in prim_names:
        prim = "MPI_" + primitive
        mpi_prim = Bs_data.find_all('func',{'name':prim})
        sendcount = 0
        for item in mpi_prim:
            if primitive in collectives:
                tmp = int( item.get("count") )
                if tmp > sendcount:
                    sendcount = tmp
            else:
                sendcount += int( item.get("count") )
        if ( sendcount > 0 ):
            print(prim,"calls = %d"%sendcount)
