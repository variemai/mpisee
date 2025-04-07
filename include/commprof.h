/*****************************************************************************\
* This source file is part of mpisee
* Copyright (C) 2024  Ioannis Vardas - vardas@par.tuwien.ac.at
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <https://www.gnu.org/licenses/>
******************************************************************************/
#ifndef COMMPROF_H_
#define COMMPROF_H_

#include "mpi.h"
#include <unordered_map>
#include "utils.h"
#include <vector>

#define MAX_DIMS 8
extern int prof_enabled;
extern prof_attrs **local_data;
// extern prof_attrs **local_comms;
extern std::vector<MPI_Comm> comms_table;
extern std::unordered_map<MPI_Request, MPI_Comm> requests_map;
extern std::unordered_map<MPI_Win, MPI_Comm> comm_map;
extern std::vector<std::pair<prof_meta_pair*, MPI_Group>> free_array;

extern int local_cid;
extern int my_coms;
extern int keys[2]; // keyval[0]  contains metadata
                      // keyval[1]  contains profiling data
extern MPI_Comm dummy_comm;

extern "C" {
int
namekey(void);
}
extern "C" {
int
namedel(MPI_Comm comm, int keyval, void *attr, void *s);
}

extern "C" {
void
profile_this(MPI_Comm comm, int64_t count,MPI_Datatype datatype,int prim,
             double t_elapsed,int v);
}

extern "C" {
int
win_namekey(void);
}
#endif
