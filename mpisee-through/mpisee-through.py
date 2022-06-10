#!/usr/bin/env python3
"""
Module documentation.
"""

# Imports
import sys
import csv
from operator import add
import argparse
#from art import *

from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE, SIG_DFL)

# Global variables
RED = "\033[1;31m"
BLUE = "\033[1;34m"
CYAN = "\033[1;36m"
GREEN = "\033[0;32m"
RESET = "\033[0;0m"
BOLD = "\033[;1m"
REVERSE = "\033[;7m"

# Class declarations

# Function declarations


def compact_proc_list(proc_list):
    proc_list.sort()
    proc_str = ""
    start_proc = proc_list[0]
    for i in range(1, len(proc_list)):
        if proc_list[i - 1] + 1 != proc_list[i]:
            proc_str += (
                f"{start_proc}-{proc_list[i-1]}"
                if start_proc != proc_list[i - 1]
                else f"{start_proc}"
            )
            proc_str += ", "
            start_proc = proc_list[i]

    proc_str += (
        f"{start_proc}-{proc_list[-1]}"
        if start_proc != proc_list[-1]
        else f"{start_proc}"
    )

    return proc_str


def print_mapping(mapping):
    mapping[0] = mapping[0].split(":")[1].strip()
    node_to_procs = {}
    for proc in mapping:
        proc_rank, node = proc.split()
        proc_rank = int(proc_rank)
        if node in node_to_procs:
            node_to_procs[node].append(proc_rank)
        else:
            node_to_procs[node] = [proc_rank]

    for node, proc_list in node_to_procs.items():

        print(f"{node}: {compact_proc_list(proc_list)}")

    print(end="\n")


def prepare_data(file_path):
    with open(file_path, "r") as file:
        csv_file = csv.reader(file)
        for i in range(7): #skip those lines for now
            next(csv_file)
        mapping = next(csv_file)
        next(csv_file)                    #skip also this line
        index_to_colname = next(csv_file)
        colname_to_index = {
            index_to_colname[i]: i for i in range(0, len(index_to_colname))
        }
        raw_data = []
        for row in csv_file:
            row_parsed = []
            for j in range(0, len(row)):
                row_parsed.append(
                    float(row[j]) if colname_to_index["Comm"] != j else row[j]
                )

            raw_data.append(row_parsed)

    data_groupBy_comm, comm_to_procs = groupBy_comm(raw_data, colname_to_index)

    table = create_table(data_groupBy_comm, colname_to_index, index_to_colname)

    return table, mapping, comm_to_procs


def groupBy_comm(data, colname_to_index):
    data_groupBy_comm = {}
    comm_to_procs = {}
    for row in data:
        comm = row[colname_to_index["Comm"]]
        row[colname_to_index["Comm"]] = ""
        if comm in data_groupBy_comm:
            data_groupBy_comm[comm].append(row)
            comm_to_procs[comm].append(int(row[colname_to_index["Rank"]]))
        else:
            data_groupBy_comm[comm] = [row]
            comm_to_procs[comm] = [int(row[colname_to_index["Rank"]])]

    return data_groupBy_comm, comm_to_procs


def create_table(data_groupBy_comm, colname_to_index, index_to_colname):
    table = []

    time_indexes = [kv[1] for kv in colname_to_index.items() if kv[0].endswith("Time")]

    for comm in data_groupBy_comm:
        maximum = list(data_groupBy_comm[comm][0])
        minimum = list(data_groupBy_comm[comm][0])
        cum_sum = list(data_groupBy_comm[comm][0])

        for i in range(1, len(data_groupBy_comm[comm])):
            row = data_groupBy_comm[comm][i]
            maximum = list(map(max, maximum, row))
            minimum = list(map(min, minimum, row))
            cum_sum = list(map(add, cum_sum, row))

        average = list(
            map(
                lambda x: x / len(data_groupBy_comm[comm])
                if not isinstance(x, str)
                else x,
                cum_sum,
            )
        )

        comm_size = int(data_groupBy_comm[comm][0][colname_to_index["Size"]])
        for i in time_indexes:
            call = "_".join(index_to_colname[i].split("_")[:-1])
            call_mean, call_min, call_max = average[i], minimum[i], maximum[i]
            call_tot_vol = int(cum_sum[colname_to_index[call + "_Volume"]])
            nb_tot_call = int(cum_sum[colname_to_index[call + "_Calls"]])

            if nb_tot_call > 0:
                table.append(
                    {
                        "call": call,
                        "call_mean": call_mean,
                        "call_min": call_min,
                        "call_max": call_max,
                        "call_tot_vol": call_tot_vol,
                        "nb_tot_call": nb_tot_call,
                        "comm": comm,
                    }
                )

    return table


def print_cct(table, comm_to_procs, comm_limit):
    comm_limit = (
        comm_limit
        if comm_limit != None and comm_limit < len(comm_to_procs)
        else len(comm_to_procs)
    )
    nb_comm_printed = 0
    for comm in comm_to_procs.keys():
        if nb_comm_printed >= comm_limit:
            break

        comm_calls = list(filter(lambda x: x["comm"] == comm, table))
        comm_calls.sort(key=lambda x: x["call_mean"], reverse=True)

        sys.stdout.write(BOLD)
        print("COMM".ljust(20) + "SIZE".ljust(10) + "PROCS")
        sys.stdout.write(RESET)
        print(
            f"{comm}".ljust(20)
            + f"{len(comm_to_procs[comm])}".ljust(10)
            + f"{compact_proc_list(comm_to_procs[comm])}",
            end="\n\n",
        )
        sys.stdout.write(RESET)
        print(
            "\t"
            + "Call".ljust(20)
            + "Mean[s]".rjust(10)
            + "Min[s]".rjust(10)
            + "Max[s]".rjust(10)
            + "Volume".rjust(15)
            + "#Calls".rjust(10)
        )

        print("\t" + "".join(["-"] * 75))

        for i in range(0, len(comm_calls)):
            call = comm_calls[i]["call"]
            call_mean = comm_calls[i]["call_mean"]
            call_min = comm_calls[i]["call_min"]
            call_max = comm_calls[i]["call_max"]
            call_tot_vol = comm_calls[i]["call_tot_vol"]
            nb_tot_call = comm_calls[i]["nb_tot_call"]
            print(
                "\t"
                + f"{call}".ljust(20)
                + f"{call_mean:.4f}".rjust(10)
                + f"{call_min:.4f}".rjust(10)
                + f"{call_max:.4f}".rjust(10)
                + f"{call_tot_vol}".rjust(15)
                + f"{nb_tot_call}".rjust(10)
            )

        print("\n")
        nb_comm_printed += 1

    return


def print_tc(table, comm_to_procs, limit):

    table.sort(key=lambda x: x["call_mean"], reverse=True)

    print("".join(["-"] * 100))
    print(
        "Call".ljust(20)
        + "Mean[s]".rjust(10)
        + "Min[s]".rjust(10)
        + "Max[s]".rjust(10)
        + "Volume".rjust(15)
        + "#Calls".rjust(10)
        + "Comm".rjust(15)
        + "Size".rjust(10)
    )

    print("".join(["-"] * 100))

    end_range = limit if limit != None and limit < len(table) else len(table)
    for i in range(0, end_range):
        call = table[i]["call"]
        call_mean = table[i]["call_mean"]
        call_min = table[i]["call_min"]
        call_max = table[i]["call_max"]
        call_tot_vol = table[i]["call_tot_vol"]
        nb_tot_call = table[i]["nb_tot_call"]
        comm = table[i]["comm"]
        comm_size = len(comm_to_procs[comm])
        print(
            f"{call}".ljust(20)
            + f"{call_mean:.4f}".rjust(10)
            + f"{call_min:.4f}".rjust(10)
            + f"{call_max:.4f}".rjust(10)
            + f"{call_tot_vol}".rjust(15)
            + f"{nb_tot_call}".rjust(10)
            + f"{comm}".rjust(15)
            + f"{comm_size}".rjust(10)
        )

    return


def print_header():
     print("\n" + "".join(["*"] * 70))
     for i in range(0, 5):
        if i == 2:
            print("*" + "mpisee".center(68) + "*")
        else:
            print("*" + f"*".rjust(69))

     print("".join(["*"] * 70), end="\n\n")
    #tprint("mpisee")


def main():
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument("-f", "--full", action="store_true", help="print full info")
    my_parser.add_argument("-cct", action="store_true", help="print cct")
    my_parser.add_argument("-tc", action="store_true", help="print tc")
    my_parser.add_argument("-tc_limit", action="store", type=int, help="print tc")
    my_parser.add_argument(
        "-cct_limit",
        action="store",
        type=int,
        help="Limit number of communicators to ptint in cct view",
    )
    my_parser.add_argument(
        "file_path", metavar="path", type=str, help="Path to the csv file"
    )
    args = my_parser.parse_args()

    table, mapping, comm_to_procs = prepare_data(args.file_path)

    print_header()
    print_mapping(mapping)

    if args.cct:
        print_cct(table, comm_to_procs, args.cct_limit)

    if args.tc:
        print_tc(table, comm_to_procs, args.tc_limit)


# Main body
if __name__ == "__main__":
    main()
