##
# Testfile
#
# @file
# @version 0.1

LIBF=../lib
LIBPATH=$(abspath $(LIBF))

all:
	$(CC) -g -L$(LIBPATH) comm_group.c -o comm_group.out -lmpisee
	$(CC) -g -L$(LIBPATH) comm_topo.c -o comm_topo.out -lmpisee
	$(CC) -g -L$(LIBPATH) comm_split.c -o comm_split.out -lmpisee
	$(CC) -g -L$(LIBPATH) comm_split_v2.c -o comm_split_v2.out -lmpisee
	$(CC) -g -L$(LIBPATH) alltoallv.c -o alltoallv.out  -lmpisee
	$(CC) -g -L$(LIBPATH) cart_sub.c -o cart_sub.out  -lmpisee
	$(CC) -g -L$(LIBPATH) tuwqsort.c -o tuwqsort.out  -lmpisee
	$(CC) -g -L$(LIBPATH) wait.c -o wait.out  -lmpisee
	$(CC) -g -L$(LIBPATH) waitall.c -o waitall.out  -lmpisee
	$(CC) -g -L$(LIBPATH) sendrecv.c -o sendrecv.out  -lmpisee
	$(CC) -g -L$(LIBPATH) barrier.c -o barrier.out  -lmpisee
	$(CC) -g -L$(LIBPATH) stress_split.c -o stress_split.out -lmpisee
	$(CC) -g -L$(LIBPATH) stress_waitall.c -o stress_waitall.out -lmpisee

# end
