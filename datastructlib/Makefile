EXECS=exec
CFLAGS = -Wall -fno-inline -O3 -flto -march=native
SRCS = arith.c except.c assert.c memchk.c stack.c queue.c list.c atom.c table.c
OBJS = $(SRCS:.c=.o)
MAIN=exec
LIBFLAGS= -shared -ldl
AR=ar rcs

.PHONY: depend clean

all: LIB
	@echo  Data Structures Library has been compiled

shared: LIB
	@echo  Shared Library has been compiled

SHLIB: $(OBJS)
	$(CC) $(CFLAGS) $(LIBFLAGS) $(OBJS) -o libciface.so

$(MAIN): $(OBJS) main.c
	mpicc $(CFLAGS) main.c -o $(MAIN) $(OBJS)

LIB: $(OBJS)
	$(AR) libciface.a $(OBJS)

.c.o:
	$(CC) -fPIC $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	$(RM) *.o *~ $(MAIN) *.a *.so

depend: $(SRCS)
	makedepend $(INCLUDES) $^

#DO NOT DELETE
