CC=gcc

# We intentionally turn off vectorization so we can compare
# no vectorization against intrinsics
OPT=-O2 -fno-tree-vectorize -ftree-vectorizer-verbose=3

BINARY = sse
OBJ=  src/timer.o \
      src/loopblockingsse.o \
      src/initialise_arrays.o \
      src/nonblocking.o \
      src/blocking.o \
      src/blocking_unroll.o \
      src/blocking_sse.o

all: ${OBJ}
	${CC} ${OPT} ${OBJ} -o ${BINARY}

clean:
	rm -f ${OBJ}
	rm -f ${BINARY}

%.o: %.c
	${CC} ${OPT} -I./include -c -o $@ $<
