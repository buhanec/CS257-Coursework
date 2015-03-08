CC=gcc

# We intentionally turn off vectorization so we can compare
# no vectorization against intrinsics
OPT=-O2 -fno-tree-vectorize -ftree-vectorizer-verbose=3

BINARY = vectorsum
OBJ=  src/timer.o \
      src/vectorsumsse.o \
      src/initialise_arrays.o \
      src/nonsse.o \
      src/unroll.o \
      src/sse.o \

all: ${OBJ}
	${CC} ${OPT} ${OBJ} -o ${BINARY}

clean:
	rm -f ${OBJ}
	rm -f ${BINARY}

%.o: %.c
	${CC} ${OPT} -I./include -c -o $@ $<
