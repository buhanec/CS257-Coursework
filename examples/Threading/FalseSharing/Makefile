CC=gcc
OPT=-O2
OPENMP=-fopenmp

BINARY = falsesharing
OBJ=  src/falsesharing.o \
      src/nofalseshare.o \
      src/initialise_arrays.o \
      src/timer.o 

all: ${OBJ}
	${CC} ${OPT} ${OBJ} ${OPENMP} -o ${BINARY}

clean:
	rm -f ${OBJ}
	rm -f ${BINARY}

%.o: %.c
	${CC} ${OPT} ${OPENMP} -I./include -c -o $@ $<
