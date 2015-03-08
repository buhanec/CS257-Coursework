CC=gcc
OPT=-O2


BINARY = loopinterchange
OBJ=  src/loopinterchange.o \
      src/initialise_arrays.o \
      src/rowfirst.o \
      src/columnfirst.o \
      src/timer.o 

all: ${OBJ}
	${CC} ${OPT} ${OBJ} -o ${BINARY}

clean:
	rm -f ${OBJ}
	rm -f ${BINARY}

%.o: %.c
	${CC} ${OPT} -I./include -c -o $@ $<
