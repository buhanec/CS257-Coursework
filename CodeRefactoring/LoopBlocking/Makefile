CC=gcc
OPT=-O2

BINARY = loopblocking
OBJ=  src/timer.o \
      src/loopblocking.o \
      src/initialise_arrays.o \
      src/nonblocking.o \
      src/blocking.o

all: ${OBJ}
	${CC} ${OPT} ${OBJ} -o ${BINARY}

clean:
	rm -f ${OBJ}
	rm -f ${BINARY}

%.o: %.c
	${CC} ${OPT} -I./include -c -o $@ $<
