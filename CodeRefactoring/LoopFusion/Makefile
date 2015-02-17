CC=gcc
OPT=-O2


BINARY = loopfusion
OBJ=  src/timer.o \
      src/gen_validation.o \
      src/initialise_arrays.o \
      src/fused.o \
      src/unfused.o \
      src/fused_peeled.o \
      src/fused_peeled_nointer.o \
      src/loopfusion.o

all: ${OBJ}
	${CC} ${OPT} ${OBJ} -o ${BINARY}

clean:
	rm -f ${OBJ}
	rm -f ${BINARY}

%.o: %.c
	${CC} ${OPT} -I./include -c -o $@ $<
