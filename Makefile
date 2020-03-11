CFLAGS= -fopenmp
INCS=
# use the following line to "Purify" the code
#CC=purify gcc
CC=gcc
SRCS=src/BaumWelch.c src/viterbi.c src/hmmutils.c \
  src/emutils.c src/nrutil.c src/esthmm.c src/hmmrand.c \
  src/logmath.c src/fwd_bwd.c src/fileutils.c

all :	TRACE

TRACE: src/esthmm.o src/BaumWelch.o src/nrutil.o src/hmmutils.o \
    src/emutils.o src/logmath.o src/fwd_bwd.o src/viterbi.o src/fileutils.o
	 $(CC) -fopenmp -lgsl -lgslcblas -o TRACE src/esthmm.o src/nrutil.o \
    src/emutils.o src/hmmutils.o src/logmath.o src/fwd_bwd.o \
    src/BaumWelch.o src/viterbi.o src/fileutils.o -lm

clean:
	rm src/*.o 
# DO NOT DELETE THIS LINE -- make depend depends on it.

