CFLAGS= -fopenmp
INCS=
# use the following line to "Purify" the code
#CC=purify gcc
CC=gcc
SRCS=BaumWelch.c viterbi.c hmmutils.c \
  sequence.c genseq.c nrutil.c esthmm.c hmmrand.c \
  logmath.c vithmm.c fwd_bwd.c

all :	esthmm viterbi

esthmm: esthmm.o BaumWelch.o nrutil.o hmmutils.o \
    sequence.o logmath.o fwd_bwd.o viterbi.o
	 $(CC) -fopenmp -lgsl -lgslcblas -o esthmm esthmm.o nrutil.o \
    sequence.o hmmutils.o logmath.o fwd_bwd.o \
    BaumWelch.o viterbi.o -lm
viterbi: vithmm.o viterbi.o nrutil.o hmmutils.o sequence.o \
    BaumWelch.o logmath.o fwd_bwd.o viterbi.o
	 $(CC) -fopenmp -lgsl -lgslcblas -o viterbi vithmm.o viterbi.o nrutil.o \
    sequence.o logmath.o hmmutils.o  BaumWelch.o fwd_bwd.o \
    -lm 
clean:
	rm *.o 
# DO NOT DELETE THIS LINE -- make depend depends on it.

