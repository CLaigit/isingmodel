
COMPILER1 =  gcc
COMPILER2 = /usr/local/cuda-5.5/bin/nvcc
OBJS1 = ising.c
LIB = -lm
ifdef gflag
 CFLAG += -g
endif

default: clean ising
all: clean ising

ising:	${OBJS}
	${COMPILER1}	${CFLAG}	$(LIB)	${OBJS1}	-O3	-o	ising
cuda:	${OBJS}
	${COMPILER2}	${CFLAG}	$(LIB)	${OBJS1}	-O3	-o	ising
clean:
	-rm	-f	*.o	ising
