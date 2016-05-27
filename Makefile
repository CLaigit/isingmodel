
COMPILER1 =  gcc
COMPILER2 = /usr/local/cuda-5.5/bin/nvcc
OBJS1 = src/ising.c
OBJS2 = src/ising2d.cu
LIB = -lm
# CFLAG = -arch=sm_20
CFLAG += -Wall

ifdef gflag
 CFLAG += -g
endif

ifdef picture
 CFLAG += -DPICORE=1
endif

default: clean ising
all: clean ising cuda

ising:	${OBJS}
	${COMPILER1}	${CFLAG}	$(LIB)	${OBJS1}	-O3	-o	src/ising.o
cuda:	${OBJS}
	${COMPILER2}	${CFLAG}	$(LIB)	${OBJS2}	-O3	-o	src/ising2d.o
clean:
	-rm	-f	src/*.o
cleandata:
	-rm	-f	data/*.txt
