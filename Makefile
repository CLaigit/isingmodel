
COMPILER =  gcc
OBJS1 = ising.c
LIB = -lm
ifdef gflag
 CFLAG += -g
endif

default: clean ising
all: clean ising

ising:	${OBJS}
	${COMPILER}	${CFLAG}	$(LIB)	${OBJS1}	-O3	-o	ising

clean:
	-rm	-f	*.o	ising
