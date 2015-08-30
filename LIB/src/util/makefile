#  For IBM-RS6000
#CC     = xlC -+
#CFLAGS = -g 

#  For SUN - g++
CC     = g++
CFLAGS = -g

INCL   = -I../include -I../template
LIBDIR  = -L../

.c.o:
	$(CC) $(CCFLAGS) -c $<

all: totex makedep

totex: totex.o
	$(CC) $(CCFLAGS) $(INCL) -o totex   totex.o   $(LIBDIR) -ldo -lm
makedep: makedep.o 
	$(CC) $(CCFLAGS) $(INCL) -o makedep makedep.o $(LIBDIR) -ldo -lm

totex.o:        totex.c 
makedep.o: 	makedep.c
