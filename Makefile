# Makefile for C version of CUSS
# Andrew Turpin
# Wed  2 Mar 2011 13:29:52 EST

CC=g++
#CFLAGS= -Wall -O9 -funroll-loops -L c:/Program\ Files\ \(x86\)/GnuWin32/lib -I c:/Program\ Files\ \(x86\)/GnuWin32/include -static -Wl,--stack,500000000
CFLAGS= -Wall -O9 -funroll-loops
#CFLAGS= -Wall -O3 -funroll-loops

OBJ= main.o smalloc.o tree.o chunk.o

all: buss printTree

buss: $(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o buss $(OBJ) -lgsl -lgslcblas

printTree: tree.o chunk.o printTree.o smalloc.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o printTree tree.o chunk.o printTree.o smalloc.o -lgsl -lgslcblas

main.o: Makefile ent.h tree.h chunk.h
log.o: Makefile log.h
smalloc.o: Makefile smalloc.h
hashtable.o: Makefile hashtable.h
tree.o: Makefile tree.h
chunk.o: Makefile chunk.h

clean:
	\rm -f *.o buss

nosource:
	\rm -f $(OBJ)

ent.h:
	R --slave < makeLogLut.r > ent.h
