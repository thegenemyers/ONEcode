# makefile for Richard's ONEcode repository

DESTDIR = ~/bin

CFLAGS = -O3
#CFLAGS = -g				# for debugging

ALL = seqconvert seqextract seqstat

all: $(ALL)

install:
	cp $(ALL) $(DESTDIR)

clean:
	$(RM) *.o *~ $(ALL) *.gch
	$(RM) A.* B.* C.* small.1seq
	$(RM) -r *.dSYM

### object files

ONElib.o: ONElib.c ONElib.h 
	$(CC) $(CFLAGS) -c $^

UTILS_OBJS=hash.o dict.o array.o utils.o
UTILS_HEADERS=utils.h array.h dict.h hash.h
$(UTILS_OBJS): utils.h $(UTILS_HEADERS)

SEQIO_OPTS = -DONEIO
SEQIO_LIBS = -lm -lz 

seqio.o: seqio.c seqio.h 
	$(CC) $(CFLAGS) $(SEQIO_OPTS) -c $^

### programs

seqconvert: seqconvert.c seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

seqextract: seqextract.c seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

seqstat: seqstat.c seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

ONElogan: ONElogan.c ONElib.o
	$(CC) -fopenmp $(CFLAGS) $^ -o $@

#TEST

test: seqstat seqconvert seqextract
	./seqstat -b small.fa
	./seqconvert -1 small.fa > small.1seq
	./seqstat -l small.1seq
	./seqconvert -o A.fa.gz small.1seq
	gunzip A.fa.gz
	diff A.fa small.fa
	./seqconvert -t -1 -H -o B.1seq small.1seq
	./seqconvert -fa B.1seq
	./seqconvert -U -o C.fa B.1seq
	diff C.fa small.fa
	./seqextract -fa -f seq6 -c 3:10-20 small.1seq

### end of file
