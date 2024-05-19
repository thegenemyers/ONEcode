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

### end of file
