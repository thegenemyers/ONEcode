# makefile for vgp-tools/src containing library and core utilities

DEST_DIR = ~/bin

#CFLAGS= -O3 -Wall -fPIC -Wextra -Wno-unused-result -fno-strict-aliasing -DNDEBUG # NDEBUG drops asserts
CFLAGS= -g -Wall -Wextra -fno-strict-aliasing  # for debugging

CCPP=g++

LIB = libONE.a
PROGS = ONEstat ONEview

all: $(LIB) $(PROGS)

clean:
	$(RM) *.o ONEstat ONEview $(LIB) ZZ* TEST/ZZ* ONEcpptest.cpp ONEcpptest
	$(RM) -r *.dSYM

install:
	cp $(PROGS) $(DEST_DIR)

package:
	make clean
	tar -zcf ONE-core.tar.gz *.c *.h Makefile

### library

LIB_OBJS = ONElib.o

ONElib.o: ONElib.h

$(LIB): $(LIB_OBJS)
	ar -cr $@ $^
	ranlib $@

### programs

ONEstat: ONEstat.c utils.o $(LIB)
	$(CC) $(CFLAGS) -o $@ $^

ONEview: ONEview.c utils.o $(LIB)
	$(CC) $(CFLAGS) -o $@ $^

### test

ONEcpptest.cpp: ONElib.hpp
	\ln -s ONElib.hpp $@

ONEcpptest: ONEcpptest.cpp ONElib.o
	$(CCPP) -D TEST_HEADER -o $@ $^

test: ONEview ONEcpptest
	./ONEview TEST/small.seq
	cd test ; ../ONEview -b -o ZZ-small.1seq small.seq
	./ONEcpptest TEST/ZZ-small.1seq
	cd TEST ; source t1.sh ; cd ..

### end of file
