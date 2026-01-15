# makefile for vgp-tools/src containing library and core utilities

DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -fPIC -Wextra -Wno-unused-result -fno-strict-aliasing -DNDEBUG # NDEBUG drops asserts
#CFLAGS = -g -Wall -Wextra -fno-strict-aliasing  # for debugging

CCPP = g++

LIB = libONE.a
PROGS = ONEstat ONEview

all: $(LIB) $(PROGS)

clean:
	$(RM) *.o ONEstat ONEview $(LIB) ZZ* TEST/ZZ* ONEcpptest.cpp ONEcpptest
	$(RM) -r *.dSYM
	$(RM) ONEcode*.so

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

ONEstat: ONEstat.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $^

ONEview: ONEview.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $^

### test

test: ONEview TEST
	./ONEview TEST/small.seq
	./ONEview -b -o TEST/ZZ-small.1seq TEST/small.seq
	bash -c "cd TEST ; source t1.sh ; source t2.sh ; cd .."

ONEcpptest.cpp: ONElib.hpp
	\ln -s ONElib.hpp $@

ONEcpptest: ONEcpptest.cpp ONElib.o
	$(CCPP) -D TEST_HEADER -o $@ $^

cpptest: ONEcpptest
	./ONEcpptest TEST/ZZ-small.1seq

### Python bindings

python: $(LIB)
	$(CCPP) -O3 -Wall -shared -std=c++11 -fPIC -undefined dynamic_lookup \
		$$(python3 -m pybind11 --includes) pyONElib.cpp \
		-o ONEcode$$(python3-config --extension-suffix) ONElib.o

python-test: python
	python3 -c "import ONEcode; print('ONEcode module loaded successfully')"

python-install: python
	pip3 install -e .

### end of file
