CFLAGS=-O3 -Wall -std=c++11
LFLAGS=-lNGS -lz
INCLUDE=-I/usr/include/NGS
GXX=g++-4.7

DensityStrand: DensityStrand.o
	$(GXX) $^ -o $@ $(LFLAGS)

%.o: %.cpp
	$(GXX) $(CFLAGS) $(INCLUDE) -c $^ -o $@

.PHONY: clean
clean: 
	rm -f *.o DensityStrand
