N ?= 8										# Make default variable
Np ?= 1

CXX=mpic++
INCDIR=.
CXXFLAGS=-O3 -I$(INCDIR)

SRC =$(INCDIR)/tma.h
SRC+=$(wildcard $(INCDIR)/tma/*.h)
SRC+=$(wildcard $(INCDIR)/tma/*/*.h)

main.cpp: $(SRC)

test: main.cpp $(SRC)
	@$(CXX) $(CXXFLAGS) -o $@ $<

run: test
	mpirun -np $(Np) ./test $(N)

clean:
	@rm -f test
