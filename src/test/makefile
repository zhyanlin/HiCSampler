SYSTEM     = 
#LIBFORMAT  = static_pic
COMPILER = g++
#INCL_DIR =  -I../include
HEADERS = 
#LIBRARY =
SOURCES =

#CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio124/cplex
#CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio124/concert
#CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
OBJECTS = poissonRegression.o testPoissonRegresson.o  




# CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread
CCFLAGS = -O3 -static -lstdc++ -std=c++11  -DARMA_DONT_USE_WRAPPER  -pthread -lm -lz -lgsl -lgslcblas
#-lm -pthread -llapack  -static -static-libgfortran -lblas
#CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 

testPoissonRegresson: $(OBJECTS)
	$(COMPILER) -o $@ $(OBJECTS) $(CCFLAGS) 

all: 
	testPoissonRegresson;

clean:
	rm -f *.o







testPoissonRegresson.o:./testPoissonRegresson.cpp
	$(COMPILER) $(CCFLAGS) -c $<

%.o:../%.cpp
	$(COMPILER) $(CCFLAGS) -c $<
