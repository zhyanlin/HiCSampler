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
OBJECTS = poissonRegression.o bspline.o curvefit.o contactMap.o readOptions.o utils.o optimize.o KDE.o MRF.o adaptMCMC.o  HiCSampler.o

OBJECTS2=contactMap.o readOptions.o utils.o  IFevaluation.o

OBJECTS_singleEntry=singleEntry.o utils.o contactMap.o readOptions.o
OBJECTS_decomp_ll=decomp_ll.o utils.o contactMap.o readOptions.o
#CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
# CFLAGS = -g 


#CONCERTINCDIR = $(CONCERTDIR)/include
#CPLEXINCDIR   = $(CPLEXDIR)/include



#

# CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread
CCFLAGS = -O3 -fpermissive -static -lstdc++ -std=c++11  -DARMA_DONT_USE_WRAPPER  -pthread -lm -lz -lgsl -lgslcblas
#-lm -pthread -llapack  -static -static-libgfortran -lblas
#CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 

all: 
	make HiCSampler

clean:
	rm -f *.o HiCSampler

HiCSampler: $(OBJECTS)
	$(COMPILER) -o $@ $(OBJECTS) $(CCFLAGS) 
IFevaluation: $(OBJECTS2)
	$(COMPILER) -o $@ $(OBJECTS2) $(CCFLAGS) 
singleEntry: $(OBJECTS_singleEntry)
	$(COMPILER) -o $@ $(OBJECTS_singleEntry) $(CCFLAGS)
decomp_ll: $(OBJECTS_decomp_ll)
	$(COMPILER) -o $@ $(OBJECTS_decomp_ll) $(CCFLAGS)

IFevaluation.o:../src/misc/IFevaluation.cpp
	$(COMPILER) $(CCFLAGS) -c $<

decomp_ll.o:../src/misc/decomp_ll.cpp
	$(COMPILER) $(CCFLAGS) -c $<
singleEntry.o:../src/misc/singleEntry.cpp
	$(COMPILER) $(CCFLAGS) -c $<
%.o:./src/%.cpp
	$(COMPILER) $(CCFLAGS) -c $<

