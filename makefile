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
OBJECTS = contactMap.o readOptions.o utils.o optimize.o KDE.o MRF.o adaptMCMC.o diffusionMCMC.o adaptiveHiC.o

OBJECTS2=contactMap.o readOptions.o utils.o optimize.o IFevaluation.o


#CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
# CFLAGS = -g 


#CONCERTINCDIR = $(CONCERTDIR)/include
#CPLEXINCDIR   = $(CPLEXDIR)/include



#

# CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread
CCFLAGS = -O3  -lstdc++ -std=c++11  -DARMA_DONT_USE_WRAPPER  -pthread -lm -lz
#-lm -pthread -llapack  -static -static-libgfortran -lblas
#CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 

all: 
	adaptiveHiC; IFevaluation

clean:
	rm -f *.o adaptiveHiC

adaptiveHiC: $(OBJECTS)
	$(COMPILER) -o $@ $(OBJECTS) $(CCFLAGS) 
IFevaluation: $(OBJECTS2)
	$(COMPILER) -o $@ $(OBJECTS2) $(CCFLAGS) 

IFevaluation.o:../src/misc/IFevaluation.cpp
	$(COMPILER) $(CCFLAGS) -c $<


%.o:../src/%.cpp
	$(COMPILER) $(CCFLAGS) -c $<

