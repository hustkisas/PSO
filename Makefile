SYSTEM     = x86-64_sles10_4.1
LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

CPLEXDIR      = /opt/ilog/cplex
CONCERTDIR    = /opt/ilog/concert
# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

CCOPT = -g   -fPIC -fexceptions  -DIL_STD
#CCOPT = -m32 -O -fPIC -fexceptions -DIL_STD

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXJARDIR   = $(CPLEXDIR)/lib/cplex.jar
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert  -lm -lpthread 
CLNFLAGS  = -L$(CPLEXLIBDIR) -lcplex  -lm -lpthread


CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

EXDIR         = $(CPLEXDIR)/examples
EXSRC         = $(EXDIR)/src
EXINC         = $(EXDIR)/include
EXDATA        = $(EXDIR)/data

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 


decomp: main.o mersenne.o DataGenerator.o DataModel.o OptModel.o CutPool.o WolseyStrengthening.o MixingSetDecomposition.o TwoConstrainSetCutGenerator.o IntersectionCutGenerator.o ResourceCutGenerator.o CutGenerator.o PortfolioInstanceGenerator.o CPXTimer.o PSCInstanceGenerator.o stoc1.o userintf.o
	$(CCC) $(CCFLAGS) main.o mersenne.o DataGenerator.o DataModel.o OptModel.o CutPool.o WolseyStrengthening.o MixingSetDecomposition.o TwoConstrainSetCutGenerator.o IntersectionCutGenerator.o ResourceCutGenerator.o CutGenerator.o PortfolioInstanceGenerator.o CPXTimer.o PSCInstanceGenerator.o stoc1.o userintf.o -o decomp $(CCLNFLAGS)
main.o: main.cpp
	$(CCC) -c $(CCFLAGS) main.cpp -o main.o
mersenne.o: mersenne.cpp
	$(CCC) -c $(CCFLAGS) mersenne.cpp -o mersenne.o
DataGenerator.o: DataGenerator.cpp
	$(CCC) -c $(CCFLAGS) DataGenerator.cpp -o DataGenerator.o
DataModel.o: DataModel.cpp
	$(CCC) -c $(CCFLAGS) DataModel.cpp -o DataModel.o
OptModel.o: OptModel.cpp
	$(CCC) -c $(CCFLAGS) OptModel.cpp -o OptModel.o
CutPool.o: CutPool.cpp
	$(CCC) -c $(CCFLAGS) CutPool.cpp -o CutPool.o
WolseyStrengthening.o: WolseyStrengthening.cpp 
	$(CCC) -c $(CCFLAGS) WolseyStrengthening.cpp -o WolseyStrengthening.o
MixingSetDecomposition.o: MixingSetDecomposition.cpp   
	$(CCC) -c $(CCFLAGS) MixingSetDecomposition.cpp -o MixingSetDecomposition.o
TwoConstrainSetCutGenerator.o: TwoConstrainSetCutGenerator.cpp   
	$(CCC) -c $(CCFLAGS) TwoConstrainSetCutGenerator.cpp -o TwoConstrainSetCutGenerator.o
IntersectionCutGenerator.o: IntersectionCutGenerator.cpp   
	$(CCC) -c $(CCFLAGS) IntersectionCutGenerator.cpp -o IntersectionCutGenerator.o  
ResourceCutGenerator.o: ResourceCutGenerator.cpp   
	$(CCC) -c $(CCFLAGS) ResourceCutGenerator.cpp -o ResourceCutGenerator.o  
CutGenerator.o: CutGenerator.cpp   
	$(CCC) -c $(CCFLAGS) CutGenerator.cpp -o CutGenerator.o  
PortfolioInstanceGenerator.o: PortfolioInstanceGenerator.cpp    
	$(CCC) -c $(CCFLAGS) PortfolioInstanceGenerator.cpp -o PortfolioInstanceGenerator.o  
CPXTimer.o: CPXTimer.cpp    
	$(CCC) -c $(CCFLAGS) CPXTimer.cpp -o CPXTimer.o  
PSCInstanceGenerator.o: PSCInstanceGenerator.cpp    
	$(CCC) -c $(CCFLAGS) PSCInstanceGenerator.cpp -o PSCInstanceGenerator.o  
stoc1.o: stoc1.cpp   
	$(CCC) -c $(CCFLAGS) stoc1.cpp -o stoc1.o  
userintf.o: userintf.cpp   
	$(CCC) -c $(CCFLAGS) userintf.cpp -o userintf.o  

