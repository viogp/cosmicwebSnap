####################################################
###          User definable stuff
#Monopole
NB_R = 84
R_MAX = 65.
#Logarithmic binning
N_LOGINT = 12
#Behavior options
DEFINEOPTIONS = -D_VERBOSE -D_LONGIDS
DEFINEOPTIONS += -D_LOGBIN
#DEFINEOPTIONS = -D_HAVE_OMP
$(info ${DEFINEOPTIONS})
### End of user-definable stuff
####################################################

#Transform to inverse scales
I_R_MAX = $(shell echo "scale=5;1./$(R_MAX)" | bc)
LOG_R_MAX = $(shell echo "scale=9;l($(R_MAX))/l(10)" | bc -l)

# DEFINES for the OpenMP version
DEFINEFLAGSCPU = $(DEFINEOPTIONS)
DEFINEFLAGSCPU += -DNB_R=$(NB_R) -DI_R_MAX=$(I_R_MAX) -DLOG_R_MAX=$(LOG_R_MAX)
DEFINEFLAGSCPU += -DN_LOGINT=$(N_LOGINT)

# COMPILER AND OPTIONS
COMPCPU = gcc
OPTCPU = -Wall -O3 -fopenmp $(DEFINEFLAGSCPU)

#INCLUDES AND LIBRARIES
INCLUDECOM = -I./src
LIBCPU = -lm -lhdf5

#.o FILES
#CUTE
DEF = src/define.o
COM = src/common.o
PM = src/pm.o
TREE = src/tree.o
IO = src/io.o
NEIGH = src/neighbors.o
CORR = src/correlator.o
MAIN = src/main.c
OFILES = $(DEF) $(COM) $(PM) $(TREE) $(IO) $(NEIGH) $(CORR) $(MAIN)

#FINAL GOAL
EXE = CUTE_box

#RULES
default : $(EXE)
#RULE TO MAKE .o's FROM .c's
$(DEF) : src/define.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(COM) : src/common.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(CORR) : src/correlator.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(TREE) : src/tree.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(IO) : src/io.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(NEIGH) : src/neighbors.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(PM) : src/pm.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)

#RULES TO MAKE THE FINAL EXECUTABLES
$(EXE) : $(OFILES)
	$(COMPCPU) $(OPTCPU) $(OFILES) -o $(EXE) $(INCLUDECOM) $(LIBCPU)

#CLEANING RULES
clean :
	rm -f ./src/*.o

cleaner :
	rm -f ./src/*.o ./src/*~ *~ CUTE_box
