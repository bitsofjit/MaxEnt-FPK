# makefile
# Purpose: Create the executable or the static library file for MAXENT
#
include makefile.inc

#LIB  = ../lib/libmaxent.a

MACPORTS = /opt/local

LIBPATH = -L$(MACPORTS)/lib/

EXT_LIBS  = $(LIBPATH) -llapack -lcblas -lf77blas -latlas

# On Linux can use: EXT_LIBS  = -lblas -llapack /usr/lib/libg2c.so.0
# (if DEFS is undefined in makefile.inc, then blas/lapack are not needed)

EXEC = ../bin/maxent
EXECP = ../bin/maxentfpk

SRC      = nrtype.f90 computeprior.f90 maxent.f90

SRCEXEC  = nrtype.f90 computeprior.f90 maxent.f90 main.f90

SRCEXECP = nrtype.f90 computeprior.f90 maxent.f90 maxent2dfpk.f90

OBJS = ${SRC:.f90=.o}

OBJSEXEC  = ${SRCEXEC:.f90=.o}
OBJSEXECP = ${SRCEXECP:.f90=.o}

MOD  = ${SRC:.f90=.mod}

maxent: 
	@echo "Compiling $(G95) $(FFLAGS) $(DEFS) $(SRCEXEC)"
	@$(G95) $(FFLAGS) $(DEFS) $(SRCEXEC)
	@echo "Linking..."
	@$(G95) -o $(EXEC) $(OBJSEXEC) $(EXT_LIBS)

maxentfpk: 
	@echo "Compiling $(G95) $(FFLAGS) $(DEFS) $(SRCEXECP)"
	@$(G95) $(FFLAGS) $(DEFS) $(SRCEXECP)
	@echo "Linking..."
	@$(G95) -o $(EXECP) $(OBJSEXECP) $(EXT_LIBS)

maxentlib: 
	@echo "Compiling $(G95) $(FFLAGS) $(DEFS) $(SRC)"
	@$(G95) $(FFLAGS) $(DEFS) $(SRC)
	@echo "Archiving ... $(AR) $(LIB) $(OBJS)"
	@$(AR) $(LIB) $(OBJS)

clean   :
	@echo "Removing *.o, *.mod, ../lib/libmaxent.a, and ../bin/maxent"
	@$(RM) $(OBJSEXEC)
	@$(RM) $(MOD)
	@$(RM) $(LIB)
	@$(RM) $(EXEC)
