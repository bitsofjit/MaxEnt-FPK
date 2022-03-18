# makefile
# Purpose: Create the EXECUTABLE or the STATIC LIBRARY FILE for MAXENT
#

SRCDIR = src
TESTDIR = tests/sample.dat
EXEC = bin/maxent
MAKE = make
RM   = rm -rf

#
# To create the executable file in the bin directory, execute "make" or 
# "make maxent" on the command line
#
maxent:	
	    cd $(SRCDIR); $(MAKE) maxent
#
# To create the library file in lib directory, execute "make maxentlib" 
# on the command line
#
maxentlib :	
	    cd $(SRCDIR); $(MAKE) maxentlib

maxentfpk : 
	        cd $(SRCDIR); $(MAKE) maxentfpk

check     :     
	   $(MAKE) maxent; $(EXEC) < $(TESTDIR)

#	   for i in $(TESTDIR)/*.dat; do 
#	     echo " ";
#	     echo "******************************************************* " ; 
#	     echo "Using input test data file $$I " ; \
#	     echo "******************************************************* " ; 
#	     echo " "; 
#	   done ;
# $(RM) check.dat

clean :	
	    cd $(SRCDIR); make clean
