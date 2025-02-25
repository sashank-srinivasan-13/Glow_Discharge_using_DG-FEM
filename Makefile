PROG = LDG
CC = mpiCC -O2 -flto
SRCS = backup.h dgsolve.h fluxes.h fluxpen.h global_var.h initialize.h limiter.h matrix.h num_int.h out.h poscheck.h poisson.h qreshape.h runmain.h updatevals.h timesteps.h backup.cpp dgsolve.cpp  fluxes.cpp fluxpen.cpp global_var.cpp initialize.cpp limiter.cpp matrix.cpp num_int.cpp out.cpp poscheck.cpp  poisson.cpp qreshape.cpp runmain.cpp updatevals.cpp timesteps.cpp

#LIBS = -I/home/sriniv35/lapack/LAPACKE/include -I/home/sriniv35/Boost -L/home/sriniv35/lapack -llapacke -llapack -lblas
LIBS = -I/home/roger/a/sriniv35/C++_LIbraries/lapack/LAPACKE/include -I/home/roger/a/sriniv35/C++_LIbraries/Boost -L/home/roger/a/sriniv35/C++_LIbraries/lapack -llapacke -llapack -lblas


$(PROG): $(SRCS)
	 $(CC) -o $(PROG) $(SRCS) $(INCLUDES) $(LIBS) -o $(PROG)


.PHONY: clean cleanresults

clean:
	rm -f $(PROG)
	rm *.dat *.plt

cleanresults:
	rm *.dat *.plt

