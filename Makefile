FC=gfortran
FFLAGS=-O3

PROGRAM=sub.f90 normal_modes.f90
OBJECT=$(PROGRAM:%.f90=%.o)

NM_calc: $(OBJECT)
	$(FC) $(FFLAGS) $^ -o NM_calc.x

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
