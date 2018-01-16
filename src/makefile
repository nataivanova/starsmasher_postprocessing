FC = $(shell which gfortran)
CPP = $(shell which g++)

FF = -O4 -ffixed-line-length-132
CF = -O4

LIB = -L/usr/lib -lstdc++ #-lc

#       macro definition
#.SUFFIXES: .c .o .C .f .for.cpp



#       implicit rules
.for.o: $*.for
	$(FC) $(LIBS) $(FF) -c  ${*}.for -o  $*.o -static

.f.o:   $*.f
	$(FC) $(LIBS) $(FF) -c  ${*}.f -o  $*.o -static

%.o:%.cpp
	$(CPP) -c $(CF) $<


process_fobj = extract_from_sph.o read_file.o ejecta.o orbital_parameter.o classification.o omega.o useeostable.o getTemperature.o opacfile.o  opac_rho_t_temp.o
process_cobj = roche_potential.o deg_of_ion.o
process_obj = $(process_fobj) $(process_cobj)
${process_fobj}: common_sph_var.h

all:	sortit enc extract

extract: $(process_obj)
	 $(FC) -o extract $(process_obj) $(LIB)

sortit: sortit.o
	$(FC) $(LIBS) $(FF) -o sortit sortit.f

enc: en_comparison.o
	$(FC) $(LIBS) $(FF) -o enc en_comparison.f



clean:
	$(shell which rm) *.o enc sortit extract
