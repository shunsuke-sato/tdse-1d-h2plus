FC = gfortran -O0 ## gfotran

LN = -llapack -lblas #other



VPATH = src:object
SRC = $(shell cd src ;ls *.f90 ;cd ..)
OBJ = $(SRC:.f90=.o)
OBJ_dir = $(addprefix object/,$(OBJ))

PROG = chain1d

$(PROG):global_variables.o $(OBJ)
	$(FC) -o $(PROG) $(OBJ_dir) $(LN)

main.o:main.f90
	$(FC) -c $< $(LN);mv $@  object 
%.o:%.f90
	$(FC) -c $< $(LN);mv $@  object 


clean:
	rm  -f  object/*.o  *.mod chain1d
clean_complete:
	rm  -f *~  */*~ */*/*~ object/*.o  */*.mod *.mod chain1d */#*
