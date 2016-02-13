
Objects=lensing.o entree_lensing.o rng.o

LIBS = -openmp -L/home/nis/gavazzi/lib -lsoftlens 
OPT = -I/home/nis/gavazzi/include

#FC=gfortran $(OPT)
FC=ifort $(OPT)

#======Cible du Makefile======#
all : lensing 

#=========Executable=========#
lensing : $(Objects)
	$(FC)  -o lensing $(Objects) $(LIBS)


#======Fichier intermedaire======#
entree_lensing.mod: entree_lensing.o entree_lensing.f90
	$(FC)  -c entree_lensing.f90

lensing.o: lensing.f90 entree_lensing.o   rng.o
	$(FC)  -c lensing.f90

entree_lensing.o: entree_lensing.f90
	$(FC)  -c entree_lensing.f90
#==========Autre cible==========#
clean:
	rm -f *.o *.mod

rng.o: rng.f90
	$(FC) -fpp -c rng.f90



#=========Executable=========#
test_rng : rng.o
	$(FC)  -o test_rng rng.o  $(LIBS)

