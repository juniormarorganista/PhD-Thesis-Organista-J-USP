# makefile for gvts
OBJ1 = main.o condictions.o case_01.o case_02.o case_03.o case_04.o case_05.o case_06.o case_07.o case_08.o case_09.o case_10.o case_11.o case_12.o case_13.o case_14.o conform.o derivs.o drvs.o write_solutions.o euler_rk.o loop.o normerror.o poisson.o solvers.o boundarie_vorticity_cavity_flow.o 

main : $(OBJ1)
	gfortran -pg $(OBJ1) -o prg

clean: 
	rm *.o *.out

cleanall: 
	rm -rf *.o prg *.dat *.pdf *.png *.vtk *.txt fort.* *.out

run: 
	./prg 

.f.o:	
	gfortran -O0 -g -Wall -fbounds-check -c $<
#	gfortran -pg -c -O3 -funroll-all-loops $<
