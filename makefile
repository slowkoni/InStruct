CC=gcc -O2 -Wall -ggdb -march=native

InStruct: InStruct.o data_interface.o mcmc.o initial.o check_converg.o result_analysis.o DPMM.o random.o quantile.o nrutil.o poly_geno.o
	$(CC) -Wall -o InStruct InStruct.o data_interface.o mcmc.o initial.o check_converg.o DPMM.o result_analysis.o quantile.o random.o nrutil.o poly_geno.o -lm

InStruct.o: InStruct.c data_interface.h mcmc.h random.h initial.h check_converg.h result_analysis.h quantile.h nrutil.h poly_geno.h
	$(CC) -Wall -c InStruct.c

data_interface.o: data_interface.c data_interface.h nrutil.h
	$(CC) -Wall -c data_interface.c

mcmc.o: mcmc.c mcmc.h data_interface.h random.h initial.h nrutil.h DPMM.h poly_geno.h
	$(CC) -Wall -c mcmc.c

poly_geno.o: poly_geno.c poly_geno.h mcmc.h data_interface.h random.h initial.h nrutil.h
	$(CC) -Wall -c poly_geno.c

DPMM.o: DPMM.c DPMM.h random.h nrutil.h
	$(CC) -Wall -c DPMM.c

initial.o: initial.c initial.h random.h nrutil.h
	$(CC) -Wall -c initial.c

result_analysis.o: result_analysis.c result_analysis.h data_interface.h mcmc.h initial.h quantile.h nrutil.h
	$(CC) -Wall -c result_analysis.c

check_converg.o: check_converg.c check_converg.h mcmc.h initial.h data_interface.h nrutil.h
	$(CC) -Wall -c check_converg.c

quantile.o: quantile.c quantile.h nrutil.h
	$(CC) -Wall -c quantile.c

random.o: random.c random.h nrutil.h
	$(CC) -Wall -c random.c

nrutil.o: nrutil.c nrutil.h
	$(CC) -Wall -c nrutil.c

