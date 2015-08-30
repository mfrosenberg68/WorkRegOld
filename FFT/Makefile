CFLAGS = -DGNU -w -O3 -Ddoubleonly
##CFLAGS = -g -DGNU -w -O3
INCL = -I$(HOME)/include -I$(HOME)/lib/include -I$(HOME)/lib/template -L$(HOME)/lib
LIBS = -ldo 
CC     = g++
#FILE   = *.C 

#OBJ = random.o qc_utilities.o
OBJ = random.o function.o datac.o qc_utilities.o
OBJ2 = random.o function.o qc_utilities.o
OBJCD = random.o function.o datacd.o qc_utilities.o
OBJCC = random.o function.o datacc.o qc_utilities.o
OBJFG = random.o function.o datafgop.o qc_utilities.o
PROGS = fft4gaussreg fft4tpgaus fft4tp2gaus fft4tpcd fft4gausscd fft4g2gaus fft4tpfgop fft4gfgop fft4tp scanfft4g2gaus scanfft4tpcd scanfft4gausscd scanfft4gfgop scanfft4tpfgop scanfft4tp2gaus scanfft4tpcc scanfft4gausscc scanfft4asvcc scanfft4g2asv scanfft4asvcd scanfft4glpcc scanfft4glpcd scanfft4asv2exp scanfft4gaus2exp scanfft4tpexp scanfft4tpdelta scanfft4gdelta scanfft4lpdelta

BACKUPDIR = $(HOME)/Shared/BACKUPS

.C.o:
	$(CC) $(CFLAGS) $(INCL) -c $< 

#all: Regfft_cos FFT4TP RegGaussFft FFT4TPGauss

all: singel scan

singel: fft4gaussreg fft4tpgaus fft4tp2gaus fft4tpcd fft4gausscd fft4g2gaus fft4tpfgop fft4gfgop fft4tp  

test: fft_test fft_cos

scan: scanfft4g2gaus scanfft4tpcd scanfft4gausscd scanfft4gfgop scanfft4tpfgop scanfft4tp2gaus scanfft4tpcc scanfft4gausscc scanfft4asvcc scanfft4g2asv scanfft4asvcd scanfft4glpcc scanfft4glpcd scanfft4asv2exp scanfft4gaus2exp scanfft4tpexp scanfft4tpdelta scanfft4gdelta scanfft4lpdelta

#OBj = gmp_first.o

fft_test:fft_test.o 
	$(CC) $(CFLAGS) -o fft_test fft_test.o 

fft_cos:fft_cos.o 
	$(CC) $(CFLAGS) -o fft_cos fft_cos.o 

fft4tp:fft4tp.o $(OBJ2)
	$(CC) $(CFLAGS) $(INCL) -o fft4tp fft4tp.o $(OBJ2) $(LIBS)

fft4gaussreg:fft4gaussreg.o $(OBJ)
	$(CC) $(CFLAGS) $(INCL) -o fft4gaussreg fft4gaussreg.o $(OBJ) $(LIBS)

fft4tpgaus:fft4tpgaus.o $(OBJ)
	$(CC) $(CFLAGS) $(INCL) -o fft4tpgaus fft4tpgaus.o $(OBJ) $(LIBS)

fft4tp2gaus:fft4tp2gaus.o $(OBJ)
	$(CC) $(CFLAGS) $(INCL) -o fft4tp2gaus fft4tp2gaus.o $(OBJ) $(LIBS)

fft4g2gaus:fft4g2gaus.o $(OBJ)
	$(CC) $(CFLAGS) $(INCL) -o fft4g2gaus fft4g2gaus.o $(OBJ) $(LIBS)

fft4tpcd:fft4tpcd.o $(OBJCD)
	$(CC) $(CFLAGS) $(INCL) -o fft4tpcd fft4tpcd.o $(OBJCD) $(LIBS)

fft4tpfgop:fft4tpfgop.o $(OBJFG)
	$(CC) $(CFLAGS) $(INCL) -o fft4tpfgop fft4tpfgop.o $(OBJFG) $(LIBS)

fft4gfgop:fft4gfgop.o $(OBJFG)
	$(CC) $(CFLAGS) $(INCL) -o fft4gfgop fft4gfgop.o $(OBJFG) $(LIBS)

fft4gausscd:fft4gausscd.o $(OBJCD)
	$(CC) $(CFLAGS) $(INCL) -o fft4gausscd fft4gausscd.o $(OBJCD) $(LIBS)

scanfft4tp2gaus:scanfft4tp2gaus.o $(OBJ)
	$(CC) $(CFLAGS) $(INCL) -o scanfft4tp2gaus scanfft4tp2gaus.o $(OBJ) $(LIBS)

scanfft4g2gaus:scanfft4g2gaus.o $(OBJ)
	$(CC) $(CFLAGS) $(INCL) -o scanfft4g2gaus scanfft4g2gaus.o $(OBJ) $(LIBS)

scanfft4tpcd:scanfft4tpcd.o $(OBJCD)
	$(CC) $(CFLAGS) $(INCL) -o scanfft4tpcd scanfft4tpcd.o $(OBJCD) $(LIBS)

scanfft4gausscd:scanfft4gausscd.o $(OBJCD)
	$(CC) $(CFLAGS) $(INCL) -o scanfft4gausscd scanfft4gausscd.o $(OBJCD) $(LIBS)

scanfft4gfgop:scanfft4gfgop.o $(OBJFG)
	$(CC) $(CFLAGS) $(INCL) -o scanfft4gfgop scanfft4gfgop.o $(OBJFG) $(LIBS)

scanfft4tpfgop:scanfft4tpfgop.o $(OBJFG)
	$(CC) $(CFLAGS) $(INCL) -o scanfft4tpfgop scanfft4tpfgop.o $(OBJFG) $(LIBS)

scanfft4tpcc:scanfft4tpcc.o $(OBJCC)
	$(CC) $(CFLAGS) $(INCL) -o scanfft4tpcc scanfft4tpcc.o $(OBJCC) $(LIBS)

scanfft4gausscc:scanfft4gausscc.o $(OBJCC)
	$(CC) $(CFLAGS) $(INCL) -o scanfft4gausscc scanfft4gausscc.o $(OBJCC) $(LIBS)

scanfft4asvcc:scanfft4asvcc.o $(OBJCC)
	$(CC) $(CFLAGS) $(INCL) -o scanfft4asvcc scanfft4asvcc.o $(OBJCC) $(LIBS)

scanfft4g2asv:scanfft4g2asv.o $(OBJ)
	$(CC) $(CFLAGS) $(INCL) -o scanfft4g2asv scanfft4g2asv.o $(OBJ) $(LIBS)

scanfft4asvcd:scanfft4asvcd.o $(OBJCD)
	$(CC) $(CFLAGS) $(INCL) -o scanfft4asvcd scanfft4asvcd.o $(OBJCD) $(LIBS)

scanfft4glpcd:scanfft4glpcd.o $(OBJCD)
	$(CC) $(CFLAGS) $(INCL) -o scanfft4glpcd scanfft4glpcd.o $(OBJCD) $(LIBS)

scanfft4glpcc:scanfft4glpcc.o $(OBJCC)
	$(CC) $(CFLAGS) $(INCL) -o scanfft4glpcc scanfft4glpcc.o $(OBJCC) $(LIBS)

scanfft4tpexp:scanfft4tpexp.o $(OBJ)
	$(CC) $(CFLAGS) $(INCL) -o scanfft4tpexp scanfft4tpexp.o $(OBJ) $(LIBS)

scanfft4asv2exp:scanfft4asv2exp.o $(OBJ)
	$(CC) $(CFLAGS) $(INCL) -o scanfft4asv2exp scanfft4asv2exp.o $(OBJ) $(LIBS)

scanfft4gaus2exp:scanfft4gaus2exp.o $(OBJ)
	$(CC) $(CFLAGS) $(INCL) -o scanfft4gaus2exp scanfft4gaus2exp.o $(OBJ) $(LIBS)

scanfft4tpdelta:scanfft4tpdelta.o $(OBJ)
	$(CC) $(CFLAGS) $(INCL) -o scanfft4tpdelta scanfft4tpdelta.o $(OBJ) $(LIBS)

scanfft4gdelta:scanfft4gdelta.o $(OBJ)
	$(CC) $(CFLAGS) $(INCL) -o scanfft4gdelta scanfft4gdelta.o $(OBJ) $(LIBS)

scanfft4lpdelta:scanfft4lpdelta.o $(OBJ)
	$(CC) $(CFLAGS) $(INCL) -o scanfft4lpdelta scanfft4lpdelta.o $(OBJ) $(LIBS)

clean:
	rm -fv *.o *~

clean_out:
	rm *.out

backup:
	tar cvf fftreg.tar *.[hc] *.[HC] Makefile
	gzip -9fv fftreg.tar
	cp -v fftreg.tar.gz $(BACKUPDIR)

rm_progs:
	rm -vf $(PROGS)
