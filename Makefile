.c.o:
	$(CC) -c $(CFLAGS) $*.c

CC = gcc -std=gnu99
CFLAGS = -Wall -I/usr/local/include -g -O

LDFLAGS = 
LIBS = -lm -llapack -lblas -larpack -lgsl -lnetcdf -lhdf5_hl -lhdf5 -lz

PROGS = global

OBJS = main.o \
       getparam.o \
       probgen.o \
       gridgen.o \
       create_matrix.o \
       setup_arpack.o \
       find_sigma.o \
       arpack_driver.o \
       batch_driver.o \
       full_driver.o \
       shearlayerkcrit_driver.o \
       couette.o \
       shearlayer.o \
       writematrix.o \
       eigensolve.o \
       eigensolve_full.o \
       wnetcdf.o \
       dataprofile.o

global : $(OBJS)
	$(CC) -o global $(OBJS) $(LDFLAGS) $(LIBS)

clean:
	$(RM) $(PROGS) gmon.{out,sum} *.o *.s *~ core core.[0-9]*[0-9]