#
# Makefile for MDTCP
#
# Byoungseon Jeon
# Dept. Applied Science, UC.Davis
# July 21, 2009
#
# suffix definition for compiling
.SUFFIXES: .o .f90
#
# select compiler
F90 = /opt/openmpi.intel/bin/mpif90
#F90 = /Network/Servers/asgard.das.ucdavis.edu/Volumes/HD_USER/Home/bjeon/../anurag/sw-local/bin/mpif90
FLAGS = -openmp -O3 -fast -ipo -no-prec-div #-g
# Object files
OBJTS = data.o main.o init.o force.o vverlet.o result.o \
	fft235.o kernel.o zfft3d.o pztrans.o 
TARGET = tcp_pbc_mpi
#
# generation of executable
${TARGET}:data.o main.o init.o force.o vverlet.o result.o \
	fft235.o kernel.o zfft3d.o pztrans.o 
	${F90} ${FLAGS} -o ${TARGET} ${OBJTS} ${LIB} ${INCLUDE}
#
# generation of object files
.f90.o:
	${F90} ${FLAGS} -c $< ${LIB}
.f.o:
	${F90} ${FLAGS} -c $< ${LIB}
#${F90} -O3 -fomit-frame-pointer -I.. -c $< ${LIB}

#
# clean object files and executable
clean:
	rm -rf *.o *.f90~ core ${TARGET}
