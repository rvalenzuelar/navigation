# Makefile for navigation corrections software

FC = gfortran
CC = gcc
FFLAGS = -g
CFLAGS = -g
LDFLAGS = 

# Libraries
NETCDF = /usr/local
LIBS = -L${NETCDF}/lib -lnetcdf -lnetcdff -lcurl -lhdf5 -lhdf5_hl
INCLUDES = -I${NETCDF}/include
CLIBS = -lgeotiff -ltiff

all: netcdf2text cns_noaa-p3_rv aster2txt

*.o:
	${FC} ${FFLAGS} ${LDFLAGS} -c $? ${INCLUDES}

cns_noaa-p3_rv: cns_noaa-p3_rv.f90 chol_inv.f
	${FC} ${FFLAGS} ${LDFLAGS} -o $@ $? ${LIBS} ${INCLUDES}

netcdf2text: netcdf2text.f90
	${FC} ${FFLAGS} ${LDFLAGS} -o $@ $? ${LIBS} ${INCLUDES}

aster2txt: aster2txt.o
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ $? ${CLIBS} ${CINCLUDES}

clean:
	rm -f core *.o *~ 

