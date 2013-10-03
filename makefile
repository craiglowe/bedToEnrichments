include /home/lowec/code/cTools/use_kent_new_gsl.mk

L += -lm -lz
MYLIBDIR = /home/lowec/kent/src/lib/${MACHTYPE}
MYLIBS =  ${MYLIBDIR}/jkweb.a

A = bedToEnrichments
H = bedLong.h
O = bedLong.o bedToEnrichments.o

bedToEnrichments: ${O} ${MYLIBS}
	${CC} ${COPT} -o ${A}${EXE} $O ${MYLIBS} $L

bedLong.o: bedLong.c bedLong.h
bedToEnrichments.o: bedToEnrichments.c bedLong.h

clean:
	rm -f ${A}${EXE} ${O}

