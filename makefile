
L=-pthread
HG_DEFS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DMACHTYPE_${MACHTYPE}

##########
#
# EDIT the lines below so that they point to your kent source
# and gsl header files and libraries 
#
HG_INC += -I/home/lowec/kent/src/hg/inc -I/home/lowec/kent/src/inc
L += /home/lowec/kent/src/lib/${MACHTYPE}/jkweb.a
HG_INC += -I/home/lowec/src/gsl/gsl-1.16_install/include
L += /home/lowec/src/gsl/gsl-1.16_install/lib/libgsl.a /home/lowec/src/gsl/gsl-1.16_install/lib/libgslcblas.a
#
# END of basic editing
#
##########

##########
#
# If you compiled your kent source with SSL support
# then you will need to edit the below lines as well
#
ifeq (${USE_SSL},1)
	L+=-lssl -lcrypto
	HG_DEFS+=-DUSE_SSL
endif
##########

##########
#
# If you compiled your kent source with sam/bam and tabix support
# then you will need to edit the below lines as well.
# If none of that sounded familiar to you, then you probably don't
# need to edit this.
#
ifeq (${USE_BAM},1)
	SAMINC = ${SAMDIR}
	SAMLIB = ${SAMDIR}/lib/libbam.a
	HG_INC += -I${SAMINC}
	L+=${SAMLIB}
	HG_DEFS+=-DUSE_BAM
	HG_DEFS+=-DKNETFILE_HOOKS
	TABIXINC = ${TABIXDIR}
	TABIXLIB = ${TABIXDIR}/libtabix.a
	HG_INC += -I${TABIXINC}
	L+=${TABIXLIB} -lz
	HG_DEFS+=-DUSE_TABIX
	HG_DEFS+=-DKNETFILE_HOOKS
endif
#
# End of sam/bam and tabix editing
#
##########

CC=gcc
ifeq (${COPT},)
    COPT=-O -g
endif
ifeq (${CFLAGS},)
    CFLAGS=
endif

SYS = $(shell uname -s)

ifeq (${HG_WARN},)
  ifeq (${SYS},Darwin)
      HG_WARN = -Wall -Wno-unused-variable
      HG_WARN_UNINIT=
  else
    ifeq (${SYS},SunOS)
      HG_WARN = -Wall -Wformat -Wimplicit -Wreturn-type
      HG_WARN_UNINIT=-Wuninitialized
    else
      HG_WARN = -Wall -Werror -Wformat -Wimplicit -Wreturn-type
      HG_WARN_UNINIT=-Wuninitialized
    endif
  endif
  # -Wuninitialized generates a warning without optimization
  ifeq ($(findstring -O,${COPT}),-O)
     HG_WARN += ${HG_WARN_UNINIT}
  endif
endif

CFLAGS += ${HG_WARN}

%.o: %.c
	${CC} ${COPT} ${CFLAGS} ${HG_DEFS} ${HG_WARN} ${HG_INC} ${XINC} -o $@ -c $<

L += -lm -lz

A = bedToEnrichments
H = bedLong.h
O = bedLong.o bedToEnrichments.o

bedToEnrichments: ${O} ${MYLIBS}
	${CC} ${COPT} -o ${A} $O ${MYLIBS} $L

bedLong.o: bedLong.c bedLong.h
bedToEnrichments.o: bedToEnrichments.c bedLong.h

clean:
	rm -f ${A} ${O}

