include ../macros.make

SHELL = /bin/sh
FFLAGS = $(PHYSFLAGS) $(R8FLAG) $(OMPFLAG) $(MICFLAGS)
LIBS = $(GPTLLIB)

BASENAME = wsm6_kernel wsm6_driver
SRCS = $(addsuffix .F90, $(BASENAME))
OBJS = $(addsuffix .o, $(BASENAME))

.SUFFIXES : .o .f90 .F90

$(WSM6EXE): $(OBJS)
	$(FCX) -o $@ $(ACCFLAGS) $(MICFLAGS) $(OBJS) $(GPTLOBJS) $(LIBS)
	rm *.o


%.o %.F90: 
	$(CPP) $(CPP_FLAGS) $*.F90 > $*.f90
	${FCX} -c ${FFLAGS} $*.f90

ifeq ($(HW) $(findstring intel,$(ARCH)),gpu intel)
  CPP_FLAGS += -DF2C -DINLINE
  LIBS += $(LIBCUDA)
#
wsm6_kernel.o: wsm6_kernel.F90
	$(CPP) $(CPP_FLAGS) $*.F90 > $*.f90
	$(F2C) --Generate=CUDA --Keep --Real8 --Constant KM=$(MAXNZ),KTS=1,KTE=41,ITS=1,ITE=10242 --Timer --ResidentOn=CPU $*.f90
	$(M4) -I $(M4_INCLUDE) $*.m4 > $*.cu
	$(NVCC) -c $(NVCC_FLAGS) $(GPU_FLAGS) -I$(F2C_HOME)/include $*.cu
endif
#
clean:
	$(RM) *.f90 *.o *.mod
	- find ./ -type l -exec rm  -f  {} \;



