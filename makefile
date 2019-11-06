# Find all source files, create a list of corresponding object files
SRCS = $(wildcard src/*.f90)
OBJS = $(patsubst %.f90,%.o,$(SRCS))

# Compiler/Linker settings
FC      = gfortran
FCFLAGS = -O3 -ffast-math -fexpensive-optimizations -flto -s
FLFLAGS = 

PROGRAM = backflow

# Compiler steps for all objects
$(OBJS) : %.o : %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

# Linker
$(PROGRAM) : $(OBJS)
	$(FC) $(FLFLAGS) -o $@ $^

#clean:
#	rm -f *.o *.mod
#	@echo "Cleanup complete!"

# Dependencies

jump_integration.o : quadpack_double.o
main.o: jump_integration.o eispack.o
