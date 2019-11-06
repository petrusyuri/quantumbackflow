# Find all source files, create a list of corresponding object files
SRCS = src/quadpack_double.f90 src/eispack.f90 src/jump_integration.f90 src/main.f90
OBJS = $(patsubst %.f90,%.o,$(SRCS))

# Compiler/Linker settings
FC      = gfortran
FCFLAGS = -Wall -O3 -ffast-math -fexpensive-optimizations -flto -s
FLFLAGS = 

PROGRAM = quantumbackflow

# Compiler steps for all objects
$(OBJS) : %.o : %.f90
    $(FC) $(FCFLAGS) -o $@ $<

# Linker
$(PROGRAM) : $(OBJS)
    $(FC) $(FLFLAGS) -o $@ $^

clean:
	rm -f *.o *.mod
	@echo "Cleanup complete!"

# Dependencies

jump_integration.o : quadpack_double.o
main.o: jump_integration.o eispack.o
