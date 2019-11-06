# Find all source files, create a list of corresponding object files
SRCS = $(wildcard src/*.f90)
OBJS = $(patsubst %.f90,%.o,$(SRCS))

# Compiler/Linker settings
FC      = gfortran
FCFLAGS = -O3 -ffast-math -fexpensive-optimizations -flto -s
FLFLAGS = 

PROGRAM = backflow

all: $(PROGRAM)

# Compiler steps for all objects
$(OBJS) : %.o : %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

# Linker
$(PROGRAM) : $(OBJS)
	$(FC) $(FLFLAGS) -o $@ $^
	
clean:
	rm -rf *.mod $(OBJS)

# Dependencies

src/jump_integration.o : src/quadpack_double.o
src/main.o: src/jump_integration.o src/eispack.o
