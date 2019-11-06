# Directories
SRCDIR = ./src/
MODDIR = ./mod/
OBJDIR = ./obj/

# Find all source files, create a list of corresponding object files
SRCS = $(wildcard $(SRCDIR)*.f90)
OBJS = $(patsubst %.f90, %.o, $(patsubst $(SRCDIR), $(OBJDIR), $(SRCS)))

# Compiler/Linker settings
FC      = gfortran
FCFLAGS = -O3 -ffast-math -fexpensive-optimizations -flto -s -JDir $(MODDIR)
FLFLAGS = 

PROGRAM = backflow

all: $(PROGRAM)

# Compiler steps for all objects
$(OBJS) : %.o : %.f90
	mkdir -p $(MODDIR)
	mkdir -p $(OBJDIR)
	$(FC) $(FCFLAGS) -c -o $@ $<

# Linker
$(PROGRAM) : $(OBJS)
	mkdir -p $(MODDIR)
	mkdir -p $(OBJDIR)
	$(FC) $(FLFLAGS) -o $@ $^
	
clean:
	rm -rf *.mod $(OBJS)

# Dependencies

src/jump_integration.o : src/quadpack_double.o
src/main.o: src/jump_integration.o src/eispack.o
