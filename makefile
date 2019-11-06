# Directories
SRCDIR = ./src/
MODDIR = ./mod/
OBJDIR = ./obj/

# Find all source files, create a list of corresponding object files
SRCS = $(wildcard $(SRCDIR)*.f90)
OBJS = $(patsubst %.f90, %.o, $(patsubst $(SRCDIR), $(OBJDIR), $(SRCS)))

# Compiler/Linker settings
FC      = gfortran
FCFLAGS = -O3 -ffast-math -fexpensive-optimizations -flto -s -J $(MODDIR)
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

$(SRCDIR)jump_integration.o : $(SRCDIR)quadpack_double.o
$(SRCDIR)main.o: $(SRCDIR)jump_integration.o $(SRCDIR)eispack.o
