# Directories
SRCDIR = ./src/
MODDIR = ./mod/
OBJDIR = ./obj/
OUTDIR = ./out/

# Find all source files, create a list of corresponding object files
SRCS = $(wildcard $(SRCDIR)*.f90)
OBJS = $(patsubst %.f90, %.o, $(patsubst $(SRCDIR)%, $(OBJDIR)%, $(SRCS)))
MODS = $(patsubst %.f90, %.mod, $(patsubst $(SRCDIR)%, $(MODDIR)%, $(SRCS)))

# Compiler/Linker settings
FC      = gfortran
FCFLAGS = -O3 -ffast-math -fexpensive-optimizations -flto -s -J $(MODDIR)
FLFLAGS = 

PROGRAM = backflow

all: prep $(PROGRAM)

# Compiler steps for all objects
$(OBJS) : $(OBJDIR)%.o : $(SRCDIR)%.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

# Linker
$(PROGRAM) : $(OBJS)
	$(FC) $(FLFLAGS) -o $@ $^

prep:
	mkdir -p $(MODDIR)
	mkdir -p $(OBJDIR)
	mkdir -p $(OUTDIR)

clean:
	rm -rf $(MODS) $(OBJS)

# Dependencies

$(OBJDIR)jump_integration.o : $(OBJDIR)quadpack_double.o
$(OBJDIR)main.o: $(OBJDIR)jump_integration.o $(OBJDIR)eispack.o
