# project name (generate executable with this name)
TARGET   = backflow

FC       = gfortran
# compiling flags here
FCFLAGS   = -Wall -O3 -ffast-math -fexpensive-optimizations -flto -s

LINKER   = gfortran
# linking flags here
LFLAGS   =

# change these to proper directories where each file should be
SRCDIR   = src

SOURCES  := $(wildcard $(SRCDIR)/*.f90)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.f90=%.o)
rm       = rm -f

$(TARGET): $(OBJECTS)
	@$(LINKER) $(OBJECTS) $(LFLAGS) -o $@
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.f90
	@$(FC) $(FCFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

clean:
	@$(rm) $(OBJECTS)
	@echo "Cleanup complete!"
