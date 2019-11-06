# project name (generate executable with this name)
TARGET   = quantumbackflow

FC       = gfortran
# compiling flags here
FCFLAGS   = -Wall -O3 -ffast-math -fexpensive-optimizations -flto -s

LINKER   = gfortran
# linking flags here
LFLAGS   =

# change these to proper directories where each file should be
SRCDIR   = src
OBJDIR   = obj

SOURCES  := $(wildcard $(SRCDIR)/*.f90)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
rm       = rm -f


$(TARGET): $(OBJECTS)
	@$(LINKER) $(OBJECTS) $(LFLAGS) -o $@
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.f90
	@$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

clean:
	@$(rm) $(OBJECTS)
	@echo "Cleanup complete!"
