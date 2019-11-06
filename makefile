# project name (generate executable with this name)
TARGET   = backflow

#compiler
FC       = gfortran
FCFLAGS   = -Wall -O3 -ffast-math -fexpensive-optimizations -flto -s

#dependencies
DEPS = src/quadpack_double.f90 src/eispack.f90 src/jump_integration.f90

#objects
OBJ = quadpack_double.o eispack.o jump_integration.o main.o

#compiling project
$(TARGET): $(OBJ)
	$(FC) -o $@ $^ $(FCFLAGS)
	@echo "Compiled "$<" successfully!"

#compiling dependencies
%.o: %.f90 $(DEPS)
	$(FC) -c -o $@ $< $(FCFLAGS)

clean:
	rm -f *.o *.mod
	@echo "Cleanup complete!"
