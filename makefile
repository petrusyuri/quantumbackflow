# compiler
FC = gfortran

# compile flags
FCFLAGS = -Wall -O3 -ffast-math -fexpensive-optimizations -flto -s

# link flags
FLFLAGS =

# source files and objects
SRCS = $(patsubst %.F, %.o, $(wildcard *.F))

# program name
PROGRAM = backflow

all: $(PROGRAM)

$(PROGRAM): $(SRCS)
    $(FC) $(FLFLAGS) -o $@ $<

%.o: %.F
    $(FC) $(FCFLAGS) -o $@ $<

clean:
    rm -f *.o *.mod
