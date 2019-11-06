CC = gfortran
CFLAGS = -Wall -O3 -ffast-math -fexpensive-optimizations -flto -s

backflow: quadpack_double.o eispack.o jump_integration.o main.o
	$(CC) -o backflow quadpack_double.o eispack.o jump_integration.o main.o
