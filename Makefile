CC = mpicc
FLAGS = -c -O2 -Wall -Wextra
LINK_FLAGS =
OBJS = main.o column.o gauss_jordan.o mpi_wrappers.o
EXECUTABLE_NAME = gj
TAR_COMPRESSED = gj.tar
GZ_TAR_COMPRESSED = gj.tar.gz
FILES=*.c *.h Makefile tests.sh

# Compile
all: $(OBJS)
	$(CC) $(OBJS) $(LINK_FLAGS) -o $(EXECUTABLE_NAME)
	rm -f $(OBJS)

main.o: main.c
	$(CC) $(FLAGS) main.c

column.o: column.c
	$(CC) $(FLAGS) column.c

gauss_jordan.o: gauss_jordan.c
	$(CC) $(FLAGS) gauss_jordan.c

mpi_wrappers.o: mpi_wrappers.c
	$(CC) $(FLAGS) mpi_wrappers.c

clean:
	rm -f $(EXECUTABLE_NAME)
	rm -f $(TAR_COMPRESSED)
	rm -f $(GZ_TAR_COMPRESSED)

run:
	mpiexec -n $(PN) ./$(EXECUTABLE_NAME) $(N) $(GROUP)

test:
	./tests.sh

zcompress:
	tar -zcvf $(GZ_TAR_COMPRESSED) $(FILES)

compress:
	tar -cvf $(TAR_COMPRESSED) $(FILES)

# Usage
help:
	@:
		$(info Run with: make run PN=NUMBER_OF_MPI_PROCESSES N=DIMENSION GROUP=FLAG(0 or 1))
