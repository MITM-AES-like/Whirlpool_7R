CC= g++
CFLAGS= -fopenmp -maes -mavx2 -O3 -std=c++11
DFLAGS= -DDoFF=32

all: F36 F40 F44 F48 F56 F64

F36:
	$(CC) -o $@ main.cpp $(CFLAGS) $(DFLAGS) -Domp_nb_threads=1 -DPARTIAL_InBits=36

F40:
	$(CC) -o $@ main.cpp $(CFLAGS) $(DFLAGS) -Domp_nb_threads=1 -DPARTIAL_InBits=40

F44:
	$(CC) -o $@ main.cpp $(CFLAGS) $(DFLAGS) -Domp_nb_threads=1 -DPARTIAL_InBits=44

F48:
	$(CC) -o $@ main.cpp $(CFLAGS) $(DFLAGS) -Domp_nb_threads=4 -DPARTIAL_InBits=48

F56:
	$(CC) -o $@ main.cpp $(CFLAGS) $(DFLAGS) -Domp_nb_threads=4 -DPARTIAL_InBits=56

F64:
	$(CC) -o $@ main.cpp $(CFLAGS) $(DFLAGS) -Domp_nb_threads=4 -DPARTIAL_InBits=64

.PHONY:clean
clean:
	rm F36 F40 F44 F48 F56 F64