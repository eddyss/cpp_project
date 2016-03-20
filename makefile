all: prog

prog: mt19937.o main.cpp storage.o
	g++ mt19937.o storage.o main.cpp -o prog

storage.o: storage.cpp
	g++ -c storage.cpp

mt19937.o: mt19937.c
	g++ -c mt19937.c

clean:
	rm -rf *.o

mrproper: clean
	rm -rf prog
