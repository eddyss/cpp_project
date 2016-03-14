all: prog

prog: mt19937.o main.cpp
	g++ mt19937.o main.cpp -o prog

mt19937.o: mt19937.c
	g++ -c mt19937.c

clean:
	rm -rf *.o

mrproper: clean
	rm -rf prog
