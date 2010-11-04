all: lu

lu: main.o lu_decomposition.o
	g++ main.o lu_decomposition.o -o lu

main.o: main.cc
	g++ -c main.cc

lu_decomposition.o: lu_decomposition.cc
	g++ -c lu_decomposition.cc

clean: 
	rm *.o lu
