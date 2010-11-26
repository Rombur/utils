all : reordering gcd

reordering : reordering.o  
	g++ -o reordering -g reordering.o 
reordering.o : reordering.cc
	g++ -c -g -Wall reordering.cc

gcd : gcd.o
	g++ -o  gcd gcd.o
gcd.o : gcd.cc

.PHONY : clean
clean :
	-rm *.o
	-rm gcd
	-rm reordering
