lipschitz: lipschitz.o
	g++ -o lipschitz -g lipschitz.o
lipschitz.o : lipschitz.cc
	g++ -c -g -Wall lipschitz.cc

.PHONY : clean
clean :
	-rm *.o
	-rm lipschitz
