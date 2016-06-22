INCLUDE  = -I/usr/include/eigen3/
LIBS     = -lgsl -lblas
CXXFLAGS = -g

all: prog clean

#################################################

prog: prog.o
	clang++ -g $(INCLUDE) prog.o $(LIBS) -o prog

prog.o: prog.cpp
	clang++ -c -g $(INCLUDE) prog.cpp -o prog.o

clean:
	rm  *.o

