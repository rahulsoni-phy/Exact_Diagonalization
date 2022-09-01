OBJS = functions.o main.o
CCC = g++
OPT = -O3
DEBUG = #-g3
MKL_LIB = -L /usr/lib/lapack/liblapack.so.3 -L /usr/lib/libblas/libblas.so.3 -llapacke 

all : $(OBJS)
	$(CCC) $(DEBUG) $(OBJS) -o main $(MKL_LIB)

functions.o : functions.cpp
	$(CCC) -c $(OPT) $(DEBUG) functions.cpp $(MKL_LIB)

main.o : main.cpp
	$(CCC) -c $(OPT) $(DEBUG) main.cpp $(MKL_LIB)


clean:
	rm -f *.o *~ main

