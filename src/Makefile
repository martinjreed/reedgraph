#CPPFLAGS = -g
#CPPFLAGS = -O2
all: rgalgorithms.so

clean:
	rm -f *.so *.o

rgalgorithms.so: rgalgorithms-igraph.cpp rgalgorithms.cpp rgalgorithms.h
	R CMD SHLIB -o rgalgorithms.so rgalgorithms-igraph.cpp rgalgorithms.cpp 

test:
	echo $(CPPFLAGS)

