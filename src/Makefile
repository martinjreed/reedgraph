#CPPFLAGS = -g
#CPPFLAGS = -O2
all: rgalgorithms.so

clean:
	rm -f *.so *.o

rgalgorithms.so: rgalgorithms-igraph.cpp rgalgorithms.cpp rgalgorithms.h igraph_modified.cpp
	rm -f rgalgorithms.so # not sure why this is need in windows
	R CMD SHLIB rgalgorithms.cpp rgalgorithms_scratch.cpp rgalgorithms-igraph.cpp igraph_modified.cpp

test:
	echo $(CPPFLAGS)

