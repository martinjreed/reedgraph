#CPPFLAGS = -g
#CPPFLAGS = -O2
all: rgalgorithms.so

clean:
	rm -f *.so *.o

rgalgorithms.so: rgalgorithms.cpp rgalgorithms.h rgalgorithms_scratch.cpp
	rm -f rgalgorithms.so # not sure why this is need in windows
	R CMD SHLIB rgalgorithms.cpp rgalgorithms_scratch.cpp 

test:
	echo $(CPPFLAGS)

