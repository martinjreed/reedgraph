CPPFLAGS = -g
#CPPFLAGS = -O2
all: rgalgorithms.so rgalgorithms.o

clean:
	rm -f rgalgorithms.so rgalgorithms.o rgalgorithms_scratch.o *.gch

rgalgorithms.so: rgalgorithms.o rgalgorithms_scratch.o
	R CMD SHLIB rgalgorithms.o rgalgorithms_scratch.o

test:
	echo $(CPPFLAGS)

#	i686-pc-linux-gnu-gcc -std=gnu99 -shared -Wl,-O1 -o rgalgorithms.so rgalgorithms.o -L/home/mjreed/R/i686-pc-linux-gnu-library/2.10/Rcpp/lib -lRcpp -Wl,-rpath,/home/mjreed/R/i686-pc-linux-gnu-library/2.10/Rcpp/lib -L/usr/lib/R/lib -lR

rgalgorithms.h.gch: rgalgorithms.h
	ccache i686-pc-linux-gnu-g++ -I/usr/lib/R/include  -I/usr/local/include -I/home/mjreed/R/i686-pc-linux-gnu-library/current/Rcpp/include -fpic -Wno-deprecated $(CPPFLAGS) -march=i686 -pipe -c rgalgorithms.h


rgalgorithms.o: rgalgorithms.cpp rgalgorithms.h rgalgorithms.h.gch
	 ccache i686-pc-linux-gnu-g++ -I/usr/lib/R/include  -I/usr/local/include -I/usr/local/include/igraph -I/home/mjreed/R/i686-pc-linux-gnu-library/current/Rcpp/include -fpic -Wno-deprecated $(CPPFLAGS) -march=i686 -pipe -c rgalgorithms.cpp -o rgalgorithms.o


rgalgorithms_scratch.o: rgalgorithms_scratch.cpp rgalgorithms.h rgalgorithms.h.gch
	ccache i686-pc-linux-gnu-g++ -I/usr/lib/R/include  -I/usr/local/include -I/home/mjreed/R/i686-pc-linux-gnu-library/current/Rcpp/include -fpic -Wno-deprecated $(CPPFLAGS) -march=i686 -pipe -c rgalgorithms_scratch.cpp -o rgalgorithms_scratch.o

#R CMD COMPILE CCFLAGS="-g" rgalgorithms.cpp

#rgalgorithms.so: rgalgorithms.o
#	i686-pc-linux-gnu-g++ -g -shared -Wl,-O1 -o rgalgorithms.so rgalgorithms.o -L/home/mjreed/R/i686-pc-linux-gnu-library/2.10/Rcpp/lib -lRcpp -Wl,-rpath,/home/mjreed/R/i686-pc-linux-gnu-library/2.10/Rcpp/lib -L/usr/lib/R/lib -lR

#rgalgorithms.so: rgalgorithms.o
#	i686-pc-linux-gnu-g++ -shared -g -o rgalgorithms.so rgalgorithms.o -L/home/mjreed/R/i686-pc-linux-gnu-library/2.10/Rcpp/lib -lRcpp -Wl,-rpath,/home/mjreed/R/i686-pc-linux-gnu-library/2.10/Rcpp/lib -L/usr/lib/R/lib -lR

#rgalgorithms.o: rgalgorithms.cpp
#	i686-pc-linux-gnu-g++ -I/usr/lib/R/include  -I/usr/local/include   `Rscript -e 'Rcpp:::CxxFlags()'` -Wno-deprecated -fpic  -march=i686 -pipe -c rgalgorithms.cpp -o rgalgorithms.o

	#i686-pc-linux-gnu-g++ -I/usr/lib/R/include  -I/usr/local/include   `Rscript -e 'Rcpp:::CxxFlags()'` -Wno-deprecated -fpic  -march=i686 -pipe -c rgalgorithms.cpp -o rgalgorithms.o
#	i686-pc-linux-gnu-g++ -I/usr/lib/R/include  -I/usr/local/include   `Rscript -e 'Rcpp:::CxxFlags()'` -Wno-deprecated -fpic  -g -c rgalgorithms.cpp -o rgalgorithms.o

#rgalgorithms.so: rgalgorithms.o
#	i686-pc-linux-gnu-g++ -g -shared -Wl,-O1 -o rgalgorithms.so rgalgorithms.o -L/home/mjreed/R/i686-pc-linux-gnu-library/2.10/Rcpp/lib -lRcpp -Wl,-rpath,/home/mjreed/R/i686-pc-linux-gnu-library/2.10/Rcpp/lib -L/usr/lib/R/lib -lR

#rgalgorithms.o: rgalgorithms.cpp
#	i686-pc-linux-gnu-g++ -I/usr/lib/R/include  -I/usr/local/include   `Rscript -e 'Rcpp:::CxxFlags()'` -g -fpic -c rgalgorithms.cpp -o rgalgorithms.o

