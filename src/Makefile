#CPPFLAGS = -g
#CPPFLAGS = -O2
all: rgalgorithms.so

clean:
	rm -f *.so *.o

rgalgorithms.so: rgalgorithms.cpp rgalgorithms.h rgalgorithms_scratch.cpp max-flow.cpp max-flow.hpp
	rm -f rgalgorithms.so # not sure why this is needed in windows
	R CMD SHLIB rgalgorithms.cpp rgalgorithms_scratch.cpp max-flow.cpp


debug:  rgalgorithms.so

rgalgorithms.o: rgalgorithms.cpp rgalgorithms.h
	g++ -arch i386 -I/Library/Frameworks/R.framework/Resources/include -I/Library/Frameworks/R.framework/Resources/include/i386 -DNDEBUG  -I/usr/local/include `Rscript -e 'Rcpp:::CxxFlags()'` -fPIC  -g   -c rgalgorithms.cpp -o rgalgorithms.o

rgalgorithms_scratch.o: rgalgorithms_scratch.cpp rgalgorithms.h
	g++ -arch i386 -I/Library/Frameworks/R.framework/Resources/include -I/Library/Frameworks/R.framework/Resources/include/i386 -DNDEBUG  -I/usr/local/include `Rscript -e 'Rcpp:::CxxFlags()'` -fPIC  -g   -c rgalgorithms_scratch.cpp -o rgalgorithms_scratch.o

max-flow.o: max-flow.cpp max-flow.hpp
	g++ -arch i386 -I/Library/Frameworks/R.framework/Resources/include -I/Library/Frameworks/R.framework/Resources/include/i386 -DNDEBUG  -I/usr/local/include `Rscript -e 'Rcpp:::CxxFlags()'` -fPIC  -g   -c max-flow.cpp -o max-flow.o

debugrgalgorithms.so: rgalgorithms.o rgalgorithms_scratch.o max-flow.o
	g++ -arch i386 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/usr/local/lib -o rgalgorithms.so rgalgorithms.o rgalgorithms_scratch.o max-flow.o /Library/Frameworks/R.framework/Versions/2.15/Resources/library/Rcpp/lib/i386/libRcpp.a -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation




