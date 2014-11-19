#CPPFLAGS = -g
CPPFLAGS = -O2
all: rgalgorithms.so

RSCRIPT= -I/Library/Frameworks/R.framework/Versions/3.0/Resources/library/Rcpp/include

clean:
	rm -f *.so *.o

rgalgorithms.so: rgalgorithms.cpp rgalgorithms.h rgalgorithms_scratch.cpp max-flow.cpp max-flow.hpp
	rm -f rgalgorithms.so # not sure why this is needed in windows
	R CMD SHLIB rgalgorithms.cpp rgalgorithms_scratch.cpp max-flow.cpp


debug:  debugrgalgorithms.so

rgalgorithms.o: rgalgorithms.cpp rgalgorithms.h
	clang++ -stdlib=libstdc++ -I/Library/Frameworks/R.framework/Resources/include -DNDEBUG  -I/usr/local/include   `Rscript -e 'Rcpp:::CxxFlags()'` -I/opt/local/include -fPIC  -mtune=core2 -g -O  -c rgalgorithms.cpp -o rgalgorithms.o

rgalgorithms_scratch.o: rgalgorithms_scratch.cpp rgalgorithms.h
	clang++ -stdlib=libstdc++ -I/Library/Frameworks/R.framework/Resources/include -DNDEBUG  -I/usr/local/include   `Rscript -e 'Rcpp:::CxxFlags()'` -I/opt/local/include -fPIC  -mtune=core2 -g -O  -c rgalgorithms_scratch.cpp -o rgalgorithms_scratch.o

max-flow.o: max-flow.cpp max-flow.hpp
	clang++ -stdlib=libstdc++ -I/Library/Frameworks/R.framework/Resources/include -DNDEBUG  -I/usr/local/include   `Rscript -e 'Rcpp:::CxxFlags()'` -I/opt/local/include -fPIC  -mtune=core2 -g -O  -c max-flow.cpp -o max-flow.o

debugrgalgorithms.so: rgalgorithms.o rgalgorithms_scratch.o max-flow.o
	clang++ -stdlib=libstdc++ -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/usr/local/lib -L/usr/local/lib -o rgalgorithms.so rgalgorithms.o rgalgorithms_scratch.o max-flow.o -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation
	touch debugrgalgorithms.so





