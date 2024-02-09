
OFILES = base/Matrix.o \
	 base/Volume.o \
	 base/Fields.o \
	 base/LinearAlgebra.o \
	 base/Camera.o \
	 base/Color.o \
	 base/ImplicitVolumeShapes.o \
	 base/ImplicitVectorShapes.o \
	 base/ImplicitColors.o \
	 base/ImplicitMatrixShapes.o \
	 base/Models.o \
	 base/ImgProc.o \
	 base/VolumeRenderer.o


ROOTDIR = .
LIB = $(ROOTDIR)/lib/libstarter.a

CXX = g++ -Wall -g -O2 -fPIC -fopenmp -std=c++14


PYTHONINCLUDE = -I/usr/include/python3.8


INCLUDES = -I ./include/ $(PYTHONINCLUDE)  -I ../build/include -I /usr/local/include/ -I/usr/include/



.C.o: $<
	$(CXX) -c $(INCLUDES) $< -o $@

base: $(OFILES)
	ar rv $(LIB) $?
	$(CXX) base/Main.C $(INCLUDES) -L./lib -lstarter -L ../build/lib -lOpenImageIO -lOpenImageIO_Util -o bin/volume
clean:
	rm -rf *.o bin/* base/*.o core ./lib/libstarter.a  *~ swig/*~ swig/*.so swig/*.o swig/*.cxx swig/*.pyc swig/bishop.py* doc/html doc/latex python/*.pyc 


