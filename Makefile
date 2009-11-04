
OBJS:= ulec.o macopt.o  MatrixAlgebrac.o Quad_plus_line.o Line_model.o Sample.o PES.o Data_generator.o
HEADERS=Array.h Line_model.h MatrixAlgebra.h Min.h Point.h Quad_plus_line.h Sample.h macopt.h nrutil.h r.h rand.h ulec.h Data_generator.h
#CXXFLAGS:=-O2 -DNDEBUG  -DUSE_BLAS -framework Accelerate
CXXFLAGS:= -O2 -DNDEBUG 

all: linemin

linemin: $(OBJS) $(HEADERS) linemin.o
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) linemin.o

%.o:%.cpp  $(HEADERS)
	$(CXX) -c  $(CXXFLAGS)  -o $@ $<

clean :
	rm -f $(OBJS) linemin.o
