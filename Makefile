
OBJS:= ulec.o macopt.o  MatrixAlgebrac.o Quad_plus_line.o Line_model.o Sample.o PES.o
HEADERS=Array.h Line_model.h MatrixAlgebra.h Min.h Point.h Quad_plus_line.h Sample.h macopt.h nrutil.h r.h rand.h ulec.h
CXXFLAGS:=-O2

all: linemin

linemin: $(OBJS) $(HEADERS) linemin.o
	$(CXX) -o $@ $(OBJS) linemin.o

%.o:%.cpp  $(HEADERS)
	$(CXX) -c  $(CXXFLAGS)  -o $@ $<

clean :
	rm -f $(OBJS) linemin.o
