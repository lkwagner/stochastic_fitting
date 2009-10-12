
OBJS:= ulec.o macopt.o  MatrixAlgebrac.o Quad_plus_line.o Line_model.o Sample.o
HEADERS=Prob_model.h ulec.h macopt.h  Min.h Quad_plus_line.h Sample.h
CXXFLAGS:=-O2

all: linemin

linemin: $(OBJS) $(HEADERS) linemin.o
	$(CXX) -o $@ $(OBJS) linemin.o

clean :
	rm -f $(OBJS) linemin.o
