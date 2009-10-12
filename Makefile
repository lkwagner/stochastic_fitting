
OBJS:= Prob_model.o ulec.o macopt.o  MatrixAlgebrac.o
HEADERS=Prob_model.h ulec.h macopt.h Point.h Min.h
CXXFLAGS:=-O2

all: linemin

linemin: $(OBJS) $(HEADERS) linemin.o
	$(CXX) -o $@ $(OBJS) linemin.o

clean :
	rm -f $(OBJS) linemin.o
