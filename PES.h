#ifndef PES_H_INCLUDED
#define PES_H_INCLUDED

#include "Array.h"
#include <vector>
using namespace std;

class Potential_energy_surface { 
public:
  void gen_pes(int n);
  int ndim() { return C.GetDim(0); } 
  double eval_pes(const vector <double> & x);
  void get_minima(Array1 <double> & min) { min=minima; } 
private:
  Array2 <double> C;
  Array1 <double> minima;
  Array1 <double> offdiag;
};

#endif //PES_H_INCLUDED