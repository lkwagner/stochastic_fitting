#ifndef DATA_GENERATOR_H_INCLUDED
#define DATA_GENERATOR_H_INCLUDED

#include <vector>
#include <fstream>

using namespace std; 

class Data_generator { 
public:
  virtual void eval(const vector <double> & x, const double desired_err, 
                    double & f, double & err)=0;
  virtual int ndim()=0;
};

//------------------------------------------------------------------

#include "PES.h"

class Random_quadratic:public Data_generator { 
public:
  virtual void eval(const vector <double> & x, const double desired_err, 
                    double & f, double & err) { 
    f=pes.eval_pes(x)+rng.gasdev()*desired_err;
    err=desired_err;
  }
  virtual int ndim() { return pes.ndim(); } 
  Random_quadratic(int n) { 
    pes.gen_pes(n);
  }
private:
  Potential_energy_surface pes;
};

//------------------------------------------------------------------
/*
class MonteCarlo_caller:public Data_generator { 
public:
  virtual void eval(const vector <double> & x, const double desired_err, 
                    double & f, double & err);
  virtual int ndim()  { return n; } 
  MonteCarlo_caller(istream & is);
private:
  int n;
};
 */

//------------------------------------------------------------------

#endif //DATA_GENERATOR_H_INCLUDED