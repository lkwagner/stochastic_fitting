#ifndef DATA_GENERATOR_H_INCLUDED
#define DATA_GENERATOR_H_INCLUDED

#include <vector>
#include <fstream>
#include "ulec.h"
#include <string>

using namespace std; 

class Data_generator { 
public:
  virtual void eval(const vector <double> & x, const double desired_err, 
                    double & f, double & err)=0;
  virtual int gradient(const vector <double> & x, vector <double> & grad, vector <double> & err) { 
    return 0;
  }

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
  virtual int gradient(const vector <double> & x, vector <double> & grad, vector <double> & err) { 
    grad.resize(x.size());
    err.resize(x.size());
    pes.eval_grad(x,grad);
  
    for(vector<double>::iterator i=grad.begin(); i!=grad.end(); i++) *i+=sigma_grad*rng.gasdev();
    for(vector<double>::iterator i=err.begin(); i!=err.end(); i++) *i=sigma_grad;
    return 1;
  }

  virtual int ndim() { return pes.ndim(); } 
  Random_quadratic(int n) { 
    pes.gen_pes(n);
    sigma_grad=0.1;
  }
private:
  Potential_energy_surface pes;
  double sigma_grad;
};

//------------------------------------------------------------------

class MonteCarlo_caller:public Data_generator { 
public:
  virtual void eval(const vector <double> & x, const double desired_err, 
                    double & f, double & err);
  virtual int ndim()  { return n; } 
  MonteCarlo_caller(istream & is);
private:
  int n;
  double estimated_rms;
  int min_blocks;
  string command; //command should take an argument that controls the number of blocks
                  //x will be written to mc_input in space-delimited format
                  //expects an output in mc_output of the function value and its error bar
};


//------------------------------------------------------------------


class Random_cubic:public Data_generator { 
public:
  virtual void eval(const vector <double> & x, const double desired_err, 
                    double & f, double & err) { 
    f=pes.eval_pes(x)+rng.gasdev()*desired_err;
    int n=cubic.GetDim(0);
    assert(pes.ndim()==n);
    assert(x.size() >= n);
    Array1 <double> minima;
    pes.get_minima(minima);
    for(int i=0; i< n; i++) {
      for(int j=0; j< n; j++) { 
        for(int k=0; k< n; k++) { 
          f+=cubic(i,j,k)*(x[i]-minima[i])*(x[j]-minima[j])*(x[k]-minima[k]);
        }
      }
    }
    err=desired_err;
  }
  virtual int ndim() { return pes.ndim(); } 
  Random_cubic(istream & is) { 
    int n;
    double cube;
    is >> n >> cube;
    pes.gen_pes(n);
    cubic.Resize(n,n,n);
    for(int i=0; i< n; i++) { 
      for(int j=0; j< n; j++) { 
        for(int k=0; k< n; k++) { 
          cubic(i,j,k)=cube*rng.gasdev();
        }
      }
    }
  }
private:
  Potential_energy_surface pes;
  Array3 <double> cubic;
};




//------------------------------------------------------------------

#endif //DATA_GENERATOR_H_INCLUDED
