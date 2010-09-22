#ifndef LINE_MODEL_H_INCLUDED
#define LINE_MODEL_H_INCLUDED
#include <vector>
#include <iostream>
#include <cassert>
#include "Array.h"
using namespace std;

struct Data_point { 
public:
  double t;
  double val;
  double inverr;  // 1/sigma, which speeds us up quite a bit
};

class Line_data { 
public:
  vector <double> direction; //the direction in which we searched (n-fold)
  vector <double> start_pos; 
  vector <Data_point> data;
  void store(ostream & os);
  void read(istream & is); 
};

struct Fix_information { 
  int enforce;
  double min;  //Location of the minimum (in units of t)
  double valmin; //Value of the objective function at the minimum
  double curve; // d^2E/dt^2 at the minimum
  Fix_information() { enforce=0; }
};




class Line_model { 
public:
  virtual const void generate_guess(const Line_data &,const Fix_information &, vector <double> & c)=0;
  virtual const double prob(const Line_data & ,const Fix_information &, const vector <double> & c);
  virtual const void convert_c(const Fix_information &, const vector <double> & c_in, 
                          vector <double> & c_out) { error("fix not implemented\n"); }
  //convert c from nonfixed to fixed c
  virtual const void downconvert_c(const Fix_information &, 
                                   const vector <double> & c_in, vector <double> & c_out){
    error("downconversion not supported for this Line_model");
  }
  virtual const int nparms(int fix)=0;
  //Note the c here is the full one, ie, without any fix. Use convert_c if you want 
  //to use this with fixes, same applies to the minimum..which you should know anyway
  //if you fix it!
  virtual const double func(const vector <double> & c, double t) const = 0;
  virtual const void minimum(const vector <double> & c, vector <double> & min)=0;
  //return the curvature at the minum;
  virtual double curve(const vector <double> & c) { error("Curvature not supported for this Line_model"); }
  virtual double funcmin(const vector <double> & c) { error("function min not supported for this Line_model"); }
};


class Morse_model: public Line_model { 
public:
  virtual const void generate_guess(const Line_data &,const Fix_information &, vector <double> & c);
  virtual const void minimum(const vector <double> & c, vector <double> & min);
  virtual const int nparms(int fix) { if(fix) return 1; else return 4; }
  virtual const double  func(const vector <double> & c,double t) const ;
  virtual const void convert_c(const Fix_information &, const vector <double> & c_in, 
                                vector <double> & c_out);
  
};

class Cubic_model: public Line_model { 
public:
  virtual const void generate_guess(const Line_data &,const Fix_information &, vector <double> & c);
  virtual const void minimum(const vector <double> & c, vector <double> & min);
  virtual const int nparms(int fix) { if(fix) return 1; else return 4; }
  virtual const double  func(const vector <double> & c,double t) const ;
  virtual const void convert_c(const Fix_information &, const vector <double> & c_in, 
                               vector <double> & c_out);
  virtual double curve(const vector <double> & c);
  virtual double funcmin(const vector <double> & c);
  virtual const void downconvert_c(const Fix_information &, 
                                   const vector <double> & c_in, vector <double> & c_out);


};


class Quadratic_model: public Line_model { 
public:
  virtual const void generate_guess(const Line_data &,const Fix_information &, vector <double> & c);
  virtual const void minimum(const vector <double> & c, vector <double> & min);
  virtual const int nparms(int fix) { if(fix) return 0; else return 3; }
  virtual const double  func(const vector <double> & c,double t) const;
  virtual const void convert_c(const Fix_information &, const vector <double> & c_in, 
                               vector <double> & c_out);

  virtual double curve(const vector <double> & c);
  virtual const void downconvert_c(const Fix_information &, 
                                   const vector <double> & c_in, vector <double> & c_out);
  
};

class Linear_model: public Line_model { 
public:
  virtual const void generate_guess(const Line_data &,const Fix_information &, vector <double> & c);
  virtual const void minimum(const vector <double> & c, vector <double> & min) { 
    min.resize(1); min[0]=0;
  }
  virtual const int nparms(int fix) { if(fix) return 0; else return 2; }
  virtual const double  func(const vector <double> & c,double t) const {
    return c[0]+c[1]*t;
  }
  virtual const void convert_c(const Fix_information &, const vector <double> & c_in, 
                               vector <double> & c_out) { 
    error("No conversion for linear model"); } 
  virtual double curve(const vector <double> & c) { return 0; } 
  virtual const void downconvert_c(const Fix_information &, 
                                   const vector <double> & c_in, vector <double> & c_out)
  {  error("No downconversion for linear model");}
  
  
};




#endif //LINE_MODEL_H_INCLUDED
