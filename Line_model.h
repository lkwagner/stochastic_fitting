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
  //vector <double> t;  //point n is start_pos+t*direction
  //vector <double> val;
  //vector <double> err; 
  vector <Data_point> data;
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
  virtual const int nparms(int fix)=0;
  //Note the c here is the full one, ie, without any fix. Use convert_c if you want 
  //to use this with fixes, same applies to the minimum..which you should know anyway
  //if you fix it!
  virtual const double func(const vector <double> & c, double t)=0 ;
  virtual const void minimum(const vector <double> & c, vector <double> & min)=0;
  
};


class Morse_model: public Line_model { 
public:
  virtual const void generate_guess(const Line_data &,const Fix_information &, vector <double> & c);
  virtual const void minimum(const vector <double> & c, vector <double> & min);
  virtual const int nparms(int fix) { if(fix) return 1; else return 4; }
  virtual const double func(const vector <double> & c,double t);
  virtual const void convert_c(const Fix_information &, const vector <double> & c_in, 
                                vector <double> & c_out);
  
};



class Quadratic_model: public Line_model { 
public:
  virtual const void generate_guess(const Line_data &,const Fix_information &, vector <double> & c);
  virtual const void minimum(const vector <double> & c, vector <double> & min);
  virtual const int nparms(int fix) { if(fix) return 0; else return 3; }
  virtual const double func(const vector <double> & c,double t);
  virtual const void convert_c(const Fix_information &, const vector <double> & c_in, 
                               vector <double> & c_out);
  
};


#endif //LINE_MODEL_H_INCLUDED