#ifndef LINE_MODEL_H_INCLUDED
#define LINE_MODEL_H_INCLUDED
#include <vector>
#include <iostream>
#include <cassert>
#include "Array.h"
using namespace std;

class Line_data { 
public:
  vector <double> direction; //the direction in which we searched (n-fold)
  vector <double> start_pos; 
  vector <double> t;  //point n is start_pos+t*direction
  vector <double> val;
  vector <double> err; 
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
  virtual void generate_guess(const Line_data &,const Fix_information &, vector <double> & c)=0;
  virtual double prob(const Line_data & ,const Fix_information &, const vector <double> & c);
  virtual void convert_c(const Fix_information &, const vector <double> & cin, 
                         const vector <double> & cout) { error("fix not implemented\n"); } 
  
  //Note the c here is the full one, ie, without any fix. Use convert_c if you want 
  //to use this with fixes, same applies to the minimum..which you should know anyway
  //if you fix it!
  virtual double func(const vector <double> & c, double t)=0 ;
  virtual void minimum(const vector <double> & c, vector <double> & min)=0;

};


class Morse_model: public Line_model { 
public:
  virtual void generate_guess(const Line_data &,const Fix_information &, vector <double> & c);
  virtual void minimum(const vector <double> & c, vector <double> & min);
  virtual double func(const vector <double> & c,double t);
  
};



#endif //LINE_MODEL_H_INCLUDED