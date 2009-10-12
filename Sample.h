#ifndef SAMPLE_H_INCLUDED 
#define SAMPLE_H_INCLUDED
#include "Line_model.h"
#include <iostream>
using namespace std;

struct Fit_info { 
  vector <double> cavg;   //the parameters
  vector <double> cer;
  vector <double> min;    //the average minimum.  May have several concatenated
  vector <double> minerr;
  void print(ostream & os);
};


void sample(Line_model & mod, const Line_data & line, Fit_info & finfo,
            int verbose=1);
  


#endif //SAMPLE_H_INCLUDED