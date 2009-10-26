#ifndef SAMPLE_H_INCLUDED 
#define SAMPLE_H_INCLUDED
#include "Line_model.h"
#include "Quad_plus_line.h"
#include <iostream>
using namespace std;

struct Fit_info { 
  vector <double> cavg;   //the parameters
  vector <double> cer;
  vector <double> min;    //the average minimum.  May have several concatenated
  vector <double> minerr;
  void print(ostream & os);
};


void find_good_guess(Line_model & mod, const Line_data & data,  Fix_information & fix,
                     vector <double> & c); 


void sample(Line_model & mod, const Line_data & line, Fit_info & finfo,
            int verbose=1);
  

void optimize_quad(Quad_plus_line & quad, vector <Line_data> & data, vector <Line_model *> & models,
                   vector <double> & c);

//startc can be a reasonable guess from a previous fit
void sample(Quad_plus_line & quad, vector <Line_data> & data, vector <Line_model *> & models, 
            Fit_info & finfo, vector <double> & startc, int verbose=1);

#endif //SAMPLE_H_INCLUDED