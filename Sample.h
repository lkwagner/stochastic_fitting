#ifndef SAMPLE_H_INCLUDED 
#define SAMPLE_H_INCLUDED
#include "Line_model.h"
#include "Quad_plus_line.h"
#include "Hess_grad.h"
#include <iostream>
using namespace std;

struct Fit_info { 
  vector <double> cavg;   //the parameters
  vector <double> cer;
  vector <double> min;    //the average minimum.  May have several concatenated
  vector <double> minerr;
  double curve, curveerr; //curvature at the first minimum
  double funcmin, funcminerr; //Error at the function minimum
  void print(ostream & os);
};


struct Walker { 
  vector <double> c;
  double prob;
};



void find_good_guess(Line_model & mod, const Line_data & data,  Fix_information & fix,
                     vector <double> & c); 


void sample(Line_model & mod, const Line_data & line, Fit_info & finfo,
            vector <Walker> & allwalkers, int verbose=1);
  

void optimize_quad(Quad_plus_line & quad, vector <Line_data> & data, vector <Line_model *> & models,
                   vector <double> & c);

void shake_quad(Quad_plus_line & quad, vector <Line_data> & data, 
                vector <Line_model *> & models, vector <Fix_information> & fixes,
                vector <double> & c);  
//startc can be a reasonable guess from a previous fit
void sample(Quad_plus_line & quad, vector <Line_data> & data, vector <Line_model *> & models, 
            Fit_info & finfo, vector <double> & startc, vector <Walker> & allwalkers, int verbose=1);

void shake_quad_grad(vector<Line_data> & lines, 
    vector<Line_model *> & models,
    vector<Gradient_data> & gradients, vector<double> & c,int restart);



#endif //SAMPLE_H_INCLUDED
