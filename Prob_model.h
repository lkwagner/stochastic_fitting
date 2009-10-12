#ifndef PROB_MODEL_H_INCLUDED
#define PROB_MODEL_H_INCLUDED

#include <vector>
#include <fstream>
#include <cmath>
#include "Point.h"
#include <iostream>
using namespace std;


struct Fit_info { 
  vector <double> cavg;   //the parameters
  vector <double> cer;
  vector <double> min;    //the average minimum.  May have several concatenated
  vector <double> minerr;
  vector <double> curv;   //curvature at the minimum
  vector <double> curverr; 
};

//--------------------------------------------------------------


class Prob_model { 
public:
  //c are the parameters(coefficients) in the model
  //p are the data points
  //returns log(P(D|M))
  virtual double probability(const vector <double> & c)=0;
  virtual void gradient(const vector <double> & c, 
                        vector<double> & grad)=0;
  virtual double gradient(const vector <double> & c, 
                        int dir)=0;
  
  virtual void read(istream & is)=0;
  virtual void generate_guess(vector<double> & c)=0;
  virtual bool is_ok(const vector <double> & c)=0;
  
  virtual void niceprint(const vector <double> & c, const vector <double> & cerr)=0;
  virtual void minimum(const vector <double> & c, vector <double> & min) { 
    cerr << "Minimum not supported (or doesn't exist) for this probability model" << endl;
    exit(1);
  }
  
  //report the curvature at the minumum of this set of parameters
  virtual void curvature(const vector <double> & c, vector <double> & curve) { 
    cerr << "Curvature not supported (or doesn't exist) for this probability model" << endl;
    exit(1);

  }
  
  virtual void set_data(const vector <Point> & dat)=0;
  virtual ~Prob_model() { }
};


//--------------------------------------------------------------

class Morse_model:public Prob_model { 
public:
  virtual double probability(const vector <double> & c);
  virtual void gradient(const vector <double> & c, 
                        vector<double> & grad);
  virtual double gradient(const vector <double> & c, 
                          int dir);
  virtual void read(istream & is);
  virtual bool is_ok(const vector <double> & c);
  virtual void niceprint(const vector <double> & c, const vector <double> & cerr);
  virtual void generate_guess(vector<double> & c);
  virtual void minimum(const vector <double> & c, vector <double>  & min);
  virtual void curvature(const vector <double> & c, vector <double> & curve);

  virtual void set_data(const vector <Point> & dat) { data=dat;}

private:
  vector <Point> data;
};


/*
 c[0]*x+c[1]
 */
class Linear_model:public Prob_model { 
public:
  virtual double probability(const vector <double> & c);
  virtual void gradient(const vector <double> & c, 
                        vector<double> & grad);
  virtual double gradient(const vector <double> & c, 
                          int dir);
  virtual void read(istream & is);
  virtual bool is_ok(const vector <double> & c);
  virtual void niceprint(const vector <double> & c, const vector <double> & cerr);
  virtual void generate_guess(vector<double> & c);
  virtual void set_data(const vector <Point> & dat) { data=dat;}

private:
    vector <Point> data;
};

//--------------------------------------------------------------

//###########################################################################


class Prob_distribution { 
public:
  Prob_distribution(int np, double mi,
                    double spac) {
    npoints=np;
    mymin=mi;
    spacing=spac;
    density.resize(npoints);
    for(vector<double>::iterator i=density.begin();
        i!=density.end(); i++) *i=0;
  }
  //-----------------------------------------------------
  
  void accumulate(double a, double prob) {
    int place=int((a-mymin+spacing/2)/spacing);
    //cout << a <<" place " << place << endl;
    if(place >0 && place<npoints){
      density[place]+=prob; 
      //cout << "prob " << prob << endl;
    }
  }
  
  //-----------------------------------------------------

  void write(const string name) { 
    ofstream out(name.c_str());
    double sum=0;
    for(vector<double>::iterator i=density.begin();
        i!=density.end(); i++) sum+=*i;
    
    for(int i=0; i< npoints; i++) {
      out << mymin+spacing*i << "   " << density[i]/sum << endl;
    }
  }
  
  
  
  //-----------------------------------------------------
  void avg(double & avg, double & err);
  //-----------------------------------------------------

  double skew();
  
  void clear() { 
    for(vector<double>::iterator i=density.begin();
        i!= density.end(); i++) *i=0.0;
  }
  
private:
  vector <double> density;
  double mymin;
  double spacing;
  int npoints;
};

//###########################################################################

int check_model_consistency(Prob_model & mod, vector <double> & c);


//c should contain the best guess..will be returned as the maximum probability point
//data should be data on a line..
//cavg will be the average min, and cerr the variance
void sample_min(vector <double> & c, Prob_model & mod, 
                Fit_info & finfo, int verbose=1 );
//------------------------------------------------------------------


#endif //PROB_MODEL_H_INCLUDED
