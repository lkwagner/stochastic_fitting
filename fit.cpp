#include <cassert>
#include <cmath>
#include <vector>
#include <iostream>
#include "ulec.h"
#include "Point.h"
#include "Prob_model.h"
#include "macopt.h"
#include "Min.h"
using namespace std;



//------------------------------------------------------------------

typedef vector<double>::iterator dit_t;



void append_number(string & str, int num){
  char strbuff[100];
  sprintf(strbuff, "%d", num);
  str+=strbuff;
}

/*
 File format:
 pos1 pos2  en err
 etc.
 */

typedef vector<double>::iterator dit_t;
int main(int argc, char ** argv) { 
  cout.precision(10);
  double min, min_err;
  
  Morse_model mod;
  mod.read(cin);
  vector <double> c;
  mod.generate_guess(c);
  vector <double> cerr=c;
  for(dit_t i=cerr.begin(); i!= cerr.end(); i++) *i=0;
  mod.niceprint(c,cerr);

  vector <double> best_c;
  double best_p=-1e99;
  
  int nit=1000;
  int nparms=c.size();
  
  for(int i=0; i< nit; i++) { 
    Least_squares_opt opt(nparms,0,.1,50,1);
    mod.generate_guess(opt.c);
    opt.mod=&mod;    
    opt.optimize();
    double p=mod.probability(opt.c);
    cout.flush();
    cout << "prob " << p << endl;
    if(p > best_p && mod.is_ok(opt.c) ) { 
      best_c=opt.c;
      best_p=p;
      cout << "\n new best ";
      for(dit_t i=opt.c.begin(); i!= opt.c.end(); i++) cout << *i << " ";
      cout << p << endl;
    }
  }
 

  
  cout << "final from maximum likelihood: \n";
  mod.niceprint(best_c, cerr);

  check_model_consistency(mod,  best_c);
  
  sample_min(best_c, mod);  
  

}




//------------------------------------------------------------------
