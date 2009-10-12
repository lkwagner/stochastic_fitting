#include <cassert>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>

#include "Array.h"
#include "MatrixAlgebra.h"
#include "Line_model.h"
#include "Sample.h"
#include "ulec.h"
using namespace std;



//------------------------------------------------------------------

typedef vector<double>::iterator dit_t;



void append_number(string & str, int num){
  char strbuff[100];
  sprintf(strbuff, "%d", num);
  str+=strbuff;
}

//------------------------------------------------------------------

void generate_fake_line(Line_data & data) { 
  Morse_model mod;
  vector <double> c;
  c.resize(4);
  c[0]=2.1;
  c[1]=3.2;
  c[2]=.05;
  c[3]=4.5;
  int npts=5;
  double range=3.0;
  double sigma=.001;
  for(int i=-npts; i < npts; i++) { 
    double t= c[3]+i*range/npts;
    data.t.push_back(t);
    data.val.push_back(mod.func(c,t)+rng.gasdev()*sigma);
    data.err.push_back(sigma);
  }
  
  for(dit_t t=data.t.begin(),v=data.val.begin(),e=data.err.begin();
      t!=data.t.end(); 
      t++,v++,e++) { 
    cout << *t << " " << *v << " " << *e << endl;
  }
}


int main(int argc, char ** argv) { 
  Line_data data;
  generate_fake_line(data);
  Fit_info finfo;
  Morse_model mod;
  sample(mod,data, finfo);
  finfo.print(cout);
}
//------------------------------------------------------------------



