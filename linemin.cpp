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
#include "Quad_plus_line.h"
#include "Min.h"
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
  data.direction.resize(1);
  data.start_pos.resize(1);
  data.direction[0]=1.0;
  data.start_pos[0]=0.0;
  for(dit_t t=data.t.begin(),v=data.val.begin(),e=data.err.begin();
      t!=data.t.end(); 
      t++,v++,e++) { 
    cout << *t << " " << *v << " " << *e << endl;
  }
}

//------------------------------------------------------------------
void test_quad_plus_line(vector <Line_data> & data, vector <Line_model *> & models) { 
  Quad_plus_line quad;
  vector <double> c;
  quad.generate_guess(data,models, c);
  
  int niterate=1000;
  vector <double> best_c=c;
  
  double best_p=-1e99;
  
  int nit=1000;
  int nparms=c.size();
  
  for(int i=0; i< nit; i++) { 
    Quad_opt opt(nparms,0,.1,50,1);
    quad.generate_guess(data,models,opt.c);
    opt.mod=&quad;    
    opt.data=&data;
    opt.models=&models;
    opt.optimize();
    double p=quad.prob(data, models,opt.c);
    cout.flush();
    //cout << "prob " << p << endl;
    if(p > best_p ) { 
      best_c=opt.c;
      best_p=p;
      if(1) { 
        cout << "\n new best ";
        for(dit_t i=opt.c.begin(); i!= opt.c.end(); i++) cout << *i << " ";
        cout << p << endl;
      }
    }
  }
  c=best_c;
  
}

//------------------------------------------------------------------



int main(int argc, char ** argv) { 
  Line_data data;
  generate_fake_line(data);
  Fit_info finfo;
  Morse_model mod;
  //sample(mod,data, finfo);
  finfo.print(cout);
  
  vector <Line_data> datas;
  datas.push_back(data);
  vector <Line_model *> models;
  models.push_back(&mod);
  test_quad_plus_line(datas,models);
  
}
//------------------------------------------------------------------



