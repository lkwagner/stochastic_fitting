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
#include "PES.h"
#include "Data_generator.h"
#include "generate_line.h"
using namespace std;



//------------------------------------------------------------------

typedef vector<double>::iterator dit_t;

//------------------------------------------------------------------

void fit_line(Line_model & mod, Line_data & data, Fit_info & finfo) { 
  vector <Walker> allwalkers;
  cout << "sampling " << endl;
  sample(mod, data, finfo, allwalkers);
  cout << "done sample " << endl;
  double rms=rms_deviation(mod, data,finfo.cavg);
  double avgerr=0;
  int npoints=0;
  for(vector<Data_point>::const_iterator d=data.data.begin(); d!=data.data.end(); d++) { 
    avgerr+= 1.0/d->inverr;
    npoints++;
  }
  avgerr/=npoints;

  cout << "#rms deviation " << rms <<  " average uncertainty " << avgerr << endl;
  
  //Check to see how many points are outside of 3 error bars..
  //It should actually happen occasionally, so we should be careful not to 
  //freak too much if it does.
  int noutside=0;
  for(vector<Data_point>::const_iterator d=data.data.begin(); d!=data.data.end(); d++) { 
    double f=mod.func(finfo.cavg,d->t);
    double deviation=fabs(f-d->val);
    double err=1.0/d->inverr;
    if(deviation > 3*err) {
      cout << "#outside case!  t= " << d->t << " data point " << d->val << " +/- " << err 
      << " fitted value " << f << endl;
    }
    noutside++;
  }
  
  
  cout << "#data " << endl;
  for(vector<Data_point>::const_iterator d=data.data.begin(); d!=data.data.end(); d++) { 
    cout << d->t << " " << d->val <<  " " << 1.0/d->inverr << endl;
  }
  
  cout << endl << "#average parameter fit " << endl;
  double begin=data.data[0].t, end=(data.data.end()-1)->t;
  double res=(end-begin)/400.0;
  for(double t=begin; t < end; t+=res) { 
    cout << t << " " << mod.func(finfo.cavg, t) << " 0.0 " << endl;
  }
}

//------------------------------------------------------------------

int main(int argc, char ** argv) { 
  cout.precision(14);
  string dummy;
  Line_data data;
  data.direction.resize(1);
  data.direction[0]=1;
  data.start_pos.resize(1);
  data.start_pos[0]=1;
  Line_model * mod=NULL;
  if(argc <=1 ) error("usage: line_fit inputfile");
  ifstream in(argv[1]);
  while(in >> dummy) { 
    if(dummy=="morse_mod" && mod == NULL) mod=new Morse_model;
    else if(dummy=="cubic_mod" && mod == NULL) mod=new Cubic_model;
    else if(dummy=="quadratic_mod" && mod == NULL) mod=new Quadratic_model;
    else if(dummy=="linear_mod" && mod==NULL) mod=new Linear_model;
    else if(dummy=="seed") { long int s1,s2;
      in >> s1 >> s2;
      rng.seed(s1,s2);
    }
    else if(dummy=="data") { 
      while(in) { 
        Data_point pt;
        in >> pt.t >> pt.val >> pt.inverr;
        pt.inverr=1.0/pt.inverr;
        data.data.push_back(pt);
      }
    }
    else { cerr << "Didn't understand " << dummy << endl;
      exit(132);
    }
  }
  if(mod==NULL) mod=new Linear_model;
  if(data.data.size()==0) { 
    cerr << "Didn't find any data section!" << endl;
    exit(143);
  }
  cout << "done reading " << endl;
  Fit_info finfo;
  fit_line(*mod, data, finfo);
  
  
  delete mod;
}
//------------------------------------------------------------------



