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
using namespace std;



//------------------------------------------------------------------

typedef vector<double>::iterator dit_t;




//------------------------------------------------------------------
/*
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
*/

void generate_line(double range, int npts, Potential_energy_surface & pes,
                   Line_data & data) { 
  assert(data.start_pos.size()==data.direction.size());
  assert(data.start_pos.size()==pes.ndim());
  int ndim=pes.ndim();
  vector <double> x(ndim);
  double sigma=.01;
  for(int i=0; i< npts; i++) { 
    double t=(i-npts/2)*range/npts;
    for(int d=0; d< ndim; d++) { 
      x[d]=data.start_pos[d]+t*data.direction[d];
    }
    double f=pes.eval_pes(x);
    Data_point pt;
    pt.t=t;
    pt.val=f+rng.gasdev()*sigma;
    pt.inverr=1.0/sigma;
    data.data.push_back(pt);
  }
  
  /*
  
  for(dit_t t=data.t.begin(),v=data.val.begin(),e=data.err.begin();
      t!=data.t.end(); 
      t++,v++,e++) { 
    cout << "gendata: " << data.start_pos[0]+(*t)*data.direction[0] << "  " 
         << *t << " " << *v << " " << *e << endl;
  }
  cout << "gendata: \n";
  cout << "end " << endl;
  */
}

//------------------------------------------------------------------

//------------------------------------------------------------------

void update_directions(vector < vector < double> > & hess, 
                       vector < vector <double> > & directions) { 
  int n=hess.size();
  assert(hess[0].size()==n);
  assert(directions.size()==n);
  assert(directions[0].size()==n);
  
  Array2 <double> H(n,n);
  Array2 <double> Ev(n,n);
  Array1 <double> evals(n);
  
  for(int i=0;i < n; i++) {
    for(int j=0; j< n; j++) { 
      H(i,j)=hess[i][j];
    }
  }
  EigenSystemSolverRealSymmetricMatrix(H,evals, Ev);

  cout << "eigenvalues: ";
  for(int i=0; i< n; i++) cout << evals(i) << " ";
  cout << endl;
  
  cout << "New directions : " << endl;
  for(int i=0; i< n; i++) { 
    for(int j=0; j< n; j++) { 
      directions[i][j]=Ev(j,i); 
      cout << Ev(i,j) << " ";
    }
    cout << endl;
  }

}

//------------------------------------------------------------------
void check_directions(Quad_plus_line & quad, vector <Line_model *> &  models,
                      vector <Line_data> & datas, vector <double> & c) { 
  int nfit=c.size();
  double baseprob=quad.prob(datas, models, c);
  double del=1e-6;
  for(int i=0; i< nfit; i++) { 
    c[i]+=del;
    double plusprob=quad.prob(datas,models,c);
    c[i]-=2.0*del;
    double minusprob=quad.prob(datas,models,c);
    c[i]+=del;
    double curve=(plusprob+minusprob-2.0*baseprob)/(del*del);
    cout << "curvature  " << i << " : " << curve << endl;
  }
}


//------------------------------------------------------------------

int main(int argc, char ** argv) { 
  Line_data data;
  Fit_info finfo;
  Quadratic_model mod;
  
  
  int n=2;
  if(argc > 1) { 
    n=atoi(argv[1]);
  }
  Potential_energy_surface pes;
  pes.gen_pes(n);
  data.start_pos.resize(n);
  data.direction.resize(n);
  for(int i=0; i< n; i++) data.start_pos[i]=.1;
  data.direction[0]=1.0;
  generate_line(1.0,20,pes,data);
  
  
  //sample(mod,data, finfo);
  //finfo.print(cout);
  vector <Line_model *> models;

  vector <Line_data> datas;
  
  Quad_plus_line quad;
  vector <double> c;
  vector <double> currmin(n);
  for(int i=0; i< n; i++) { currmin[i]=.1; } 
  int nit=1; 
  vector < vector < double> > directions(n);
  for(int i=0; i < n; i++) directions[i].resize(n);
  for(int i=0; i< n; i++) 
    for(int j=0; j< n; j++) directions[i][j]= (i==j)?1.0:0.0;
  
  for(int it=0; it < nit; it++) { 
    for(int d=0; d< n; d++) { 
      Line_data tdata;
      tdata.start_pos.resize(n);
      tdata.direction.resize(n);
      tdata.start_pos=currmin;
      tdata.direction=directions[d];
      generate_line(1.0,10,pes,tdata);
      sample(mod, tdata, finfo);
      for(int i=0; i< n; i++) { 
        currmin[i]+=finfo.min[0]*directions[d][i];
      }
      cout << "currmin_from_line ";
      for(int i=0; i< n; i++) { cout << currmin[i] << "  "; } 
      cout << endl;      
      
      datas.push_back(tdata);
      models.push_back(&mod);
    }
    //
    int use_quad=1;
    if(use_quad) { 
      //optimize_quad(quad, datas,models,c);
      if(it==0) quad.generate_guess(datas, models, c);
      sample(quad, datas, models, finfo, c);
      c=finfo.cavg;
      vector <double> quadmin(n);
      quad.get_minimum(c,n,quadmin);
      vector < vector <double> > hess;
      quad.get_hessian(c,n,hess);
      
      cout << "currmin_from_quad: ";
      for(int i=0; i< n; i++) { cout << quadmin[i] << "  "; } 
      cout << endl;
      cout << "hessian: \n";
      for(int i=0; i< n; i++) {
        for(int j=0; j< n; j++) { cout << hess[i][j] << " "; }
        cout << endl;
      }
      check_directions(quad, models, datas, c);
      update_directions(hess,directions);
    }
  }
  
}
//------------------------------------------------------------------



