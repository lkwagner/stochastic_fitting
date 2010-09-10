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



//------------------------------------------------------------------



void update_directions(vector < vector < double> > & hess, 
                       vector <double> & evals_return,
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
  evals_return.resize(n);
  for(int i=0; i< n;i++) evals_return[i]=evals[i];
  for(int i=0; i< n; i++) { 
    for(int j=0; j< n; j++) { 
      directions[i][j]=Ev(j,i); 
    }
  }
  

}

//------------------------------------------------------------------


void average_directions(Quad_plus_line & quad,vector <Walker> &  allwalkers,
                        vector <double> & evals,vector <vector <double> > & directions) { 
  int n=directions.size();
  assert(directions.size()==directions[0].size());

  vector < vector <double> > hess;
  vector <double> tevals(n);
  vector <double> posevals(n);
  vector < vector <double> > tdirections(n);
  vector <vector <double> > posdirections(n);
  for(vector< vector <double> >::iterator i=tdirections.begin(); 
      i!= tdirections.end(); i++) i->resize(n);
  for(vector< vector <double> >::iterator i=posdirections.begin(); 
      i!= posdirections.end(); i++) i->resize(n);
  
  vector <ofstream * > eigenout(n);
  for(int i=0; i< n; i++) { 
    string nm="eigen";
    append_number_fixed(nm,i);
    eigenout[i]= new ofstream(nm.c_str());
  }
  
  int nwalkers=allwalkers.size();
  for(int i=0; i< n; i++) { 
    evals[i]=0;
    posevals[i]=0;
    for(int j=0; j< n; j++) { 
      directions[i][j]=0;
      posdirections[i][j]=0;
    }
  }
  
  int npos=0;
  
  double bestprob=allwalkers[0].prob;
  double bestposprob=-1e99;
  for(vector<Walker>::iterator w=allwalkers.begin(); w!=allwalkers.end(); w++) {
    quad.get_hessian(w->c,n,hess);
    update_directions(hess,tevals, tdirections);
    if(w->prob > bestprob) { 
      directions=tdirections;
      evals=tevals;
      bestprob=w->prob;
    }
    
    
    //cout << "eval ";
    
    //for(int i=0; i< n; i++) { 
      //evals[i]+=tevals[i]/nwalkers;
      //cout << tevals[i] << " ";
      //for(int j=0; j< n; j++) { 
        //directions[i][j]+=tdirections[i][j]/nwalkers;
       // if(tevals[n-1] > 0) *eigenout[i] << j << " " << tdirections[i][j] << endl;
      //}
      
    //}
    //cout << endl;
    if(tevals[n-1]> 0) { //our eigenvector routine sorts the eigenvalues
      npos++;
      if(w->prob > bestposprob) { 
        posevals=tevals;
        posdirections=tdirections;
        bestposprob=w->prob;
      }
      //for(int i=0; i< n; i++) { 
      //  posevals[i]+=tevals[i];
      //  for(int j=0; j< n; j++) { 
       //   posdirections[i][j]+=tdirections[i][j];
       // }
      //}
    }
  }
  
  for(int i=0; i< n; i++) { 
    delete eigenout[i];
  }

  cout << "besteigen ";
  for(int i=0; i< n; i++) cout << evals[i] << " ";
  cout << endl;
  
  
  if(npos >0) { 
    cout << "Found positive definite matrices! prob " 
    << bestposprob << " best nonposdef " << bestprob << endl;
     //for(int i=0; i< n; i++) { 
     //  posevals[i]/=npos;
     //  for(int j=0; j< n; j++) { 
     //    posdirections[i][j]/=npos;
     //  }
     //}
    cout << "bestpos ";
    for(int i=0; i< n; i++) cout << posevals[i] << " ";
    cout << endl;
    evals=posevals;
    directions=posdirections;
  }
  
}



int main(int argc, char ** argv) { 
  
  int nit=15; 
  
  string dummy;
  Data_generator * pes=NULL;
  Line_model * mod=NULL;
  vector <double> currmin;
  double trust_rad=0.4;
  if(argc <=1 ) error("usage: linemin inputfile");
  ifstream in(argv[1]);
  int nsweeps_keep=0;
  int use_quad=1;
  int powell_update=0;
  double start_sigma=0.01;
  while(in >> dummy) { 
    if(dummy=="iterations") in >> nit;
    else if(dummy=="random_quad") { 
      int ndim; in >> ndim;
      pes=new Random_quadratic(ndim);
    }
    else if(dummy=="montecarlo_caller") { 
      pes=new MonteCarlo_caller(in);
    }
    else if(dummy=="random_cubic") { 
      pes=new Random_cubic(in);
    }
    else if(dummy=="start_sigma") in >> start_sigma;
    else if(dummy=="powell") { powell_update=1; use_quad=0; }
    else if(dummy=="morse_mod" && mod == NULL) mod=new Morse_model;
    else if(dummy=="cubic_mod" && mod == NULL) mod=new Cubic_model;
    else if(dummy=="trust_radius") in >> trust_rad;
    else if(dummy=="no_hessian") use_quad=0;
    else if(dummy=="seed") { long int s1,s2;
      in >> s1 >> s2;
      rng.seed(s1,s2);
    }
    else if(dummy=="currmin") {
      if(pes==NULL) error("PES not defined before minimum");
      int ndim=pes->ndim();
      double dum=0;
      for(int i=0; i< ndim; i++) {
        in >> dum; currmin.push_back(dum);
      }
    }
    else if(dummy=="nsweeps_keep") { in >> nsweeps_keep; } 
    else { cerr << "Didn't understand " << dummy << endl;
      exit(132);
    }
  }
  
  
  if(pes==NULL) pes=new Random_quadratic(2);
  if(mod==NULL) mod=new Quadratic_model;
  int n=pes->ndim();

  if(!nsweeps_keep) nsweeps_keep=pes->ndim();

  Fit_info finfo;

  
  Quad_plus_line quad;
  vector <double> c;
  //vector <double> currmin(n);
  vector <double> sigma(n);
  //for(int i=0; i< n; i++) { currmin[i]=10; } 
  vector < vector < double> > directions(n);
  for(int i=0; i < n; i++) directions[i].resize(n);
  for(int i=0; i< n; i++) sigma[i]=start_sigma;
  for(int i=0; i< n; i++) 
    for(int j=0; j< n; j++) directions[i][j]= (i==j)?1.0:0.0;
  
  cout << "begin " << endl;
  vector <Line_model *> models;
  
  vector <Line_data> datas;
  Data_point firstpoint=gen_data(*pes, directions[0],currmin,0,sigma[0]);
  double last_energy=firstpoint.val;
  cout << "energy " << last_energy << " +/- " << 1.0/firstpoint.inverr << endl;
  for(int it=0; it < nit; it++) { 
    vector <double> min_start=currmin;
    double max_decrease=0;
    int max_decrease_index=0;
    for(int d=0; d< n; d++) { 
      Line_data tdata;
      tdata.start_pos.resize(n);
      tdata.direction.resize(n);
      tdata.start_pos=currmin;
      tdata.direction=directions[d];
      find_minimum(trust_rad,10, *pes, *mod, tdata, finfo,sigma);
      
      //assuming that c[0] is E_0..dirty, but works for now..
      double val0=finfo.cavg[0];
      cout << "energy " << val0 << " +/- " << finfo.cer[0] 
           <<  " decrease " << last_energy-val0 << endl;
      if(last_energy-val0 > max_decrease) {
        max_decrease=last_energy-val0;
        max_decrease_index=d;
      }
      last_energy=val0;
        
      
      for(int i=0; i< n; i++) { 
        currmin[i]+=finfo.min[0]*directions[d][i];
      }
      cout << "currmin_from_line ";
      for(int i=0; i< n; i++) { cout << currmin[i] << "  "; } 
      cout << endl;      
      
      datas.push_back(tdata);
      models.push_back(mod);
      if(it >= nsweeps_keep) {
        datas.erase(datas.begin());
        models.erase(models.begin());
      }
    }
    //
    if(use_quad) { 
      quad.generate_guess(datas,models,c);
      vector <Fix_information> fixes;
      shake_quad(quad,datas,models,fixes,c);
      
      //--Just for making a "pretty" graph
      int totlines=models.size();
      cout << "------------------- All lines" << endl;
      for(int d=0; d< totlines; d++) { 
        vector <Walker> allwalkers;
        sample(*(models[d]), datas[d], finfo, allwalkers,1);
        cout << "line " << d  << endl;
        for(int i=0; i< min(10,int(allwalkers.size())); i++) { 
          vector <double> c=allwalkers[i].c;
          cout << datas[d].start_pos[0] << " + (u+"<< c[1] << ")*" << datas[d].direction[0] 
          << ", " << datas[d].start_pos[1] << " + (u+"<< c[1] << ")*" << datas[d].direction[1]
          << " , " <<  c[0] << " + " << c[2] << "*u*u+" << c[3]<< " *u*u*u " << endl;
          
        }
      }
      
      //-----done with the graph
      
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
      vector <double> evals;
      update_directions(hess,evals, directions);
      
      cout << "eigenvalues: ";
      for(int i=0; i< n; i++) cout << evals[i] << " ";
      cout << endl;
      
      cout << "New directions : " << endl;
      for(int i=0; i< n; i++) { 
        for(int j=0; j< n; j++) { 
          cout << directions[i][j] << " ";
        }
        cout << endl;
      }
      /*
      average_directions(quad, allwalkers,evals, directions);
      cout << "************Averaging over directions \n";
       
      cout << "eigenvalues: ";
      for(int i=0; i< n; i++) cout << evals[i] << " ";
      cout << endl;
      
      cout << "New directions : " << endl;
      for(int i=0; i< n; i++) { 
        for(int j=0; j< n; j++) { 
          cout << directions[i][j] << " ";
        }
        cout << endl;
      }
       */
      
    }
    if(powell_update) { 
      vector <double> totaldir(n);
      double norm=0;
      for(int d=0; d< n; d++) { 
        totaldir[d]=currmin[d]-min_start[d];
        norm+=totaldir[d]*totaldir[d];
      }
      norm=sqrt(norm);
      for(int d=0; d< n; d++) totaldir[d]/=norm;
      directions[max_decrease_index]=totaldir;
      Line_data tdata;
      tdata.start_pos.resize(n);
      tdata.direction.resize(n);
      tdata.start_pos=currmin;
      tdata.direction=totaldir;      
      find_minimum(trust_rad,10,*pes,*mod,tdata,finfo,sigma);
      
      //reset the directions to make sure we don't get stuck..
      if((it+1)%n == 0) { 
        for(int i=0; i< n; i++) 
          for(int j=0; j< n; j++) directions[i][j]= (i==j)?1.0:0.0;
      }
      
      cout << "energy " << finfo.cavg[0] << " +/- " << finfo.cer[0] 
      <<  " decrease " << last_energy-finfo.cavg[0] << " (finishing step) " <<  endl;
      last_energy=finfo.cavg[0];

      for(int i=0; i< n; i++) { 
        cout << "powell_dir ";
        for(int j=0; j< n; j++) cout << directions[i][j] << " ";
        cout << endl;
      }
      cout << "powell_dir \n";
    }
  }
  
  delete pes;
  delete mod;
}
//------------------------------------------------------------------



