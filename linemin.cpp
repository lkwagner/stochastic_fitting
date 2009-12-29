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
using namespace std;



//------------------------------------------------------------------

typedef vector<double>::iterator dit_t;

//------------------------------------------------------------------
inline void append_number_fixed(string & str, int num){
  char strbuff[100];
  sprintf(strbuff, "%03d", num);
  str+=strbuff;
}


//------------------------------------------------------------------
Data_point gen_data(Data_generator & pes, const vector <double> & direction, 
                    const vector <double> & start_pos, const double t, double desired_sigma) { 
  int ndim=pes.ndim();
  assert(ndim==start_pos.size());
  assert(ndim==direction.size());
  Data_point pt;
  vector <double> x(ndim);
  for(int d=0; d< ndim; d++) { 
    x[d]=start_pos[d]+t*direction[d];
  }
  //double f=pes.eval_pes(x);
  double f,err;
  pes.eval(x,desired_sigma, f,err);
  pt.t=t;
  pt.val=f;
  pt.inverr=1.0/err;
  cout << "data: " << pt.t << " " << pt.val << " " << err << endl;
  return pt;
}

//------------------------------------------------------------------

void generate_line(double range, int npts, Data_generator & pes,
                   Line_data & data) { 
  assert(data.start_pos.size()==data.direction.size());
  assert(data.start_pos.size()==pes.ndim());
  int ndim=pes.ndim();
  vector <double> x(ndim);
  double sigma=.1;
  for(int i=0; i< npts; i++) { 
    double t=(i-npts/2)*range/npts;
    data.data.push_back(gen_data(pes,data.direction,data.start_pos,t,sigma));
  }
}


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



//------------------------------------------------------------------
//Using a set of trial positions, estimate the RMS deviation.
double rms_deviation(const Line_model & mod, const Line_data & data, const vector <double> & c) { 
  double rms=0;
  int npoints=0;
  for(vector<Data_point>::const_iterator d=data.data.begin(); d!=data.data.end(); d++) { 
    double f=mod.func(c,d->t);
    rms+=(d->val-f)*(d->val-f);
    npoints++;
  }
  return sqrt(rms/npoints);
  
}

//------------------------------------------------------------------
//This returns true if the values are significantly different, fairly arbitrary, 
//but it should be several sigmas.
bool sig_diff(Data_point & p1, Data_point & p2) { 
  return fabs(p1.val-p2.val) > 4.0*sqrt(1.0/(p1.inverr*p1.inverr)+1.0/(p2.inverr*p2.inverr));
}

void find_minimum(double range, int maxpts, Data_generator & pes,
                  Line_model & mod, Line_data & data, Fit_info & finfo, 
                  vector <double> & sigma) { 
  assert(data.start_pos.size()==data.direction.size());
  assert(data.start_pos.size()==pes.ndim());
  assert(data.data.size()==0);
  assert(sigma.size()==data.direction.size());
  int ndim=data.direction.size();
  
  static int call_num=0;
  string logname="find_minimum";
  append_number_fixed(logname, call_num);
  ofstream log(logname.c_str());
  log.precision(15);
  double start_sigma=0;;
  for(int d=0; d< ndim; d++) start_sigma+=data.direction[d]*data.direction[d]*sigma[d];
  double start_t=-range;
  double end_t=range;
  Line_data tempdata=data;
  Data_point p1,p2,p3,p4;
  double golden=1.0; //(1+sqrt(5.0))/2.0;
  double base_t=0;
  //Let's do a search to find the minimum
  p1=gen_data(pes,data.direction,data.start_pos,start_t,start_sigma);
  p2=gen_data(pes,data.direction,data.start_pos,start_t+(end_t-start_t)/(1.0+golden),start_sigma);
  p3=gen_data(pes,data.direction,data.start_pos,end_t,start_sigma);
  vector <Data_point> allpts;
  allpts.push_back(p1); allpts.push_back(p2); allpts.push_back(p3);
  log << "#starting sigma " << start_sigma << endl;
  //First bracket the minimum.
  while(1) { 
    log << "#step! points at  " << p1.t << " " << p2.t << " " << p3.t << endl;
    if(sig_diff(p1,p2) && sig_diff(p2,p3) && p2.val < p1.val && p2.val < p3.val) { 
      log << "#found set of bracketing points " << p1.t << " " << p2.t << " " << p3.t << endl;
      base_t=p2.t;
      break;
    }
    else if(!sig_diff(p1,p2) || !sig_diff(p2,p3)) { 
      //First check to see if we're just on opposite sides of the minimum and recenter if so:
      int shifted_minimum=0;
      if(!sig_diff(p1,p2) && sig_diff(p2,p3)) { 
        log << "#checking for minimum between " << p1.t << " and " << p2.t << endl;

        p4=gen_data(pes,data.direction,data.start_pos,(p1.t+p2.t)/2.0,start_sigma);
        allpts.push_back(p4);
        if(sig_diff(p4,p1) && sig_diff(p4,p2)) {
          start_t=p4.t-range;
          end_t=p4.t+range;
          p2=p4;
          shifted_minimum=1;
        }
      }
      else if(sig_diff(p1,p2) && !sig_diff(p2,p3)) { 
        log << "#checking for minimum between " << p2.t << " and " << p3.t << endl;
        p4=gen_data(pes,data.direction,data.start_pos,(p2.t+p3.t)/2.0,start_sigma);
        allpts.push_back(p4);
        if(sig_diff(p4,p3) && sig_diff(p4,p2)) {
          start_t=p4.t-range;
          end_t=p4.t+range;
          p2=p4;
          shifted_minimum=1;
        }        
      }
      if(!shifted_minimum) { 
        start_sigma*=0.5;
        log << "#reducing sigma to " << start_sigma << " in bracketing routine \n";
        p2=gen_data(pes,data.direction,data.start_pos,start_t+(end_t-start_t)/(1.0+golden),start_sigma);
        allpts.push_back(p2);
      }
      p1=gen_data(pes,data.direction,data.start_pos,start_t,start_sigma);
      allpts.push_back(p1);
      p3=gen_data(pes,data.direction,data.start_pos,end_t,start_sigma);
      allpts.push_back(p3);
    }
    else {  
      if(p1.val < p2.val) {
        start_t-=range;
        end_t-=range;
        p4=gen_data(pes,data.direction,data.start_pos,start_t,start_sigma);
        allpts.push_back(p4);
        p3=p2;
        p2=p1;
        p1=p4;
      }
      if( p3.val < p2.val) {
        end_t+=range;
        start_t+=range;
        p4=gen_data(pes,data.direction,data.start_pos,end_t,start_sigma);
        allpts.push_back(p4);
        p1=p2;
        p2=p3;
        p3=p4;
      }
      log << "#changing range to " << start_t << " to " << end_t << endl;
    }
  }
  
  log << "#all points used in the bracketing " << endl;
  for(vector<Data_point>::iterator d=allpts.begin(); d!= allpts.end(); d++) { 
    log << d->t << " " << d->val << "  " << 1.0/d->inverr << endl;
  }
  log << endl;

  log << "#minimum value point at t= " << p2.t << " val= " << p2.val << " +/- " << 1.0/p2.inverr << endl;
  //Now that we've bracketed the minimum, we can fill in the data a little bit 
  //for the fit.  Note that we may throw away a little data potentially, but 
  //including it seems more trouble than it's worth and will slow down the Monte
  //Carlo later.
  data.data.clear();
  data.data.push_back(p1); data.data.push_back(p2); data.data.push_back(p3);
  data.data.push_back(gen_data(pes,data.direction,data.start_pos,p2.t+0.6*range,start_sigma));
  data.data.push_back(gen_data(pes,data.direction,data.start_pos,p2.t-0.6*range,start_sigma));
  data.data.push_back(gen_data(pes,data.direction,data.start_pos,p2.t+0.8*range,start_sigma));
  data.data.push_back(gen_data(pes,data.direction,data.start_pos,p2.t-0.8*range,start_sigma));
  data.data.push_back(gen_data(pes,data.direction,data.start_pos,p2.t+0.3*range,start_sigma));
  data.data.push_back(gen_data(pes,data.direction,data.start_pos,p2.t-0.3*range,start_sigma));

  sample(mod, data, finfo);
  base_t=finfo.min[0];
  
  double rms=rms_deviation(mod, data,finfo.cavg);
  double avgerr=0;
  int npoints=0;
  for(vector<Data_point>::const_iterator d=data.data.begin(); d!=data.data.end(); d++) { 
    avgerr+= 1.0/d->inverr;
    npoints++;
  }
  avgerr/=npoints;

  log << "#rms deviation " << rms <<  " average uncertainty " << avgerr << endl;
  
  //Check to see how many points are outside of 3 error bars..
  //It should actually happen occasionally, so we should be careful not to 
  //freak too much if it does.
  int noutside=0;
  for(vector<Data_point>::const_iterator d=data.data.begin(); d!=data.data.end(); d++) { 
    double f=mod.func(finfo.cavg,d->t);
    double deviation=fabs(f-d->val);
    double err=1.0/d->inverr;
    if(deviation > 3*err) {
      log << "#outside case!  t= " << d->t << " data point " << d->val << " +/- " << err 
      << " fitted value " << f << endl;
    }
    noutside++;
  }
  
  
  log << "#data " << endl;
  for(vector<Data_point>::const_iterator d=data.data.begin(); d!=data.data.end(); d++) { 
    log << d->t << " " << d->val <<  " " << 1.0/d->inverr << endl;
  }
  
  log << endl << "#average parameter fit " << endl;
  base_t=finfo.min[0];
  double begin=base_t-range*1.2, end=base_t+range*1.2;
  double res=(end-begin)/400.0;
  for(double t=begin; t < end; t+=res) { 
    log << t << " " << mod.func(finfo.cavg, t) << " 0.0 " << endl;
  }
  
  

  call_num++;
}

//------------------------------------------------------------------

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
        
      
      //generate_line(1.0,10,*pes,tdata);
      //sample(mod, tdata, finfo);
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
      //optimize_quad(quad, datas,models,c);
      /*
      vector <Walker> allwalkers;
      if(it==0) quad.generate_guess(datas, models, c);
      sample(quad, datas, models, finfo, c, allwalkers);
      cout << "*****************Average Hessian direction choosing\n";
      c=finfo.cavg;
       */
      quad.generate_guess(datas,models,c);
      vector <Fix_information> fixes;
      shake_quad(quad,datas,models,fixes,c);
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



