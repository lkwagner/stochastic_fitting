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
inline void append_number_fixed(string & str, int num){
  char strbuff[100];
  sprintf(strbuff, "%03d", num);
  str+=strbuff;
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
  int nsweeps_keep=2;
  
  while(in >> dummy) { 
    if(dummy=="iterations") in >> nit;
    if(dummy=="random_quad") { 
      int ndim; in >> ndim;
      pes=new Random_quadratic(ndim);
    }
    if(dummy=="montecarlo_caller") { 
      pes=new MonteCarlo_caller(in);
    }
    if(dummy=="morse_mod" && mod == NULL) mod=new Morse_model;
    if(dummy=="cubic_mod" && mod == NULL) mod=new Cubic_model;
    if(dummy=="trust_radius") in >> trust_rad;
    if(dummy=="currmin") {
      if(pes==NULL) error("PES not defined before minimum");
      int ndim=pes->ndim();
      double dum=0;
      for(int i=0; i< ndim; i++) {
        in >> dum; currmin.push_back(dum);
      }
    }
    if(dummy=="nsweeps_keep") { in >> nsweeps_keep; } 
  }
  
  if(pes==NULL) pes=new Random_quadratic(2);
  if(mod==NULL) mod=new Quadratic_model;
  int n=pes->ndim();
  
  Fit_info finfo;

  
  Quad_plus_line quad;
  vector <double> c;
  //vector <double> currmin(n);
  vector <double> sigma(n);
  //for(int i=0; i< n; i++) { currmin[i]=10; } 
  vector < vector < double> > directions(n);
  for(int i=0; i < n; i++) directions[i].resize(n);
  for(int i=0; i< n; i++) sigma[i]=0.01;
  for(int i=0; i< n; i++) 
    for(int j=0; j< n; j++) directions[i][j]= (i==j)?1.0:0.0;
  
  cout << "begin " << endl;
  vector <Line_model *> models;
  
  vector <Line_data> datas;
  
  for(int it=0; it < nit; it++) { 
    
    for(int d=0; d< n; d++) { 
      Line_data tdata;
      tdata.start_pos.resize(n);
      tdata.direction.resize(n);
      tdata.start_pos=currmin;
      tdata.direction=directions[d];
      find_minimum(trust_rad,10, *pes, *mod, tdata, finfo,sigma);

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
      update_directions(hess,directions);

    }
  }
  
  //Try removing the history to see what happens..
  /*
  for(int it=0; it < nit; it++) { 
    cout << "####################Removing the first " << it << " iterations " << endl;
    for(int d=0; d< n; d++) {
      datas.erase(datas.begin());
      models.erase(models.begin());
    }
    quad.generate_guess(datas, models, c);
    sample(quad, datas, models, finfo, c);
    cout << "finished sample " << endl;
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
    update_directions(hess,directions);
  }
  */
  delete pes;
  delete mod;
}
//------------------------------------------------------------------



