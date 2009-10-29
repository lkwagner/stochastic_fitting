#include "Sample.h"
#include <cmath>
#include "ulec.h"
#include <iostream> 
#include "Min.h"

typedef vector<double>::iterator dit_t;

inline void append_number(string & str, int num){
  char strbuff[100];
  sprintf(strbuff, "%d", num);
  str+=strbuff;
}

void Fit_info::print(ostream & os) { 
   assert(cavg.size()==cer.size());
   for(dit_t a=cavg.begin(),e=cer.begin(); a!= cavg.end() && e!=cavg.end();
       a++,e++) { 
       os << *a << " +/- " << *e << endl;
   }
}


typedef vector<double>::iterator dit_t;

struct Walker { 
  vector <double> c;
  double prob;
};

double gradient(Line_model & mod, const Line_data & data, const Fix_information & fix,
    const vector <double> & c, int d) { 
  double prob=mod.prob(data, fix, c);
  double del=1e-12;
  vector <double> tmpc=c;
  tmpc[d]+=del;
  double px=mod.prob(data, fix, tmpc);
  return (px-prob)/del;
}

void find_good_guess(Line_model & mod, const Line_data & data,  Fix_information & fix,
      vector <double> & c) { 
  mod.generate_guess(data, fix,c);
  vector <double> best_c=c;
  
  double best_p=-1e99;
  
  int nit=1000;
  int nparms=c.size();
  
  
  
  for(int i=0; i< nit; i++) { 
    Least_squares_opt opt(nparms,0,.1,50,1);
    mod.generate_guess(data,fix,opt.c);
    opt.mod=&mod;    
    opt.fix=&fix;
    opt.data=&data;
    opt.optimize();
    double p=mod.prob(data, fix,opt.c);
    cout.flush();
//cout << "prob " << p << endl;
    if(p > best_p ) { 
      best_c=opt.c;
      best_p=p;
      if(0) { 
        cout << "\n new best ";
        for(dit_t i=opt.c.begin(); i!= opt.c.end(); i++) cout << *i << " ";
        cout << p << endl;
      }
    }
  }
  c=best_c;
}

void sample(Line_model & mod, const Line_data & data, Fit_info & finfo,
            int verbose) { 
  Fix_information fakefix;
  Walker walker;
  find_good_guess(mod,data, fakefix, walker.c);
  walker.prob=mod.prob(data, fakefix, walker.c);
  
  int nstep=100000;
  int warmup=nstep/3;
  Walker nw=walker;
  int nparms=walker.c.size();
  vector <double> tmpc;
  vector <double> tstep(nparms);
  for(int p=0; p < nparms; p++) tstep[p]=2e-12;
  double acc=0;
  //cout << "initial probability " << walker.prob << endl;
  vector <double> avgs(nparms);
  vector <double> vars(nparms);
  int navgpts=0;
  for(dit_t i=avgs.begin(); i!= avgs.end(); i++) *i=0.0;
  for(dit_t i=vars.begin(); i!= vars.end(); i++) *i=0.0;
  vector <double> min, curve;
  mod.minimum(walker.c,min);
  int nminima=min.size();
  finfo.min.resize(min.size());
  finfo.minerr.resize(min.size());
  for(dit_t i=finfo.min.begin(); i!= finfo.min.end(); i++) *i=0.0;
  for(dit_t i=finfo.minerr.begin(); i!= finfo.minerr.end(); i++) *i=0.0;
  
  
  for(int step=0; step < nstep; step++) { 
    nw=walker;
    for(int p=0; p < nparms; p++) { 
      double ststep=sqrt(tstep[p]);
      double ts=tstep[p];
      double delta=1e-9;
      double grad=gradient(mod, data, fakefix, walker.c, p);
      nw.c[p]=walker.c[p]+ststep*rng.gasdev()+ts*grad;
      double ngrad=gradient(mod, data, fakefix, nw.c,p);
      nw.prob=mod.prob(data, fakefix, nw.c);
      
      double diff=nw.c[p]-walker.c[p];
      double num=-(-diff-ts*ngrad)*(-diff-ts*ngrad);
      double den=-(diff-ts*grad)*(diff-ts*grad);
      double accprob=exp(nw.prob-walker.prob+(num-den)/(2*ts));
      if(accprob+rng.ulec()> 1.0) {
        acc++;
        walker=nw;
        if(tstep[p] < 1.0) 
          tstep[p]*=1.1;
      }
      else { 
        nw.c=walker.c;
        tstep[p]*=.9;
      }
    }
    
    if(step > warmup) { 
      for(int p=0; p < nparms; p++) { 
        double oldavg=avgs[p];
        double oldvar=vars[p];
        avgs[p]=oldavg+(walker.c[p]-oldavg)/(navgpts+1);
        if(navgpts >0) { 
          vars[p]=(1-1.0/navgpts)*oldvar+(navgpts+1)*(avgs[p]-oldavg)*(avgs[p]-oldavg);
        }
      }
      mod.minimum(walker.c,min);
      for(int m=0; m < nminima; m++) { 
        double oldmin=finfo.min[m];
        double oldvar=finfo.minerr[m];
        finfo.min[m]=oldmin+(min[m]-oldmin)/(navgpts+1);
        if(navgpts > 0) { 
          finfo.minerr[m]=(1-1.0/navgpts)*oldvar+(navgpts+1)*(finfo.min[m]-oldmin)*(finfo.min[m]-oldmin);
        }
      }
      
      
      navgpts++;
    }
  }
  
  for(int p=0; p < nparms; p++) 
    vars[p]=sqrt(vars[p]);
  
  if(verbose) { 
    cout << "acceptance " << acc/(nstep*nparms) << endl;
    cout << "timesteps   avg   err" << endl;
    for(int p=0; p < nparms; p++) {
      cout << tstep[p]  << "  " << avgs[p] << "  " << vars[p] << endl;
    }
    cout << endl;
    
  }
  
  finfo.cavg=avgs;
  finfo.cer=vars;

}


//#############################################################################


double curvature(Quad_plus_line & quad, vector <Line_model *> &  models,
                      vector <Line_data> & datas, vector <double> & c) { 
  int nfit=c.size();
  double baseprob=quad.prob(datas, models, c);
  double del=1e-6;
  double totalcurve=0.0;
  for(int i=0; i< nfit; i++) { 
    c[i]+=del;
    double plusprob=quad.prob(datas,models,c);
    c[i]-=2.0*del;
    double minusprob=quad.prob(datas,models,c);
    c[i]+=del;
    double curve=(plusprob+minusprob-2.0*baseprob)/(del*del);
    totalcurve+=curve;
    //cout << "curvature  " << i << " : " << curve << endl;
  }
  return totalcurve;
}

void optimize_quad(Quad_plus_line & quad, vector <Line_data> & data, 
                   vector <Line_model *> & models,
                   vector <double> & c) { 
  int nparms=c.size();
  Quad_opt opt(nparms,0,.01,50,1);
  opt.c=c;
  opt.mod=&quad;    
  opt.data=&data;
  opt.models=&models;
  opt.optimize();
  c=opt.c;
}

void shake_quad(Quad_plus_line & quad, vector <Line_data> & data, 
                    vector <Line_model *> & models,
                         vector <double> & c) { 
  quad.generate_guess(data,models, c);
  vector <double> best_c=c;
  optimize_quad(quad, data, models, c);
  double best_p=quad.prob(data,models, c);
  cout << "initial probability " << best_p << endl;
  int nit=1000;
  //int nit=10;
  int nparms=c.size();
  
  for(int i=0; i< nit; i++) { 
    quad.generate_guess(data,models,c);
    optimize_quad(quad, data, models, c);
    double p=quad.prob(data, models,c);
    cout.flush();
    //cout << "prob " << p << endl;
    if(p > best_p  ) { 
      best_c=c;
      best_p=p;
      if(1) { 
        cout << " it " << i << " ";
        for(dit_t i=c.begin(); i!= c.end(); i++) cout << *i << " ";
        cout << p  << endl;
      }
    }
  }
  c=best_c;
  
}

//------------------------------------------------------------------------------

double gradient(Quad_plus_line & quad, vector <Line_data> & data, vector <Line_model *> & models,
    const vector <double> & c, int d) { 
  double prob=quad.prob(data, models, c);
  double del=1e-12;
  vector <double> tmpc=c;
  tmpc[d]+=del;
  double px=quad.prob(data,models, tmpc);
  return (px-prob)/del;
}

//------------------------------------------------------------------------------

int metropolis_step(Quad_plus_line & quad, vector <Line_data> & data, 
                     vector <Line_model * > & models, vector <double> & tstep, 
                     Walker & walker, int p) { 
  Walker nw=walker;
  double ststep=sqrt(tstep[p]);
  double ts=tstep[p];
  double delta=1e-9;
  //double grad=0, ngrad=0;
  double grad=gradient(quad, data, models, walker.c, p);
  nw.c[p]=walker.c[p]+ststep*rng.gasdev()+ts*grad;
  double ngrad=gradient(quad, data, models, nw.c,p);
  nw.prob=quad.prob(data, models, nw.c);
  
  double diff=nw.c[p]-walker.c[p];
  double num=-(-diff-ts*ngrad)*(-diff-ts*ngrad);
  double den=-(diff-ts*grad)*(diff-ts*grad);
  double accprob=exp(nw.prob-walker.prob+(num-den)/(2*ts));
  if(accprob+rng.ulec()> 1.0) {
    walker=nw;
    if(tstep[p] < 1.0) 
      tstep[p]*=1.1;
    return 1;
  }
  else { 
    nw.c=walker.c;
    tstep[p]*=.9;
    return 0;
  }
  
}

//------------------------------------------------------------------------------



void sample(Quad_plus_line & quad, vector <Line_data> & data, vector <Line_model *> & models, 
            Fit_info & finfo, vector<double> & startc, int verbose) { 
  int nconfig=100;
  vector <Walker> configs(nconfig);
  Walker walker;
  quad.generate_guess(data, models, walker.c);
  assert(startc.size()==walker.c.size());
  walker.c=startc;
  shake_quad(quad,data, models, walker.c);
  walker.prob=quad.prob(data, models, walker.c);
  
  for(vector <Walker>::iterator i=configs.begin(); i!= configs.end(); i++)
    *i=walker;
  
  int nstep=10000;
  //int nstep=  10000;
  int decorr=100;
  int warmup=200;
  //Walker nw=walker;
  int nparms=walker.c.size();
  vector <double> tmpc;
  vector <double> tstep(nparms);
  for(int p=0; p < nparms; p++) tstep[p]=2e-2;
  double acc=0;
  //cout << "initial probability " << walker.prob << endl;
  vector <double> avgs(nparms);
  vector <double> vars(nparms);
  int navgpts=0;
  for(dit_t i=avgs.begin(); i!= avgs.end(); i++) *i=0.0;
  for(dit_t i=vars.begin(); i!= vars.end(); i++) *i=0.0;
  vector <double> min, curve;
  //quad.minimum(walker.c,min);
  //int nminima=min.size();
  finfo.min.resize(min.size());
  finfo.minerr.resize(min.size());
  for(dit_t i=finfo.min.begin(); i!= finfo.min.end(); i++) *i=0.0;
  for(dit_t i=finfo.minerr.begin(); i!= finfo.minerr.end(); i++) *i=0.0;
  vector <Walker> allwalkers;
  
  for(int step=0; step < nstep; step++) { 
    //nw=walker;
    for(vector<Walker>::iterator w=configs.begin(); w!=configs.end(); w++) { 
      for(int p=0; p < nparms; p++) { 
        if(metropolis_step(quad, data, models, tstep, *w, p)) acc++;
      }
    
      if(step > warmup && step%decorr==0) { 
        for(int p=0; p < nparms; p++) { 
          double oldavg=avgs[p];
          double oldvar=vars[p];
          avgs[p]=oldavg+(w->c[p]-oldavg)/(navgpts+1);
          if(navgpts >0) { 
            vars[p]=(1-1.0/navgpts)*oldvar+(navgpts+1)*(avgs[p]-oldavg)*(avgs[p]-oldavg);
          }
        }
        allwalkers.push_back(*w);
        navgpts++;
      }
    }
  }
  
  for(int p=0; p < nparms; p++) 
    vars[p]=sqrt(vars[p]);
  
  if(verbose) { 
    cout << "acceptance " << acc/(nstep*nparms) << endl;
    cout << "timesteps   avg   err" << endl;
    for(int p=0; p < nparms; p++) {
      cout << tstep[p]  << "  " << avgs[p] << "  " << vars[p] << endl;
    }
    cout << endl;
    
  }
  
  finfo.cavg=avgs;
  finfo.cer=vars;
  
  
  vector <double> min_c(nparms), max_c(nparms);
  for(int i=0; i< nparms; i++) { 
    min_c[i]=max_c[i]=allwalkers[0].c[i];
  }
  int npoints=50;
  for(vector<Walker>::iterator w=allwalkers.begin(); w!=allwalkers.end();
      w++) { 
    for(int i=0; i< nparms; i++) { 
      if(w->c[i] < min_c[i]) min_c[i]=w->c[i];
      if(w->c[i] > max_c[i]) max_c[i]=w->c[i];
    }
  }
  cout << "parameter information:: \n";
  cout << "number   min   max \n";
  for(int i=0; i< nparms; i++) { 
    cout << i << " " << min_c[i] << " " << max_c[i] << endl;
  }
  
  for(int i=0; i< nparms; i++) { 
    
    vector <double> counts(npoints);
    double tiny=1e-06;
    double invsize=1.0/(max_c[i]-min_c[i]+tiny);
    for(vector<double>::iterator j=counts.begin(); j!=counts.end(); j++) *j=0.0;
    
    string pathname="path";
    append_number(pathname,i);
    ofstream pathout(pathname.c_str());
    for(vector<Walker>::iterator w=allwalkers.begin(); w!=allwalkers.end();
        w++) { 
      int pos=int(npoints*(w->c[i]-min_c[i])*invsize);
      pathout << w->c[i] << endl;
      if(pos >=0 && pos < counts.size() ) {
        counts[pos]++;
      }
      else { cout << "wierd pos " << pos << " c " << w->c[i]  << endl; }
    }
    pathout.close();
     
     
    
    int nwalkers=allwalkers.size();
    string outname="dist";
    append_number(outname,i);
    ofstream out(outname.c_str());
    for(int j=0; j< npoints; j++) { 
      out << min_c[i]+j/(invsize*npoints) << "  " << counts[j]/nwalkers << endl;
    }
    out.close();
     
  }
}


