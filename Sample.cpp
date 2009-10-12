#include "Sample.h"
#include <cmath>
#include "ulec.h"
#include <iostream> 
#include "Min.h"

typedef vector<double>::iterator dit_t;


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
      if(1) { 
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
