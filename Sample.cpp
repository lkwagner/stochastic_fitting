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
    if(p > best_p  && mod.curve(opt.c) > 0) { 
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
  
  /*
  double sig=0.2;
  cout << "shaking with a Gaussian " << endl;
  for(int i=0; i< nit; i++) { 
    for(int j=0; j< best_c.size(); j++) c[j]=best_c[j]*(1+sig*rng.gasdev());
    Least_squares_opt opt(nparms,0,.1,50,1);
    mod.generate_guess(data,fix,opt.c);
    opt.mod=&mod;    
    opt.fix=&fix;
    opt.data=&data;
    opt.optimize();
    double p=mod.prob(data, fix,opt.c);    
    cout.flush();
    //cout << "prob " << p << endl;
    if(p > best_p || isnan(best_p) ) { 
      best_c=opt.c;
      best_p=p;
      if(1) { 
        cout << " it " << i << " ";
        for(dit_t i=opt.c.begin(); i!= opt.c.end(); i++) cout << *i << " ";
        cout << p  << endl;
      }
    }
  }
   */
  c=best_c;
  
}

void sample(Line_model & mod, const Line_data & data, Fit_info & finfo,
             vector <Walker> & allwalkers, int verbose) { 
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
  finfo.curve=finfo.curveerr=finfo.funcmin=finfo.funcminerr=0.0;
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
      double curve=mod.curve(walker.c);
      double funcmin=mod.funcmin(walker.c);
      double oldavg=finfo.curve;
      finfo.curve=finfo.curve+(curve-finfo.curve)/(navgpts+1);
      if(navgpts>0) 
        finfo.curveerr=(1-1.0/navgpts)*finfo.curveerr+
          (navgpts+1)*(finfo.curve-oldavg)*(finfo.curve-oldavg);
      oldavg=finfo.funcmin;
      finfo.funcmin=finfo.funcmin+(funcmin-finfo.funcmin)/(navgpts+1);
      if(navgpts>0) 
        finfo.funcminerr=(1-1.0/navgpts)*finfo.funcminerr+
          (navgpts+1)*(finfo.funcmin-oldavg)*(finfo.funcmin-oldavg);

      for(int m=0; m < nminima; m++) { 
        double oldmin=finfo.min[m];
        double oldvar=finfo.minerr[m];
        finfo.min[m]=oldmin+(min[m]-oldmin)/(navgpts+1);
        if(navgpts > 0) { 
          finfo.minerr[m]=(1-1.0/navgpts)*oldvar+(navgpts+1)*(finfo.min[m]-oldmin)*(finfo.min[m]-oldmin);
        }
      }
      if(step%500==0) { 
        allwalkers.push_back(walker);
      }
      
      
      navgpts++;
    }
  }
  
  for(int p=0; p < nparms; p++) 
    vars[p]=sqrt(vars[p]);
  
  for(int m=0; m < nminima; m++) { 
    finfo.minerr[m]=sqrt(finfo.minerr[m]);
  }
  finfo.curveerr=sqrt(finfo.curveerr);
  finfo.funcminerr=sqrt(finfo.funcminerr);
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



void optimize_quad(Quad_plus_line & quad, vector <Line_data> & data, 
                   vector <Line_model *> & models, vector <Fix_information> & fixes,
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
                    vector <Line_model *> & models, vector <Fix_information> & fixes,
                         vector <double> & c) { 
  quad.generate_guess(data,models, c);
  vector <double> best_c=c;
  double sig=0.6;

  optimize_quad(quad, data, models, fixes, c);
  best_c=c;
  double best_p=quad.prob(data,models, c,fixes);
  cout << "initial probability ";
  for(dit_t i=c.begin(); i!= c.end(); i++) cout << *i << " ";
  cout << best_p  << endl;
  
  //int nit=1000;
  int nit=10;
  int nparms=c.size();
  int nsubit=nparms;
  int ndim=data[0].direction.size();
  cout << "individual shake " << endl;
  for(int i=0; i< nit*nsubit; i++) { 
    int j=int(best_c.size()*rng.ulec());
    c=best_c;
    c[j]=best_c[j]*(1+sig*rng.gasdev());
    optimize_quad(quad, data, models,fixes,c);
    double p=quad.prob(data, models, c,fixes);
    if( (p > best_p && quad.has_minimum(c,ndim)) || isnan(best_p) ) { 
      best_c=c;
      best_p=p;
      if(1) { 
        cout << " it " << i << " ";
        cout << p  << endl;
      }
    }
  }
  
  c=best_c;
  
  
}

//------------------------------------------------------------------------------

double optimize_quad_grad(vector<Line_data> & lines, 
    vector<Line_model *> & models, 
    vector<Gradient_data> & gradients, vector<double> & c) { 
    Quad_opt_with_grad optgrad(c.size(),0,1e-8,30,1);
    optgrad.data=&lines;
    optgrad.gdata=&gradients;
    optgrad.models=&models;
    optgrad.c=c;
    double p=optgrad.optimize();
    c=optgrad.c;
    return p;
}

void shake_quad_grad(vector<Line_data> & lines, 
    vector<Line_model *> & models,
    vector<Gradient_data> & gradients, vector<double> & c, int restart) {

  Quad_plus_line quad;
  Hess_grad gquad;
  vector <double> c_tmp;
  quad.generate_guess(lines,models,c_tmp);
  double p_tmp=optimize_quad_grad(lines,models,gradients,c_tmp);
  double best_p;vector<double> best_c;

  best_p=p_tmp;
  best_c=c_tmp;
  int nstarts=20;
  for(int i=0; i< nstarts; i++) { 
     quad.generate_guess(lines,models,c_tmp);
     p_tmp=optimize_quad_grad(lines,models,gradients,c_tmp);
     if(p_tmp > best_p && !isnan(p_tmp)) {
       best_p=p_tmp;
       best_c=c_tmp;
     }
  }

  if(restart) { 
    double restart_p=optimize_quad_grad(lines,models,gradients,c);
    cout << "restart_p " << restart_p << " from intial guess " << best_p << endl;
    if( (restart_p > best_p) &&  !isnan(restart_p)) {
      best_p=restart_p;
      best_c=c;
    }
  }
  int ndim=lines[0].direction.size(); 
  int nfit=3*lines.size()+ndim*gradients.size();

  cout << "initial prob " << best_p/nfit << endl;

 
  double sig=0.6;
  int nit=10;
  int nparms=c.size();
  int nsubit=nparms;
  cout << "individual shake " << endl;
    for(int i=0; i< nit*nsubit; i++) { 
      int j=int(best_c.size()*rng.ulec());
      c=best_c;
      c[j]=best_c[j]*(1+sig*rng.gasdev());
      double p=optimize_quad_grad(lines,models,gradients,c);
      if( (p > best_p && quad.has_minimum(c,ndim)) && !isnan(p) ) { 
        best_c=c;
        best_p=p;
        if(1) { 
          cout << " it " << i << " " << p/nfit << endl;
        }
      }
    }
  
  
  c=best_c;

}


//------------------------------------------------------------------------------

int metropolis_step(Quad_plus_line & quad, vector <Line_data> & data, 
                     vector <Line_model * > & models, vector <Fix_information> & fixes,vector <double> & tstep, 
                     Walker & walker, int p) { 
  Walker nw=walker;
  vector <Fix_information> savefixes=fixes;
  double ststep=sqrt(tstep[p]);
  double ts=tstep[p];
  double del=1e-12;
  vector <double> c=walker.c;
  double del1=quad.delta_prob(data,models,fixes,c,p,walker.c[p]+del);
  double grad=(del1-walker.prob)/del;
  
  //double grad=gradient(quad, data, models,fixes, walker.c, p,walker.prob);
  nw.c[p]=walker.c[p]+ststep*rng.gasdev()+ts*grad;
  //do it in this order so that the fixes are up to date on acceptance
  double del2=quad.delta_prob(data,models,fixes,c,p,nw.c[p]+del);
  nw.prob=quad.delta_prob(data,models,fixes,c,p,nw.c[p]);
  double ngrad=(del2-nw.prob)/del;
  //double ngrad=gradient(quad, data, models,fixes, nw.c,p,nw.prob);
  //nw.prob=quad.prob(data, models, nw.c,fixes);
  
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
    //nw.c=walker.c;
    fixes=savefixes;
    tstep[p]*=.9;
    return 0;
  }
  
}

//------------------------------------------------------------------------------

void write_paths(vector <Walker> & allwalkers) { 
  int nparms=allwalkers[0].c.size();

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


//------------------------------------------------------------------------------

void check_quad_model(Quad_plus_line & quad,vector <Line_data> &  data, 
                      vector <Line_model *> & models, vector <Fix_information> & fixes,
                      vector <double> & c) { 
  int ndim=data[0].direction.size();
  vector <double> m; vector <vector <double> > H;
  quad.get_minimum(c,ndim,m);
  
  quad.get_hessian(c,ndim,H);
  cout << "minima " << endl;
  for(int i=0; i< ndim; i++) cout << m[i] << " ";
  cout << endl;
  cout << "Hessian " << endl;
  for(int i=0; i < ndim; i++) {
    for(int j=0;j < ndim; j++) cout << H[i][j] << " ";
    cout << endl;
  }
  
  int nparms=c.size();
  double delx=0.01;
  vector <double> real_dels(nparms);
  double base=quad.prob(data,models,c,fixes);
  for(int i=0; i< nparms; i++) { 
    double tmp=c[i];
    c[i]+=delx;
    real_dels[i]=quad.prob(data,models, c,fixes);
    c[i]=tmp;
  }
  quad.prob(data,models,c,fixes);

  for(int i=1; i< nparms; i++) { 
    double firstdel=quad.delta_prob(data, models, fixes,c,i,c[i]+delx);
    double seconddel=quad.delta_prob(data,models, fixes, c,i,c[i]-delx);
    cout << "real delta " << real_dels[i] << " delta " << firstdel << " diff " << real_dels[i]-firstdel
    << endl;
  }
}


//------------------------------------------------------------------------------

void sample(Quad_plus_line & quad, vector <Line_data> & data, vector <Line_model *> & models, 
            Fit_info & finfo, vector<double> & startc, vector <Walker> & allwalkers, int verbose) { 
  int nconfig=100;
  vector <Walker> configs(nconfig);
  vector <Fix_information> fixes;
  
  //check_quad_model(quad, data, models, fixes,startc);
  //exit(1);
  
  Walker walker;
  quad.generate_guess(data, models, walker.c);
  //assert(startc.size()==walker.c.size());
  if(startc.size() < walker.c.size()) { 
    int orig=startc.size();
    startc.resize(walker.c.size());
    for(int i=orig; i < walker.c.size(); i++) 
      startc[i]=walker.c[i];
  }
  
  walker.c=startc;
  shake_quad(quad,data, models, fixes, walker.c);
  walker.prob=quad.prob(data, models, walker.c,fixes);
  //exit(0);
  for(vector <Walker>::iterator i=configs.begin(); i!= configs.end(); i++)
    *i=walker;
  //int nstep=10000;
  int nstep=  16;
  int decorr=100;
  //int warmup=200;
  int warmup=10;
  int nparms=walker.c.size();
  vector <double> tmpc;
  vector <double> tstep(nparms);
  for(int p=0; p < nparms; p++) tstep[p]=2e-2;
  double acc=0;
  vector <double> avgs(nparms);
  vector <double> vars(nparms);
  int navgpts=0;
  for(dit_t i=avgs.begin(); i!= avgs.end(); i++) *i=0.0;
  for(dit_t i=vars.begin(); i!= vars.end(); i++) *i=0.0;
  vector <double> min, curve;
  finfo.min.resize(min.size());
  finfo.minerr.resize(min.size());
  for(dit_t i=finfo.min.begin(); i!= finfo.min.end(); i++) *i=0.0;
  for(dit_t i=finfo.minerr.begin(); i!= finfo.minerr.end(); i++) *i=0.0;
  //vector <Walker> allwalkers;

  for(int step=0; step < nstep; step++) { 
    for(vector<Walker>::iterator w=configs.begin(); w!=configs.end(); w++) {
      //update the fixes for this walker.
      quad.prob(data,models,w->c,fixes);
      for(int d=0; d< decorr; d++) { 
        for(int p=0; p < nparms; p++) { 
          if(metropolis_step(quad, data, models, fixes, tstep, *w, p)) acc++;
        }
      }
      if(step > warmup) { 
        
        for(int p=0; p < nparms; p++) { 
          double oldavg=avgs[p];
          double oldvar=vars[p];
          //cout << p; cout.flush();
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
    cout << "acceptance " << acc/(nstep*decorr*nparms) << endl;
    cout << "timesteps   avg   err" << endl;
    for(int p=0; p < nparms; p++) {
      cout << tstep[p]  << "  " << avgs[p] << "  " << vars[p] << endl;
    }
    cout << endl;
    
  }
  
  finfo.cavg=avgs;
  finfo.cer=vars;
  //write_paths(allwalkers);
  
}


