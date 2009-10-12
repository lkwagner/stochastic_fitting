
#include "Prob_model.h"
#include <iostream>
#include "ulec.h"
#include <cassert>
#include <string>

typedef vector<double>::iterator dit_t;


void split(string & text, string & separators, vector<string> & words) {
  int n = text.length();
  int start, stop;
  
  start = text.find_first_not_of(separators);
  while ((start >= 0) && (start < n)) {
    stop = text.find_first_of(separators, start);
    if ((stop < 0) || (stop > n)) stop = n;
    words.push_back(text.substr(start, stop - start));
    start = text.find_first_not_of(separators, stop+1);
  }
}

//###########################################################################

//---------------------------------------------------------------------

void read_into_data(istream & is, vector <Point> & data) { 
  Point pt;
  string line; string space=" ";
  while(getline(cin, line)) { 
    vector <string> words;
    split(line, space, words);
    if(words.size() > 0) { 
      int ndim=words.size()-2;
      pt.x.resize(ndim);
      for(int d=0; d< ndim; d++) {
        pt.x[d]=atof(words[d].c_str());
      }
      pt.en=atof(words[ndim].c_str());
      pt.err=atof(words[ndim+1].c_str());
      //if(pt.x[0] > 0.81 && pt.x[0] < 1.11
      //   && pt.x[1] > 1.01 && pt.x[1] < 2.39) 
      data.push_back(pt);
    }
  }
}

//---------------------------------------------------------------------

void Morse_model::read(istream & is) { 
  read_into_data(is, data);
}
//------------------------------------------------------------------------

void Morse_model::generate_guess(vector <double> & c) { 
  int ndim=data[0].x.size();
  int nparms=1+2*ndim+ndim*(ndim+1)/2;
  c.resize(nparms);
  //cout << " ndim " << ndim << " nparms " << nparms << endl;
  
  vector <double> maxd(ndim), mind(ndim);
  for(int d=0; d< ndim; d++) { 
    maxd[d]=-1e99; mind[d]=1e99;
    for(vector<Point>::iterator i=data.begin(); 
        i!= data.end(); i++) { 
      if(maxd[d] < i->x[d]) maxd[d]=i->x[d];
      if(mind[d] > i->x[d]) mind[d]=i->x[d];
    }
  }

  for(dit_t i=c.begin(); i!=c.end(); i++) *i=4*rng.ulec();
  c[0]=data[0].en;
  for(int d=0; d< ndim; d++) { 
    c[2*d+1]=mind[d]+(maxd[d]-mind[d])*(0.5+rng.gasdev()*.1);
  }
  
}

//------------------------------------------------------------------------
//does the matrix multiplication for a general Morse function
inline double fval_quad(const vector <double> & c,const vector <double> & x) { 
  double f=c[0];
  int n=x.size();
  assert(c.size() >= 1+2*n+n*(n+1)/2);
  int count=1;
  vector <double> xshift;
  for(int i=0; i< n; i++) { 
    xshift.push_back(1-exp(-c[count+1]*(x[i]-c[count])));
    //xshift.push_back(x[i]-c[count]);
    count++; count++;
  }
  
  for(int i=0; i< n; i++) { 
    for(int j=i; j< n; j++) { 
      double fac=1.0;
      if(i!=j) fac=2.0;
      f+=fac*xshift[i]*c[count]*xshift[j];
      count++;
    }
  }
  return f;
}

//---------------------------------------------------------------------
void Morse_model::minimum(const vector <double> & c, vector <double>  & min) { 
  int n=data[0].x.size();
  min.resize(n);
  assert(c.size() >= 1+2*n+n*(n+1)/2);
  int count=1;
  for(int i=0; i< n; i++) {
    min[i]= c[count];
    count+=2;
  }
}

//---------------------------------------------------------------------


void Morse_model::curvature(const vector <double> & c, vector <double>  & curve) { 
  int n=data[0].x.size();
  curve.resize(n);
  assert(n==1);
  assert(c.size() >= 1+2*n+n*(n+1)/2);
  curve[0]=2*c[2]*c[2]*c[3];
}

//---------------------------------------------------------------------
double Morse_model::probability(const vector <double> & c) { 

  double f=0;
  
  int ndatapts=data.size();
  for(vector<Point>::const_iterator i=data.begin(); i != data.end();
      i++) { 
    double fp=fval_quad(c,i->x);
    f-=(i->en-fp)*(i->en-fp)/(2*i->err*i->err);
  }
  return f;
}

//---------------------------------------------------------------------

void Morse_model::gradient(const vector <double> & c,
                               vector <double> & grad) { 
  double delta=1e-12;
  int n=c.size();
  grad.resize(n);
  vector <double> cp;
  double base=probability(c);
  for(int i=0; i< n; i++) { 
    cp=c;
    cp[i]+=delta;
    double px=probability(cp);
    grad[i]=(px-base)/delta;
  }
}
//---------------------------------------------------------------------

double Morse_model::gradient(const vector <double> & c,
                               int dir) { 
  double delta=1e-9;
  int n=c.size();
  vector <double> cp;
  double base=probability(c);
  cp=c;
  cp[dir]+=delta;
  double px=probability(cp);
  return (px-base)/delta;
}

//---------------------------------------------------------------------
bool Morse_model::is_ok(const vector <double> & c) {
  int count=1;
  int n=data[0].x.size();
  
  for(int i=0; i< n; i++) { 
    double min=1e99, max=-1e99;
    for(vector<Point>::const_iterator j=data.begin(); j!= data.end(); j++) { 
      if(j->x[i] < min) min=j->x[i];
      if(j->x[i] > max) max=j->x[i];
    }
    //if(c[count] < 0) { return false; }
    if (c[count] < min || c[count] > max) {  
      //cout << "rejecting an maximal point because it's out of range of the data set.  You may want"
      //" to check if you're really bracketing the minimum." 
      //<< " direction " << i << " value " << c[count] << endl;
      return false; 
    }
    //if(c[count+1] < 0) {  return false; }
    count+=2;
  }
  
  for(int i=0; i< n; i++) { 
    for(int j=i; j< n; j++) { 
      if(i==j && c[count] < 0) {  
        cout << "rejection because of negative correlation matrix \n";
        return false;}
      count++;
    }
  }
  return true;
  
}

//---------------------------------------------------------------------
void Morse_model::niceprint(const vector <double> & c, const vector <double> & cerr) {
  int n=data[0].x.size();
  
  cout << "min energy " << c[0] << " +/- " << cerr[0] << endl;
  assert(c.size() >= 1+2*n+n*(n+1)/2);
  int count=1;
  vector <double> xshift;
  for(int i=0; i< n; i++) { 
    cout << "gmin " << c[count] << " +/- " << cerr[count] << " b " << c[count+1] 
    << " +/- " << cerr[count+1] << endl;
    count++; count++;
  }
  
  cout << "correlation matrix " << endl;
  for(int i=0; i< n; i++) { 
    int counta=count;
    for(int j=i; j< n; j++) { 
      cout << c[count] << " ";
      count++;
    }
    
    cout << " +/- ";
    for(int j=i; j< n; j++) { 
      cout << cerr[counta] << " ";
      counta++;
    }
    cout << endl;
  }
  
  if(n==1) { 
    cout << "gplot function: f(x)="
    << c[0]<< " + " << c[3] << "*(1-exp(-"<< c[2] 
    << "*(x-" << c[1] << ")))**2\n";
  }
  
}


//###########################################################################
//---------------------------------------------------------------------

void Linear_model::read(istream & is) { 
  read_into_data(is, data);
}
//------------------------------------------------------------------------

void Linear_model::generate_guess(vector <double> & c) { 
  int ndim=data[0].x.size();
  assert(ndim==1);
  int nparms=2;
  c.resize(nparms);
  //cout << " ndim " << ndim << " nparms " << nparms << endl;
  //this is probably not terribly optimal, but it shouldn't be too hard to brute-force our way through it.
  for(int d=0; d< 1; d++) { 
    c[d]=30*(0.5+rng.gasdev()*.1);
  }
  c[1]=data[0].x[0]+0.5*rng.gasdev();
  
}

//---------------------------------------------------------------------

double Linear_model::probability(const vector <double> & c) { 
  
  double f=0;
  
  int ndatapts=data.size();
  for(vector<Point>::const_iterator i=data.begin(); i != data.end();
      i++) { 
    double x=i->x[0];
    double fp=c[0]*x+c[1];
    f-=(i->en-fp)*(i->en-fp)/(2*i->err*i->err);
  }
  return f;
}

//---------------------------------------------------------------------

void Linear_model::gradient(const vector <double> & c,
                           vector <double> & grad) { 
  double delta=1e-12;
  int n=c.size();
  grad.resize(n);
  vector <double> cp;
  double base=probability(c);
  for(int i=0; i< n; i++) { 
    cp=c;
    cp[i]+=delta;
    double px=probability(cp);
    grad[i]=(px-base)/delta;
  }
}
//---------------------------------------------------------------------

double Linear_model::gradient(const vector <double> & c,
                             int dir) { 
  double delta=1e-12;
  int n=c.size();
  vector <double> cp;
  double base=probability(c);
  cp=c;
  cp[dir]+=delta;
  double px=probability(cp);
  return (px-base)/delta;
}

//---------------------------------------------------------------------
bool Linear_model::is_ok(const vector <double> & c) {
  return true;  
}

//---------------------------------------------------------------------
void Linear_model::niceprint(const vector <double> & c, const vector <double> & cerr) {
  int n=data[0].x.size();
  cout << "Function: c[0]*x+c[1] \n";
  for(int i=0; i< 2; i++) { 
    cout << "c[" <<  i << "] = " << c[i] << " +/- " << cerr[i] << endl;
  }
    
}


//---------------------------------------------------------------------
//###########################################################################
void Prob_distribution::avg(double & avg, double & err) { 
  double sum=0;
  for(vector<double>::iterator i=density.begin();
      i!=density.end(); i++) sum+=*i;
  
  avg=0;
  double a=mymin;
  for(vector <double>::iterator i=density.begin();
      i!= density.end(); i++) {
    a+=spacing;
    avg+=a*(*i)/sum;
  }
  
  a=mymin;
  err=0;
  for(vector <double>::iterator i=density.begin();
      i!= density.end(); i++) {
    a+=spacing;
    err+=(a-avg)*(a-avg)*(*i)/sum;
  }
  err=sqrt(err);
}

//---------------------------------------------------------------------

double Prob_distribution::skew() { 
  double mu,sigma;
  avg(mu,sigma);
  
  double sum=0.0;
  for(vector<double>::iterator i=density.begin();
      i!=density.end(); i++) sum+=*i;
  
  double mu3=0.0;
  
  double a=mymin;
  for(vector <double>::iterator i=density.begin(); 
      i!= density.end(); i++) { 
    a+=spacing;
    mu3+=(a-mu)*(a-mu)*(a-mu)*(*i)/sum;
  }
  return mu3/(sigma*sigma*sigma);
}




//###########################################################################


struct Walker { 
  vector <double> c;
  double prob;
  vector <double> grad;
  Walker(int size) { 
    c.resize(size);
    grad.resize(size);
  }
  void update(Prob_model * mod) { 
    prob=mod->probability(c);
    mod->gradient(c,grad);
  }
  
  void updatedir(Prob_model * mod, int d) { 
    prob=mod->probability(c);
    vector <double> tmpc=c;
    double del=1e-12;
    tmpc[d]+=del;
    double px=mod->probability(tmpc);
    grad[d]=(px-prob)/del;
  }
};

//---------------------------------------------------------------------


int check_model_consistency(Prob_model & mod,vector <double> & c) { 
  double del=1e-12;
  int nparms=c.size();
  Walker walk(nparms);
  walk.c=c;
  walk.update(&mod);
  Walker nw=walk;
  for(int p=0; p < nparms; p++) { 
    nw=walk;
    nw.c[p]+=del;
    nw.update(&mod);
    double ratio=nw.prob-walk.prob;
    double deriv=(ratio)/del;
    cout << "numerical derivative " << deriv 
      << " analytic " << walk.grad[p] <<    endl;
    
  }
  
}


//---------------------------------------------------------------------


void sample_min(vector <double> & c, Prob_model & mod, 
                Fit_info & finfo, int verbose ) { 
  
  Walker walker(c.size());
  walker.c=c;

  walker.update(&mod);
  Walker best_walker=walker;
  
  int nstep=100000;
  int warmup=nstep/3;
  Walker nw=walker;
  int nparms=c.size();
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
  mod.minimum(c,min);
  mod.curvature(c,curve);
  int nminima=min.size();
  finfo.min.resize(min.size());
  finfo.minerr.resize(min.size());
  finfo.curv.resize(curve.size());
  finfo.curverr.resize(curve.size());
  for(dit_t i=finfo.min.begin(); i!= finfo.min.end(); i++) *i=0.0;
  for(dit_t i=finfo.minerr.begin(); i!= finfo.minerr.end(); i++) *i=0.0;
  for(dit_t i=finfo.curv.begin(); i!= finfo.curv.end(); i++) *i=0.0;
  for(dit_t i=finfo.curverr.begin(); i!= finfo.curverr.end(); i++) *i=0.0;
  
  
  double avg_min=0;
  for(int step=0; step < nstep; step++) { 
    nw=walker;
    for(int p=0; p < nparms; p++) { 
      double ststep=sqrt(tstep[p]);
      double ts=tstep[p];
      double delta=1e-9;
      walker.updatedir(&mod,  p);
      nw.c[p]=walker.c[p]+ststep*rng.gasdev()+ts*walker.grad[p];
      nw.updatedir(&mod, p);
      
      
      double diff=nw.c[p]-walker.c[p];
      double num=-(-diff-ts*nw.grad[p])*(-diff-ts*nw.grad[p]);
      double den=-(diff-ts*walker.grad[p])*(diff-ts*walker.grad[p]);
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
    
    if(walker.prob > best_walker.prob) best_walker=walker;
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
      mod.curvature(walker.c,curve);
      for(int m=0; m < nminima; m++) { 
        double oldmin=finfo.min[m];
        double oldvar=finfo.minerr[m];
        finfo.min[m]=oldmin+(min[m]-oldmin)/(navgpts+1);
        if(navgpts > 0) { 
          finfo.minerr[m]=(1-1.0/navgpts)*oldvar+(navgpts+1)*(finfo.min[m]-oldmin)*(finfo.min[m]-oldmin);
        }
        
        oldmin=finfo.curv[m];
        oldvar=finfo.curverr[m];
        finfo.curv[m]=oldmin+(curve[m]-oldmin)/(navgpts+1);
        if(navgpts > 0) { 
          finfo.curverr[m]=(1-1.0/navgpts)*oldvar+(navgpts+1)*(finfo.curv[m]-oldmin)*(finfo.curv[m]-oldmin);
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
    
  /*
    string nm="minimum";
    min1.write(nm);
    cout << "minimum skew " << min1.skew() << endl;
    string nm2="minimum2";
    min2.write(nm2);
   */
  }
  
  c=best_walker.c;
  finfo.cavg=avgs;
  finfo.cer=vars;
}

//---------------------------------------------------------------------

