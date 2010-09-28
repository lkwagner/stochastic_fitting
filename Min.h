#ifndef MIN_H_INCLUDED
#define MIN_H_INCLUDED

#include "macopt.h"
#include "Point.h"
#include "Line_model.h"
#include <vector>
#include <iostream>
#include "Quad_plus_line.h"
#include "Hess_grad.h"
using namespace std;




//------------------------------------------------------------------------

class Least_squares_opt:public Macopt { 
public:
  Least_squares_opt(int _n,int _verbose, double _tol,
                    int _itmax,int _rich):Macopt(_n,_verbose,_tol,_itmax,_rich) {
                    } 
  double func(double * _p) { 
    
    double f=0;
    for(int i=0; i< c.size(); i++) c[i]=_p[i+1];
    f=-mod->prob(*data, *fix, c);
    return f;
     
  }  
  
  double dfunc(double * _p, double * _g) {
    int m=c.size()+1;
    double base=func(_p);
    double step=1e-8;
    for(int i=1; i< m; i++) {
      _p[i]+=step;
      double nwfunc=func(_p);
      _p[i]-=step;
      _g[i]=(nwfunc-base)/step;
    }
    _g[0]=base;
    return base;
  }
  void iteration_print(double f, double gg, double tol,  int itn) {
    //cout <<"macopt: ";
    //for(int i=0; i< c.size(); i++) cout << c[i] << " ";
    //cout <<  endl;
  }
  
  void optimize() { 
    double p[c.size()+1];
    for(int i=0; i< c.size(); i++) p[i+1]=c[i];
    //maccheckgrad(p,c.size(),.0001, c.size());
    macoptII(p,c.size());
  }
  
  Line_model * mod;
  const Line_data * data;
  const Fix_information * fix;
  vector <double>  c; 
  
};


//------------------------------------------------------------------

class Quad_opt:public Macopt { 
public:
  Quad_opt(int _n,int _verbose, double _tol,
                    int _itmax,int _rich):Macopt(_n,_verbose,_tol,_itmax,_rich) {
  } 
  double func(double * _p) { 
    
    double f=0;
    for(int i=0; i< c.size(); i++) c[i]=_p[i+1];
    f=-mod->prob(*data, *models, c, fixes);
    return f;
    
  }
  
  double dfunc(double * _p, double * _g) {
    int m=c.size();
    for(int i=0; i< c.size(); i++) c[i]=_p[i+1];
    
    double base=mod->prob(*data,*models,c,fixes);
    
    double step=1e-8;
    
    for(int i=0; i< m; i++) { 
      vector <Fix_information> tmpfixes=fixes;
      double save=c[i];
      double nwfunc=mod->delta_prob(*data,*models,tmpfixes,c,i,c[i]+step);
      double backfunc=mod->delta_prob(*data,*models,tmpfixes,c,i,c[i]-2.0*step);
      _g[i+1]=-(nwfunc-backfunc)/(2.0*step);
      fixes=tmpfixes;
      c[i]=save;
    }      
    _g[0]=base;
    return base;
  }
  void iteration_print(double f, double gg, double tol,  int itn) {
    //cout <<"macopt: ";
    //for(int i=0; i< c.size(); i++) cout << c[i] << " ";
    //cout <<  endl;
  }
  
  void optimize() { 
    double p[c.size()+1];
    for(int i=0; i< c.size(); i++) p[i+1]=c[i];
    //maccheckgrad(p,c.size(),.0001, c.size());
    macoptII(p,c.size());
  }
  
  Quad_plus_line * mod;
  const vector <Line_data> * data;
  const vector <Line_model * > * models;
  vector <Fix_information> fixes;
  vector <double>  c; 
  
};


//------------------------------------------------------------------


class Quad_opt_approx:public Macopt { 
public:
  Quad_opt_approx(int _n,int _verbose, double _tol,
                    int _itmax,int _rich):Macopt(_n,_verbose,_tol,_itmax,_rich) {
  } 


  double func(double * _p) { 
    
    double f=0;
    for(int i=0; i< c.size(); i++) c[i]=_p[i+1];
    mod.set_fixes(*data, c,fixes);
    int nlines=fixes.size();

    vector<Fit_info>::iterator fit=fits->begin();
    vector<Fix_information>::iterator fix=fixes.begin();
    for(int l=0; l< nlines; l++) { 
      f+=(fix->curve-fit->curve)*(fix->curve-fit->curve)/(2*fit->curveerr*fit->curveerr);
      f+=(fix->min-fit->min[0])*(fix->min-fit->min[0])/(2*fit->minerr[0]*fit->minerr[0]);
      f+=(fix->valmin-fit->funcmin)*(fix->valmin-fit->funcmin)/(2*fit->funcminerr*fit->funcminerr);
      fit++;
      fix++;
    }
    return f;
  }
  
  double dfunc(double * _p, double * _g) {
    int m=c.size();
    double base=func(_p);
    double step=1e-7;
    for(int i=1; i< m+1; i++) {
      double savep=_p[i];
      _p[i]+=step;
      double nwfunc=func(_p);
      _p[i]-=2.0*step;
      double backfunc=func(_p);
      _g[i]=(nwfunc-backfunc)/(2.0*step);
      _p[i]=savep;
                                     
    }
     
    _g[0]=base;
    return base;
  }
  void iteration_print(double f, double gg, double tol,  int itn) {
    //cout <<"macopt: ";
    //for(int i=0; i< c.size(); i++) cout << c[i] << " ";
    //cout <<  endl;
  }
  
  double optimize() { 
    double p[c.size()+1];
    fixes.resize(data->size());
    for(int i=0; i< c.size(); i++) p[i+1]=c[i];
    //maccheckgrad(p,c.size(),.0001, c.size());
    macoptII(p,c.size());
    for(int i=0; i < c.size(); i++) c[i]=p[i+1];
    return func(p);
  }
  
  Quad_plus_line  mod;
  vector <Line_data> * data;
  vector <Fit_info> * fits;
  vector <Fix_information> fixes;
  vector <double>  c; 
  
};



class Quad_opt_with_grad:public Macopt { 
public:
  Quad_opt_with_grad(int _n,int _verbose, double _tol,
                    int _itmax,int _rich):Macopt(_n,_verbose,_tol,_itmax,_rich) {
  } 
  //------------------------------------
  double func(double * _p) { 
    
    double f=0;
    for(int i=0; i< c.size(); i++) c[i]=_p[i+1];
    f=-mod.prob(*data, *models, c, fixes);
    vector<double>cg=c;
    cg.erase(cg.begin());
    f+=-hessg.prob(*gdata,cg);
    return f;
    
  }
  //--------------------------------------------
  
  double dfunc(double * _p, double * _g) {
    int m=c.size();
    for(int i=0; i< c.size(); i++) c[i]=_p[i+1];
    
    double base=mod.prob(*data,*models,c,fixes);
    for(int i=0; i< m; i++) _g[i+1]=0.0; 
    double step=1e-8;
   
    for(int i=0; i< m; i++) { 
      vector <Fix_information> tmpfixes=fixes;
      double save=c[i];
      double nwfunc=mod.delta_prob(*data,*models,tmpfixes,c,i,c[i]+step);
      double backfunc=mod.delta_prob(*data,*models,tmpfixes,c,i,c[i]-2.0*step);
      _g[i+1]=-(nwfunc-backfunc)/(2.0*step);
      fixes=tmpfixes;
      c[i]=save;
    }
   
    vector <double> cg=c;
    cg.erase(cg.begin());
    vector <double> hg_grad;
    hessg.grad_prob(*gdata,cg,hg_grad);
    for(int i=0; i< m-1; i++) { 
      _g[i+2]-=hg_grad[i];
    }
    


    _g[0]=-(base+hessg.prob(*gdata,cg));
    return base;
  }
  //---------------------------------------------------------
  void iteration_print(double f, double gg, double tol,  int itn) {
    //cout <<"opt_val_grad " << f << endl;
    //for(int i=0; i< c.size(); i++) cout << c[i] << " ";
    //cout <<  endl;
  }
  
  double optimize() { 
    double p[c.size()+1];
    for(int i=0; i< c.size(); i++) p[i+1]=c[i];
    //maccheckgrad(p,c.size(),1e-8, c.size());
    macoptII(p,c.size());
    for(int i=0; i< c.size(); i++) c[i]=p[i+1];
    return -func(p);
  }
  
  Quad_plus_line  mod;
  Hess_grad hessg;
  const vector <Line_data> * data;
  const vector <Gradient_data> * gdata;
  const vector <Line_model * > * models;
  vector <Fix_information> fixes;
  vector <double>  c; 
  
};



#endif //MIN_H_INCLUDED

