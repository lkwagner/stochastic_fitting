#include "Quad_plus_line.h"
#include "ulec.h"
#include "Sample.h"
#include "macopt.h"
#include "MatrixAlgebra.h"
//The format of c will be:
// E_0 , vector of minima, then the unique values of Hessian in the following order
// for(i=0; i< n; i++)  for(j=i; j< n; j++) H(i,j)
// Then followed by the extra variables for the line models in order, so it's 
// important to keep the order of the models and data the same between calls.
// Near the minimum, E is approximated as E_0 + (x-m)^T H (x-m)



class Minimize_curve_err:public Macopt { 
public:
  Minimize_curve_err(int _n,int _verbose, double _tol,
                    int _itmax,int _rich):Macopt(_n,_verbose,_tol,_itmax,_rich) {
  } 
  double func(double * _p) { 
    
    double f=0;
    Array2 <double> H(ndim, ndim);
    int count=1;
    for(int i=0;i< ndim; i++) {
      for(int j=i; j< ndim; j++)  { 
        H(i,j)=_p[count++];
        H(j,i)=H(i,j);
      }
    }
    int line=0;
    for(vector<Line_data>::iterator d=data.begin(); d!=data.end(); d++) { 
      double vHv=0;
      for(int i=0; i< ndim; i++) { 
        for(int j=0; j< ndim; j++) { 
          vHv+=d->direction[i]*H(i,j)*d->direction[j];
        }
      }
      f+=fabs(vHv-0.5*curve[line++]);
    }
    return f;
  }  
  
  double dfunc(double * _p, double * _g) {
    int m=a_n;
    double base=func(_p);
    double step=1e-8;
    for(int i=1; i<= m; i++) {
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
  
  void optimize(vector<double> & c) { 
    int n=(ndim+1)*ndim/2;
    double p[n+1];
    for(int i=1; i<=n; i++) p[i]=10*(rng.ulec()-0.5);
    //maccheckgrad(p,c.size(),.0001, c.size());
    macoptII(p,n);
    c.resize(n);
    for(int i=0; i < n; i++) {
      c[i]=p[i+1];
    }
  }

  vector <Line_data> data;
  vector <double>  curve; 
  int ndim;
  
};

//----------------------------------------------------------------------------

void Quad_plus_line::generate_guess(const vector <Line_data> & data, 
                                    const vector <Line_model * > & models,
                                    vector <double> & c) {
  int ndim=data[0].direction.size();
  int nlines=data.size();
  for(vector<Line_data>::const_iterator i=data.begin(); i!=data.end(); i++) {
    assert(i->direction.size()==ndim);
    assert(i->start_pos.size()==ndim);
  }
  assert(models.size()==nlines);
  
  vector <Fix_information> fixes(nlines);
  
  int nparms_quad=1+ndim+(ndim+1)*ndim/2;
  c.resize(nparms_quad);
  int ndata=data.size();
  double min_en=1e99;
  vector <double> min_en_pos(ndim);
  for(vector <Line_data>::const_iterator i=data.begin(); i!= data.end(); i++) { 
    int nvals=i->data.size();
    for(int j=0; j< nvals; j++) { 
      if(i->data[j].val < min_en) {
        min_en=i->data[j].val;
        for(int d=0; d<ndim; d++) { 
          min_en_pos[d]=i->start_pos[d]+i->data[j].t*i->direction[d];
        }
      }
    }
  }
  
  //cout << "guessing minimum at: en " << min_en << " and pos ";
  int count=0;
  c[count++]=min_en;
  for(int i=0; i< ndim; i++) { 
    c[count++]=min_en_pos[i];
  }
  
  
  vector <Fit_info> finfo(nlines);
  vector <vector <double> > guesses;
  vector <double> curves;
  for(int i=0; i< nlines; i++) { 
    Fix_information falsefix;
    vector <double> c_tmp;
    models[i]->generate_guess(data[i],falsefix,c_tmp);
    find_good_guess(*models[i],data[i],falsefix,c_tmp);
    guesses.push_back(c_tmp);
    curves.push_back(models[i]->curve(c_tmp));
  }


  Array2 <double> H(ndim,ndim),Q(ndim,ndim);
  for(int i=0; i < ndim; i++) { 
     for(int j=0; j< ndim; j++) { 
       Q(i,j)=data[nlines-ndim+i].direction[j];
     }
  }
  Array2 <double> Qinv(ndim, ndim);
  InvertMatrix(Q,Qinv,ndim);
  H=0;
  for(int i=0; i< ndim; i++) 
    for(int j=0; j< ndim; j++) 
      for(int k=0; k < ndim; k++) 
        H(i,k)+=Q(i,j)*curves[j+nlines-ndim]*Qinv(j,k);
  
  assert(count==ndim+1);
  for(int i=0; i< ndim; i++) { 
    for(int j=i; j< ndim; j++) { 
      c[count++]=H(i,j);
    }
  }
   

  set_fixes(data, c, fixes);
  for(int i=0; i< nlines; i++) { 
    vector <double> c_tmp;
    models[i]->downconvert_c(fixes[i],guesses[i],c_tmp);
    c.insert(c.end(),c_tmp.begin(), c_tmp.end());
  }
}
//----------------------------------------------------------------------------
void Quad_plus_line::generate_guess(const vector <Line_data> & data, const vector <Line_model *> & models, 
    const double e0, const vector <double> & guess_min, const vector < vector <double> > & guess_hess,
    vector <double> & c) { 
  int ndim=data[0].direction.size();
  int nlines=data.size();
  for(vector<Line_data>::const_iterator i=data.begin(); i!=data.end(); i++) {
    assert(i->direction.size()==ndim);
    assert(i->start_pos.size()==ndim);
  }
  assert(models.size()==nlines);
  assert(guess_hess.size()>=ndim);
  assert(guess_hess[0].size() >=ndim);
  assert(guess_min.size() >=ndim);
  vector <Fix_information> fixes(nlines);
  int nparms_quad=1+ndim+(ndim+1)*ndim/2;
  c.resize(nparms_quad);
  c[0]=e0;
  for(int i=0; i< ndim; i++) c[i+1]=guess_min[i];
  int count=ndim+1;
  for(int i=0; i< ndim; i++) { 
    for(int j=i; j< ndim; j++) c[count++]=guess_hess[i][j];
  }
  set_fixes(data, c, fixes);
  for(int i=0; i < nlines; i++) { 
    vector <double> c_tmp;
    models[i]->generate_guess(data[i],fixes[i],c_tmp);
    c.insert(c.end(), c_tmp.begin(), c_tmp.end());
  }
}

//-----------------------------------------------------------------------------

void Quad_plus_line::set_fixes(const vector <Line_data> & data,
                               const vector <double> & c,
                               vector <Fix_information> & fixes)  {
  int nlines=data.size();
  int ndim=data[0].direction.size();
  
  assert(data.size()==fixes.size());
  //static is thread-unsafe, but saves large amounts of time
  static vector <vector <double> > H; 
  static vector <double> m;
  get_minimum(c,ndim,m);

  get_hessian(c,ndim,H);
  for(int line=0; line< nlines; line++) {
    double vHv=0, xmHv=0,xmHxm=0;
    

    double tmp1,tmp2;
    for(int i=0; i< ndim; i++) { 
      tmp1=tmp2=0.0;
      for(int j=0; j< ndim; j++) { 
        tmp1+=H[i][j]*data[line].direction[j];
        tmp2+=H[i][j]*(data[line].start_pos[j]-m[j]);
      }
      vHv+=data[line].direction[i]*tmp1;
      xmHv+=(data[line].start_pos[i]-m[i])*tmp1;
      xmHxm+=(data[line].start_pos[i]-m[i])*tmp2;
    }
     
    fixes[line].enforce=1;
    fixes[line].curve=2.0*vHv;
    double tmin=-xmHv/vHv;
    fixes[line].valmin=c[0]+xmHxm+tmin*tmin*vHv+2.0*tmin*xmHv;
    fixes[line].min=tmin;
    //cout << "tmin " << tmin << " xmHxm " << xmHxm << " xmHv " << xmHv << endl;
  }
}


//-----------------------------------------------------------------------------

void Quad_plus_line::update_fixes(const vector <Line_data> & data,
                                vector <double> & c,
                               vector <Fix_information> & fixes, int p, 
                                  double cp)  {
  int nlines=data.size();
  int ndim=data[0].direction.size();
  static vector <vector <double> > H; 
  static vector <double> m;
  
  if(p==0) { 
    for(int line=0; line < nlines; line++) { 
      fixes[line].valmin+=(cp-c[0]);
    }
  }
  else if(p < ndim+1) { 
    //cout <<"-----------------------------------------\n";
    get_minimum(c,ndim,m);
    
    get_hessian(c,ndim,H);

/*
   vector <Fix_information> testfix(nlines);

    cout << "minima " << endl;
    for(int i=0; i< ndim; i++) cout << m[i] << " ";
    cout << endl;
    cout << "Hessian " << endl;
    for(int i=0; i < ndim; i++) {
      for(int j=0;j < ndim; j++) cout << H[i][j] << " ";
      cout << endl;
    }
    
    cout << "minimum " << p-1 << " changed to " << cp << endl;
    
    set_fixes(data, c, testfix);
    for(int line=0; line < nlines; line++) { 
      //cout << "premin " << fixes[line].min << " " << testfix[line].min << endl;
      //cout << "precurve " << fixes[line].curve << " " << testfix[line].curve << endl;
      cout << "prevalmin " << fixes[line].valmin << " " << testfix[line].valmin << endl;
    }
    */
    for(int line=0; line < nlines; line++) { 
      double vHv=fixes[line].curve*0.5;
      int i=p-1;
      double tmp1=0;
      double tmp2=0;
      double deltaxm=c[p]-cp;
      for(int j=0; j< ndim; j++) { 
        tmp1+=H[i][j]*data[line].direction[j];
        if(i!=j)
          tmp2+=H[i][j]*(data[line].start_pos[j]-m[j]);
      }
      double xm=data[line].start_pos[i]-c[p];
      double xmp=data[line].start_pos[i]-cp;
      double delxmHv=deltaxm*tmp1;
      double oldmin=fixes[line].min;
      double oldxmHv= -vHv*oldmin;
      double xmHv=oldxmHv+delxmHv;
      double delxmHxm=2.0*deltaxm*tmp2+H[i][i]*(xmp*xmp-xm*xm);
      fixes[line].enforce=1;
      fixes[line].min+= -delxmHv/vHv;
      double tmin=fixes[line].min;
      //cout << "deltaxmHxm " << delxmHxm << " xmHv " << xmHv << " oldxmHv " << oldxmHv << endl;
      //cout << "valmin del " << delxmHxm+2.0*tmin*xmHv-2.0*oldmin*oldxmHv 
      //     <<  " old " << fixes[line].valmin << endl;
      //cout << "delxmHxm " << delxmHxm << "  delsec " << 2.0*tmin*xmHv-2.0*oldmin*oldxmHv << endl;
      fixes[line].valmin+=delxmHxm+2.0*tmin*xmHv-2.0*oldmin*oldxmHv+vHv*(tmin*tmin-oldmin*oldmin);
    }
    
    
    /*
    c[p]=cp;
    
    set_fixes(data, c, testfix);
    for(int line=0; line < nlines; line++) { 
      //cout << "min " << fixes[line].min << " " << testfix[line].min << endl;
      //cout << "curve " << fixes[line].curve << " " << testfix[line].curve << endl;
      cout << "valmin " << fixes[line].valmin << " " << testfix[line].valmin << endl;
    }
     */
    //exit(0);
  }
  else if(p < ndim+1+ndim*(ndim+1)/2) { 
    //cout << "----------------------------------------\n";
    int index=p-ndim-1;
    int i,j, count=0;
    for(i=0; i< ndim; i++) {
      if(count+ndim-i > index) { 
        j=index-count+i;
        break;
      }
      else count+= ndim-i;
    }
    
    //cout << "i= " << i << " j= " << j << endl;
    get_minimum(c,ndim,m);    
    //get_hessian(c,ndim,H);
    double fac=1.0;
    if(i!=j) fac=2.0;
    for(int line=0; line < nlines; line++) { 
      double oldvHv=fixes[line].curve*0.5;
      double oldxmHv=-fixes[line].min*oldvHv;
      double oldmin=fixes[line].min;
      double delH=cp-c[p];
      double delvHv=fac*data[line].direction[i]*delH*data[line].direction[j];
      double delxmHv;
      if(i==j) 
        delxmHv=(data[line].start_pos[i]-m[i])*delH*data[line].direction[j];
      else 
        delxmHv=(data[line].start_pos[i]-m[i])*delH*data[line].direction[j]
                        +data[line].direction[i]*delH*(data[line].start_pos[j]-m[j]);
                    

      double delxmHxm=fac*(data[line].start_pos[i]-m[i])*delH*(data[line].start_pos[j]-m[j]);
      double vHv=oldvHv+delvHv;
      double xmHv=oldxmHv+delxmHv;
      //cout << "xmHv " << xmHv << " delxmHv " << delxmHv << endl;
      fixes[line].curve+=2.0*delvHv;
      fixes[line].min=-xmHv/vHv;
      double tmin=fixes[line].min;
      fixes[line].valmin+=delxmHxm+tmin*tmin*vHv+2.0*tmin*xmHv-oldmin*oldmin*oldvHv-2.0*oldmin*oldxmHv;
    }
    /*
    c[p]=cp;
    vector <Fix_information> testfix(nlines);
    set_fixes(data, c, testfix);
    for(int line=0; line < nlines; line++) { 
      cout << "min " << fixes[line].min << " " << testfix[line].min << " diff "
           << fixes[line].min-testfix[line].min << endl;
      cout << "curve " << fixes[line].curve << " " << testfix[line].curve << " diff "
      << fixes[line].curve-testfix[line].curve << endl;
      cout << "valmin " << fixes[line].valmin << " " << testfix[line].valmin << " diff "
      << fixes[line].valmin-testfix[line].valmin << endl;
    }
     */
    //fixes=testfix;
    
  }
  
  c[p]=cp;
  //for other parameters, do nothing..
  
  /*
  assert(data.size()==fixes.size());
  //static is thread-unsafe, but saves large amounts of time
  for(int line=0; line< nlines; line++) {
    double vHv=0, xmHv=0,xmHxm=0;
    
    
    double tmp1,tmp2;
    for(int i=0; i< ndim; i++) { 
      tmp1=tmp2=0.0;
      for(int j=0; j< ndim; j++) { 
        tmp1+=H[i][j]*data[line].direction[j];
        tmp2+=H[i][j]*(data[line].start_pos[j]-m[j]);
      }
      vHv+=data[line].direction[i]*tmp1;
      xmHv+=(data[line].start_pos[i]-m[i])*tmp1;
      xmHxm+=(data[line].start_pos[i]-m[i])*tmp2;
    }
    
    fixes[line].enforce=1;
    fixes[line].curve=2.0*vHv;
    double tmin=-xmHv/vHv;
    fixes[line].valmin=c[0]+xmHxm+tmin*tmin*vHv+2.0*tmin*xmHv;
    fixes[line].min=tmin;
  }
   */
}


//-----------------------------------------------------------------------------


void Quad_plus_line::get_minimum(const vector <double> & c, int ndim,
                                 vector <double> & m) {
  m.resize(ndim);
  int count=1;
  for(int i=0; i< ndim; i++) m[i]=c[count++];
}

//-----------------------------------------------------------------------------


void Quad_plus_line::get_hessian(const vector <double> & c, int ndim,
                                 vector <vector <double> > & H) { 
  assert(c.size() >= 1+ndim+ndim*(ndim+1)/2);
  H.resize(ndim);
  int count=ndim+1;
  for(int i=0; i< ndim; i++) H[i].resize(ndim);

  for(int i=0; i<ndim; i++) { 
    for(int j=i; j< ndim; j++) { 
      H[i][j]=c[count++];
      H[j][i]=H[i][j];
    }
  }
  
}


//-----------------------------------------------------------------------------

bool Quad_plus_line::has_minimum(const vector <double> & c, int ndim) { 
  assert(c.size() >= 1+ndim+ndim*(ndim+1)/2);
  Array2 <double> H(ndim,ndim);
  int count=ndim+1;  
  for(int i=0; i<ndim; i++) { 
    for(int j=i; j< ndim; j++) { 
      H(i,j)=c[count++];
      H(j,i)=H(i,j);
    }
  }
  return is_posdef(ndim, H);
}


//-----------------------------------------------------------------------------

double Quad_plus_line::prob(const vector <Line_data> & data, 
                            const vector <Line_model * > & models,      
                            const vector <double> & c, vector <Fix_information> & fixes) {
  //vector <Fix_information> fixes(data.size());
  assert(data.size()==models.size());
  fixes.resize(data.size());
  int nlines=data.size();
  set_fixes(data,c,fixes);
  double f=0;
  int ndim=data[0].direction.size();
  vector<double>::const_iterator i=c.begin()+1+ndim+(ndim+1)*ndim/2;
  for(int line=0; line < nlines; line++) {
    int n=models[line]->nparms(1);
    vector <double> ctmp(i,i+n);
    i+=n;
    f+=models[line]->prob(data[line],fixes[line],ctmp );
  }
  //cout << "prob " << f << endl;
  //exit(0);
  //if(isnan(f)) { 
  //  cout << "probnan! ";
   // for(vector<double>::const_iterator i=c.begin(); i!= c.end(); i++) 
   //   cout << *i << " ";
  // cout << endl;
  //}
  return f;
}

//p is the parameter number to change and cp is the value to change it to. 
//returns the new probability, not the change..
double Quad_plus_line::delta_prob(const vector <Line_data> & data, 
                                   const vector <Line_model * > & models,   
                                   vector <Fix_information> & fixes,
                                   vector <double> & c, int p, double cp) {
  assert(fixes.size()==data.size());
  assert(data.size()==models.size());
  int nlines=data.size();
  double f=0;
  int ndim=data[0].direction.size();
  /*
  for(int line=0; line < nlines; line++) {
    int n=models[line]->nparms(1);
    vector <double> ctmp(i,i+n);
    i+=n;
    f-=models[line]->prob(data[line],fixes[line],ctmp );
  }*/
  
  update_fixes(data,c,fixes, p, cp);
  vector<double>::const_iterator i=c.begin()+1+ndim+(ndim+1)*ndim/2;

  //i=c.begin()+1+ndim+(ndim+1)*ndim/2;
  for(int line=0; line < nlines; line++) {
    int n=models[line]->nparms(1);
    vector <double> ctmp(i,i+n);
    i+=n;
    f+=models[line]->prob(data[line],fixes[line],ctmp );
  }
  
  //if(isnan(f)) { 
  //  cout << "deltanan! ";
  //  for(vector<double>::const_iterator i=c.begin(); i!= c.end(); i++) 
  //    cout << *i << " ";
  //  cout << endl;
  //}
  
  return f;
}


//-----------------------------------------------------------------------------

#include "MatrixAlgebra.h"

void generate_posdef_matrix(int n, Array2 <double> & C) { 
  C.Resize(n,n);
  Array2 <double> A(n,n), B(n,n), evecs(n,n);
  Array1 <double> evals(n);

  bool posdef=false;
  do { 
    for(int i=0; i< n; i++) { 
      for(int j=i+1; j < n; j++) { 
        A(i,j)=2.0*(rng.ulec()-0.5);
        A(j,i)=0.0;
      }
      A(i,i)=1.0*rng.ulec()+0.5;
    }
    B=A;
    TransposeMatrix(B,n);
    MultiplyMatrices(A,B,C,n);
    EigenSystemSolverRealSymmetricMatrix(C,evals,evecs);
    posdef=true;
    for(int i=0; i< n; i++) { 
      //cout << evals(i) << " ";
      if(evals(i) <= 0) {
        posdef=false;
        break;
      }
    }
    //cout << endl;
  }
  while(!posdef) ;
}
//-----------------------------------------------------------------------------


bool is_posdef(int n, Array2 <double> & C) { 
  bool posdef=false;
  Array2 <double> evecs(n,n); 
  Array1 <double> evals(n);
  EigenSystemSolverRealSymmetricMatrix(C,evals,evecs);
  posdef=true;
  for(int i=0; i< n; i++) { 
    //cout << evals(i) << " ";
    if(evals(i) <= 0) {
      posdef=false;
      break;
    }
  }
  return posdef;
  
}
