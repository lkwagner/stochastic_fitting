#include "Quad_plus_line.h"
#include "ulec.h"

//The format of c will be:
// E_0 , vector of minima, then the unique values of Hessian in the following order
// for(i=0; i< n; i++)  for(j=i; j< n; j++) H(i,j)
// Then followed by the extra variables for the line models in order, so it's 
// important to keep the order of the models and data the same between calls.
// Near the minimum, E is approximated as E_0 + (x-m)^T H (x-m)
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
    int nvals=i->val.size();
    for(int j=0; j< nvals; j++) { 
      if(i->val[j] < min_en) {
        min_en=i->val[j];
        for(int d=0; d<ndim; d++) { 
          min_en_pos[d]=i->start_pos[d]+i->t[j]*i->direction[d];
        }
      }
    }
  }
  
  //cout << "guessing minimum at: en " << min_en << " and pos ";
  int count=0;
  c[count++]=min_en;
  for(int i=0; i< ndim; i++) { 
    //c[i+1]= data[ndata-1].start_pos[i]*(1+.5*rng.gasdev());
    c[count++]=min_en_pos[i];
    //cout << min_en_pos[i] << " ";
  }
  //cout << endl;
  Array2 <double> H(ndim, ndim);
  generate_posdef_matrix(ndim,H);
  assert(count==ndim+1);
  for(int i=0; i< ndim; i++) { 
    for(int j=i; j< ndim; j++) { 
      c[count++]=H(i,j);
    }
  }

  set_fixes(data, c, fixes);
  for(int i=0; i< nlines; i++) { 
    vector <double> c_tmp;
    models[i]->generate_guess(data[i],fixes[i],c_tmp);
    c.insert(c.end(),c_tmp.begin(), c_tmp.end());
  }
  
}

//-----------------------------------------------------------------------------

void Quad_plus_line::set_fixes(const vector <Line_data> & data,
                               const vector <double> & c,
                               vector <Fix_information> & fixes)  {
  int nlines=data.size();
  int ndim=data[0].direction.size();
  
  assert(data.size()==fixes.size());
  vector <vector <double> > H; vector <double> m;
  get_minimum(c,ndim,m);

  get_hessian(c,ndim,H);
  for(int line=0; line< nlines; line++) {
    double vHv=0, xmHv=0,xmHxm=0;
    for(int i=0; i< ndim; i++) { 
      for(int j=0; j< ndim; j++) { 
        vHv+=data[line].direction[i]*H[i][j]*data[line].direction[j];
        xmHv+=(data[line].start_pos[i]-m[i])*H[i][j]*data[line].direction[j];
        xmHxm+=(data[line].start_pos[i]-m[i])*H[i][j]*(data[line].start_pos[j]-m[j]);
      }
    }
    fixes[line].enforce=1;
    fixes[line].curve=2.0*vHv;
    double tmin=-xmHv/vHv;
    fixes[line].valmin=c[0]+xmHxm+tmin*tmin*vHv+2.0*tmin*xmHv;
    fixes[line].min=tmin;
  }
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

double Quad_plus_line::prob(const vector <Line_data> & data, 
                            const vector <Line_model * > & models,      
                            const vector <double> & c) {
  vector <Fix_information> fixes(data.size());
  assert(data.size()==models.size());
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
      A(i,i)=5.0*rng.ulec()+0.5;
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
