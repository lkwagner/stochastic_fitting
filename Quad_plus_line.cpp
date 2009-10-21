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
  c[0]=data[0].val[0];
  for(int i=0; i< ndim; i++) { 
    c[i+1]=data[0].start_pos[i]*(1+.05*rng.gasdev());
  }
  for(int i=ndim; i < nparms_quad; i++) c[i]=30*(rng.ulec()-0.5);

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
  int ndim=data[0].direction[0];
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

