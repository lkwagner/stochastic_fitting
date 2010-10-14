
#include "Hess_grad.h"
#include <cassert>
#include "Array.h"
typedef vector<double>::iterator dit_t;

//----------------------------------------------------------------------

void Gradient_data::store(ostream  & os) { 
  int ndim=grad.size();
  os << "ndim " << ndim << endl;
  os << "x "; for(vector<double>::iterator i=x.begin(); i!=x.end(); i++) os << *i << " ";
  os << endl;
  os << "grad "; for(vector<double>::iterator i=grad.begin(); i!=grad.end(); i++) os << *i << " ";
  os << endl;
  os << "sigma "; for(vector<double>::iterator i=sigma.begin(); i!=sigma.end(); i++) os << *i << " ";
  os << endl;

}
//----------------------------------------------------------------------

void Gradient_data::read(istream & is) { 
  int ndim; string dummy;
  is >> dummy; assert(dummy=="ndim");
  x.resize(ndim);
  grad.resize(ndim);
  sigma.resize(ndim);
  is >> dummy; assert(dummy=="x");
  for(dit_t i=x.begin(); i!=x.end(); i++) is >> *i;
  is >> dummy; assert(dummy=="grad");
  for(dit_t i=grad.begin(); i!=grad.end(); i++) is >> *i;
  is >> dummy; assert(dummy=="sigma");
  for(dit_t i=sigma.begin(); i!=sigma.end(); i++) is >> *i;
  

}

//----------------------------------------------------------------------
void Hess_grad::generate_guess(const vector <Gradient_data> & data, vector <double> & c) { 
  int ndata=data.size();
  assert(ndata >0);
  int n=data[0].grad.size();
  assert(data[0].grad.size()==data[0].x.size()==data[0].sigma.size());
  int nparms=n*(n-1)/2+n;
  c.resize(nparms);
  for(int i=0; i< n; i++) c[i]=0.0;
  int count=n;
  for(int i=0; i< n; i++) { 
    for(int j=i; j< n; j++) { 
      if(i==j) c[count++]=1.0;
      else c[count++]=0.0;
    }
  }
}

//----------------------------------------------------------------------
double Hess_grad::prob(const vector <Gradient_data> & data, const vector <double> & c)  { 
   assert(data.size() > 0);
   int n=data[0].grad.size();

   double prob=0;
   vector <double> mod_grad; mod_grad.resize(n);
  
   Array2 <double> hess(n,n);
   int count=n;
   for(int i=0; i<n; i++) { 
     for(int j=i; j< n; j++) { 
       hess(i,j)=hess(j,i)=c[count++];
     }
   }
   
   for(vector<Gradient_data>::const_iterator g=data.begin(); g!=data.end(); g++) {
     for(vector<double>::iterator i=mod_grad.begin(); i!=mod_grad.end(); i++) *i=0.0;
     
     for(int i=0;i < n; i++) { 
       for(int j=0; j < n; j++) { 
          mod_grad[i]+=2.0*hess(i,j)*(g->x[j]-c[j]);
          count++;
       }
     }
     for(int i=0; i< n; i++) { 
        prob-=(mod_grad[i]-g->grad[i])*(mod_grad[i]-g->grad[i])/(2*g->sigma[i]*g->sigma[i]);
     }
   }
   return prob;
}
//----------------------------------------------------------------------

double Hess_grad::grad_prob(const vector <Gradient_data> & data, const vector <double> & c,
    vector <double> & grad) { 
  grad.resize(c.size());
  assert(data.size() > 0);
  int n=data[0].grad.size();
  for(vector<double>::iterator i=grad.begin(); i!=grad.end(); i++) *i=0;

  Array2 <double> hess(n,n);
  int count=n;
  for(int i=0; i<n; i++) { 
    for(int j=i; j< n; j++) { 
      hess(i,j)=hess(j,i)=c[count++];
    }
  }

  //for(int i=0; i< n; i++) { 
  //  cout << "hess ";
  //  for(int j=0; j< n; j++) { 
  //    cout << hess(i,j) << " ";
  //  }
  //  cout << endl;
  //}
  vector <double> tmpvec(n);
  for(vector<Gradient_data>::const_iterator g=data.begin(); g!=data.end(); g++) { 
    for(int i=0; i< n; i++) { 
      tmpvec[i]=0.0;
      for(int j=0; j< n; j++) { 
        tmpvec[i]+=2.0*hess(i,j)*(g->x[j]-c[j]);
      }
      tmpvec[i]-=g->grad[i];
      tmpvec[i]/=g->sigma[i]*g->sigma[i];
    }

    for(int i=0; i< n; i++) { 
      for(int j=0; j< n; j++) { 
        grad[j]+=2.0*tmpvec[i]*hess(i,j);
      }
    }

    count=n;
    for(int i=0; i< n; i++) { 
      for(int j=i; j< n; j++) { 
        grad[count]-=2.0*tmpvec[i]*(g->x[j]-c[j]);
        if(i!=j) grad[count]-=2.0*tmpvec[j]*(g->x[i]-c[i]);
        count++;
      }
    }
  }
}

//----------------------------------------------------------------------


void test_hess_gradients(const vector <Gradient_data> & data, const vector<double> & c) { 
  cout << "testing gradients " << endl;
  vector <double> grad,cp=c;
  Hess_grad grad_hess;

  for(vector<Gradient_data>::const_iterator i=data.begin(); i!=data.end(); i++) {
    cout << "x " << i->x[0] << " g " << i->grad[0] << " sigma " << i->sigma[0] << endl;
  }
  cout << "coefficients ";
  for(vector <double>::iterator i=cp.begin(); i!=cp.end(); i++) cout << *i << " ";
  cout << endl;

  grad_hess.grad_prob(data,c,grad);
  double base=grad_hess.prob(data,c);
  cout << "calculated probability " << base << endl;
  double del=1e-9;
  int nparms=grad.size();
  for(int i=0; i< nparms; i++) { 
    double sv=c[i];
    cp[i]+=del;
    double d1=grad_hess.prob(data,cp);
    cp[i]-=2.0*del;
    double d2=grad_hess.prob(data,cp);
    double g_finite=(d1-d2)/(2*del);
    cp[i]=sv;
    cout << "checking gradient: finite diff " << g_finite << "  calc " << grad[i] << " diff " << grad[i]-g_finite << endl;
  }

}

//----------------------------------------------------------------------

