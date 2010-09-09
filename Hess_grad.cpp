
#include "Hess_grad.h"
#include <cassert>
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
   for(vector<Gradient_data>::const_iterator g=data.begin(); g!=data.end(); g++) {
     for(vector<double>::iterator i=mod_grad.begin(); i!=mod_grad.end(); i++) *i=0.0;
     int count=n;
     for(int i=0;i < n; i++) { 
       for(int j=i; j < n; j++) { 
          mod_grad[i]+=c[count]*(g->x[j]-c[j]);
          mod_grad[j]+=c[count]*(g->x[i]-c[i]);
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
