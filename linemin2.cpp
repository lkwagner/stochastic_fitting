#include "Minimum_inference.h"
using namespace std;


int main(int argc, char ** argv) { 
  
  int nit=15; 
  
  string dummy;
  int ndim=10;
  Data_generator * pes=new Random_quadratic(ndim);
  Line_model * mod=new Cubic_model;
  vector <double> currmin;
 
  Minimum_inference min_infer;


  Search_data search(ndim);
  search.trust_rad=.4;
  search.sigma=1.0;
 
  


  for(int it=0; it < nit; it++) { 
    for(int d=0; d< ndim; d++) { 
      min_infer.addLine(*pes,*mod,search,d);
      cout << "currx "; for(int i=0; i< ndim; i++) cout << search.currx[i] << " ";
      cout << endl;
    }

    vector <vector <double> > hess; double e0; vector <double> hess_min;
    min_infer.calcHess(search.trust_rad,e0,hess_min,hess);
    for(int i=0; i< ndim; i++) { 
      cout << "hess: " ;
      for(int j=0; j< ndim; j++) { 
        cout << hess[i][j] << " ";
      }
      cout << endl;
    }
  }


  delete pes;
  delete mod;
}
//------------------------------------------------------------------



