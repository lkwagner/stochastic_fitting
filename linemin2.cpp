#include "Minimum_inference.h"
using namespace std;





int main(int argc, char ** argv) { 
  
  int nit=15; 
  
  string dummy;
  int ndim=6;
  Data_generator * pes=new Random_quadratic(ndim);
  Line_model * mod=new Cubic_model;
  vector <double> currmin;
 
  Minimum_inference min_infer;


  Search_data search(ndim);
  search.trust_rad=.4;
  search.sigma=1.0;
 
  

  vector <vector <double> > hess; double e0; vector <double> hess_min;

  for(int it=0; it < nit; it++) { 
    for(int d=0; d< ndim; d++) { 
      min_infer.addLine(*pes,*mod,search,d);
      cout << "currx "; for(int i=0; i< ndim; i++) cout << search.currx[i] << " ";
      cout << endl;
      Gradient_data tmpgrad; tmpgrad.x=search.currx;
      if(pes->gradient(tmpgrad.x,tmpgrad.grad,tmpgrad.sigma)) { 
        cout << "adding gradient!" << endl;
        min_infer.addGradient(tmpgrad);
      }
    }

    int restart=1; if(it==0) restart=0;
    min_infer.calcHess(search.trust_rad,e0,hess_min,hess,restart);
    cout << "hess: " << endl;
    for(int i=0; i< ndim; i++) { 
      cout << "hess: " ;
      for(int j=0; j< ndim; j++) { 
        cout << hess[i][j] << " ";
      }
      cout << endl;
    }
    vector <double> evals;
    search.update_directions(hess,evals);
    cout << "evals ";
    for(int n=0; n < ndim; n++)  cout << evals[n] << " ";
    cout << endl;
  }


  delete pes;
  delete mod;
}
//------------------------------------------------------------------



