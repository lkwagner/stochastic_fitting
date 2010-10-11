#include "Minimum_inference.h"
using namespace std;





int main(int argc, char ** argv) { 
  
  int nit=15; 
  
  string dummy;
  Data_generator * pes;
  Line_model * mod=new Cubic_model;
  vector <double> currmin;
  double start_sigma=-1,trust_rad=-1;
  while(cin >> dummy) { 
    if(dummy=="iterations") cin >> nit;
    else if(dummy=="random_quad") { 
      int ndim; cin >> ndim;
      pes=new Random_quadratic(ndim);
    }
    else if(dummy=="montecarlo_caller") { 
      pes=new MonteCarlo_caller(cin);
    }
    else if(dummy=="random_cubic") { 
      pes=new Random_cubic(cin);
    }
    else if(dummy=="start_sigma") cin >> start_sigma;
    else if(dummy=="trust_radius") cin >> trust_rad;
    else if(dummy=="seed") { long int s1,s2;
      cin >> s1 >> s2;
      rng.seed(s1,s2);
    }
    else if(dummy=="currmin") {
      if(pes==NULL) error("PES not defined before minimum");
      int ndim=pes->ndim();
      double dum=0;
      for(int i=0; i< ndim; i++) {
        cin >> dum; currmin.push_back(dum);
      }
    }
    else { cerr << "Didn't understand " << dummy << endl;
      exit(132);
    }
  }
  int ndim=pes->ndim();
  if(trust_rad <0) error("invalid value for trust_rad");
  if(start_sigma<0) error("invalid value for start_sigma");
  if(currmin.size()!= ndim) error("currmin has wrong number of dimensions");

  Minimum_inference min_infer;
  Search_data search(ndim);
  search.trust_rad=trust_rad;
  search.sigma=start_sigma;
  search.currx=currmin;
 
  

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



