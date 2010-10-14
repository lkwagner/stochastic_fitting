#include "Minimum_inference.h"
using namespace std;


void proj(vector <double> & u, vector<double> & v, vector <double> & p) { 
  int n=u.size();
  assert(u.size()==v.size());
  p.resize(n);
  double vdotu=0,udotu=0;
  for(int i=0; i< n; i++) {
    vdotu+=v[i]*u[i];
    udotu+=u[i]*u[i];
  }
  for(int i=0; i< n; i++) { 
    p[i]=u[i]*vdotu/udotu;
  }

}

//orthogonalize vector k with respect to vectors 0-k-1
void gram_schmidt(vector <vector <double> > & vecs, int k) { 
  int ndim=vecs.size();
  vector <double> u(ndim),p(ndim);
  u=vecs[k];
  for(int d=0; d< k; d++) {
    proj(vecs[d],u,p);
    for(int d1=0; d1 < ndim; d1++) { 
      u[d1]-=p[d1];
    }
  }   
  double norm=0;
  for(int d=0; d< ndim; d++) norm+=u[d]*u[d];
  norm=sqrt(norm);
  for(int d=0; d< ndim; d++) vecs[k][d]=u[d]/norm;
}


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


  //first iteration..let's do kind of a conjugate gradient approach
  cout << "getting first directions " << endl; 
  for(int d=0; d< ndim; d++) { 
    Gradient_data tmpgrad; tmpgrad.x=search.currx;
    if(pes->gradient(tmpgrad.x,tmpgrad.grad,tmpgrad.sigma)) { 
      cout << " have gradient " << endl;
      min_infer.addGradient(tmpgrad);
      double norm=0;
      for(int d1=0; d1 < ndim; d1++) norm+=tmpgrad.grad[d1]*tmpgrad.grad[d1];
      norm=sqrt(norm);
      for(int d1=0; d1 < ndim; d1++) search.directions[d][d1]=tmpgrad.grad[d1]/norm;
      gram_schmidt(search.directions,d);
      min_infer.addLine(*pes,*mod,search,d);
      cout << "currx "; for(int i=0; i< ndim; i++) cout << search.currx[i] << " ";
      cout << endl;
      cout << "new direction "; for(int i=0; i< ndim; i++) cout << search.directions[d][i] << " ";
      cout << endl;
      ofstream statesave("state");
      min_infer.saveState(statesave);

    }
    else { 
      for(int d1=0;d1 < ndim; d1++) { search.directions[d][d1]=d==d1?1:0; }
    }
  }
  cout << "==============moving to normal iterations " << endl;

  //regular iterations after this

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
      ofstream statesave("state");
      min_infer.saveState(statesave);
      
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



