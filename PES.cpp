#include "MatrixAlgebra.h"
#include "ulec.h"
#include "PES.h"
#include <iostream>
#include "Quad_plus_line.h"

using namespace std;
void Potential_energy_surface::gen_pes(int n) { 
  generate_posdef_matrix(n,C);
  /*
  C.Resize(n,n);
  for(int i=0; i< n; i++) { 
    for(int j=0; j< n; j++) { 
      if(i==j) C(i,j)=1;
      else if(fabs(i-j)==1) C(i,j)=-0.5;
      else C(i,j)=0.0;
    }
  }
   */
            
  Array2 <double> evecs(n,n);
  Array1 <double> evals(n);
  EigenSystemSolverRealSymmetricMatrix(C,evals,evecs);

  minima.Resize(n);
  for(int i=0; i< n; i++) minima[i]=2.0*(rng.ulec()-0.5);
  
  for(int i=0; i< n; i++) { 
    cout << "Hgen: ";
    for(int j=0; j< n; j++) { 
      cout << C(i,j) << " ";
    }
    cout << endl;
  }
  cout << endl;
  cout << "minima ";
  for(int i=0; i< n; i++) cout << minima[i] << " ";
  cout << endl;
  
  cout << "Evecs " << endl;
  for(int i=0; i< n; i++) { 
    cout << "d" << i << " ";
    for(int j=0;j< n; j++) cout << evecs(i,j) << " ";
    cout << endl;
  }
  
  cout << "Evals " << endl;
  for(int i=0; i< n; i++) cout << evals(i) << " ";
  cout << endl;
}


double Potential_energy_surface::eval_pes(const vector <double> & x) { 
  //cout << "sizes " << x.size() << " " << C.GetDim(0) << endl;
  assert(x.size()==C.GetDim(0));
  
  int n=C.GetDim(0);
  double f=0;
  for(int i=0; i< n; i++) { 
    for(int j=0; j< n; j++) { 
      f+=(x[i]-minima[i])*C(i,j)*(x[j]-minima[j]);
    }
  }
  return f;
}


void Potential_energy_surface::eval_grad(const vector <double> & x,vector <double> & grad) { 
 assert(x.size()==C.GetDim(0));
 int n=C.GetDim(0);
 grad.resize(n);
 for(vector<double>::iterator g=grad.begin(); g!=grad.end(); g++) *g=0.0;
 for(int i=0; i< n; i++) { 
   for(int j=0; j< n; j++) { 
     grad[i]+=2.0*C(i,j)*(x[j]-minima[j]);
   }
 }

}

