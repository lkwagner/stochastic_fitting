#include <cassert>
#include <cmath>
#include <vector>
#include <iostream>
#include "ulec.h"
#include "Point.h"
#include "Prob_model.h"
#include "macopt.h"
#include "Min.h"
#include <iomanip>

#include "Array.h"
#include "MatrixAlgebra.h"

using namespace std;



//------------------------------------------------------------------

typedef vector<double>::iterator dit_t;



void append_number(string & str, int num){
  char strbuff[100];
  sprintf(strbuff, "%d", num);
  str+=strbuff;
}

/*
 File format:
 pos1 pos2  en err
 etc.
 */


//------------------------------------------------------------------

void fit_line(Prob_model & mod, Fit_info & finfo, int verbose=1) { 
  vector <double> c;
  mod.generate_guess(c);
  finfo.cer=c;
  vector <double> best_c=c;
  for(dit_t i=finfo.cer.begin(); i!= finfo.cer.end(); i++) *i=0;
  
  if(verbose) 
    mod.niceprint(c,finfo.cer);
  
  double best_p=-1e99;
  
  int nit=1000;
  int nparms=c.size();
  
  for(int i=0; i< nit; i++) { 
    Least_squares_opt opt(nparms,0,.1,50,1);
    mod.generate_guess(opt.c);
    opt.mod=&mod;    
    opt.optimize();
    double p=mod.probability(opt.c);
    cout.flush();
//cout << "prob " << p << endl;
    if(p > best_p && mod.is_ok(opt.c) ) { 
      best_c=opt.c;
      best_p=p;
      if(verbose) { 
        cout << "\n new best ";
        for(dit_t i=opt.c.begin(); i!= opt.c.end(); i++) cout << *i << " ";
        cout << p << endl;
      }
    }
  }
  
  if(verbose) { 
    cout << "final from maximum likelihood: \n";
    mod.niceprint(best_c, finfo.cer);
  }
  //check_model_consistency(mod,  best_c);
  
  sample_min(best_c, mod,finfo, verbose);    
  
}

//------------------------------------------------------------------

//given a set of positions in data, set the energies
void fill_hilbert(vector <Point> & data) { 
  static bool first_run=true;
  assert(data.size() >=1);
  int n=data[0].x.size();
  
  vector <vector <double> > H;
  H.resize(n);
  for(int i=0; i< n; i++) H[i].resize(n);
  assert(n==3);
  //H[0][0]=1.69943; H[1][0]=H[0][1]=0.480298; H[2][0]=H[0][2]=0.726374;
  //H[1][1]=0.384995; H[2][1]=H[1][2]=0.0218349;
  //H[2][2]=0.621747;  
  H[0][0]=1.69943; H[1][0]=H[0][1]=0.480298; H[2][0]=H[0][2]=0.726374;
  H[1][1]=0.384995; H[2][1]=H[1][2]=0.0218349;
  H[2][2]=0.621747;  
  
  if(first_run) { 
    first_run=false;
    Array2 <double> A(n,n);
    for(int i=0; i< n; i++) { 
      for(int j=0; j < n; j++) { 
        A(i,j)=H[i][j];
      }
    }
    Array1 <double> evals(n);
    Array2 <double> evecs(n,n);
    EigenSystemSolverRealSymmetricMatrix(A, evals, evecs);
    cout << "eigenvalues " ;
    for(int i=0; i< n; i++) cout << evals(i) << " ";
    cout << endl;
    cout << "eigenvectors: " << endl;
    for(int i=0; i< n; i++) { 
      for(int j=0; j< n; j++) { 
        cout << setw(15) << evecs(i,j);
      }
      cout << endl;
    }
  }
  
  for(vector<Point>::iterator d=data.begin(); d!= data.end(); d++) { 
    d->en=0;
    d->err=0;
    for(int i=0; i< n; i++) {
      for(int j=0; j< n; j++) { 
        d->en+=d->x[i]*H[i][j]*d->x[j];
      }
    }
    //cout << "data: ";
    //for(int i=0; i< n; i++) cout << d->x[i] << " ";
    //cout << " en " << d->en << endl;
    /*
    d->en+=d->x[0]*d->x[0]-2.0*d->x[0];
    for(int i=1; i< n; i++) { 
      d->en+=d->x[i]*(d->x[i]*2.0-2.0*d->x[i-1]);
      //for(int j=0; j< n; j++) { 
      //  d->en+=d->x[i]*d->x[j]/(i+j+1); //The plus 1 is because of indices
      //}
    }
     */
  }
}

void add_noise(double noise, vector <Point> & data ) { 
  for(vector<Point>::iterator d=data.begin(); d!= data.end(); d++) { 
    d->en+=noise*rng.gasdev();
    d->err=noise;
  }
}


void generate_points(const vector <double> & direction, 
                     const vector <double> & center,
                     double trust, int npoints, vector <Point> & data,
                     vector <Point> & linear_data) { 
  data.resize(npoints);
  int n=center.size();
  assert(center.size()==direction.size());
  linear_data.resize(npoints);
  for(int i=0; i< npoints; i++) { 
    data[i].x=center;
    double beta=trust*(double(i)/npoints-0.5);
    linear_data[i].x.resize(1);
    linear_data[i].x[0]=beta;
    for(int j=0; j< n; j++) { 
      data[i].x[j]=center[j]+direction[j]*beta;
    }
  }
}
//------------------------------------------------------------------
void move(const vector <double> & direction, vector <double> & currmin, 
              double trust, int npoints, Prob_model & mod, Fit_info & finfo) { 
  int n=currmin.size();
  
  vector <Point> data,linear_data;
  generate_points(direction,currmin,trust,npoints, data,linear_data);
  fill_hilbert(data);
  add_noise(0.001,data);
  for(int i=0; i< npoints; i++) { 
    linear_data[i].en=data[i].en;
    linear_data[i].err=data[i].err;
  }
  mod.set_data(linear_data);
  fit_line(mod, finfo,0);
  for(int i=0; i < n; i++) { 
    currmin[i]+=direction[i]*finfo.min[0];
  }
  cout << "new point in subiteration " << " : " ;
  for(int d=0; d< n; d++) cout <<  currmin[d] << "  " ;
  cout << endl;
  
}
//------------------------------------------------------------------

void find_new_directions(vector <Fit_info> & finfo, 
                         vector < vector <double> > & directions) {
  int n=finfo.size();
  assert(n==directions.size());
  assert(n==directions[0].size());
  Array2 <double> U(n,n,0.0);
  Array2 <double> D(n,n,0.0);
  for(int i=0; i< n; i++) { 
    for(int j=0; j< n; j++) { 
      U(i,j)=directions[j][i];
    }
    D(i,i)=1.0/finfo[i].curv[0];
  }
  
  cout << "U : " << endl;
  for(int i=0; i< n; i++) { 
    for(int j=0; j< n; j++) { 
      cout << setw(15) << U(i,j);
    }
    cout << endl;

  }
  
  
  Array2 <double> tmp(n,n,0.0);
  MultiplyMatrices(U,D,tmp,n);
  Array2 <double> Ainv(n,n,0.0), A(n,n,0.0);
  TransposeMatrix(U,n);
  MultiplyMatrices(tmp,U,Ainv,n);
  InvertMatrix(Ainv,A,n);
  
  cout << "curvatures " << endl;
  for(int i=0; i< n; i++) { 
    cout << D(i,i) << "  ";
  }
  cout << endl;
  
  cout << "Guess at the Hessian matrix: " << endl;
  for(int i=0; i< n; i++) { 
    for(int j=0; j< n; j++) { 
      cout << setw(15) << A(i,j) ;
    }
    cout << endl;
  }
  
  Array1 <double> evals(n);
  Array2 <double> evecs(n,n);
  EigenSystemSolverRealSymmetricMatrix(A, evals, evecs);
  
  cout << "eigenvalues " ;
  for(int i=0; i< n; i++) cout << evals(i) << " ";
  cout << endl;
  cout << "eigenvectors: " << endl;
  for(int i=0; i< n; i++) { 
    for(int j=0; j< n; j++) { 
      cout << setw(15) << evecs(i,j);
    }
    cout << endl;
  }
  
  
  
}

//------------------------------------------------------------------

int main(int argc, char ** argv) { 
  cout.precision(10);
  
  Morse_model mod;
  //Fit_info finfo;
  int n=3;
  vector <double> currmin;
  currmin.resize(n);
  for(dit_t i=currmin.begin(); i!= currmin.end(); i++) *i=1.0;
  vector <vector <double> > directions;
  directions.resize(n);
  for(int i=0; i< n; i++) { 
    directions[i].resize(n);
    for(int j=0; j< n; j++) directions[i][j]=0.0;
    directions[i][i]=1.0;
  }
  
  double trust=3.0;
  int npoints=8;
  
  for(int it=0; it < 10; it++) { 
    vector <Fit_info> finfo(n);

    for(int douter=0; douter < n; douter++) { 

      vector <double> oldmin=currmin;
      //Do a round of minimizations
      for(int d=0; d< n; d++) { 
        move(directions[d],currmin, trust, npoints, mod, finfo[d]);
      }
      
      
      //Find the total move and minimize along that direction.
      
      vector <double> n_direction(n);
      double norm=0;
      for(int d=0; d < n; d++) { 
        n_direction[d]=currmin[d]-oldmin[d];
        norm+=n_direction[d]*n_direction[d];
      }
      norm=sqrt(norm);
      for(int d=0; d< n; d++) n_direction[d]/=norm;
      
      Fit_info nfinfo;
      move(n_direction,currmin, trust, npoints, mod, nfinfo);
    
      cout << "grep_point " ;
      for(int d=0; d< n; d++) cout <<  currmin[d] << "  " ;
      cout << endl;
      
      for(int i=0; i < n-1; i++) { 
        directions[i]=directions[i+1];
        finfo[i]=finfo[i+1];
      }
      directions[n-1]=n_direction;
      finfo[n-1]=nfinfo;
      
      cout << "new directions " << endl;
      for(int i=0; i < n; i++) { 
        for(int j=0; j< n; j++) { 
          cout << setw(10) << directions[i][j] << " ";
        }
        cout << endl;
      }
      
    }
    find_new_directions(finfo, directions);
 
  }
  
  /*
  mod.read(cin);
  Fit_info finfo;
  fit_line(mod, finfo,1);
  mod.niceprint(finfo.cavg, finfo.cer);
  
  cout << "minimum " << finfo.min[0] << " +/- " << finfo.minerr[0] << endl;
  */

}

//------------------------------------------------------------------



