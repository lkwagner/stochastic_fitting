#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;


int main() { 
  int n=3;
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
  
  vector <double> delx(n);
  for(int d=0; d< n; d++) delx[d]=1;
  
  for(int it=0; it < 100; it++)  { 
    for(int d=0; d< n; d++) { 
      delx[d]=0;
      for(int i=0; i< n; i==d-1?i+=2:i++)  {
        delx[d]+=H[i][d]*delx[i]/H[d][d];
      }
      cout << "it " << it;
      for(int i=0; i< n; i++)  cout << setw(15) << delx[i];
      cout << endl;
    }
  }
  return 0;
}