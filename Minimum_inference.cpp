#include "Minimum_inference.h"
#include "generate_line.h"
#include "Min.h"

int Minimum_inference::addLine(Data_generator & pes, Line_model & mod, Search_data & search, int d) {  
  cout << "addline " << endl;
   assert(search.currx.size()==search.directions[d].size());
   int ndim=search.currx.size();
   Line_data nwdata;
   nwdata.direction=search.directions[d];
   nwdata.start_pos=search.currx;
   Fit_info  finfo;
   int maxpts=10;
   vector <double> sig_line(ndim);
   for(vector<double>::iterator i=sig_line.begin(); i!=sig_line.end(); i++) *i=search.sigma;
   find_minimum(search.trust_rad,maxpts, pes, mod, nwdata,finfo, sig_line);
   lines.push_back(nwdata);
   models.push_back(&mod);
   fits.push_back(finfo);
   for(int i=0; i< ndim; i++) { 
     search.currx[i]+=finfo.min[0]*search.directions[d][i];
   }


   cout << "Minimum " << finfo.min[0] << " +/- " << finfo.minerr[0]  << endl;
   cout << "Curvature " << finfo.curve << " +/- " << finfo.curveerr << endl;
   cout << "Funcmin " << finfo.funcmin << " +/- " << finfo.funcminerr << endl;
}

//----------------------------------------------------------------------

int Minimum_inference::addForce(Data_generator & pes, double sigma, vector <double> & x) { } 
void Minimum_inference::readState(istream & is) { } 



void Minimum_inference::calcHess(double trust_rad, double & E0, vector <double> & min, 
    vector < vector < double > > & hess) {
  vector <double> c;
  double best_err=1e18;
  Quad_plus_line quad;
  vector <Fix_information> fixes;
  shake_quad(quad,lines,models,fixes,c);
  int ndim=lines[0].direction.size();
  quad.get_hessian(c,ndim,hess);
  quad.get_minimum(c,ndim,min);

} 
void Minimum_inference::saveState(ostream & os) { } 

