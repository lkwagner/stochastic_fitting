#ifndef MINIMUM_INFERENCE_H_INCLUDED
#define MINIMUM_INFERENCE_H_INCLUDED
using namespace std;
#include "Line_model.h"
#include "Hess_grad.h"
#include "Data_generator.h"
#include "Sample.h"



struct Search_data { 
  Search_data();
  Search_data(int ndim) { 
    currx.resize(ndim);
    directions.resize(ndim);
    for(int i=0; i< ndim; i++) { 
      currx[i]=0.0;
      directions[i].resize(ndim);
      for(int j=0; j< ndim; j++) directions[i][j]= i==j?1:0;
    }
  }
  double trust_rad,sigma;
  vector <double> currx;
  vector < vector < double> > directions;
 
};


struct Minimum { 
  vector <double> min;
  vector <double> err; //standard deviation
};

class Minimum_inference { 
  public:
    //modify the internal state of the object, don't modify the parameters
    int addLine(Data_generator & pes, Line_model & mod,  Search_data & search,int d);
    int addForce(Gradient_data & grad);
    void readState(istream & is);

    //Don't modify the internal state, only return information
    void calcHess(double trust_rad, double & E0, vector <double> & min, 
        vector < vector < double > > & hess);
    void saveState(ostream & os);
  private:
    vector <Line_data> lines;
    vector <Fit_info> fits;
    vector <Gradient_data> gradients;
    vector <Line_model *> models;
};


#endif //MINIMUM_INFERENCE_H_INClUDED
