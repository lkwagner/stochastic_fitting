#ifndef QUAD_PLUS_LINE_H_INCLUDED
#define QUAD_PLUS_LINE_H_INCLUDED
#include <vector>
#include <cmath>
#include "Line_model.h"

using namespace std;


class Quad_plus_line { 
public:
  
  void generate_guess(const vector <Line_data> & data, const vector <Line_model *> & models,
                      vector <double> & c);
  void generate_guess(const vector <Line_data> & data, const vector <Line_model *> & models, 
      const double e0, const vector <double> & guess_min, const vector < vector <double> > & guess_hess,
      vector <double> & c);
  double prob(const vector <Line_data> & data, const vector <Line_model * > & models,
              const vector <double> & c, vector <Fix_information> & fixes);
  double delta_prob(const vector <Line_data> & data, const vector <Line_model * > & models,   
                    vector <Fix_information> & fixes,vector <double> & c, int p, double cp);
  
  
  void set_fixes(const vector <Line_data> & data, 
                 const vector <double> & c,
                 vector <Fix_information> & fixes);
  void update_fixes(const vector <Line_data> & data, 
                     vector <double> & c,
                    vector <Fix_information> & fixes, int p, double cp);
  void get_hessian(const vector <double> & c, int ndim, 
                   vector < vector < double> > & H);
  void get_minimum(const vector <double> & c, int ndim, 
                   vector < double>  & m);
    
  bool has_minimum(const vector <double> & c, int ndim);
  
};

#include "Array.h"
bool is_posdef(int n, Array2 <double> & C); 

void generate_posdef_matrix(int n, Array2 <double> & );

#endif //QUAD_PLUS_LINE_H_INCLUDED
