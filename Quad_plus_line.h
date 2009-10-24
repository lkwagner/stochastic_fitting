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
  double prob(const vector <Line_data> & data, const vector <Line_model * > & models,
              const vector <double> & c);
  
  
  void set_fixes(const vector <Line_data> & data, 
                 const vector <double> & c,
                 vector <Fix_information> & fixes);
  void get_hessian(const vector <double> & c, int ndim, 
                   vector < vector < double> > & H);
  void get_minimum(const vector <double> & c, int ndim, 
                   vector < double>  & m);

  
};

#include "Array.h"

void generate_posdef_matrix(int n, Array2 <double> & );

#endif //QUAD_PLUS_LINE_H_INCLUDED