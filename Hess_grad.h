#ifndef HESS_GRAD_H_INCLUDED
#define HESS_GRAD_H_INCLUDED
#include <vector>
using namespace std;
class Gradient_data { 
  public:
    vector <double> grad;
    vector <double> x;  //position at which the gradient was evaluated
    vector <double> sigma; //uncertainty..
};


class Hess_grad { 
  public:
    void generate_guess(const vector <Gradient_data> & data, vector <double> & c);
    double prob(const vector <Gradient_data> & data, const vector <double> & c);
};

#endif //HESS_GRAD_H_INCLUDED
