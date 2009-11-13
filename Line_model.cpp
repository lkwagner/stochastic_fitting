#include "Line_model.h"
#include "ulec.h"

typedef vector<double>::const_iterator cdit_t;



const double Line_model::prob(const Line_data & data ,const Fix_information & fix, 
                        const vector <double> & c) { 
  static vector <double> realc;  //unsafe in multi-threaded operation, but makes us faster
  if(fix.enforce) { 
    convert_c(fix,c,realc);
  }
  else { 
    realc=c;
  }
  
  double f=0;
  double fp, t1;
  for(vector <Data_point>::const_iterator i=data.data.begin(); 
      i!=data.data.end(); i++) {
    fp=func(realc,i->t);
    t1=i->val-fp;
    f-=t1*t1*0.5*(i->inverr)*(i->inverr);
  }

  return f;
}

//###############################################################################

//------------------------------------------------------------------------------

const void Morse_model::minimum(const vector <double> & c, vector <double> & min) { 
  assert(c.size() >=4);
  min.resize(1);
  min[0]=c[3];
}

//------------------------------------------------------------------------------

const double Morse_model::func(const vector <double> & c, double t) const { 
  assert(c.size() >=4 );
  double tmp=1-exp(-c[2]*(t-c[3]));
  return c[0]+c[1]*tmp*tmp;
}


//------------------------------------------------------------------------------

const void Morse_model::generate_guess(const Line_data & data, const Fix_information & fix,
                                   vector <double> & c) { 
  if(fix.enforce==1) { 
    c.resize(1);
    c[0]=50*(rng.ulec()-0.5);
  }
  else { 
    c.resize(4);
    c[0]=data.data[0].val;
    
    c[1]=1.0*(rng.ulec()-0.5);
    c[2]=1.0*(rng.ulec()-0.5);
    double avg=0;
    
    for(vector<Data_point>::const_iterator t=data.data.begin(); t!= data.data.end(); t++) {
      avg+=t->t;
      if(c[0] < t->val) c[0]=t->val;
    }
    avg/=data.data.size();
    c[3]=avg;
    //cout << "guessing " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << endl;

  }
  
  
}

//------------------------------------------------------------------------------

const void Morse_model::convert_c(const Fix_information & fix, const vector <double> & c_in, 
                                   vector <double> & c_out) { 
  assert(c_in.size() >=1);
  c_out.resize(4);
  c_out[2]=c_in[0];
  c_out[1]=fix.curve/(2.0*c_out[2]*c_out[2]);
  c_out[0]=fix.valmin;
  c_out[3]=fix.min;
}
//###############################################################################

//------------------------------------------------------------------------------

const void Cubic_model::minimum(const vector <double> & c, vector <double> & min) { 
  assert(c.size() >=4);
  min.resize(1);
  min[0]=c[1];
}

//------------------------------------------------------------------------------

const double Cubic_model::func(const vector <double> & c, double t) const { 
  assert(c.size() >=4 );
  return c[0]+c[2]*(t-c[1])*(t-c[1])+c[3]*(t-c[1])*(t-c[1])*(t-c[1]);
}


//------------------------------------------------------------------------------

const void Cubic_model::generate_guess(const Line_data & data, const Fix_information & fix,
                                       vector <double> & c) { 
  if(fix.enforce==1) { 
    c.resize(1);
    c[0]=1.0*(rng.ulec()-0.5);
  }
  else { 
    c.resize(4);
    //c[0]=data.data[0].val;
    
    c[2]=1.0*(rng.ulec()-0.5);
    c[3]=1.0*(rng.ulec()-0.5);
    double avg=0;
    
    for(vector<Data_point>::const_iterator t=data.data.begin(); t!= data.data.end(); t++) {
      avg+=t->t;
      if(c[0] < t->val) c[0]=t->val;
    }
    avg/=data.data.size();
    c[1]=avg;
    //cout << "guessing " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << endl;
    
  }
  
  
}

//------------------------------------------------------------------------------

const void Cubic_model::convert_c(const Fix_information & fix, const vector <double> & c_in, 
                                  vector <double> & c_out) { 
  assert(c_in.size() >=1);
  c_out.resize(4);
  c_out[3]=c_in[0];
  c_out[2]=fix.curve*0.5;
  c_out[0]=fix.valmin;
  c_out[1]=fix.min;
}

//------------------------------------------------------------------------------


double Cubic_model::curve(const vector <double> & c) { 
  assert(c.size() >=4);
  return c[2]*2.0;
}

const void Cubic_model::downconvert_c(const Fix_information & fix, 
                         const vector <double> & c_in, vector <double> & c_out){
  c_out.resize(1);
  assert(c_in.size() >= 4);
  c_out[0]=c_in[3];
}


//###############################################################################

const void Quadratic_model::minimum(const vector <double> & c, vector <double> & min) { 
  assert(c.size() >=3);
  min.resize(1);
  min[0]=c[1];
}

//------------------------------------------------------------------------------

const double Quadratic_model::func(const vector <double> & c, double t) const { 
  assert(c.size() >=3 );
  return c[0]+(t-c[1])*(t-c[1])*c[2];
}
//------------------------------------------------------------------------------

const void Quadratic_model::generate_guess(const Line_data & data, const Fix_information & fix,
                                       vector <double> & c) { 
  if(fix.enforce==1) { 
    c.resize(0);
  }
  else { 
    c.resize(3);
    c[0]=data.data[0].val;
    
    c[2]=10.0*(rng.ulec()-0.5);
    double avg=0;
    
    for(vector<Data_point>::const_iterator t=data.data.begin(); t!= data.data.end(); t++) {
      avg+=t->t;
      if(t->val < c[0]) c[0]=t->val;
    }
    avg/=data.data.size();
    c[1]=avg;
    //cout << "guessing " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << endl;
    
  }
  
  
}

//------------------------------------------------------------------------------

const void Quadratic_model::convert_c(const Fix_information & fix, const vector <double> & c_in, 
                                  vector <double> & c_out) { 
  assert(c_in.size() >=0);
  c_out.resize(3);
  c_out[2]=0.5*fix.curve;
  c_out[0]=fix.valmin;
  c_out[1]=fix.min;
}

