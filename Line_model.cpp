#include "Line_model.h"
#include "ulec.h"

typedef vector<double>::const_iterator cdit_t;



double Line_model::prob(const Line_data & data ,const Fix_information & fix, 
                        const vector <double> & c) { 
  vector <double> realc;
  if(fix.enforce) { 
    convert_c(fix,c,realc);
  }
  else { 
    realc=c;
  }
  
  double f=0;
  assert(data.t.size()==data.val.size());
  assert(data.t.size()==data.err.size());
  
  for(cdit_t t=data.t.begin(), v=data.val.begin(), e=data.err.begin(); 
      t != data.t.end(); t++,v++,e++) { 
    double fp=func(c,*t);
    f-=(*v-fp)*(*v-fp)/(2*(*e)*(*e));
  }
  return f;
}


//------------------------------------------------------------------------------

void Morse_model::minimum(const vector <double> & c, vector <double> & min) { 
  assert(c.size() >=4);
  min.resize(1);
  min[0]=c[3];
}


double Morse_model::func(const vector <double> & c, double t) { 
  assert(c.size() >=4 );
  double tmp=1-exp(-c[2]*(t-c[3]));
  return c[0]+c[1]*tmp*tmp;
}

void Morse_model::generate_guess(const Line_data & data, const Fix_information & fix,
                                   vector <double> & c) { 
  if(fix.enforce==1) { 
    c.resize(1);
    c[0]=50*(rng.ulec()-0.5);
  }
  else { 
    c.resize(4);
    c[0]=data.val[0];
    
    c[1]=10.0*(rng.ulec()-0.5);
    c[2]=10.0*(rng.ulec()-0.5);
    double avg=0;
    
    for(cdit_t t=data.t.begin(); t!= data.t.end(); t++) avg+=*t;
    avg/=data.t.size();
    c[3]=avg;
    //cout << "guessing " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << endl;

  }
  
  
}
