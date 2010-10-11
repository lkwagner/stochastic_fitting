#include "Data_generator.h"

inline void append_number(string & str, int num){
  char strbuff[100];
  sprintf(strbuff, "%d", num);
  str+=strbuff;
}


MonteCarlo_caller::MonteCarlo_caller(istream & is) { 
  is >> n >> estimated_rms >> command >> force_command;
  min_blocks=10;
}

void MonteCarlo_caller::eval(const vector <double> & x, const double desired_err, 
                             double & f, double & err) { 
  double tol=0.5; //accept this fractional error in the error..

  ofstream out("mc_input");
  for(vector<double>::const_iterator i=x.begin(); i!=x.end(); i++) out << *i << " ";
  out << endl;
  out.close();
  int it=0;
  int maxit=3;
  do { 
    string run_command=command+" ";
    int nblocks=max(min_blocks, int(estimated_rms*estimated_rms/(desired_err*desired_err)+0.9999));
    cout << "run try " << it << " desired error " << desired_err 
    <<  " estimated rms " << estimated_rms << " nblocks " << nblocks;
    cout.flush();
    append_number(run_command,nblocks);
    if(system(run_command.c_str())) { 
      cout << "Error evaluating function!" << endl;
      exit(1);
    }
    
    ifstream in("mc_output");
    in >> f >> err;
    in.close();
  
    cout << " actual error " << err << endl;
    estimated_rms=err*sqrt(double(nblocks));
    it++;
  } while(it < maxit && (err-desired_err)/desired_err > tol ); 
  //note that we don't care if the err is smaller than desired, only larger
}

//----------------------------------------------------------------------
int MonteCarlo_caller::gradient(const vector <double>&x, 
    vector<double>&grad,vector<double>&err) { 
  if(force_command=="noforce") return 0;
  ofstream out("mc_input");
  for(vector<double>::const_iterator i=x.begin(); i!=x.end(); i++) out << *i << " ";
  out << endl;
  out.close();

  if(system(force_command.c_str())) { 
    cout << "Error evaluating gradient!" << endl;
    exit(1);
  }
  ifstream in("mc_output");
  grad.resize(x.size()); err.resize(x.size());
  for(vector<double>::iterator i=grad.begin(); i!=grad.end(); i++) in >> *i;
  for(vector<double>::iterator i=err.begin(); i!= err.end(); i++) in >> *i;
  in.close();


}


//#########################################################################
