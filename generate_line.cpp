#include "generate_line.h"
//----------------------------------------------------------------------


//------------------------------------------------------------------
Data_point gen_data(Data_generator & pes, const vector <double> & direction, 
                    const vector <double> & start_pos, const double t, double desired_sigma) { 
  int ndim=pes.ndim();
  assert(ndim==start_pos.size());
  assert(ndim==direction.size());
  Data_point pt;
  vector <double> x(ndim);
  for(int d=0; d< ndim; d++) { 
    x[d]=start_pos[d]+t*direction[d];
  }
  //double f=pes.eval_pes(x);
  double f,err;
  pes.eval(x,desired_sigma, f,err);
  pt.t=t;
  pt.val=f;
  pt.inverr=1.0/err;
  //cout << "data: " << pt.t << " " << pt.val << " " << err << endl;
  return pt;
}


//------------------------------------------------------------------
//Using a set of trial positions, estimate the RMS deviation.
double rms_deviation(const Line_model & mod, const Line_data & data, const vector <double> & c) { 
  double rms=0;
  int npoints=0;
  for(vector<Data_point>::const_iterator d=data.data.begin(); d!=data.data.end(); d++) { 
    double f=mod.func(c,d->t);
    rms+=(d->val-f)*(d->val-f);
    npoints++;
  }
  return sqrt(rms/npoints);
  
}

//------------------------------------------------------------------
//This returns true if the values are significantly different, fairly arbitrary, 
//but it should be several sigmas.
bool sig_diff(Data_point & p1, Data_point & p2) { 
  return fabs(p1.val-p2.val) > 4.0*sqrt(1.0/(p1.inverr*p1.inverr)+1.0/(p2.inverr*p2.inverr));
}

void find_minimum(double range, int maxpts, Data_generator & pes,
                  Line_model & mod, Line_data & data, Fit_info & finfo, 
                  vector <double> & sigma) { 
  assert(data.start_pos.size()==data.direction.size());
  assert(data.start_pos.size()==pes.ndim());
  assert(data.data.size()==0);
  assert(sigma.size()==data.direction.size());
  int ndim=data.direction.size();
  
  static int call_num=0;
  string logname="find_minimum";
  append_number_fixed(logname, call_num);
  ofstream log(logname.c_str());
  log.precision(15);
  double start_sigma=0;;
  for(int d=0; d< ndim; d++) start_sigma+=data.direction[d]*data.direction[d]*sigma[d];
  double start_t=-range;
  double end_t=range;
  Line_data tempdata=data;
  Data_point p1,p2,p3,p4;
  double golden=1.0; //(1+sqrt(5.0))/2.0;
  double base_t=0;
  //Let's do a search to find the minimum
  p1=gen_data(pes,data.direction,data.start_pos,start_t,start_sigma);
  p2=gen_data(pes,data.direction,data.start_pos,start_t+(end_t-start_t)/(1.0+golden),start_sigma);
  p3=gen_data(pes,data.direction,data.start_pos,end_t,start_sigma);
  vector <Data_point> allpts;
  allpts.push_back(p1); allpts.push_back(p2); allpts.push_back(p3);
  log << "#starting sigma " << start_sigma << endl;
  //First bracket the minimum.
  while(1) { 
    log << "#step! points at  " << p1.t << " " << p2.t << " " << p3.t << endl;
    if(sig_diff(p1,p2) && sig_diff(p2,p3) && p2.val < p1.val && p2.val < p3.val) { 
      log << "#found set of bracketing points " << p1.t << " " << p2.t << " " << p3.t << endl;
      base_t=p2.t;
      break;
    }
    else if(!sig_diff(p1,p2) || !sig_diff(p2,p3)) { 
      //First check to see if we're just on opposite sides of the minimum and recenter if so:
      int shifted_minimum=0;
      if(!sig_diff(p1,p2) && sig_diff(p2,p3)) { 
        log << "#checking for minimum between " << p1.t << " and " << p2.t << endl;

        p4=gen_data(pes,data.direction,data.start_pos,(p1.t+p2.t)/2.0,start_sigma);
        allpts.push_back(p4);
        if(sig_diff(p4,p1) && sig_diff(p4,p2)) {
          start_t=p4.t-range;
          end_t=p4.t+range;
          p2=p4;
          shifted_minimum=1;
        }
      }
      else if(sig_diff(p1,p2) && !sig_diff(p2,p3)) { 
        log << "#checking for minimum between " << p2.t << " and " << p3.t << endl;
        p4=gen_data(pes,data.direction,data.start_pos,(p2.t+p3.t)/2.0,start_sigma);
        allpts.push_back(p4);
        if(sig_diff(p4,p3) && sig_diff(p4,p2)) {
          start_t=p4.t-range;
          end_t=p4.t+range;
          p2=p4;
          shifted_minimum=1;
        }        
      }
      if(!shifted_minimum) { 
        start_sigma*=0.5;
        log << "#reducing sigma to " << start_sigma << " in bracketing routine \n";
        p2=gen_data(pes,data.direction,data.start_pos,start_t+(end_t-start_t)/(1.0+golden),start_sigma);
        allpts.push_back(p2);
      }
      p1=gen_data(pes,data.direction,data.start_pos,start_t,start_sigma);
      allpts.push_back(p1);
      p3=gen_data(pes,data.direction,data.start_pos,end_t,start_sigma);
      allpts.push_back(p3);
    }
    else {  
      if(p1.val < p2.val) {
        start_t-=range;
        end_t-=range;
        p4=gen_data(pes,data.direction,data.start_pos,start_t,start_sigma);
        allpts.push_back(p4);
        p3=p2;
        p2=p1;
        p1=p4;
      }
      if( p3.val < p2.val) {
        end_t+=range;
        start_t+=range;
        p4=gen_data(pes,data.direction,data.start_pos,end_t,start_sigma);
        allpts.push_back(p4);
        p1=p2;
        p2=p3;
        p3=p4;
      }
      log << "#changing range to " << start_t << " to " << end_t << endl;
    }
  }
  
  log << "#all points used in the bracketing " << endl;
  for(vector<Data_point>::iterator d=allpts.begin(); d!= allpts.end(); d++) { 
    log << d->t << " " << d->val << "  " << 1.0/d->inverr << endl;
  }
  log << endl;

  log << "#minimum value point at t= " << p2.t << " val= " << p2.val << " +/- " << 1.0/p2.inverr << endl;
  //Now that we've bracketed the minimum, we can fill in the data a little bit 
  //for the fit.  Note that we may throw away a little data potentially, but 
  //including it seems more trouble than it's worth and will slow down the Monte
  //Carlo later.
  data.data.clear();
  data.data.push_back(p1); data.data.push_back(p2); data.data.push_back(p3);
  data.data.push_back(gen_data(pes,data.direction,data.start_pos,p2.t+0.6*range,start_sigma));
  data.data.push_back(gen_data(pes,data.direction,data.start_pos,p2.t-0.6*range,start_sigma));
  data.data.push_back(gen_data(pes,data.direction,data.start_pos,p2.t+0.8*range,start_sigma));
  data.data.push_back(gen_data(pes,data.direction,data.start_pos,p2.t-0.8*range,start_sigma));
  data.data.push_back(gen_data(pes,data.direction,data.start_pos,p2.t+0.3*range,start_sigma));
  data.data.push_back(gen_data(pes,data.direction,data.start_pos,p2.t-0.3*range,start_sigma));
  vector <Walker> allwalkers;
  sample(mod, data, finfo, allwalkers);
  base_t=finfo.min[0];
  
  double rms=rms_deviation(mod, data,finfo.cavg);
  double avgerr=0;
  int npoints=0;
  for(vector<Data_point>::const_iterator d=data.data.begin(); d!=data.data.end(); d++) { 
    avgerr+= 1.0/d->inverr;
    npoints++;
  }
  avgerr/=npoints;

  log << "#rms deviation " << rms <<  " average uncertainty " << avgerr << endl;
  
  //Check to see how many points are outside of 3 error bars..
  //It should actually happen occasionally, so we should be careful not to 
  //freak too much if it does.
  int noutside=0;
  for(vector<Data_point>::const_iterator d=data.data.begin(); d!=data.data.end(); d++) { 
    double f=mod.func(finfo.cavg,d->t);
    double deviation=fabs(f-d->val);
    double err=1.0/d->inverr;
    if(deviation > 3*err) {
      log << "#outside case!  t= " << d->t << " data point " << d->val << " +/- " << err 
      << " fitted value " << f << endl;
    }
    noutside++;
  }
  
  
  log << "#data " << endl;
  for(vector<Data_point>::const_iterator d=data.data.begin(); d!=data.data.end(); d++) { 
    log << d->t << " " << d->val <<  " " << 1.0/d->inverr << endl;
  }
  
  log << endl << "#average parameter fit " << endl;
  base_t=finfo.min[0];
  double begin=base_t-range*1.2, end=base_t+range*1.2;
  double res=(end-begin)/400.0;
  for(double t=begin; t < end; t+=res) { 
    log << t << " " << mod.func(finfo.cavg, t) << " 0.0 " << endl;
  }
  
  /*
  int nrandom=10;
  for(int rand=0; rand < nrandom; rand++) { 
    log << endl << "#random sampling of fits " << endl;
    vector <double> ctmp=finfo.cavg;
    int nparm=finfo.cavg.size();
    for(int i=0; i< nparm; i++) ctmp[i]+=rng.gasdev()*finfo.cer[i];
    for(double t=begin; t < end; t+=res) { 
      log << t << " " << mod.func(ctmp, t) << " 0.0 " << endl;
    }
  }*/
  

  call_num++;
}

//------------------------------------------------------------------


