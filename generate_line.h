#ifndef GENERATE_LINE_H_INCLUDED
#define GENERATE_LINE_H_INCLUDED
#include "Data_generator.h"
#include "Line_model.h"
#include "Sample.h"
#include "PES.h"
//-------------------------------------------------------
void find_minimum(double range, int maxpts, Data_generator & pes,
                  Line_model & mod, Line_data & data, Fit_info & finfo, 
                  vector <double> & sigma); 
Data_point gen_data(Data_generator & pes, const vector <double> & direction, 
                    const vector <double> & start_pos, const double t, double desired_sigma);

inline void append_number_fixed(string & str, int num){
  char strbuff[100];
  sprintf(strbuff, "%03d", num);
  str+=strbuff;
}

double rms_deviation(const Line_model & mod, const Line_data & data, const vector <double> & c);

#endif //GENERATE_LINE_H_INCLUDED

