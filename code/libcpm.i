%module libcpm
%{
#include <iostream>
#include <vector>

using namespace std;

extern vector< vector<double> > linear_least_squares_wrap(vector< vector<double> >);
%}

%include "std_vector.i"
namespace std {
  %template(VecDouble) vector<double>;
  %template(VecVecdouble) vector< vector<double> >;
}

extern std::vector< std::vector<double> > linear_least_squares_wrap (std::vector< std::vector<double> >);
