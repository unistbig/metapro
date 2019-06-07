#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h> /* exit, EXIT_FAILURE*/
#include <stdio.h>
#include <string>
#include <string.h>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <math.h> /*erf*/
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <boost/math/distributions/normal.hpp> // for normal_distribution
#include <boost/math/distributions.hpp>
//#include <boost/math/distributions/hypergeometric.hpp>
//#include <boost/math/distributions/binomial.hpp>
//#include <boost/math/policies/policy.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <ginac/ginac.h> // Symbolic integration
#include <boost/algorithm/string/replace.hpp>
#include <limits>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
//#include <gsl/gsl_errno.h>
#include <map>

using namespace std;
using std::setw;
using std::setprecision;
using std::numeric_limits;
using namespace GiNaC;
using namespace Rcpp;

void split(const string& str, vector<string>& tokens, const string& delimiters = " ")
{
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  string::size_type pos = str.find_first_of(delimiters, lastPos);
  while(string::npos!=pos || string::npos != lastPos){
    tokens.push_back(str.substr(lastPos, pos-lastPos));
    lastPos = str.find_first_not_of(delimiters, pos);
    pos = str.find_first_of(delimiters, lastPos);
  }
}

int eCount(vector<double> a){
  // a is sorted p-value vectors divided by p-value cutoff
  int count,k;
  k=a.size();
  count = 0;
  for(int i=0; i<k; i++)
  {
    if(a[i]<=1)
      count++;
    else
      return count;
  }
  return count;
}

template <typename T>
std::string to_string(T const& value) {
  stringstream sstr;
  sstr << value;
  return sstr.str();
}

template <typename type>
type signn(type value){
  return type((value>0)-(value<0));
}

template<typename SequenceT, typename Range1T, typename Range2T>
void replace_all(SequenceT & Input, const Range1T & Search,
                 const Range2T & Format);

bool myfunction(string i, string j) {
  return(i==j);
}

// [[Rcpp::export]]
List ordmeta_mono(NumericVector p)
{
  int nStudy;
  int j;
  double minP, P_temp, baseP, LB;
  vector<double> P;
  ex templete, formula;
  symbol x("x"), y("y");
  vector<int> index_na;
  P.clear();
  for(int i=0; i<p.size(); i++)
  {
    if(NumericVector::is_na(p[i]))
    {
      index_na.push_back(i+1);
      continue;
    }else
      P.push_back(p[i]);
  }
  int psize = P.size();
  NumericVector Psort(psize);
  copy(P.begin(), P.end(), Psort.begin());
  sort(Psort.begin(), Psort.end());

  nStudy = Psort.size();

  minP = 1;
  for(j=0; j<nStudy; j++)
  {
    P_temp = boost::math::ibeta(j+1, nStudy-j,Psort[j]);
    if(P_temp<minP){minP=P_temp;baseP=Psort[j];}
  }
  formula = 1;

  vector<int> numCohort;
  for(j=0;j<p.size();j++){
    if(p[j]<=baseP)
    {
      numCohort.push_back(j+1);
    }
  }
  for(j=0; j<(nStudy-1); j++)
  {
    LB = gsl_cdf_beta_Pinv(minP, j+1, nStudy-j); //Fi_inv(nStudy, j+1, minP);
    formula = formula.subs(y==x);
    templete=integral(x, LB, y, formula);
    formula = templete.eval_integ();
    formula = formula*(j+1);
  }
  LB = gsl_cdf_beta_Pinv(minP, nStudy,1); //Fi_inv(nStudy, nStudy, minP);

  formula = formula.subs(y==x);
  templete=integral(x, LB, 1, formula);
  formula = templete.eval_integ();
  formula = formula*nStudy;
  double D = ex_to<numeric>(formula).to_double();
  double RES=1.0-D;

  NumericVector NumCohort(numCohort.size());
  copy(numCohort.begin(),numCohort.end(),NumCohort.begin());

  NumericVector Index_NA(index_na.size());
  copy(index_na.begin(), index_na.end(), Index_NA.begin());
  List DATAFRAME;
  DATAFRAME["p"]=RES;
  DATAFRAME["eff.p.index"]= NumCohort;
  DATAFRAME["index.NA"] = Index_NA;
  return DATAFRAME;
}

// [[Rcpp::export]]
List ordmeta_signed(NumericVector p, NumericVector effect_size, bool sort_by_decreasing=true)
{
  NumericVector newP(p.size());
  for(int i=0; i<p.size(); i++)
  {
    if(NumericVector::is_na(p[i]))
      newP[i] = p[i];
    else
    {
      if(sort_by_decreasing)
      {
        if(effect_size[i]>=0){newP[i] = p[i]/2;}
        if(effect_size[i]<0){newP[i] = 1-p[i]/2;}
      }else{
        if(effect_size[i]>=0){newP[i] = 1-p[i]/2;}
        if(effect_size[i]<0){newP[i] = p[i]/2;}
      }
    }
  }

  List a=ordmeta_mono(newP);
  return a;
}
