#include <Rcpp.h>
using namespace Rcpp;
const double fuzz = __DBL_EPSILON__;

//Implement quantile(type=8) from R.

// [[Rcpp::export]]
NumericVector q8(std::vector<double> x, NumericVector probs) {
    int np = probs.size();
    int n = x.size();
    NumericVector nppm(np);
    IntegerVector j(np);
    NumericVector qs(np);
    NumericVector h(np);
    
    nppm = 1.0/3.0 + probs * ((double)n + 1.0/3.0);
    j = floor(nppm + fuzz * 4.0);
    
    for(int i = 0; i < np; i++) {
        h[i] = nppm[i] - (double)j[i];
        if(std::abs(h[i]) < fuzz) {
            h[i] = 0.0;
        }
    }
    
    std::sort(x.begin(), x.end());
    double minx = x[0];
    double maxx = x[x.size()-1];
    x.insert(x.begin(),minx);
    x.insert(x.begin(),minx);
    x.push_back(maxx);
    x.push_back(maxx);
    for(int i = 0; i < np; i++) {
        qs[i] = x[j[i] + 1];
        if(std::abs(h[i] - 1) < fuzz) qs[i] = x[j[i] + 2];
        if(h[i] > 0 & h[i] < 1) {
            qs[i] = (1 - h[i]) * x[j[i] + 1] + h[i] * x[j[i] + 2];
        }
    }
    return qs;    
}

// Fast version of mq
// [[Rcpp::export]]
NumericMatrix mqC(NumericVector y, NumericVector q, IntegerVector l, NumericVector d) {
    int ny = y.size();
    int nq = q.size();
    int nd = d.size();
    int nl = l.size();
    int st = max(l);
    NumericVector qtmp(nq);
    
    NumericMatrix out(ny, nq);
    NumericVector tmp(nd);
    for(int i = 0; i < st; i++) {
        for(int ii = 0; ii< nq; ii ++) {
            out(i,ii) = NA_REAL;
        }
    }
    for(int i=st; i<ny; i++) {
        for(int j=0; j<nd; j++) {
            tmp[j] = d[j] * y[i - l[j]];
        }
        qtmp = q8(Rcpp::as<std::vector<double> >(tmp), q);
        for(int ii = 0; ii <qtmp.size(); ii++) {
           out(i,ii) = qtmp[ii]; 
        }
    }
    
    return out;
}



