#include <Rcpp.h>
using namespace Rcpp;

//Implement quantile(type=8) from R.

// [[Rcpp::export]]
NumericVector q8(NumericVector x, NumericVector probs) {
    int np = probs.size();
    int n = x.size();
    NumericVector nppm(np);
    IntegerVector j(np);
    NumericVector qs(np);
    NumericVector h(np);

    std::vector<double> xx = Rcpp::as<std::vector<double> >(x);
    
    double fuzz = 4*(std::numeric_limits<double>::denorm_min());
        
    for(int i = 0; i < np; i++) {
        nppm[i] = 1.0/3.0 + probs[i] * ((double)n + 1.0/3.0);
        j[i] = std::floor(nppm[i] + fuzz);
        h[i] = nppm[i] - (double)j[i];
        if(std::abs(h[i]) < fuzz) {
            h[i] = 0;
        }
    }
    
    std::sort(xx.begin(), xx.end());
    double minx = xx[0];
    double maxx = xx[xx.size()-1];
    xx.insert(xx.begin(),minx);
    xx.insert(xx.begin(),minx);
    xx.push_back(maxx);
    xx.push_back(maxx);
    for(int i = 0; i < np; i++) {
        qs[i] = xx[j[i]+1];
        if(abs(h[i]-1)<fuzz/4) qs[i] = xx[j[i] + 2];
        if(h[i]>0 & h[i]<1) {
            qs[i] = (1 - h[i]) * xx[j[i] + 1] + h[i] * xx[j[i] + 2];
        }
    }
    //return List::create(_["j"] = j, _["h"] = h, _["xx"] = xx, _["nppm"] = nppm, _["fuzz"] =fuzz , _["qs"] = qs); 
    return qs;    
}
// @ImportFrom Rcpp sourceCpp
// @useDynlib mqr
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
        qtmp = q8(tmp, q);
        for(int ii = 0; ii <qtmp.size(); ii++) {
           out(i,ii) = qtmp[ii]; 
        }
    }
    
    return out;
}



