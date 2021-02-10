#include <RcppArmadillo.h>

//' @name logLik
//' @title Log likelihood function for normal distribution with spatial
//' dependency.
//'
//' @param mDist distance matris;
//' @param s2 variance from the covariance model;
//' @param phi variance from the covariance model;
//' @param vDiff column vector of data (subtracted the mean vector)
//' @export

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::colvec logLik(const arma::mat& mDist, const arma::colvec& vDiff, double phi, double s2) {
  arma::mat mCov = s2 * exp(- mDist / phi);
  arma::colvec output = arma::zeros(1);

  output = 0.5 * log(det(mCov)) + 0.5 * trans(vDiff) * inv(mCov) * vDiff;

  return output;
}
