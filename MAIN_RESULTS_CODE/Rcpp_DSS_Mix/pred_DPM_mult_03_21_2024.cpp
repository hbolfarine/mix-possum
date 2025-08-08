#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

double dmvnorm_cpp(const arma::rowvec &x,
                   const arma::rowvec &mu,
                   const arma::mat &Sigma) {
  // Dimension
  int k = x.n_elem;
  
  // Determinant and inverse of Sigma
  double detSigma = arma::det(Sigma);
  arma::mat invSigma = arma::inv(Sigma);
  
  if (detSigma <= 0) {
    stop("Covariance matrix is not positive definite");
  }
  
  // Mahalanobis distance
  arma::rowvec diff = x - mu;
  double mahalanobis = arma::as_scalar(diff * invSigma * diff.t());
  
  // Multivariate normal density
  double density = std::exp(-0.5 * mahalanobis) /
    std::sqrt(std::pow(2.0 * M_PI, (double) k) * detSigma);
  
  return density;
}

// density.temp <- numeric(kmax)
//   expec.pred.post = matrix(0,sample_size,sample_size)
//   for(j in 1:sample_size){
//     for(i in 1:sample_size){
//       for(k in 1:kmax){
//         density.temp[k] <- Eta[j,k]*dmvnorm(y.pred[i,],mean = Mu[j,,k],sigma = Sigma[j,,,k])
//       }
//       expec.pred.post[j,i] = sum(density.temp)
//     }
// # cat(round(j/R,3),"\n")
// # cat(j)
//   }

// [[Rcpp::export]]
arma::mat expect_prob_SFM(double sample_size, 
                    int kmax,
                    const int p,
                const arma::mat &y_pred,
                const arma::mat &Eta, 
                const arma::cube &Mu,
                const arma::cube &S){
  
  arma::mat expec_pred_post(sample_size, sample_size, arma::fill::zeros);
  arma::vec dens(kmax, arma::fill::zeros);
  arma::mat Sigma;
  double prob;
  
  for(int j = 0; j < sample_size; j++){
    for(int i = 0 ; i < sample_size; i++){
      double sum_density = 0.0;
      for(int k = 0; k < kmax; k++){
        prob = Eta(j,k);
        arma::rowvec mu = Mu.slice(k).row(j);
        arma::rowvec vtemp = S.slice(k).row(j);
        Sigma = arma::reshape(vtemp,p,p);
        double val = prob*dmvnorm_cpp(y_pred.row(i), mu, Sigma);
        // std::cout << "\n prob:" << prob;
        // std::cout << "\n k:" << k;
        // std::cout << "\n MAT:" << Sigma;
        sum_density += val;
      }
      expec_pred_post(j,i) = sum_density;
    }
  }
  return expec_pred_post;
}