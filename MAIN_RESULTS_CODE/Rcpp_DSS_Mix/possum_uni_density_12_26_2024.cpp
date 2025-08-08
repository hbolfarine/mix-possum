#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List computeDensityMclustArma(
    const arma::mat& y_fit, // N x D matrix
    const int K_sel          // Number of mixture components
) {
  int N = y_fit.n_rows;
  int D = y_fit.n_cols;
  
  // Initialize output matrices
  // Assuming the number of mixture components is K_sel
  arma::mat prob_post(N, K_sel, arma::fill::zeros);
  arma::mat mean_post(N, K_sel, arma::fill::zeros);
  arma::mat sigmasq_post(N, K_sel, arma::fill::zeros);
  arma::mat pred_fit(N, D, arma::fill::zeros); // Adjust columns as needed
  
  // Obtain the R densityMclust function from the mclust package
  Function densityMclust("mclust::densityMclust");
  
  for(int i = 0; i < N; ++i){
    // Extract the i-th row as a NumericVector
    NumericVector y_i = wrap(y_fit.row(i));
    
    // Call densityMclust with the specified parameters
    List dens = densityMclust(y_i, Named("G") = K_sel, 
                              Named("plot") = false, 
                              Named("verbose") = false, 
                              Named("modelNames") = "V");
    
    // Extract the density estimate
    NumericVector dens1 = dens["density"];
    
    // Extract parameters
    List parameters = dens["parameters"];
    NumericVector pro = parameters["pro"];         // Mixing proportions
    NumericVector mean_ = parameters["mean"];      // Means of the components
    List variance = parameters["variance"];
    NumericVector sigmasq = variance["sigmasq"];  // Variances (sigmas squared)
    
    // Assign extracted values to the output matrices
    prob_post.row(i) = arma::vec(pro.begin(), pro.size(), false);
    mean_post.row(i) = arma::vec(mean_.begin(), mean_.size(), false);
    sigmasq_post.row(i) = arma::vec(sigmasq.begin(), sigmasq.size(), false);
    pred_fit.row(i) = arma::vec(dens1.begin(), dens1.size(), false);
  }
  
  // Return the results as a list
  return List::create(
    Named("prob_post") = prob_post,
    Named("mean_post") = mean_post,
    Named("sigmasq_post") = sigmasq_post,
    Named("pred_fit") = pred_fit
  );
}
