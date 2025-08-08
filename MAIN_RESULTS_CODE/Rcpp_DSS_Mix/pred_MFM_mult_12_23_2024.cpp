#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

/**
 * A helper function to compute the multivariate normal density at x,
 * with mean vector mu and covariance matrix Sigma (assumed positive-definite).
 */
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

/**
 * Replicates the R code:
 *
 *   expec.pred.post = matrix(0, sample_size, sample_size)
 *   for(j in 1:sample_size) {
 *     comp = which(matrix_weights[, j] != 0)
 *     for(i in 1:sample_size) {
 *       density.temp = numeric(length(comp))
 *       t = 1
 *       for(k in comp) {
 *         mu_vec    = matrix_mu[k, (p.dim*j - (p.dim-1)) : (p.dim*j)]
 *         sigma_row = matrix_sigma[k, ((p.dim^2)*j - (p.dim^2 - 1)) : ((p.dim^2)*j)]
 *         sigma_mat = matrix( sigma_row, p.dim, p.dim )
 *         prob      = matrix_weights[k, j]
 *         density.temp[t] = prob * mvtnorm::dmvnorm(y.pred[i,], mean=mu_vec, sigma=sigma_mat)
 *         t = t+1
 *       }
 *       expec.pred.post[j, i] = sum(density.temp)
 *     }
 *   }
 *
 * We assume everything is 0-based in C++ but the logic matches the R code.
 */
// [[Rcpp::export]]
arma::mat expec_pred_post_cpp(const arma::mat &matrix_weights,
                              const arma::mat &matrix_mu,
                              const arma::mat &matrix_sigma,
                              const arma::mat &y_pred,
                              const int sample_size,
                              const int p_dim)
{
  // Output: sample_size x sample_size
  arma::mat expec_pred_post(sample_size, sample_size, arma::fill::zeros);
  
  // n_components is the # of rows in matrix_weights
  // (assuming matrix_weights is n_components x sample_size)
  int n_components = matrix_weights.n_rows;
  
  // Loop over j in [0..sample_size-1]
  for(int j = 0; j < sample_size; j++)
  {
    // Find which components are nonzero for column j
    // (logical indexing: matrix_weights.col(j) != 0)
    arma::uvec comp = arma::find(matrix_weights.col(j) != 0);
    
    // For each i in [0..sample_size-1]
    for(int i = 0; i < sample_size; i++)
    {
      // We'll accumulate the sum of densities in a local variable
      double sum_density = 0.0;
      
      // For each k in comp
      for(arma::uword idx = 0; idx < comp.n_elem; idx++)
      {
        int k = comp[idx];
        
        // 1) Extract mu_vec:
        //    In R: mu_vec = matrix_mu[k, (p.dim*j - (p.dim-1)) : (p.dim*j)]
        //    In 0-based C++:
        //       start_col = p_dim*j
        //       end_col   = p_dim*j + (p_dim - 1)
        int start_col_mu = p_dim * j;
        int end_col_mu   = start_col_mu + p_dim - 1;
        
        // row k => row(k), subvec of columns [start_col_mu..end_col_mu]
        arma::rowvec mu_vec = matrix_mu.row(k).subvec(start_col_mu, end_col_mu);
        
        // 2) Extract sigma_mat:
        //    In R: sigma_mat = matrix( matrix_sigma[k, ((p.dim^2)*j - (p.dim^2-1)) : ((p.dim^2)*j)],
        //                              p.dim, p.dim )
        //    0-based C++:
        //       start_col_sig = p_dim^2 * j
        //       end_col_sig   = start_col_sig + (p_dim^2 - 1)
        int start_col_sig = p_dim * p_dim * j;
        int end_col_sig   = start_col_sig + (p_dim*p_dim - 1);
        
        arma::rowvec sigma_row = matrix_sigma.row(k).subvec(start_col_sig, end_col_sig);
        
        // Convert sigma_row into a p_dim x p_dim matrix, in column-major order
        arma::mat sigma_mat(p_dim, p_dim);
        // Fill sigma_mat column by column
        for(int col = 0; col < p_dim; col++) {
          for(int row = 0; row < p_dim; row++) {
            sigma_mat(row, col) = sigma_row[col * p_dim + row];
          }
        }
        
        // 3) Probability weight
        double prob = matrix_weights(k, j);
        
        // 4) Evaluate the density at y_pred.row(i)
        double val = prob * dmvnorm_cpp(y_pred.row(i), mu_vec, sigma_mat);
        
        // Accumulate
        sum_density += val;
      }
      
      expec_pred_post(j, i) = sum_density;
    }
  }
  
  return expec_pred_post;
}
