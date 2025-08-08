#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double dmvNorm(const arma::vec &x, const arma::vec &mu, const arma::mat &Sigma){
  // Dimensions
  int k = x.n_elem;
  
  // Compute the determinant and inverse of Sigma
  double detSigma = arma::det(Sigma);
  arma::mat invSigma = arma::inv(Sigma);
  
  // Ensure Sigma is positive definite
  if (detSigma <= 0) {
    Rcpp::stop("Covariance matrix is not positive definite");
  }
  
  // Compute the Mahalanobis distance
  arma::vec diff = x - mu;
  double mahalanobis = arma::as_scalar(diff.t() * invSigma * diff);
  
  // Compute the density
  double density = std::exp(-0.5 * mahalanobis) /
    std::sqrt(std::pow(2 * M_PI, k) * detSigma);
  
  return density;
}

// [[Rcpp::export]]
// arma::mat pred_func(const int &samp, const int &p, const arma::mat &y, const arma::mat &W, const arma::mat &M, const arma::mat &S){
int pred_func(const int &samp, const int &p, const arma::mat &y, const arma::mat &W, const arma::mat &M, const arma::mat &S){  
  
  arma::mat P(samp,samp,arma::fill::none);
  double sum_v;
  arma::mat Sigma;
  double prob;
  double dens_temp = 0.0;
  //arma::vec mu = {2,3};

  for(int j = 0; j < (samp-1); ++j){
    arma::uvec comp = arma::find(W.col(j) != 0);
    std::cout << "\n comp:" << comp;
    for(int i = 0; i < (samp-1); ++i){
      arma::vec dens(comp.size(),arma::fill::none);
      int t = 0;
      for(int k:comp){
        arma::rowvec mu = M.row(k).cols((p*(j+1)-(p-1))-1,(p*(j+1))-1);
        arma::rowvec vtemp = S.row(k).cols((p*p*(j+1)-(p*p)),p*p*(j+1)-1);
        Sigma = arma::reshape(vtemp,p,p);
        prob = W(k,j);
    
        dens_temp = prob*dmvNorm(y.row(j).t(), mu.t(), Sigma);
        //dens = dmvNorm(y_temp, mu, Sigma);
        //std::cout << "Mat:" << Sigma;
        std::cout << "\n k:" << k;
        //std::cout << "y:" << y.row(j);

        dens(t) = dens_temp;
        t = t+1;
      }
      std::cout << "\n END" << std::endl;;
      sum_v = arma::sum(dens);
      //std::cout << "dens:\n" << dens;
      //std::cout << "Sum_dens:\n" << sum_v;
      //double sum_v = arma::sum(dens);
      P(j,i) = sum_v;
    }
  } 
  std::cout << "P_matrix:" << P;
  return 0;
}

// [[Rcpp::export]]
arma::rowvec sub_mat(const int row, const int low, const int up, const arma::mat&M){
  
  arma::rowvec subRow = M.row(row).cols(low, up);
  
  return subRow;
}

// [[Rcpp::export]]
int sub_mat_test(const int &k, const int &j, const int &p, const arma::mat &S, const arma::mat &W){
  arma::uvec comp = arma::find(W.col(j) != 0);
  arma::rowvec vtemp = S.row(k).cols((p*p*(j+1)-(p*p)),p*p*(j+1)-1);
  arma::mat Sigma = arma::reshape(vtemp,p,p);
  double sum_v = arma::sum(vtemp);
  std::cout << "Vector:" << vtemp;
  //std::cout << "k:" << k;
  std::cout << "sum:" << sum_v;
  std::cout << "lower:" << (p*p*(j+1)-(p*p));
  std::cout << "upper:" << p*p*(j+1)-1;
  std::cout << "upper:" << Sigma;
  return 0;
}

// [[Rcpp::export]]
double gen_norm(double prob, const arma::vec &y, const arma::vec &mu, const arma::mat Sigma){
  double dens;
  dens = prob*dmvNorm(y,mu,Sigma);
  return dens;
}


  //arma::vec mu = {2,3};
  //Sigma = {{1,0},{0,1}};
  //arma::vec mu = {2,3};        // mean vector (2D)
  //arma::mat Sigma = {{1,0},{0,1}};  // covariance matrix (2×2 identity)
  //arma::vec y_temp = {0,2};  

