#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::VectorXd getEigenValues(Eigen::Map<Eigen::MatrixXd> M) {
    // computeEigenvectors=false: faster, only eigenvalues needed
    Eigen::EigenSolver<Eigen::MatrixXd> es(M, false);
    // eigenvalues() returns VectorXcd (complex); .real() extracts real parts as VectorXd
    // For IPM kernels (non-negative matrices), dominant eigenvalue is real by Perron-Frobenius
    return es.eigenvalues().real();
}
