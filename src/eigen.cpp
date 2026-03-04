#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;               // 'maps' rather than copies
using Eigen::MatrixXd;          // variable size matrix, double precision
using Eigen::VectorXd;          // variable size vector, double precision
using Eigen::EigenSolver;       // general solver -- works on non-symmetric matrices

// [[Rcpp::export]]
VectorXd getEigenValues(Map<MatrixXd> M) {
    // computeEigenvectors=false: faster, only eigenvalues needed
    EigenSolver<MatrixXd> es(M, false);
    // eigenvalues() returns VectorXcd (complex); .real() extracts real parts as VectorXd
    // For IPM kernels (non-negative matrices), dominant eigenvalue is real by Perron-Frobenius
    return es.eigenvalues().real();
}
